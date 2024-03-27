# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Removes interfaces and functions we're not using in C-PAC
#     * Removes Python 2 imports
#     * Adjusts imports for C-PAC
#     * Modifies logic for S3 with global `RETRY` and `RETRY_WAIT`
#     * Logs a debugging message instead of crashing if src == dst for copyfile
#     * Explicitly lowercases "s3"
#     * Handles empty file lists
#     * Docstrings updated accordingly
#     * Style modifications

# ORIGINAL WORK'S ATTRIBUTION NOTICE:
#     Copyright (c) 2009-2016, Nipype developers

#     Licensed under the Apache License, Version 2.0 (the "License");
#     you may not use this file except in compliance with the License.
#     You may obtain a copy of the License at

#         http://www.apache.org/licenses/LICENSE-2.0

#     Unless required by applicable law or agreed to in writing, software
#     distributed under the License is distributed on an "AS IS" BASIS,
#     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#     See the License for the specific language governing permissions and
#     limitations under the License.

#     Prior to release 0.12, Nipype was licensed under a BSD license.

# Modifications Copyright (C) 2019-2024  C-PAC Developers

# This file is part of C-PAC.
"""Interface that allow interaction with data.

Currently only available interface is:

    DataSink: Generic named output from interfaces to data store

Modified from https://github.com/nipy/nipype/blob/f64bf33/nipype/interfaces/io.py
"""

import os
import shutil
from shutil import SameFileError
import time

from nipype import config
from nipype.interfaces.base import isdefined
from nipype.interfaces.io import (
    copytree,
    DataSink as NipypeDataSink,
    ProgressPercentage,
)
from nipype.utils.filemanip import copyfile, ensure_list
from nipype.utils.misc import str2bool

from CPAC.utils.monitoring import FMLOGGER, IFLOGGER

RETRY = 5
RETRY_WAIT = 5


def _get_head_bucket(s3_resource, bucket_name):
    """Try to get the header info of a bucket to check if it exists & its permissions."""
    import botocore

    # Try fetch the bucket with the name argument
    err_msg = None
    for _ in range(RETRY):
        try:
            s3_resource.meta.client.head_bucket(Bucket=bucket_name)
            return

        except botocore.exceptions.ClientError as exc:
            error_code = int(exc.response["Error"]["Code"])
            if error_code == 403:  # noqa: PLR2004
                err_msg = (
                    "Access to bucket: %s is denied; check credentials" % bucket_name
                )
                break
            if error_code == 404:  # noqa: PLR2004
                err_msg = (
                    "Bucket: %s does not exist; check spelling and try "
                    "again" % bucket_name
                )
                break
            err_msg = "Unable to connect to bucket: %s. Error message:\n%s" % (
                bucket_name,
                exc,
            )

        except Exception as exc:
            err_msg = "Unable to connect to bucket: %s. Error message:\n%s" % (
                bucket_name,
                exc,
            )

        time.sleep(RETRY_WAIT)

    if err_msg is not None:
        raise Exception(err_msg)


class DataSink(NipypeDataSink):  # noqa: D101
    def _check_s3_base_dir(self):
        # Init variables
        s3_str = "s3://"
        bucket_name = "<N/A>"
        base_directory = self.inputs.base_directory

        if not isdefined(base_directory):
            s3_flag = False
            return s3_flag, bucket_name

        # Explicitly lower-case the "s3"
        if base_directory.lower().startswith(s3_str):
            base_dir_sp = base_directory.split("/")
            base_dir_sp[0] = base_dir_sp[0].lower()
            base_directory = "/".join(base_dir_sp)

        # Check if 's3://' in base dir
        if base_directory.startswith(s3_str):
            # Expects bucket name to be 's3://bucket_name/base_dir/..'
            bucket_name = base_directory.split(s3_str)[1].split("/")[0]
            s3_flag = True
        # Otherwise it's just a normal datasink
        else:
            s3_flag = False

        # Return s3_flag
        return s3_flag, bucket_name

    _check_s3_base_dir.__doc__ == NipypeDataSink._check_s3_base_dir.__doc__

    def _fetch_bucket(self, bucket_name):
        # Import packages
        try:
            import boto3
            import botocore
        except ImportError:
            err_msg = "Boto3 package is not installed - install boto3 and try again."
            raise Exception(err_msg)

        # Init variables
        creds_path = self.inputs.creds_path

        # Get AWS credentials
        try:
            aws_access_key_id, aws_secret_access_key = self._return_aws_keys()
        except Exception as exc:
            err_msg = (
                "There was a problem extracting the AWS credentials "
                "from the credentials file provided: %s. Error:\n%s" % (creds_path, exc)
            )
            raise Exception(err_msg)

        # Try and get AWS credentials if a creds_path is specified
        if aws_access_key_id and aws_secret_access_key:
            # Init connection
            IFLOGGER.info(
                "Connecting to S3 bucket: %s with credentials...", bucket_name
            )
            # Use individual session for each instance of DataSink
            # Better when datasinks are being used in multi-threading, see:
            # http://boto3.readthedocs.org/en/latest/guide/resources.html#multithreading
            session = boto3.session.Session(
                aws_access_key_id=aws_access_key_id,
                aws_secret_access_key=aws_secret_access_key,
            )

        else:
            IFLOGGER.info("Connecting to S3 bucket: %s with IAM role...", bucket_name)

            # Lean on AWS environment / IAM role authentication and authorization
            session = boto3.session.Session()

        s3_resource = session.resource("s3", use_ssl=True)

        # And try fetch the bucket with the name argument
        try:
            _get_head_bucket(s3_resource, bucket_name)
        except Exception:
            # Try to connect anonymously
            s3_resource.meta.client.meta.events.register(
                "choose-signer.s3.*", botocore.handlers.disable_signing
            )

            IFLOGGER.info("Connecting to AWS: %s anonymously...", bucket_name)
            _get_head_bucket(s3_resource, bucket_name)

        # Explicitly declare a secure SSL connection for bucket object
        # Return the bucket
        return s3_resource.Bucket(bucket_name)

    _fetch_bucket.__doc__ = NipypeDataSink._fetch_bucket.__doc__

    def _upload_to_s3(self, bucket, src, dst):
        """Upload outputs to S3 bucket instead of on local disk."""
        # Import packages
        import hashlib
        import os

        from botocore.exceptions import ClientError

        # Init variables
        s3_str = "s3://"
        s3_prefix = s3_str + bucket.name

        # Explicitly lower-case the "s3"
        if dst[: len(s3_str)].lower() == s3_str:
            dst = s3_str + dst[len(s3_str) :]

        # If src is a directory, collect files (this assumes dst is a dir too)
        if os.path.isdir(src):
            src_files = []
            for root, dirs, files in os.walk(src):
                src_files.extend([os.path.join(root, fil) for fil in files])
            # Make the dst files have the dst folder as base dir
            dst_files = [os.path.join(dst, src_f.split(src)[1]) for src_f in src_files]
        else:
            src_files = [src]
            dst_files = [dst]

        # Iterate over src and copy to dst
        for src_idx, src_f in enumerate(src_files):
            # Get destination filename/keyname
            dst_f = dst_files[src_idx]
            dst_k = dst_f.replace(s3_prefix, "").lstrip("/")

            # See if same file is already up there
            try:
                dst_obj = bucket.Object(key=dst_k)
                dst_md5 = dst_obj.e_tag.strip('"')

                # See if same file is already there
                src_read = open(src_f, "rb").read()
                src_md5 = hashlib.md5(src_read).hexdigest()
                # Move to next loop iteration
                if dst_md5 == src_md5:
                    FMLOGGER.info("File %s already exists on S3, skipping...", dst_f)
                    continue
                FMLOGGER.info("Overwriting previous S3 file...")

            except ClientError:
                FMLOGGER.info("New file to S3")

            # Copy file up to S3 (either encrypted or not)
            FMLOGGER.info(
                "Uploading %s to S3 bucket, %s, as %s...", src_f, bucket.name, dst_f
            )
            if self.inputs.encrypt_bucket_keys:
                extra_args = {"ServerSideEncryption": "AES256"}
            else:
                extra_args = {}

            retry_exc = None
            for _ in range(RETRY):
                try:
                    bucket.upload_file(
                        src_f,
                        dst_k,
                        ExtraArgs=extra_args,
                        Callback=ProgressPercentage(src_f),
                    )
                    break
                except Exception as exc:
                    time.sleep(RETRY_WAIT)
                    retry_exc = exc

            if retry_exc is not None:
                raise retry_exc

    # List outputs, main run routine
    def _list_outputs(self):
        # Init variables
        outputs = self.output_spec().get()
        out_files = []
        # Use hardlink
        use_hardlink = str2bool(config.get("execution", "try_hard_link_datasink"))

        # Set local output directory if specified
        if isdefined(self.inputs.local_copy):
            outdir = self.inputs.local_copy
        else:
            outdir = self.inputs.base_directory
            # If base directory isn't given, assume current directory
            if not isdefined(outdir):
                outdir = "."

        # Check if base directory reflects S3 bucket upload
        s3_flag, bucket_name = self._check_s3_base_dir()
        if s3_flag:
            s3dir = self.inputs.base_directory
            # If user overrides bucket object, use that
            if self.inputs.bucket:
                bucket = self.inputs.bucket
            # Otherwise fetch bucket object using name
            else:
                try:
                    bucket = self._fetch_bucket(bucket_name)
                # If encountering an exception during bucket access, set output
                # base directory to a local folder
                except Exception as exc:
                    s3dir = "<N/A>"
                    if not isdefined(self.inputs.local_copy):
                        local_out_exception = os.path.join(
                            os.path.expanduser("~"), "s3_datasink_" + bucket_name
                        )
                        outdir = local_out_exception
                    # Log local copying directory
                    FMLOGGER.info(
                        "Access to S3 failed! Storing outputs locally at: "
                        "%s\nError: %s",
                        outdir,
                        exc,
                    )
        else:
            s3dir = "<N/A>"

        # If container input is given, append that to outdir
        if isdefined(self.inputs.container):
            outdir = os.path.join(outdir, self.inputs.container)
            s3dir = os.path.join(s3dir, self.inputs.container)

        # If sinking to local folder
        if outdir != s3dir:
            outdir = os.path.abspath(outdir)
            # Create the directory if it doesn't exist
            if not os.path.exists(outdir):
                try:
                    os.makedirs(outdir)
                except OSError as inst:
                    if "File exists" in inst.strerror:
                        pass
                    else:
                        raise (inst)

        # Iterate through outputs attributes {key : path(s)}
        for key, _files in list(self.inputs._outputs.items()):
            files = _files
            if not isdefined(files):
                continue
            IFLOGGER.debug("key: %s files: %s", key, str(files))
            files = ensure_list(files if files else [])
            tempoutdir = outdir
            if s3_flag:
                s3tempoutdir = s3dir
            for d in key.split("."):
                if d[0] == "@":
                    continue
                tempoutdir = os.path.join(tempoutdir, d)
                if s3_flag:
                    s3tempoutdir = os.path.join(s3tempoutdir, d)

            # flattening list
            if files and isinstance(files, list):
                if isinstance(files[0], list):
                    files = [item for sublist in files for item in sublist]

            # Iterate through passed-in source files
            for _src in ensure_list(files):
                src = _src
                # Format src and dst files
                src = os.path.abspath(src)
                if not os.path.isfile(src):
                    src = os.path.join(src, "")
                dst = self._get_dst(src)
                if s3_flag:
                    s3dst = os.path.join(s3tempoutdir, dst)
                    s3dst = self._substitute(s3dst)
                dst = os.path.join(tempoutdir, dst)
                dst = self._substitute(dst)
                path, _ = os.path.split(dst)

                # If we're uploading to S3
                if s3_flag:
                    self._upload_to_s3(bucket, src, s3dst)
                    out_files.append(s3dst)
                # Otherwise, copy locally src -> dst
                if not s3_flag or isdefined(self.inputs.local_copy):
                    # Create output directory if it doesn't exist
                    if not os.path.exists(path):
                        try:
                            os.makedirs(path)
                        except OSError as inst:
                            if "File exists" in inst.strerror:
                                pass
                            else:
                                raise (inst)
                    try:
                        # If src == dst, it's already home
                        if (not os.path.exists(dst)) or (os.stat(src) != os.stat(dst)):
                            # If src is a file, copy it to dst
                            if os.path.isfile(src):
                                FMLOGGER.debug(f"copyfile: {src} {dst}")
                                copyfile(
                                    src,
                                    dst,
                                    copy=True,
                                    hashmethod="content",
                                    use_hardlink=use_hardlink,
                                )
                            # If src is a directory, copy
                            # entire contents to dst dir
                            elif os.path.isdir(src):
                                if os.path.exists(dst) and self.inputs.remove_dest_dir:
                                    FMLOGGER.debug("removing: %s", dst)
                                    shutil.rmtree(dst)
                                FMLOGGER.debug("copydir: %s %s", src, dst)
                                copytree(src, dst)
                                out_files.append(dst)
                    except SameFileError:
                        FMLOGGER.debug(f"copyfile (same file): {src} {dst}")

        # Return outputs dictionary
        outputs["out_file"] = out_files

        return outputs

    _list_outputs.__doc__ = NipypeDataSink._list_outputs.__doc__


DataSink.__doc__ = NipypeDataSink.__doc__
__all__ = ["DataSink"]
