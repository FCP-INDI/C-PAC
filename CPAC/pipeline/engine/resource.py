# Copyright (C) 2021-2024  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""Resources and ResourcePools for C-PAC."""

import ast
from collections.abc import KeysView
import copy
from itertools import chain
import os
from pathlib import Path
import re
from types import NoneType
from typing import Any, Literal, NamedTuple, Optional, overload
import warnings

from nipype.interfaces import utility as util  # type: ignore [import-untyped]
from nipype.interfaces.utility import Rename  # type: ignore [import-untyped]

from CPAC.image_utils.spatial_smoothing import spatial_smoothing
from CPAC.image_utils.statistical_transforms import (
    fisher_z_score_standardize,
    z_score_standardize,
)
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.check_outputs import ExpectedOutputs
from CPAC.pipeline.utils import MOVEMENT_FILTER_KEYS, name_fork, source_set
from CPAC.registration.registration import transform_derivative
from CPAC.resources.templates.lookup_table import lookup_identifier
from CPAC.utils.bids_utils import res_in_filename
from CPAC.utils.configuration.configuration import Configuration, Preconfiguration
from CPAC.utils.datasource import (
    calc_delta_te_and_asym_ratio,
    check_for_s3,
    check_func_scan,
    create_anat_datasource,
    create_fmap_datasource,
    create_general_datasource,
    gather_echo_times,
    get_fmap_phasediff_metadata,
    get_rest,
    resolve_resolution,
)
from CPAC.utils.interfaces.datasink import DataSink
from CPAC.utils.interfaces.function import Function
from CPAC.utils.monitoring import (
    getLogger,
    WFLOGGER,
)
from CPAC.utils.outputs import Outputs
from CPAC.utils.typing import LIST_OF_LIST_OF_STR, PIPE_IDX
from CPAC.utils.utils import (
    check_prov_for_regtool,
    create_id_string,
    get_last_prov_entry,
    get_scan_params,
    read_json,
    write_output_json,
)

EXTS = [".nii", ".gz", ".mat", ".1D", ".txt", ".csv", ".rms", ".tsv"]
POOL_DICT = dict[str | tuple, "STRAT_DICT"]
STRAT_DICT = dict[str | tuple, "dict | Resource"]


class DataPaths:
    """Store subject-session specific data paths."""

    def __init__(
        self, *, data_paths: Optional[dict] = None, part_id: Optional[str] = ""
    ) -> None:
        """Initialize a ``DataPaths`` instance."""
        if not data_paths:
            data_paths = {}
        if part_id and "part_id" in data_paths and part_id != data_paths["part_id"]:
            WFLOGGER.warning(
                "both 'part_id' (%s) and data_paths['part_id'] (%s) provided. "
                "Using '%s'.",
                part_id,
                data_paths["part_id"],
                part_id,
            )
        anat: dict[str, str] | str = data_paths.get("anat", {})
        if isinstance(anat, str):
            anat = {"T1": anat}
        self.anat: dict[str, str] = anat
        self.creds_path: Optional[str] = data_paths.get("creds_path")
        self.fmap: Optional[dict] = data_paths.get("fmap")
        self.func: dict[str, dict[str, str | dict]] = data_paths.get("func", {})
        self.part_id: str = data_paths.get("subject_id", "")
        self.site_id: str = data_paths.get("site_id", "")
        self.ses_id: str = data_paths.get("unique_id", "")
        self.unique_id: str = "_".join([self.part_id, self.ses_id])
        self.derivatives_dir: Optional[str] = data_paths.get("derivatives_dir")

    def __repr__(self) -> str:
        """Return reproducible string representation of ``DataPaths`` instance."""
        return f"DataPaths(data_paths={self.as_dict()})"

    def __str__(self) -> str:
        """Return string representation of a ``DataPaths`` instance."""
        return f"<DataPaths({self.unique_id})>"

    def as_dict(self) -> dict:
        """Return ``data_paths`` dictionary.

        data_paths format::

           {"anat": {"T1w": "{T1w path}", "T2w": "{T2w path}"},
            "creds_path": {None OR path to credentials CSV},
            "func": {
                "{scan ID}": {
                    "scan": "{path to BOLD}",
                    "scan_parameters": {scan parameter dictionary},
                }
            },
            "site_id": "site-ID",
            "subject_id": "sub-01",
            "unique_id": "ses-1",
            "derivatives_dir": "{derivatives_dir path}",}
        """
        return {
            k: v
            for k, v in {
                key: getattr(self, key)
                for key in [
                    "anat",
                    "creds_path",
                    "func",
                    "site_id",
                    "subject_id",
                    "unique_id",
                    "derivatives_dir",
                ]
            }.items()
            if v
        }


@Function.sig_imports(["from typing import Optional"])
def set_iterables(
    scan: str,
    mask_paths: Optional[list[str]] = None,
    ts_paths: Optional[list[str]] = None,
) -> tuple[str, str, str]:
    """Match scan with filepath to get filepath."""
    if mask_paths is None:
        mask_paths = []
    if ts_paths is None:
        ts_paths = []
    mask_path = [path for path in mask_paths if scan in path]
    ts_path = [path for path in ts_paths if scan in path]

    return (scan, mask_path[0], ts_path[0])


def strip_template(data_label: str) -> tuple[str, dict[str, str]]:
    """Strip a template name from a data label to use as a Resource key."""
    json = {}
    # rename to template
    for prefix in ["space-", "from-", "to-"]:
        for bidstag in data_label.split("_"):
            if bidstag.startswith(prefix):
                _template_key, template_val = bidstag.split("-")
                template_name, _template_desc = lookup_identifier(template_val)
                if template_name:
                    json["Template"] = template_val
                    data_label = data_label.replace(template_val, "template")
            elif bidstag.startswith("res-"):
                _res_key, res_val = bidstag.split("-")
                json["Resolution"] = res_val
                data_label = data_label.replace(bidstag, "")
    if data_label.find("__"):
        data_label = data_label.replace("__", "_")
    return data_label, json


class ResourceData(NamedTuple):
    """Attribute and tuple access for ResourceData."""

    node: pe.Node
    """Resource Node."""
    out: str
    """Output key."""


class Resource:
    """A single Resource and its methods."""

    def __init__(self, data: tuple[pe.Node, str], json: dict | list) -> None:
        """Initialize a Resource."""
        self.data = ResourceData(*data)
        """Tuple of source Node and output key."""
        self.json = json
        """Metadata."""
        self._keys = {"data", "json"}
        """Dictionary-style subscriptable keys."""

    def keys(self) -> list[str]:
        """Return list of subscriptable keys."""
        return list(self._keys)

    def __getitem__(self, name: str) -> tuple[pe.Node, str | tuple[str]]:
        """Provide legacy dict-style get access."""
        if name in self.keys():
            return getattr(self, name)
        msg = f"Key '{name}' not set in {self}."
        raise KeyError(msg)

    def __setitem__(self, name: str, value: Any) -> None:
        """Provide legacy dict-style set access."""
        setattr(self, name, value)
        if name not in self.keys():
            self._keys.add(name)

    def __str__(self) -> str:
        """Return string representation of Resource."""
        return f"{self.data[0]}"


class _Pool:
    """All Resources."""

    def __init__(self, name: str = "") -> None:
        """Initialize a ResourcePool or StratPool."""
        self.ants_interp: str
        self.cfg: Configuration
        self.creds_paths: Optional[str]
        self.data_paths: DataPaths
        self.fsl_interp: str
        self.func_reg: bool
        self.fwhm: list[int]
        self.info: dict = {}
        self.logdir: Optional[str]
        self.name = name
        self.num_ants_cores: int
        self.num_cpus = int
        self.part_id: str
        self.pipe_list: list
        self.ses_id: str
        self.smoothing_bool: bool
        self.smooth_opts: list[str]
        self.regressors: dict | list
        self.rpool: dict
        self.run_smoothing: bool
        self.run_zscoring: bool
        self.unique_id: str
        self.zscoring_bool: bool
        self.wf: pe.Workflow
        self._regressor_dct: dict

    def __repr__(self) -> str:
        """Return reproducible ResourcePool string."""
        params = [
            f"{param}={getattr(self, param)}"
            for param in ["rpool", "name", "cfg", "pipe_list"]
            if getattr(self, param, None)
        ]
        return f'{self.__class__.__name__}({", ".join(params)})'

    def __str__(self) -> str:
        """Return string representation of ResourcePool."""
        if self.name:
            return f"ResourcePool({self.name}): {list(self.rpool)}"
        return f"ResourcePool: {list(self.rpool)}"

    def initialize_nipype_wf(self, name: str = "") -> None:
        """Initialize a new nipype workflow."""
        if name:
            name = f"_{name}"
        workflow_name = f"cpac{name}_{self.unique_id}"
        self.wf = pe.Workflow(name=workflow_name)
        self.wf.base_dir = self.cfg.pipeline_setup["working_directory"]["path"]  # type: ignore[attr-defined]
        self.wf.config["execution"] = {
            "hash_method": "timestamp",
            "crashdump_dir": os.path.abspath(
                self.cfg.pipeline_setup["log_directory"]["path"]  # type: ignore[attr-defined]
            ),
        }

    def json_outdir_ingress(
        self, filepath: Path | str, data_label: str, json: dict
    ) -> tuple[dict, tuple[str, str], str, str]:
        """Ingress sidecars from a BIDS derivatives directory."""
        desc_val = None
        for tag in data_label.split("_"):
            if "desc-" in tag:
                desc_val = tag
                break
        jsonpath = str(filepath)
        for ext in EXTS:
            jsonpath = jsonpath.replace(ext, "")
        jsonpath = f"{jsonpath}.json"

        if not os.path.exists(jsonpath):
            WFLOGGER.info(
                "\n\n[!] No JSON found for file %s.\nCreating %s..\n\n",
                filepath,
                jsonpath,
            )
            json_info = {
                "Description": "This data was generated elsewhere and "
                "supplied by the user into this C-PAC run's "
                "output directory. This JSON file was "
                "automatically generated by C-PAC because a "
                "JSON file was not supplied with the data."
            }
            json_info = {**json_info, **json}
            write_output_json(json_info, jsonpath)
        else:
            json_info = read_json(jsonpath)
            json_info = {**json_info, **json}
        if "CpacProvenance" in json_info:
            if desc_val:
                # it's a C-PAC output, let's check for pipe_idx/strat integer
                # suffixes in the desc- entries.
                only_desc = str(desc_val)

                if only_desc[-1].isdigit():
                    for _strat_idx in range(0, 3):
                        # let's stop at 3, please don't run >999 strategies okay?
                        if only_desc[-1].isdigit():
                            only_desc = only_desc[:-1]

                    if only_desc[-1] == "-":
                        only_desc = only_desc.rstrip("-")
                    else:
                        msg = (
                            "\n[!] Something went wrong with either "
                            "reading in the output directory or when "
                            "it was written out previously.\n\nGive "
                            "this to your friendly local C-PAC "
                            f"developer:\n\n{data_label!s}\n"
                        )
                        raise IOError(msg)

                # remove the integer at the end of the desc-* variant, we will
                # get the unique pipe_idx from the CpacProvenance below
                data_label = data_label.replace(desc_val, only_desc)

            # preserve cpac provenance/pipe_idx
            pipe_idx = self.generate_prov_string(json_info["CpacProvenance"])
            node_name = ""
        else:
            json_info["CpacProvenance"] = [f"{data_label}:Non-C-PAC Origin: {filepath}"]  # type: ignore [assignment]
            if "Description" not in json_info:
                json_info["Description"] = (
                    "This data was generated elsewhere and "
                    "supplied by the user into this C-PAC run's "
                    "output directory. This JSON file was "
                    "automatically generated by C-PAC because a "
                    "JSON file was not supplied with the data."
                )
            pipe_idx = self.generate_prov_string(json_info["CpacProvenance"])
            node_name = f"{data_label}_ingress"

        return json_info, pipe_idx, node_name, data_label

    def _config_lookup(self, keylist, fallback_type: type = NoneType) -> Any:
        """Lookup a config key, return None if not found."""
        try:
            return self.cfg[keylist]
        except (AttributeError, KeyError):
            return fallback_type()

    @staticmethod
    def generate_prov_string(prov: LIST_OF_LIST_OF_STR | tuple) -> tuple[str, str]:
        """Generate a string from a SINGLE RESOURCE'S dictionary of MULTIPLE PRECEDING RESOURCES (or single, if just one).

        NOTE: this DOES NOT merge multiple resources!!! (i.e. for merging-strat pipe_idx generation).
        """
        if not isinstance(prov, list):
            msg = (
                "\n[!] Developer info: the CpacProvenance "
                f"entry for {prov} has to be a list.\n"
            )
            raise TypeError(msg)
        last_entry = get_last_prov_entry(prov)
        resource = last_entry.split(":")[0]
        return (resource, str(prov))

    def get_name(self):
        return self.name

    def check_rpool(self, resource):
        if not isinstance(resource, list):
            resource = [resource]
        for name in resource:
            if name in self.rpool:
                return True
        return False

    def create_func_datasource(
        self, rest_dict: dict, wf_name="func_datasource"
    ) -> pe.Workflow:
        """Create a workflow to gather timeseries data.

        Return the functional timeseries-related file paths for each series/scan from the
        dictionary of functional files described in the data configuration (sublist) YAML
        file.

        Scan input (from inputnode) is an iterable.
        """
        wf = pe.Workflow(name=wf_name)

        inputnode = pe.Node(
            util.IdentityInterface(
                fields=["subject", "scan", "creds_path", "dl_dir"],
                mandatory_inputs=True,
            ),
            name="inputnode",
        )

        outputnode = pe.Node(
            util.IdentityInterface(
                fields=[
                    "subject",
                    "rest",
                    "scan",
                    "scan_params",
                    "phase_diff",
                    "magnitude",
                ]
            ),
            name="outputspec",
        )

        # have this here for now because of the big change in the data
        # configuration format
        # (Not necessary with ingress - format does not comply)
        if not self.check_rpool("derivatives-dir"):
            check_scan = pe.Node(
                Function(
                    input_names=["func_scan_dct", "scan"],
                    output_names=[],
                    function=check_func_scan,
                    as_module=True,
                ),
                name="check_func_scan",
            )

            check_scan.inputs.func_scan_dct = rest_dict
            wf.connect(inputnode, "scan", check_scan, "scan")

        # get the functional scan itself
        selectrest = pe.Node(
            Function(
                input_names=["scan", "rest_dict", "resource"],
                output_names=["file_path"],
                function=get_rest,
                as_module=True,
            ),
            name="selectrest",
        )
        selectrest.inputs.rest_dict = rest_dict
        selectrest.inputs.resource = "scan"
        wf.connect(inputnode, "scan", selectrest, "scan")

        # check to see if it's on an Amazon AWS S3 bucket, and download it, if it
        # is - otherwise, just return the local file path
        check_s3_node = pe.Node(
            Function(
                input_names=["file_path", "creds_path", "dl_dir", "img_type"],
                output_names=["local_path"],
                function=check_for_s3,
                as_module=True,
            ),
            name="check_for_s3",
        )

        wf.connect(selectrest, "file_path", check_s3_node, "file_path")
        wf.connect(inputnode, "creds_path", check_s3_node, "creds_path")
        wf.connect(inputnode, "dl_dir", check_s3_node, "dl_dir")
        check_s3_node.inputs.img_type = "func"

        wf.connect(inputnode, "subject", outputnode, "subject")
        wf.connect(check_s3_node, "local_path", outputnode, "rest")
        wf.connect(inputnode, "scan", outputnode, "scan")

        # scan parameters CSV
        select_scan_params = pe.Node(
            Function(
                input_names=["scan", "rest_dict", "resource"],
                output_names=["file_path"],
                function=get_rest,
                as_module=True,
            ),
            name="select_scan_params",
        )
        select_scan_params.inputs.rest_dict = rest_dict
        select_scan_params.inputs.resource = "scan_parameters"
        wf.connect(inputnode, "scan", select_scan_params, "scan")

        # if the scan parameters file is on AWS S3, download it
        s3_scan_params = pe.Node(
            Function(
                input_names=["file_path", "creds_path", "dl_dir", "img_type"],
                output_names=["local_path"],
                function=check_for_s3,
                as_module=True,
            ),
            name="s3_scan_params",
        )

        wf.connect(select_scan_params, "file_path", s3_scan_params, "file_path")
        wf.connect(inputnode, "creds_path", s3_scan_params, "creds_path")
        wf.connect(inputnode, "dl_dir", s3_scan_params, "dl_dir")
        wf.connect(s3_scan_params, "local_path", outputnode, "scan_params")

        return wf

    def get_pipe_number(self, pipe_idx):
        return self.pipe_list.index(pipe_idx)

    def get_pool_info(self):
        return self.info

    def set_pool_info(self, info_dct):
        self.info.update(info_dct)

    def get_entire_rpool(self):
        return self.rpool

    def keys(self) -> KeysView:
        """Return rpool's keys."""
        return self.rpool.keys()

    @staticmethod
    def get_raw_label(resource: str) -> str:
        """Remove ``desc-*`` label."""
        for tag in resource.split("_"):
            if "desc-" in tag:
                resource = resource.replace(f"{tag}_", "")
                break
        return resource

    def get_strat_info(self, prov, label=None, logdir=None):
        strat_info = {}
        for entry in prov:
            if isinstance(entry, list):
                strat_info[entry[-1].split(":")[0]] = entry
            elif isinstance(entry, str):
                strat_info[entry.split(":")[0]] = entry.split(":")[1]
        if label:
            if not logdir:
                logdir = self.logdir
            WFLOGGER.info(
                "\n\nPrinting out strategy info for %s in %s\n", label, logdir
            )
            write_output_json(
                strat_info, f"{label}_strat_info", indent=4, basedir=logdir
            )

    def set_json_info(self, resource, pipe_idx, key, val):
        # TODO: actually should probably be able to inititialize resource/pipe_idx
        if pipe_idx not in self.rpool[resource]:
            msg = (
                "\n[!] DEV: The pipeline/strat ID does not exist "
                f"in the resource pool.\nResource: {resource}"
                f"Pipe idx: {pipe_idx}\nKey: {key}\nVal: {val}\n"
            )
            raise Exception(msg)
        if "json" not in self.rpool[resource][pipe_idx]:
            self.rpool[resource][pipe_idx]["json"] = {}
        self.rpool[resource][pipe_idx]["json"][key] = val

    def get_json_info(self, resource, pipe_idx, key):
        # TODO: key checks
        if not pipe_idx:
            for pipe_idx, val in self.rpool[resource].items():
                return val["json"][key]
        return self.rpool[resource][pipe_idx][key]

    @staticmethod
    def get_resource_from_prov(prov: LIST_OF_LIST_OF_STR) -> Optional[str]:
        """Return the last item in the provenance list.

        Each resource (i.e. "desc-cleaned_bold" AKA nuisance-regressed BOLD
        data) has its own provenance list. the name of the resource, and
        the node that produced it, is always the last item in the provenance
        list, with the two separated by a colon :
        """
        if not len(prov):
            return None
        if isinstance(prov[-1], list):
            last_item_in_list = prov[-1][-1]
            assert isinstance(last_item_in_list, str)
            return last_item_in_list.split(":")[0]
        if isinstance(prov[-1], str):
            return prov[-1].split(":")[0]
        return None

    def regressor_dct(self, cfg) -> dict:
        """Return the regressor dictionary for the current strategy if one exists.

        Raises KeyError otherwise.
        """
        # pylint: disable=attribute-defined-outside-init
        if hasattr(self, "_regressor_dct"):  # memoized
            # pylint: disable=access-member-before-definition
            return self._regressor_dct
        key_error = KeyError(
            "[!] No regressors in resource pool. \n\n"
            "Try turning on create_regressors or "
            "ingress_regressors."
        )
        _nr = cfg["nuisance_corrections", "2-nuisance_regression"]
        if not hasattr(self, "timeseries"):
            if _nr["Regressors"]:
                self.regressors = {reg["Name"]: reg for reg in _nr["Regressors"]}
            else:
                self.regressors = []
        if self.check_rpool("parsed_regressors"):  # ingressed regressor
            # name regressor workflow without regressor_prov
            strat_name = _nr["ingress_regressors"]["Regressors"]["Name"]
            if strat_name in self.regressors:
                self._regressor_dct = self.regressors[strat_name]
                return self._regressor_dct
            self._regressor_dct = _nr["ingress_regressors"]["Regressors"]
            return self._regressor_dct
        prov = self.get_cpac_provenance("desc-confounds_timeseries")
        strat_name_components = prov[-1].split("_")
        for _ in list(range(prov[-1].count("_"))):
            reg_name = "_".join(strat_name_components[-_:])
            if isinstance(self.regressors, dict) and reg_name in self.regressors:
                self._regressor_dct = self.regressors[reg_name]
                return self._regressor_dct
        raise key_error

    def set_data(
        self,
        resource: str,
        node: pe.Node | pe.Workflow,
        output: str,
        json_info: dict,
        pipe_idx: PIPE_IDX,
        node_name: str,
        fork: bool = False,
        inject: bool = False,
    ) -> None:
        """Plug a Resource into a _Pool."""
        json_info = json_info.copy()
        cpac_prov: LIST_OF_LIST_OF_STR = []
        if "CpacProvenance" in json_info:
            cpac_prov = json_info["CpacProvenance"]
        current_prov_list = list(cpac_prov)
        new_prov_list = list(cpac_prov)  # <---- making a copy, it was already a list
        if not inject:
            new_prov_list.append(f"{resource}:{node_name}")
        try:
            _resource, new_pipe_idx = self.generate_prov_string(new_prov_list)
        except IndexError:
            msg = (
                f"\n\nThe set_data() call for {resource} has no "
                "provenance information and should not be an "
                "injection."
            )
            raise IndexError(msg)
        if not json_info:
            json_info = {
                "RawSources": [
                    resource  # <---- this will be repopulated to the full file path at the end of the pipeline building, in gather_pipes()
                ]
            }
        json_info["CpacProvenance"] = new_prov_list

        if resource not in self.keys():
            self.rpool[resource] = {}
        elif not fork:  # <--- in the event of multiple strategies/options, this will run for every option; just keep in mind
            search = False
            if self.get_resource_from_prov(current_prov_list) == resource:
                # CHANGING PIPE_IDX, BE CAREFUL DOWNSTREAM IN THIS FUNCTION
                pipe_idx = self.generate_prov_string(current_prov_list)[1]
                if pipe_idx not in self.rpool[resource].keys():
                    search = True
            else:
                search = True
            if search:
                for idx in current_prov_list:
                    if self.get_resource_from_prov(idx) == resource:
                        if isinstance(idx, list):
                            # CHANGING PIPE_IDX, BE CAREFUL DOWNSTREAM IN THIS FUNCTION
                            pipe_idx = self.generate_prov_string(idx)[1]
                        elif isinstance(idx, str):
                            pipe_idx = idx
                        break
            if pipe_idx in self.rpool[resource].keys():
                # in case the resource name is now new, and not the original
                # remove old keys so we don't end up with a new strat for every new node unit (unless we fork)
                del self.rpool[resource][pipe_idx]
        if new_pipe_idx not in self.rpool[resource]:
            self.rpool[resource][new_pipe_idx] = Resource(
                data=ResourceData(node, output), json=json_info
            )
        if new_pipe_idx not in self.pipe_list:
            self.pipe_list.append(new_pipe_idx)

    # @overload
    # def get(
    #     self: "ResourcePool",
    #     resource: list[str] | str,
    #     pipe_idx: Optional[PIPE_IDX],
    #     report_fetched: bool,
    #     optional: bool,
    # ) -> Optional[dict[dict]] | tuple[Optional[dict[dict]], Optional[str]]: ...
    # @overload
    # def get(
    #     self: "StratPool",
    #     resource: list[str] | str,
    #     pipe_idx: Optional[PIPE_IDX],
    #     report_fetched: bool,
    #     optional: bool,
    # ) -> Optional[Resource] | tuple[Optional[Resource], Optional[str]]: ...
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[PIPE_IDX],
        report_fetched: bool,
        optional: bool,
    ) -> (
        Optional[Resource | STRAT_DICT]
        | tuple[Optional[Resource | STRAT_DICT], Optional[str]]
    ):
        """Return a dictionary of strats or a single Resource."""
        if not isinstance(resource, list):
            resource = [resource]
        # if a list of potential inputs are given, pick the first one found
        for label in resource:
            if label in self.keys():
                _found = self.rpool[label]
                if pipe_idx:
                    _found = _found[pipe_idx]
                if report_fetched:
                    return _found, label
                return _found
        if optional:
            if report_fetched:
                return (None, None)
            return None
        msg = (
            "\n\n[!] C-PAC says: None of the listed resources are in "
            f"the resource pool:\n\n  {resource}\n\nOptions:\n- You "
            "can enable a node block earlier in the pipeline which "
            "produces these resources. Check the 'outputs:' field in "
            "a node block's documentation.\n- You can directly "
            "provide this required data by pulling it from another "
            "BIDS directory using 'source_outputs_dir:' in the "
            "pipeline configuration, or by placing it directly in "
            "your C-PAC output directory.\n- If you have done these, "
            "and you still get this message, please let us know "
            "through any of our support channels at: "
            "https://fcp-indi.github.io/\n"
        )
        raise LookupError(msg)

    def copy_resource(self, resource, new_name):
        try:
            self.rpool[new_name] = self.rpool[resource]
        except KeyError:
            msg = f"[!] {resource} not in the resource pool."
            raise Exception(msg)

    def update_resource(self, resource, new_name):
        # move over any new pipe_idx's
        self.rpool[new_name].update(self.rpool[resource])

    def get_pipe_idxs(self, resource):
        return self.rpool[resource].keys()

    def get_json(self, resource, strat=None):
        # NOTE: resource_strat_dct has to be entered properly by the developer
        # it has to either be rpool[resource][strat] or strat_pool[resource]
        if strat:
            resource_strat_dct = self.rpool[resource][strat]
        else:
            # for strat_pools mainly, where there is no 'strat' key level
            resource_strat_dct = self.rpool[resource]

        # TODO: the below hits the exception if you use get_cpac_provenance on
        # TODO: the main rpool (i.e. if strat=None)
        if "json" in resource_strat_dct:
            strat_json = resource_strat_dct["json"]
        else:
            msg = (
                "\n[!] Developer info: the JSON "
                f"information for {resource} and {strat} "
                f"is incomplete.\n"
            )
            raise Exception(msg)
        return strat_json

    def get_cpac_provenance(self, resource, strat=None):
        # NOTE: resource_strat_dct has to be entered properly by the developer
        # it has to either be rpool[resource][strat] or strat_pool[resource]
        if isinstance(resource, list):
            for _resource in resource:
                try:
                    return self.get_cpac_provenance(_resource, strat)
                except KeyError:
                    continue
        json_data = self.get_json(resource, strat)
        return json_data["CpacProvenance"]

    @staticmethod
    def generate_prov_list(prov_str):
        if not isinstance(prov_str, str):
            msg = (
                "\n[!] Developer info: the CpacProvenance "
                f"entry for {prov_str!s} has to be a string.\n"
            )
            raise TypeError(msg)
        return ast.literal_eval(prov_str)

    @staticmethod
    def get_resource_strats_from_prov(prov):
        # if you provide the provenance of a resource pool output, this will
        # return a dictionary of all the preceding resource pool entries that
        # led to that one specific output:
        #   {rpool entry}: {that entry's provenance}
        #   {rpool entry}: {that entry's provenance}
        resource_strat_dct = {}
        if isinstance(prov, str):
            resource = prov.split(":")[0]
            resource_strat_dct[resource] = prov
        else:
            for spot, entry in enumerate(prov):
                if isinstance(entry, list):
                    resource = entry[-1].split(":")[0]
                    resource_strat_dct[resource] = entry
                elif isinstance(entry, str):
                    resource = entry.split(":")[0]
                    resource_strat_dct[resource] = entry
        return resource_strat_dct

    def flatten_prov(self, prov):
        if isinstance(prov, str):
            return [prov]
        if isinstance(prov, list):
            flat_prov = []
            for entry in prov:
                if isinstance(entry, list):
                    flat_prov += self.flatten_prov(entry)
                else:
                    flat_prov.append(entry)
            return flat_prov
        return None

    def derivative_xfm(self, wf, label, connection, json_info, pipe_idx, pipe_x):
        if label in self.xfm:
            json_info = dict(json_info)

            # get the bold-to-template transform from the current strat_pool info
            xfm_idx = None
            xfm_label = "from-bold_to-template_mode-image_xfm"
            for entry in json_info["CpacProvenance"]:
                if isinstance(entry, list):
                    if entry[-1].split(":")[0] == xfm_label:
                        xfm_prov = entry
                        xfm_idx = self.generate_prov_string(xfm_prov)[1]
                        break

            # but if the resource doesn't have the bold-to-template transform
            # in its provenance/strategy, find the appropriate one for this
            # current pipe_idx/strat
            if not xfm_idx:
                xfm_info = []
                for pipe_idx, entry in self.get(xfm_label).items():
                    xfm_info.append((pipe_idx, entry["json"]["CpacProvenance"]))
            else:
                xfm_info = [(xfm_idx, xfm_prov)]

            for num, xfm_entry in enumerate(xfm_info):
                xfm_idx, xfm_prov = xfm_entry
                reg_tool = check_prov_for_regtool(xfm_prov)

                xfm = transform_derivative(
                    f"{label}_xfm_{pipe_x}_{num}",
                    label,
                    reg_tool,
                    self.num_cpus,
                    self.num_ants_cores,
                    ants_interp=self.ants_interp,
                    fsl_interp=self.fsl_interp,
                    opt=None,
                )
                wf.connect(connection[0], connection[1], xfm, "inputspec.in_file")

                node, out = self.get_data("T1w-brain-template-deriv", quick_single=True)
                wf.connect(node, out, xfm, "inputspec.reference")

                node, out = self.get_data(
                    "from-bold_to-template_mode-image_xfm", pipe_idx=xfm_idx
                )
                wf.connect(node, out, xfm, "inputspec.transform")

                label = f"space-template_{label}"
                json_info["Template"] = self.get_json_info(
                    "T1w-brain-template-deriv", None, "Description"
                )
                new_prov = json_info["CpacProvenance"] + xfm_prov
                json_info["CpacProvenance"] = new_prov
                new_pipe_idx = self.generate_prov_string(new_prov)
                self.set_data(
                    label,
                    xfm,
                    "outputspec.out_file",
                    json_info,
                    new_pipe_idx,
                    f"{label}_xfm_{num}",
                    fork=True,
                )

        return wf

    @property
    def filtered_movement(self) -> bool:
        """
        Check if the movement parameters have been filtered in this strat_pool.

        Returns
        -------
        bool
        """
        try:
            return "motion_estimate_filter" in str(
                self.get_cpac_provenance("desc-movementParameters_motion")
            )
        except KeyError:
            # not a strat_pool or no movement parameters in strat_pool
            return False

    def filter_name(self, cfg: Configuration) -> str:
        """
        Return the name of the filter for this strategy.

        In a strat_pool with filtered movement parameters.
        """
        motion_filters = cfg[
            "functional_preproc",
            "motion_estimates_and_correction",
            "motion_estimate_filter",
            "filters",
        ]
        if len(motion_filters) == 1 and cfg.switch_is_on(
            [
                "functional_preproc",
                "motion_estimates_and_correction",
                "motion_estimate_filter",
                "run",
            ],
            exclusive=True,
        ):
            return motion_filters[0]["Name"]
        try:
            key = "motion"
            sidecar = self.get_json("desc-movementParameters_motion")
        except KeyError:
            sidecar = None
        if sidecar is not None and "CpacVariant" in sidecar:
            if sidecar["CpacVariant"][key]:
                return sidecar["CpacVariant"][key][0][::-1].split("_", 1)[0][::-1]
        return "none"

    def post_process(self, wf, label, connection, json_info, pipe_idx, pipe_x, outs):
        input_type = "func_derivative"

        post_labels = [(label, connection[0], connection[1])]

        if re.match(r"(.*_)?[ed]c[bw]$", label) or re.match(r"(.*_)?lfcd[bw]$", label):
            # suffix: [eigenvector or degree] centrality [binarized or weighted]
            # or lfcd [binarized or weighted]
            mask = "template-specification-file"
        elif "space-template" in label:
            if "space-template_res-derivative_desc-bold_mask" in self.keys():
                mask = "space-template_res-derivative_desc-bold_mask"
            else:
                mask = "space-template_desc-bold_mask"
        else:
            mask = "space-bold_desc-brain_mask"

        mask_idx = None
        for entry in json_info["CpacProvenance"]:
            if isinstance(entry, list):
                if entry[-1].split(":")[0] == mask:
                    mask_prov = entry
                    mask_idx = self.generate_prov_string(mask_prov)[1]
                    break

        if self.smoothing_bool:
            if label in Outputs.to_smooth:
                for smooth_opt in self.smooth_opts:
                    sm = spatial_smoothing(
                        f"{label}_smooth_{smooth_opt}_{pipe_x}",
                        self.fwhm,
                        input_type,
                        smooth_opt,
                    )
                    wf.connect(connection[0], connection[1], sm, "inputspec.in_file")
                    node, out = self.get_data(
                        mask, pipe_idx=mask_idx, quick_single=mask_idx is None
                    )
                    wf.connect(node, out, sm, "inputspec.mask")

                    if "desc-" not in label:
                        if "space-" in label:
                            for tag in label.split("_"):
                                if "space-" in tag:
                                    smlabel = label.replace(tag, f"{tag}_desc-sm")
                                    break
                        else:
                            smlabel = f"desc-sm_{label}"
                    else:
                        for tag in label.split("_"):
                            if "desc-" in tag:
                                newtag = f"{tag}-sm"
                                smlabel = label.replace(tag, newtag)
                                break

                    post_labels.append((smlabel, sm, "outputspec.out_file"))

                    self.set_data(
                        smlabel,
                        sm,
                        "outputspec.out_file",
                        json_info,
                        pipe_idx,
                        f"spatial_smoothing_{smooth_opt}",
                        fork=True,
                    )
                    self.set_data(
                        "fwhm",
                        sm,
                        "outputspec.fwhm",
                        json_info,
                        pipe_idx,
                        f"spatial_smoothing_{smooth_opt}",
                        fork=True,
                    )

        if self.zscoring_bool:
            for label_con_tpl in post_labels:
                label = label_con_tpl[0]
                connection = (label_con_tpl[1], label_con_tpl[2])
                if label in Outputs.to_zstd:
                    zstd = z_score_standardize(f"{label}_zstd_{pipe_x}", input_type)

                    wf.connect(connection[0], connection[1], zstd, "inputspec.in_file")

                    node, out = self.get_data(mask, pipe_idx=mask_idx)
                    wf.connect(node, out, zstd, "inputspec.mask")

                    if "desc-" not in label:
                        if "space-template" in label:
                            new_label = label.replace(
                                "space-template", "space-template_desc-zstd"
                            )
                        else:
                            new_label = f"desc-zstd_{label}"
                    else:
                        for tag in label.split("_"):
                            if "desc-" in tag:
                                newtag = f"{tag}-zstd"
                                new_label = label.replace(tag, newtag)
                                break

                    post_labels.append((new_label, zstd, "outputspec.out_file"))

                    self.set_data(
                        new_label,
                        zstd,
                        "outputspec.out_file",
                        json_info,
                        pipe_idx,
                        "zscore_standardize",
                        fork=True,
                    )

                elif label in Outputs.to_fisherz:
                    zstd = fisher_z_score_standardize(
                        f"{label}_zstd_{pipe_x}", label, input_type
                    )

                    wf.connect(
                        connection[0], connection[1], zstd, "inputspec.correlation_file"
                    )

                    # if the output is 'space-template_desc-MeanSCA_correlations', we want 'desc-MeanSCA_timeseries'
                    oned = label.replace("correlations", "timeseries")

                    node, out = outs[oned]
                    wf.connect(node, out, zstd, "inputspec.timeseries_oned")

                    post_labels.append((new_label, zstd, "outputspec.out_file"))

                    self.set_data(
                        new_label,
                        zstd,
                        "outputspec.out_file",
                        json_info,
                        pipe_idx,
                        "fisher_zscore_standardize",
                        fork=True,
                    )

        return (wf, post_labels)

    def node_data(self, resource: str | tuple[str], **kwargs) -> ResourceData:
        """Create ResourceData objects."""
        return ResourceData(*self.get_data(resource, **kwargs))  # type: ignore[attr-defined]


class ResourcePool(_Pool):
    """A pool of Resources."""

    def __init__(
        self,
        name: str = "",
        cfg: Optional[Configuration] = None,
        pipe_list: Optional[list] = None,
        *,
        data_paths: Optional[DataPaths | dict] = None,
        part_id: Optional[str] = None,
        pipeline_name: str = "",
        wf: Optional[pe.Workflow] = None,
    ) -> None:
        """Initialize a ResourcePool."""
        super().__init__(name=name)
        if isinstance(data_paths, dict):
            data_paths = DataPaths(data_paths=data_paths)
        elif not data_paths:
            data_paths = DataPaths(part_id=part_id)
        self.data_paths = data_paths
        # pass-through for convenient access
        self.creds_path = self.data_paths.creds_path
        self.part_id = self.data_paths.part_id
        self.ses_id = self.data_paths.ses_id
        self.unique_id = self.data_paths.unique_id
        self.rpool: POOL_DICT = {}

        if not pipe_list:
            self.pipe_list = []
        else:
            self.pipe_list = pipe_list

        if cfg:
            self.cfg = cfg
        else:
            self.cfg = Preconfiguration("blank")

        self.logdir = self._config_lookup(["pipeline_setup", "log_directory", "path"])
        self.num_cpus = self._config_lookup(
            ["pipeline_setup", "system_config", "max_cores_per_participant"]
        )
        self.num_ants_cores = self._config_lookup(
            ["pipeline_setup", "system_config", "num_ants_threads"]
        )

        self.ants_interp = self._config_lookup(
            [
                "registration_workflows",
                "functional_registration",
                "func_registration_to_template",
                "ANTs_pipelines",
                "interpolation",
            ]
        )
        self.fsl_interp = self._config_lookup(
            [
                "registration_workflows",
                "functional_registration",
                "func_registration_to_template",
                "FNIRT_pipelines",
                "interpolation",
            ]
        )
        self.func_reg = self._config_lookup(
            [
                "registration_workflows",
                "functional_registration",
                "func_registration_to_template",
                "run",
            ]
        )

        self.run_smoothing = "smoothed" in self._config_lookup(
            ["post_processing", "spatial_smoothing", "output"], list
        )
        self.smoothing_bool = self._config_lookup(
            ["post_processing", "spatial_smoothing", "run"]
        )
        self.run_zscoring = "z-scored" in self._config_lookup(
            ["post_processing", "z-scoring", "output"], list
        )
        self.zscoring_bool = self._config_lookup(
            ["post_processing", "z-scoring", "run"]
        )
        self.fwhm = self._config_lookup(
            ["post_processing", "spatial_smoothing", "fwhm"]
        )
        self.smooth_opts = self._config_lookup(
            ["post_processing", "spatial_smoothing", "smoothing_method"]
        )

        if wf:
            self.wf = wf
        else:
            self.initialize_nipype_wf(pipeline_name)

        self.xfm = [
            "alff",
            "desc-sm_alff",
            "desc-zstd_alff",
            "desc-sm-zstd_alff",
            "falff",
            "desc-sm_falff",
            "desc-zstd_falff",
            "desc-sm-zstd_falff",
            "reho",
            "desc-sm_reho",
            "desc-zstd_reho",
            "desc-sm-zstd_reho",
        ]
        ingress_derivatives = False
        try:
            if self.data_paths.derivatives_dir and self._config_lookup(
                ["pipeline_setup", "outdir_ingress", "run"], bool
            ):
                ingress_derivatives = True
        except (AttributeError, KeyError, TypeError):
            pass
        if ingress_derivatives:
            self.ingress_output_dir()
        else:
            self.ingress_raw_anat_data()
            if data_paths.func:
                self.ingress_raw_func_data()
        self.ingress_pipeconfig_paths()

    def back_propogate_template_name(
        self, resource_idx: str, json_info: dict, id_string: "pe.Node"
    ) -> None:
        """Find and apply the template name from a resource's provenance."""
        if "template" in resource_idx and self.check_rpool("derivatives-dir"):
            if self.check_rpool("template"):
                node, out = self.get_data("template")
                self.wf.connect(node, out, id_string, "template_desc")
        elif "Template" in json_info:
            id_string.inputs.template_desc = json_info["Template"]
        elif (
            "template" in resource_idx and len(json_info.get("CpacProvenance", [])) > 1
        ):
            for resource in source_set(json_info["CpacProvenance"]):
                source, value = resource.split(":", 1)
                if value.startswith("template_") and source != "FSL-AFNI-bold-ref":
                    # 'FSL-AFNI-bold-ref' is currently allowed to be in
                    # a different space, so don't use it as the space for
                    # descendents
                    try:
                        ancestors = self.rpool.get(source)
                        assert ancestors is not None
                        ancestor: dict = next(iter(ancestors.items()))[1]
                        anscestor_json: dict = ancestor.get("json", {})
                        if "Description" in anscestor_json:
                            id_string.inputs.template_desc = anscestor_json[
                                "Description"
                            ]
                            return
                    except (IndexError, KeyError):
                        pass
        return

    def gather_pipes(self, wf, cfg, all=False, add_incl=None, add_excl=None):
        excl = []
        substring_excl = []
        outputs_logger = getLogger(f"{self.part_id}_expectedOutputs")
        expected_outputs = ExpectedOutputs()

        if add_excl:
            excl += add_excl

        if "nonsmoothed" not in cfg.post_processing["spatial_smoothing"]["output"]:
            excl += Outputs.native_nonsmooth
            excl += Outputs.template_nonsmooth

        if "raw" not in cfg.post_processing["z-scoring"]["output"]:
            excl += Outputs.native_raw
            excl += Outputs.template_raw

        if not cfg.pipeline_setup["output_directory"]["write_debugging_outputs"]:
            # substring_excl.append(['bold'])
            excl += Outputs.debugging

        for resource in self.keys():
            if resource not in Outputs.any:
                continue

            if resource in excl:
                continue

            drop = False
            for substring_list in substring_excl:
                bool_list = []
                for substring in substring_list:
                    if substring in resource:
                        bool_list.append(True)
                    else:
                        bool_list.append(False)
                for item in bool_list:
                    if not item:
                        break
                else:
                    drop = True
                if drop:
                    break
            if drop:
                continue

            subdir = "other"
            if resource in Outputs.anat:
                subdir = "anat"
                # TODO: get acq- etc.
            elif resource in Outputs.func:
                subdir = "func"
                # TODO: other stuff like acq- etc.

            for pipe_idx in self.rpool[resource]:
                unique_id = self.unique_id
                part_id = self.part_id
                ses_id = self.ses_id

                if "ses-" not in ses_id:
                    ses_id = f"ses-{ses_id}"

                out_dir = cfg.pipeline_setup["output_directory"]["path"]
                pipe_name = cfg.pipeline_setup["pipeline_name"]
                container = os.path.join(f"pipeline_{pipe_name}", part_id, ses_id)
                filename = f"{unique_id}_{res_in_filename(self.cfg, resource)}"

                out_path = os.path.join(out_dir, container, subdir, filename)

                out_dct = {
                    "unique_id": unique_id,
                    "out_dir": out_dir,
                    "container": container,
                    "subdir": subdir,
                    "filename": filename,
                    "out_path": out_path,
                }
                self.rpool[resource][pipe_idx]["out"] = out_dct

                # TODO: have to link the pipe_idx's here. and call up 'desc-preproc_T1w' from a Sources in a json and replace. here.
                # TODO: can do the pipeline_description.json variants here too!

        for resource in self.keys():
            if resource not in Outputs.any:
                continue

            if resource in excl:
                continue

            drop = False
            for substring_list in substring_excl:
                bool_list = []
                for substring in substring_list:
                    if substring in resource:
                        bool_list.append(True)
                    else:
                        bool_list.append(False)
                for item in bool_list:
                    if not item:
                        break
                else:
                    drop = True
                if drop:
                    break
            if drop:
                continue

            num_variant = 0
            if len(self.rpool[resource]) == 1:
                num_variant = ""
            all_jsons = [
                self.rpool[resource][pipe_idx]["json"]
                for pipe_idx in self.rpool[resource]
            ]
            unlabelled = {
                key
                for json_info in all_jsons
                for key in json_info.get("CpacVariant", {}).keys()
                if key not in (*MOVEMENT_FILTER_KEYS, "regressors")
            }
            if "bold" in unlabelled:
                all_bolds = list(
                    chain.from_iterable(
                        json_info["CpacVariant"]["bold"]
                        for json_info in all_jsons
                        if "CpacVariant" in json_info
                        and "bold" in json_info["CpacVariant"]
                    )
                )
                # not any(not) because all is overloaded as a parameter here
                if not any(
                    not re.match(
                        r"apply_(phasediff|blip)_to_timeseries_separately_.*", _bold
                    )
                    for _bold in all_bolds
                ):
                    # this fork point should only result in 0 or 1 forks
                    unlabelled.remove("bold")
                del all_bolds
            all_forks = {
                key: set(
                    chain.from_iterable(
                        json_info["CpacVariant"][key]
                        for json_info in all_jsons
                        if "CpacVariant" in json_info
                        and key in json_info["CpacVariant"]
                    )
                )
                for key in unlabelled
            }
            # del all_jsons
            for key, forks in all_forks.items():
                if len(forks) < 2:  # noqa: PLR2004
                    # no int suffix needed if only one fork
                    unlabelled.remove(key)
            # del all_forks
            for pipe_idx in self.rpool[resource]:
                pipe_x = self.get_pipe_number(pipe_idx)
                json_info = self.rpool[resource][pipe_idx]["json"]
                out_dct = self.rpool[resource][pipe_idx]["out"]

                try:
                    if unlabelled:
                        num_variant += 1
                except TypeError:
                    pass

                try:
                    del json_info["subjson"]
                except KeyError:
                    pass

                if out_dct["subdir"] == "other" and not all:
                    continue

                unique_id = out_dct["unique_id"]
                resource_idx = resource

                if isinstance(num_variant, int):
                    resource_idx, out_dct = name_fork(
                        resource_idx, cfg, json_info, out_dct
                    )
                    if unlabelled:
                        if "desc-" in out_dct["filename"]:
                            for key in out_dct["filename"].split("_")[::-1]:
                                # final `desc` entity
                                if key.startswith("desc-"):
                                    out_dct["filename"] = out_dct["filename"].replace(
                                        key, f"{key}-{num_variant}"
                                    )
                                    resource_idx = resource_idx.replace(
                                        key, f"{key}-{num_variant}"
                                    )
                                    break
                        else:
                            suff = resource.split("_")[-1]
                            newdesc_suff = f"desc-{num_variant}_{suff}"
                            resource_idx = resource_idx.replace(suff, newdesc_suff)
                id_string = pe.Node(
                    Function(
                        input_names=[
                            "cfg",
                            "unique_id",
                            "resource",
                            "scan_id",
                            "template_desc",
                            "atlas_id",
                            "fwhm",
                            "subdir",
                            "extension",
                        ],
                        output_names=["out_filename"],
                        function=create_id_string,
                    ),
                    name=f"id_string_{resource_idx}_{pipe_x}",
                )
                id_string.inputs.cfg = self.cfg
                id_string.inputs.unique_id = unique_id
                id_string.inputs.resource = resource_idx
                id_string.inputs.subdir = out_dct["subdir"]

                # grab the iterable scan ID
                if out_dct["subdir"] == "func":
                    node, out = self.rpool["scan"]["['scan:func_ingress']"]["data"]
                    wf.connect(node, out, id_string, "scan_id")

                self.back_propogate_template_name(resource_idx, json_info, id_string)
                # grab the FWHM if smoothed
                for tag in resource.split("_"):
                    if "desc-" in tag and "-sm" in tag:
                        fwhm_idx = pipe_idx.replace(f"{resource}:", "fwhm:")
                        try:
                            node, out = self.rpool["fwhm"][fwhm_idx]["data"]
                            wf.connect(node, out, id_string, "fwhm")
                        except KeyError:
                            # smoothing was not done for this resource in the
                            # engine.py smoothing
                            pass
                        break
                atlas_suffixes = ["timeseries", "correlations", "statmap"]
                # grab the iterable atlas ID
                atlas_id = None
                if not resource.endswith("desc-confounds_timeseries"):
                    if resource.split("_")[-1] in atlas_suffixes:
                        atlas_idx = pipe_idx.replace(resource, "atlas_name")
                        # need the single quote and the colon inside the double
                        # quotes - it's the encoded pipe_idx
                        # atlas_idx = new_idx.replace(f"'{temp_rsc}:",
                        #                             "'atlas_name:")
                        if atlas_idx in self.rpool["atlas_name"]:
                            node, out = self.rpool["atlas_name"][atlas_idx]["data"]
                            wf.connect(node, out, id_string, "atlas_id")
                        elif "atlas-" in resource:
                            for tag in resource.split("_"):
                                if "atlas-" in tag:
                                    atlas_id = tag.replace("atlas-", "")
                            id_string.inputs.atlas_id = atlas_id
                        else:
                            warnings.warn(
                                str(
                                    LookupError(
                                        "\n[!] No atlas ID found for "
                                        f"{out_dct['filename']}.\n"
                                    )
                                )
                            )
                nii_name = pe.Node(Rename(), name=f"nii_{resource_idx}_{pipe_x}")
                nii_name.inputs.keep_ext = True

                if resource in Outputs.ciftis:
                    nii_name.inputs.keep_ext = False
                    id_string.inputs.extension = Outputs.ciftis[resource]
                else:
                    nii_name.inputs.keep_ext = True

                if resource in Outputs.giftis:
                    nii_name.inputs.keep_ext = False
                    id_string.inputs.extension = f"{Outputs.giftis[resource]}.gii"

                else:
                    nii_name.inputs.keep_ext = True

                wf.connect(id_string, "out_filename", nii_name, "format_string")

                node, out = self.rpool[resource][pipe_idx]["data"]
                try:
                    wf.connect(node, out, nii_name, "in_file")
                except OSError as os_error:
                    WFLOGGER.warning(os_error)
                    continue

                write_json_imports = ["import os", "import json"]
                write_json = pe.Node(
                    Function(
                        input_names=["json_data", "filename"],
                        output_names=["json_file"],
                        function=write_output_json,
                        imports=write_json_imports,
                    ),
                    name=f"json_{resource_idx}_{pipe_x}",
                )
                write_json.inputs.json_data = json_info

                wf.connect(id_string, "out_filename", write_json, "filename")
                ds = pe.Node(DataSink(), name=f"sinker_{resource_idx}_{pipe_x}")
                ds.inputs.parameterization = False
                ds.inputs.base_directory = out_dct["out_dir"]
                ds.inputs.encrypt_bucket_keys = cfg.pipeline_setup["Amazon-AWS"][
                    "s3_encryption"
                ]
                ds.inputs.container = out_dct["container"]

                if cfg.pipeline_setup["Amazon-AWS"]["aws_output_bucket_credentials"]:
                    ds.inputs.creds_path = cfg.pipeline_setup["Amazon-AWS"][
                        "aws_output_bucket_credentials"
                    ]
                expected_outputs += (
                    out_dct["subdir"],
                    create_id_string(
                        self.cfg,
                        unique_id,
                        resource_idx,
                        template_desc=id_string.inputs.template_desc,
                        atlas_id=atlas_id,
                        subdir=out_dct["subdir"],
                    ),
                )
                wf.connect(nii_name, "out_file", ds, f'{out_dct["subdir"]}.@data')
                wf.connect(write_json, "json_file", ds, f'{out_dct["subdir"]}.@json')
        outputs_logger.info(expected_outputs)

    @overload
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: None = None,
        report_fetched: Literal[False] = False,
        optional: bool = False,
    ) -> Optional[STRAT_DICT]: ...
    @overload
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: PIPE_IDX,
        report_fetched: Literal[False] = False,
        optional: bool = False,
    ) -> Optional[Resource]: ...
    @overload
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: None = None,
        *,
        report_fetched: Literal[True],
        optional: bool = False,
    ) -> tuple[Optional[STRAT_DICT], Optional[str]]: ...
    @overload
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: PIPE_IDX,
        report_fetched: Literal[True],
        optional: Literal[False],
    ) -> tuple[Optional[Resource], Optional[str]]: ...
    @overload
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[PIPE_IDX] = None,
        report_fetched: bool = False,
        optional: bool = False,
    ) -> (
        Optional[Resource | STRAT_DICT]
        | tuple[Optional[Resource | STRAT_DICT], Optional[str]]
    ): ...
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[PIPE_IDX] = None,
        report_fetched: bool = False,
        optional: bool = False,
    ):
        """Return a dictionary of strats.

        Inside those, are dictionaries like ``{'data': (node, out), 'json': info}``.
        """
        return super().get(resource, pipe_idx, report_fetched, optional)

    @overload
    def get_data(
        self,
        resource: list[str] | str,
        pipe_idx: None = None,
        report_fetched: bool = False,
        quick_single: bool = False,
    ) -> ResourceData: ...
    @overload
    def get_data(
        self,
        resource: list[str] | str,
        pipe_idx: PIPE_IDX,
        report_fetched: Literal[True],
        quick_single: Literal[False] = False,
    ) -> tuple[ResourceData, str]: ...
    @overload
    def get_data(
        self,
        resource: list[str] | str,
        pipe_idx: PIPE_IDX,
        report_fetched: Literal[False] = False,
        quick_single: bool = False,
    ) -> ResourceData: ...
    @overload
    def get_data(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[PIPE_IDX],
        report_fetched: bool,
        quick_single: Literal[True],
    ) -> ResourceData: ...
    def get_data(
        self,
        resource,
        pipe_idx,
        report_fetched,
        quick_single,
    ):
        """Get ResourceData from ResourcePool."""
        _resource = self.get(resource, pipe_idx=pipe_idx, report_fetched=report_fetched)
        if report_fetched:
            if pipe_idx:
                connect, fetched = _resource
                assert isinstance(connect, Resource) and isinstance(fetched, str)
                return connect.data, fetched
        if quick_single or len(_resource) == 1:
            assert isinstance(_resource, dict)
            for value in _resource.values():
                return value.data
        assert isinstance(_resource, Resource)
        return _resource.data

    def get_strats(self, resources, debug=False) -> dict[str | tuple, "StratPool"]:
        # TODO: NOTE: NOT COMPATIBLE WITH SUB-RPOOL/STRAT_POOLS
        # TODO: (and it doesn't have to be)
        import itertools

        linked_resources = []
        resource_list = []
        if debug:
            verbose_logger = getLogger("CPAC.engine")
            verbose_logger.debug("\nresources: %s", resources)
        for resource in resources:
            # grab the linked-input tuples
            if isinstance(resource, tuple):
                linked = []
                for label in list(resource):
                    rp_dct, fetched_resource = self.get(
                        label, report_fetched=True, optional=True
                    )
                    if not rp_dct:
                        continue
                    linked.append(fetched_resource)
                resource_list += linked
                if len(linked) < 2:  # noqa: PLR2004
                    continue
                linked_resources.append(linked)
            else:
                resource_list.append(resource)

        total_pool = []
        variant_pool: dict = {}
        len_inputs = len(resource_list)
        if debug:
            verbose_logger = getLogger("CPAC.engine")
            verbose_logger.debug("linked_resources: %s", linked_resources)
            verbose_logger.debug("resource_list: %s", resource_list)
        for resource in resource_list:
            (
                rp_dct,  # <---- rp_dct has the strats/pipe_idxs as the keys on first level, then 'data' and 'json' on each strat level underneath
                fetched_resource,
            ) = self.get(
                resource,
                report_fetched=True,
                optional=True,  # oh, and we make the resource fetching in get_strats optional so we can have optional inputs, but they won't be optional in the node block unless we want them to be
            )
            if not rp_dct:
                len_inputs -= 1
                continue
            assert isinstance(rp_dct, dict)
            sub_pool = []
            if debug:
                verbose_logger.debug("len(rp_dct): %s\n", len(rp_dct))
            for strat in rp_dct.keys():
                json_info = self.get_json(fetched_resource, strat)
                cpac_prov = json_info["CpacProvenance"]
                sub_pool.append(cpac_prov)
                if fetched_resource not in variant_pool:
                    variant_pool[fetched_resource] = []
                if "CpacVariant" in json_info:
                    for key, val in json_info["CpacVariant"].items():
                        if val not in variant_pool[fetched_resource]:
                            variant_pool[fetched_resource] += val
                            variant_pool[fetched_resource].append(f"NO-{val[0]}")

            if debug:
                verbose_logger = getLogger("CPAC.engine")
                verbose_logger.debug("%s sub_pool: %s\n", resource, sub_pool)
            total_pool.append(sub_pool)

        if not total_pool:
            raise LookupError(
                "\n\n[!] C-PAC says: None of the listed "
                "resources in the node block being connected "
                "exist in the resource pool.\n\nResources:\n"
                "%s\n\n" % resource_list
            )

        # TODO: right now total_pool is:
        # TODO:    [[[T1w:anat_ingress, desc-preproc_T1w:anatomical_init, desc-preproc_T1w:acpc_alignment], [T1w:anat_ingress,desc-preproc_T1w:anatomical_init]],
        # TODO:     [[T1w:anat_ingress, desc-preproc_T1w:anatomical_init, desc-preproc_T1w:acpc_alignment, desc-brain_mask:brain_mask_afni], [T1w:anat_ingress, desc-preproc_T1w:anatomical_init, desc-brain_mask:brain_mask_afni]]]

        # TODO: and the code below thinks total_pool is a list of lists, like [[pipe_idx, pipe_idx], [pipe_idx, pipe_idx, pipe_idx], etc.]
        # TODO: and the actual resource is encoded in the tag: of the last item, every time!
        # keying the strategies to the resources, inverting it
        if len_inputs > 1:
            strats = itertools.product(*total_pool)

            # we now currently have "strats", the combined permutations of all the strategies, as a list of tuples, each tuple combining one version of input each, being one of the permutations.
            # OF ALL THE DIFFERENT INPUTS. and they are tagged by their fetched inputs with {name}:{strat}.
            # so, each tuple has ONE STRAT FOR EACH INPUT, so if there are three inputs, each tuple will have 3 items.
            new_strats: dict[str | tuple, StratPool] = {}

            # get rid of duplicates - TODO: refactor .product
            strat_str_list = []
            strat_list_list = []
            for strat_tuple in strats:
                strat_list = list(copy.deepcopy(strat_tuple))
                strat_str = str(strat_list)
                if strat_str not in strat_str_list:
                    strat_str_list.append(strat_str)
                    strat_list_list.append(strat_list)

            if debug:
                verbose_logger = getLogger("CPAC.engine")
                verbose_logger.debug("len(strat_list_list): %s\n", len(strat_list_list))
            for strat_list in strat_list_list:
                json_dct = {}
                for strat in strat_list:
                    # strat is a prov list for a single resource/input
                    strat_resource, strat_idx = self.generate_prov_string(strat)
                    strat_json = self.get_json(strat_resource, strat=strat_idx)
                    json_dct[strat_resource] = strat_json

                drop = False
                if linked_resources:
                    for linked in linked_resources:  # <--- 'linked' is each tuple
                        if drop:
                            break
                        for xlabel in linked:
                            if drop or xlabel is None:
                                break
                            xjson = copy.deepcopy(json_dct[xlabel])
                            for ylabel in linked:
                                if xlabel == ylabel or ylabel is None:
                                    continue
                                yjson = copy.deepcopy(json_dct[ylabel])

                                if "CpacVariant" not in xjson:
                                    xjson["CpacVariant"] = {}
                                if "CpacVariant" not in yjson:
                                    yjson["CpacVariant"] = {}

                                current_strat = []
                                for key, val in xjson["CpacVariant"].items():
                                    if isinstance(val, list):
                                        current_strat.append(val[0])
                                    else:
                                        current_strat.append(val)
                                current_spread = list(set(variant_pool[xlabel]))
                                for spread_label in current_spread:
                                    if "NO-" in spread_label:
                                        continue
                                    if spread_label not in current_strat:
                                        current_strat.append(f"NO-{spread_label}")

                                other_strat = []
                                for key, val in yjson["CpacVariant"].items():
                                    if isinstance(val, list):
                                        other_strat.append(val[0])
                                    else:
                                        other_strat.append(val)
                                other_spread = list(set(variant_pool[ylabel]))
                                for spread_label in other_spread:
                                    if "NO-" in spread_label:
                                        continue
                                    if spread_label not in other_strat:
                                        other_strat.append(f"NO-{spread_label}")

                                for variant in current_spread:
                                    in_current_strat = False
                                    in_other_strat = False
                                    in_other_spread = False

                                    if variant is None:
                                        in_current_strat = True
                                        if None in other_spread:
                                            in_other_strat = True
                                    if variant in current_strat:
                                        in_current_strat = True
                                    if variant in other_strat:
                                        in_other_strat = True
                                    if variant in other_spread:
                                        in_other_spread = True

                                    if not in_other_strat:
                                        if in_other_spread:
                                            if in_current_strat:
                                                drop = True
                                                break

                                    if in_other_strat:
                                        if in_other_spread:
                                            if not in_current_strat:
                                                drop = True
                                                break
                                if drop:
                                    break
                if drop:
                    continue

                # make the merged strat label from the multiple inputs
                # strat_list is actually the merged CpacProvenance lists
                pipe_idx = str(strat_list)
                new_strats[pipe_idx] = StratPool(name=pipe_idx)
                # new_strats is A DICTIONARY OF RESOURCEPOOL OBJECTS!
                # placing JSON info at one level higher only for copy convenience
                new_strats[pipe_idx].rpool["json"] = {}
                new_strats[pipe_idx].rpool["json"]["subjson"] = {}
                new_strats[pipe_idx].rpool["json"]["CpacProvenance"] = strat_list

                # now just invert resource:strat to strat:resource for each resource:strat
                for cpac_prov in strat_list:
                    resource, strat = self.generate_prov_string(cpac_prov)
                    resource_strat_dct = self.rpool[resource][strat]
                    # remember, `resource_strat_dct` is a Resource.
                    new_strats[pipe_idx].rpool[resource] = resource_strat_dct
                    # `new_strats` is A DICTIONARY OF RESOURCEPOOL OBJECTS! each one is a new slice of the resource pool combined together.
                    self.pipe_list.append(pipe_idx)
                    if "CpacVariant" in resource_strat_dct["json"]:
                        if "CpacVariant" not in new_strats[pipe_idx].rpool["json"]:
                            new_strats[pipe_idx].rpool["json"]["CpacVariant"] = {}
                        _variant = new_strats[pipe_idx].rpool["json"]["CpacVariant"]
                        assert isinstance(_variant, dict)
                        for younger_resource, variant_list in _variant.items():
                            if (
                                younger_resource
                                not in new_strats[pipe_idx].rpool["json"]["CpacVariant"]
                            ):
                                new_strats[pipe_idx].rpool["json"]["CpacVariant"][
                                    younger_resource
                                ] = variant_list
                    # preserve each input's JSON info also
                    data_type = resource.split("_")[-1]
                    if data_type not in new_strats[pipe_idx].rpool["json"]["subjson"]:
                        new_strats[pipe_idx].rpool["json"]["subjson"][data_type] = {}
                    new_strats[pipe_idx].rpool["json"]["subjson"][data_type].update(
                        copy.deepcopy(resource_strat_dct["json"])
                    )
        else:
            new_strats = {}
            for resource_strat_list in total_pool:
                # total_pool will have only one list of strats, for the one input
                for cpac_prov in resource_strat_list:  # <------- cpac_prov here doesn't need to be modified, because it's not merging with other inputs
                    resource, pipe_idx = self.generate_prov_string(cpac_prov)
                    resource_strat_dct = self.rpool[resource][pipe_idx]
                    # remember, `resource_strat_dct` is a Resource.
                    new_strats[pipe_idx] = StratPool(
                        rpool={resource: resource_strat_dct}, name=pipe_idx
                    )  # <----- again, new_strats is A DICTIONARY OF RESOURCEPOOL OBJECTS!
                    # placing JSON info at one level higher only for copy convenience
                    new_strats[pipe_idx].rpool["json"] = resource_strat_dct["json"]
                    # TODO: WARNING- THIS IS A LEVEL HIGHER THAN THE ORIGINAL 'JSON' FOR EASE OF ACCESS IN CONNECT_BLOCK WITH THE .GET(JSON)
                    new_strats[pipe_idx].rpool["json"]["subjson"] = {}
                    new_strats[pipe_idx].rpool["json"]["CpacProvenance"] = cpac_prov
                    # preserve each input's JSON info also
                    data_type = resource.split("_")[-1]
                    if data_type not in new_strats[pipe_idx].rpool["json"]["subjson"]:
                        new_strats[pipe_idx].rpool["json"]["subjson"][data_type] = {}
                    _json = new_strats[pipe_idx].rpool["json"]
                    new_strats[pipe_idx].rpool["json"]["subjson"][data_type].update(
                        copy.deepcopy(resource_strat_dct["json"])
                    )
        return new_strats

    def ingress_freesurfer(self) -> None:
        """Ingress FreeSurfer data."""
        try:
            fs_path = os.path.join(
                self.cfg.pipeline_setup["freesurfer_dir"],  # type: ignore[attr-defined]
                self.part_id,
            )
        except KeyError:
            WFLOGGER.warning("No FreeSurfer data present.")
            return

        # fs_path = os.path.join(cfg.pipeline_setup['freesurfer_dir'], part_id)
        if not os.path.exists(fs_path):
            if "sub" in self.part_id:
                fs_path = os.path.join(
                    self.cfg.pipeline_setup["freesurfer_dir"],  # type: ignore[attr-defined]
                    self.part_id.replace("sub-", ""),
                )
            else:
                fs_path = os.path.join(
                    self.cfg.pipeline_setup["freesurfer_dir"],  # type: ignore[attr-defined]
                    ("sub-" + self.part_id),
                )

            # patch for flo-specific data
            if not os.path.exists(fs_path):
                subj_ses = f"{self.part_id}-{self.ses_id}"
                fs_path = os.path.join(
                    self.cfg.pipeline_setup["freesurfer_dir"],  # type: ignore[attr-defined]
                    subj_ses,
                )
                if not os.path.exists(fs_path):
                    WFLOGGER.info(
                        "No FreeSurfer data found for subject %s", self.part_id
                    )
                    return

        # Check for double nested subj names
        if os.path.exists(os.path.join(fs_path, os.path.basename(fs_path))):
            fs_path = os.path.join(fs_path, self.part_id)

        fs_ingress = create_general_datasource("gather_freesurfer_dir")
        fs_ingress.inputs.inputnode.set(
            unique_id=self.unique_id,
            data=fs_path,
            creds_path=self.creds_path,
            dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],  # type: ignore[attr-defined]
        )
        self.set_data(
            "freesurfer-subject-dir",
            fs_ingress,
            "outputspec.data",
            {},
            "",
            "freesurfer_config_ingress",
        )

        recon_outs = {
            "pipeline-fs_raw-average": "mri/rawavg.mgz",
            "pipeline-fs_subcortical-seg": "mri/aseg.mgz",
            "pipeline-fs_brainmask": "mri/brainmask.mgz",
            "pipeline-fs_wmparc": "mri/wmparc.mgz",
            "pipeline-fs_T1": "mri/T1.mgz",
            "pipeline-fs_hemi-L_desc-surface_curv": "surf/lh.curv",
            "pipeline-fs_hemi-R_desc-surface_curv": "surf/rh.curv",
            "pipeline-fs_hemi-L_desc-surfaceMesh_pial": "surf/lh.pial",
            "pipeline-fs_hemi-R_desc-surfaceMesh_pial": "surf/rh.pial",
            "pipeline-fs_hemi-L_desc-surfaceMesh_smoothwm": "surf/lh.smoothwm",
            "pipeline-fs_hemi-R_desc-surfaceMesh_smoothwm": "surf/rh.smoothwm",
            "pipeline-fs_hemi-L_desc-surfaceMesh_sphere": "surf/lh.sphere",
            "pipeline-fs_hemi-R_desc-surfaceMesh_sphere": "surf/rh.sphere",
            "pipeline-fs_hemi-L_desc-surfaceMap_sulc": "surf/lh.sulc",
            "pipeline-fs_hemi-R_desc-surfaceMap_sulc": "surf/rh.sulc",
            "pipeline-fs_hemi-L_desc-surfaceMap_thickness": "surf/lh.thickness",
            "pipeline-fs_hemi-R_desc-surfaceMap_thickness": "surf/rh.thickness",
            "pipeline-fs_hemi-L_desc-surfaceMap_volume": "surf/lh.volume",
            "pipeline-fs_hemi-R_desc-surfaceMap_volume": "surf/rh.volume",
            "pipeline-fs_hemi-L_desc-surfaceMesh_white": "surf/lh.white",
            "pipeline-fs_hemi-R_desc-surfaceMesh_white": "surf/rh.white",
            "pipeline-fs_xfm": "mri/transforms/talairach.lta",
        }

        for key, outfile in recon_outs.items():
            fullpath = os.path.join(fs_path, outfile)
            if os.path.exists(fullpath):
                fs_ingress = create_general_datasource(f"gather_fs_{key}_dir")
                fs_ingress.inputs.inputnode.set(
                    unique_id=self.unique_id,
                    data=fullpath,
                    creds_path=self.creds_path,
                    dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],  # type: ignore[attr-defined]
                )
                self.set_data(
                    key, fs_ingress, "outputspec.data", {}, "", f"fs_{key}_ingress"
                )
            else:
                warnings.warn(
                    str(LookupError(f"\n[!] Path does not exist for {fullpath}.\n"))
                )

        return

    def ingress_output_dir(self) -> None:
        """Ingress an output directory into a ResourcePool."""
        dir_path = self.data_paths.derivatives_dir
        assert dir_path is not None
        WFLOGGER.info("\nPulling outputs from %s.\n", dir_path)

        anat = os.path.join(dir_path, "anat")
        func = os.path.join(dir_path, "func")

        outdir_anat: list[str] = []
        outdir_func: list[str] = []
        func_paths: dict = {}
        func_dict: dict = {}
        func_key = ""

        for subdir in [anat, func]:
            if os.path.isdir(subdir):
                for filename in os.listdir(subdir):
                    for ext in EXTS:
                        if ext in filename:
                            if subdir == anat:
                                outdir_anat.append(os.path.join(subdir, filename))
                            else:
                                outdir_func.append(os.path.join(subdir, filename))

        # Add derivatives directory to rpool
        ingress = create_general_datasource("gather_derivatives_dir")
        ingress.inputs.inputnode.set(
            unique_id=self.unique_id,
            data=dir_path,
            creds_path=self.creds_path,
            dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],  # type: ignore[attr-defined]
        )
        self.set_data(
            "derivatives-dir",
            ingress,
            "outputspec.data",
            {},
            "",
            "outdir_config_ingress",
        )

        for subdirs in [outdir_anat, outdir_func]:
            for filepath in subdirs:
                filename = str(filepath)
                for ext in EXTS:
                    filename = filename.split("/")[-1].replace(ext, "")

                data_label = filename.split(self.unique_id)[1].lstrip("_")

                if len(filename) == len(data_label):
                    msg = (
                        "\n\n[!] Possibly wrong participant or "
                        "session in this directory?\n\n"
                        f"Filepath: {filepath}\n\n"
                    )
                    raise Exception(msg)

                bidstag = ""
                for tag in data_label.split("_"):
                    for prefix in ["task-", "run-", "acq-", "rec"]:
                        if tag.startswith(prefix):
                            bidstag += f"{tag}_"
                            data_label = data_label.replace(f"{tag}_", "")
                data_label, json = strip_template(data_label)

                json_info, pipe_idx, node_name, data_label = self.json_outdir_ingress(
                    filepath, data_label, json
                )

                if (
                    "template" in data_label
                    and not json_info["Template"]
                    == self.cfg.pipeline_setup["outdir_ingress"]["Template"]  # type: ignore[attr-defined]
                ):
                    continue
                # Rename confounds to avoid confusion in nuisance regression
                if data_label.endswith("desc-confounds_timeseries"):
                    data_label = "pipeline-ingress_desc-confounds_timeseries"

                if len(bidstag) > 1:
                    # Remove tail symbol
                    bidstag = bidstag[:-1]
                    if bidstag.startswith("task-"):
                        bidstag = bidstag.replace("task-", "")

                # Rename bold mask for CPAC naming convention
                # and to avoid collision with anat brain mask
                if data_label.endswith("desc-brain_mask") and filepath in outdir_func:
                    data_label = data_label.replace("brain_mask", "bold_mask")

                try:
                    pipe_x = self.get_pipe_number(pipe_idx)
                except ValueError:
                    pipe_x = len(self.pipe_list)
                if filepath in outdir_anat:
                    ingress = create_general_datasource(
                        f"gather_anat_outdir_{data_label!s}_{pipe_x}"
                    )
                    ingress.inputs.inputnode.set(
                        unique_id=self.unique_id,
                        data=filepath,
                        creds_path=self.creds_path,
                        dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],  # type: ignore[attr-defined]
                    )
                    self.set_data(
                        data_label,
                        ingress,
                        "outputspec.data",
                        json_info,
                        pipe_idx,
                        node_name=f"outdir_{data_label}_ingress",
                        inject=True,
                    )
                else:
                    if data_label.endswith("desc-preproc_bold"):
                        func_key = data_label
                        func_dict[bidstag] = {}
                        func_dict[bidstag]["scan"] = str(filepath)
                        func_dict[bidstag]["scan_parameters"] = json_info
                        func_dict[bidstag]["pipe_idx"] = pipe_idx
                    if data_label.endswith("desc-brain_mask"):
                        data_label = data_label.replace("brain_mask", "bold_mask")
                    try:
                        func_paths[data_label].append(filepath)
                    except (AttributeError, KeyError, TypeError):
                        func_paths[data_label] = []
                        func_paths[data_label].append(filepath)

        if func_dict:
            self.func_outdir_ingress(
                func_dict,
                func_key,
                func_paths,
            )

        if self.cfg.surface_analysis["freesurfer"]["ingress_reconall"]:  # type: ignore[attr-defined]
            self.ingress_freesurfer()

    def ingress_func_metadata(
        self,
        num_strat=None,
    ) -> tuple[bool, bool, list[str]]:
        """Ingress metadata for functional scans."""
        name_suffix = ""
        for suffix_part in (self.unique_id, num_strat):
            if suffix_part is not None:
                name_suffix += f"_{suffix_part}"
        # Grab field maps
        diff = False
        blip = False
        fmap_rp_list = []
        fmap_TE_list = []
        if self.data_paths.fmap:
            second = False
            for orig_key in self.data_paths.fmap:
                gather_fmap = create_fmap_datasource(
                    self.data_paths.fmap, f"fmap_gather_{orig_key}_{self.part_id}"
                )
                gather_fmap.inputs.inputnode.set(
                    subject=self.part_id,
                    creds_path=self.creds_path,
                    dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],  # type: ignore[attr-defined]
                )
                gather_fmap.inputs.inputnode.scan = orig_key

                key = orig_key
                if "epi" in key and not second:
                    key = "epi-1"
                    second = True
                elif "epi" in key and second:
                    key = "epi-2"

                self.set_data(
                    key, gather_fmap, "outputspec.rest", {}, "", "fmap_ingress"
                )
                self.set_data(
                    f"{key}-scan-params",
                    gather_fmap,
                    "outputspec.scan_params",
                    {},
                    "",
                    "fmap_params_ingress",
                )

                fmap_rp_list.append(key)

                get_fmap_metadata_imports = ["import json"]
                get_fmap_metadata = pe.Node(
                    Function(
                        input_names=["data_config_scan_params"],
                        output_names=[
                            "dwell_time",
                            "pe_direction",
                            "total_readout",
                            "echo_time",
                            "echo_time_one",
                            "echo_time_two",
                        ],
                        function=get_fmap_phasediff_metadata,
                        imports=get_fmap_metadata_imports,
                    ),
                    name=f"{key}_get_metadata{name_suffix}",
                )

                self.wf.connect(
                    gather_fmap,
                    "outputspec.scan_params",
                    get_fmap_metadata,
                    "data_config_scan_params",
                )

                if "phase" in key:
                    # leave it open to all three options, in case there is a
                    # phasediff image with either a single EchoTime field (which
                    # usually matches one of the magnitude EchoTimes), OR
                    # a phasediff with an EchoTime1 and EchoTime2

                    # at least one of these rpool keys will have a None value,
                    # which will be sorted out in gather_echo_times below
                    self.set_data(
                        f"{key}-TE",
                        get_fmap_metadata,
                        "echo_time",
                        {},
                        "",
                        "fmap_TE_ingress",
                    )
                    fmap_TE_list.append(f"{key}-TE")

                    self.set_data(
                        f"{key}-TE1",
                        get_fmap_metadata,
                        "echo_time_one",
                        {},
                        "",
                        "fmap_TE1_ingress",
                    )
                    fmap_TE_list.append(f"{key}-TE1")

                    self.set_data(
                        f"{key}-TE2",
                        get_fmap_metadata,
                        "echo_time_two",
                        {},
                        "",
                        "fmap_TE2_ingress",
                    )
                    fmap_TE_list.append(f"{key}-TE2")

                elif "magnitude" in key:
                    self.set_data(
                        f"{key}-TE",
                        get_fmap_metadata,
                        "echo_time",
                        {},
                        "",
                        "fmap_TE_ingress",
                    )
                    fmap_TE_list.append(f"{key}-TE")

                self.set_data(
                    f"{key}-dwell",
                    get_fmap_metadata,
                    "dwell_time",
                    {},
                    "",
                    "fmap_dwell_ingress",
                )
                self.set_data(
                    f"{key}-pedir",
                    get_fmap_metadata,
                    "pe_direction",
                    {},
                    "",
                    "fmap_pedir_ingress",
                )
                self.set_data(
                    f"{key}-total-readout",
                    get_fmap_metadata,
                    "total_readout",
                    {},
                    "",
                    "fmap_readout_ingress",
                )

                if "phase" in key or "mag" in key:
                    diff = True

                if re.match("epi_[AP]{2}", orig_key):
                    blip = True

            if diff:
                calc_delta_ratio = pe.Node(
                    Function(
                        input_names=["effective_echo_spacing", "echo_times"],
                        output_names=["deltaTE", "ees_asym_ratio"],
                        function=calc_delta_te_and_asym_ratio,
                        imports=["from typing import Optional"],
                    ),
                    name=f"diff_distcor_calc_delta{name_suffix}",
                )

                gather_echoes = pe.Node(
                    Function(
                        input_names=[
                            "echotime_1",
                            "echotime_2",
                            "echotime_3",
                            "echotime_4",
                        ],
                        output_names=["echotime_list"],
                        function=gather_echo_times,
                    ),
                    name="fugue_gather_echo_times",
                )

                for idx, fmap_file in enumerate(fmap_TE_list, start=1):
                    try:
                        node, out_file = self.get_data(
                            fmap_file, f"['{fmap_file}:fmap_TE_ingress']"
                        )
                        self.wf.connect(
                            node, out_file, gather_echoes, f"echotime_{idx}"
                        )
                    except KeyError:
                        pass

                self.wf.connect(
                    gather_echoes, "echotime_list", calc_delta_ratio, "echo_times"
                )

        # Add in nodes to get parameters from configuration file
        # a node which checks if scan_parameters are present for each scan
        scan_params = pe.Node(
            Function(
                input_names=[
                    "data_config_scan_params",
                    "subject_id",
                    "scan",
                    "pipeconfig_tr",
                    "pipeconfig_tpattern",
                    "pipeconfig_start_indx",
                    "pipeconfig_stop_indx",
                ],
                output_names=[
                    "tr",
                    "tpattern",
                    "template",
                    "ref_slice",
                    "start_indx",
                    "stop_indx",
                    "pe_direction",
                    "effective_echo_spacing",
                ],
                function=get_scan_params,
                imports=["from CPAC.utils.utils import check, try_fetch_parameter"],
            ),
            name=f"bold_scan_params_{self.part_id}{name_suffix}",
        )
        scan_params.inputs.subject_id = self.part_id
        scan_params.inputs.set(
            pipeconfig_start_indx=self.cfg.functional_preproc["truncation"]["start_tr"],  # type: ignore[attr-defined]
            pipeconfig_stop_indx=self.cfg.functional_preproc["truncation"]["stop_tr"],  # type: ignore[attr-defined]
        )

        node, out = self.get_data("scan", "['scan:func_ingress']")
        self.wf.connect(node, out, scan_params, "scan")

        # Workaround for extracting metadata with ingress
        if self.check_rpool("derivatives-dir"):
            selectrest_json = pe.Node(
                Function(
                    input_names=["scan", "rest_dict", "resource"],
                    output_names=["file_path"],
                    function=get_rest,
                    as_module=True,
                ),
                name="selectrest_json",
            )
            selectrest_json.inputs.rest_dict = self.data_paths.as_dict()
            selectrest_json.inputs.resource = "scan_parameters"
            self.wf.connect(node, out, selectrest_json, "scan")
            self.wf.connect(
                selectrest_json, "file_path", scan_params, "data_config_scan_params"
            )

        else:
            # wire in the scan parameter workflow
            node, out = self.get_data(
                "scan-params", "['scan-params:scan_params_ingress']"
            )
            self.wf.connect(node, out, scan_params, "data_config_scan_params")

        self.set_data("TR", scan_params, "tr", {}, "", "func_metadata_ingress")
        self.set_data(
            "tpattern", scan_params, "tpattern", {}, "", "func_metadata_ingress"
        )
        self.set_data(
            "template", scan_params, "template", {}, "", "func_metadata_ingress"
        )
        self.set_data(
            "start-tr", scan_params, "start_indx", {}, "", "func_metadata_ingress"
        )
        self.set_data(
            "stop-tr", scan_params, "stop_indx", {}, "", "func_metadata_ingress"
        )
        self.set_data(
            "pe-direction", scan_params, "pe_direction", {}, "", "func_metadata_ingress"
        )

        if diff:
            # Connect EffectiveEchoSpacing from functional metadata
            self.set_data(
                "effectiveEchoSpacing",
                scan_params,
                "effective_echo_spacing",
                {},
                "",
                "func_metadata_ingress",
            )
            node, out_file = self.get_data(
                "effectiveEchoSpacing", "['effectiveEchoSpacing:func_metadata_ingress']"
            )
            self.wf.connect(node, out_file, calc_delta_ratio, "effective_echo_spacing")
            self.set_data(
                "deltaTE", calc_delta_ratio, "deltaTE", {}, "", "deltaTE_ingress"
            )
            self.set_data(
                "ees-asym-ratio",
                calc_delta_ratio,
                "ees_asym_ratio",
                {},
                "",
                "ees_asym_ratio_ingress",
            )

        return diff, blip, fmap_rp_list

    def ingress_pipeconfig_paths(self):
        """Ingress config file paths."""
        # TODO: may want to change the resource keys for each to include one level up in the YAML as well

        import pandas as pd
        import pkg_resources as p

        template_csv = p.resource_filename("CPAC", "resources/cpac_templates.csv")
        template_df = pd.read_csv(template_csv, keep_default_na=False)

        for row in template_df.itertuples():
            key = row.Key
            val = row.Pipeline_Config_Entry
            val = self.cfg.get_nested(self.cfg, [x.lstrip() for x in val.split(",")])
            resolution = row.Intended_Resolution_Config_Entry
            desc = row.Description

            if not val:
                continue

            if resolution:
                res_keys = [x.lstrip() for x in resolution.split(",")]
                tag = res_keys[-1]
            json_info = {}

            if "$FSLDIR" in val:
                val = val.replace(
                    "$FSLDIR", self.cfg.pipeline_setup["system_config"]["FSLDIR"]
                )
            if "$priors_path" in val:
                priors_path = (
                    self.cfg.segmentation["tissue_segmentation"]["FSL-FAST"][
                        "use_priors"
                    ]["priors_path"]
                    or ""
                )
                if "$FSLDIR" in priors_path:
                    priors_path = priors_path.replace(
                        "$FSLDIR", self.cfg.pipeline_setup["system_config"]["FSLDIR"]
                    )
                val = val.replace("$priors_path", priors_path)
            if "${resolution_for_anat}" in val:
                val = val.replace(
                    "${resolution_for_anat}",
                    self.cfg.registration_workflows["anatomical_registration"][
                        "resolution_for_anat"
                    ],
                )
            if "${func_resolution}" in val:
                val = val.replace(
                    "${func_resolution}",
                    self.cfg.registration_workflows["functional_registration"][
                        "func_registration_to_template"
                    ]["output_resolution"][tag],
                )

            if desc:
                template_name, _template_desc = lookup_identifier(val)
                if template_name:
                    desc = f"{template_name} - {desc}"
                json_info["Description"] = f"{desc} - {val}"
            if resolution:
                resolution = self.cfg.get_nested(self.cfg, res_keys)
                json_info["Resolution"] = resolution

                resampled_template = pe.Node(
                    Function(
                        input_names=["resolution", "template", "template_name", "tag"],
                        output_names=["resampled_template"],
                        function=resolve_resolution,
                        as_module=True,
                    ),
                    name="resampled_" + key,
                )

                resampled_template.inputs.resolution = resolution
                resampled_template.inputs.template = val
                resampled_template.inputs.template_name = key
                resampled_template.inputs.tag = tag

                # the set_data below is set up a little differently, because we are
                # injecting and also over-writing already-existing entries
                # other alternative would have been to ingress into the
                # resampled_template node from the already existing entries, but we
                # didn't do that here
                self.set_data(
                    key,
                    resampled_template,
                    "resampled_template",
                    json_info,
                    "",
                    "template_resample",
                )  # pipe_idx (after the blank json {}) should be the previous strat that you want deleted! because you're not connecting this the regular way, you have to do it manually

            elif val:
                config_ingress = create_general_datasource(f"gather_{key}")
                config_ingress.inputs.inputnode.set(
                    unique_id=self.unique_id,
                    data=val,
                    creds_path=self.creds_path,
                    dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],
                )
                self.set_data(
                    key,
                    config_ingress,
                    "outputspec.data",
                    json_info,
                    "",
                    f"{key}_config_ingress",
                )

    def ingress_raw_func_data(self):
        """Ingress raw functional data."""
        func_paths_dct = self.data_paths.func

        func_wf = self.create_func_datasource(
            func_paths_dct, f"func_ingress_{self.part_id}_{self.ses_id}"
        )
        func_wf.inputs.inputnode.set(
            subject=self.part_id,
            creds_path=self.creds_path,
            dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],
        )
        func_wf.get_node("inputnode").iterables = ("scan", list(func_paths_dct.keys()))

        self.set_data("subject", func_wf, "outputspec.subject", {}, "", "func_ingress")
        self.set_data("bold", func_wf, "outputspec.rest", {}, "", "func_ingress")
        self.set_data("scan", func_wf, "outputspec.scan", {}, "", "func_ingress")
        self.set_data(
            "scan-params",
            func_wf,
            "outputspec.scan_params",
            {},
            "",
            "scan_params_ingress",
        )

        # TODO: CHECK FOR PARAMETERS

        diff, blip, fmap_rp_list = self.ingress_func_metadata()

        # Memoize list of local functional scans
        # TODO: handle S3 files
        # Skip S3 files for now

        local_func_scans = [
            func_paths_dct[scan]["scan"]
            for scan in func_paths_dct.keys()
            if not func_paths_dct[scan]["scan"].startswith("s3://")
        ]
        if local_func_scans:
            # pylint: disable=protected-access
            self.wf._local_func_scans = local_func_scans
            if self.cfg.pipeline_setup["Debugging"]["verbose"]:
                verbose_logger = getLogger("CPAC.engine")
                verbose_logger.debug("local_func_scans: %s", local_func_scans)
        del local_func_scans

        return diff, blip, fmap_rp_list

    def func_outdir_ingress(self, func_dict: dict, key: str, func_paths: dict) -> None:
        """Ingress a functional output directory."""
        pipe_x = len(self.pipe_list)
        ingress = self.create_func_datasource(
            func_dict, f"gather_func_outdir_{key}_{pipe_x}"
        )
        ingress.inputs.inputnode.set(
            subject=self.unique_id,
            creds_path=self.creds_path,
            dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],  # type: ignore[attr-defined]
        )
        self.set_data("subject", ingress, "outputspec.subject", {}, "", "func_ingress")
        ingress.get_node("inputnode").iterables = ("scan", list(func_dict.keys()))
        self.set_data(key, ingress, "outputspec.rest", {}, "", "func_ingress")

        self.set_data("scan", ingress, "outputspec.scan", {}, "", "func_ingress")
        self.set_data(
            "scan-params",
            ingress,
            "outputspec.scan_params",
            {},
            "",
            "scan_params_ingress",
        )
        self.ingress_func_metadata()

        # Have to do it this weird way to save the parsed BIDS tag & filepath
        mask_paths_key = (
            "desc-bold_mask"
            if "desc-bold_mask" in func_paths
            else "space-template_desc-bold_mask"
        )
        ts_paths_key = "pipeline-ingress_desc-confounds_timeseries"

        # Connect func data with approproate scan name
        iterables = pe.Node(
            Function(
                input_names=["scan", "mask_paths", "ts_paths"],
                output_names=["out_scan", "mask", "confounds"],
                function=set_iterables,
            ),
            name=f"set_iterables_{pipe_x}",
        )
        iterables.inputs.mask_paths = func_paths[mask_paths_key]
        iterables.inputs.ts_paths = func_paths[ts_paths_key]
        self.wf.connect(ingress, "outputspec.scan", iterables, "scan")

        for key in func_paths:
            if key in (mask_paths_key, ts_paths_key):
                ingress_func = create_general_datasource(f"ingress_func_data_{key}")
                ingress_func.inputs.inputnode.set(
                    unique_id=self.unique_id,
                    creds_path=self.creds_path,
                    dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],  # type: ignore[attr-defined]
                )
                self.wf.connect(iterables, "out_scan", ingress_func, "inputnode.scan")
                if key == mask_paths_key:
                    self.wf.connect(iterables, "mask", ingress_func, "inputnode.data")
                    self.set_data(
                        key,
                        ingress_func,
                        "inputnode.data",
                        {},
                        "",
                        f"outdir_{key}_ingress",
                    )
                elif key == ts_paths_key:
                    self.wf.connect(
                        iterables, "confounds", ingress_func, "inputnode.data"
                    )
                    self.set_data(
                        key,
                        ingress_func,
                        "inputnode.data",
                        {},
                        "",
                        f"outdir_{key}_ingress",
                    )

    def ingress_raw_anat_data(self) -> None:
        """Ingress raw anatomical data."""
        if not self.data_paths.anat:
            WFLOGGER.warning("No anatomical data present.")
            return

        anat_flow = create_anat_datasource(f"anat_T1w_gather_{self.unique_id}")

        anat = {}
        if "T1w" in self.data_paths.anat:
            anat["T1"] = self.data_paths.anat["T1w"]

        if "T1" in anat:
            anat_flow.inputs.inputnode.set(
                subject=self.part_id,
                anat=anat["T1"],
                creds_path=self.creds_path,
                dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],  # type: ignore[attr-defined]
                img_type="anat",
            )
            self.set_data("T1w", anat_flow, "outputspec.anat", {}, "", "anat_ingress")

        if "T2w" in self.data_paths.anat:
            anat_flow_T2 = create_anat_datasource(
                f"anat_T2w_gather_{self.part_id}_{self.ses_id}"
            )
            anat_flow_T2.inputs.inputnode.set(
                subject=self.part_id,
                anat=self.data_paths.anat["T2w"],
                creds_path=self.creds_path,
                dl_dir=self.cfg.pipeline_setup["working_directory"]["path"],  # type: ignore[attr-defined]
                img_type="anat",
            )
            self.set_data(
                "T2w", anat_flow_T2, "outputspec.anat", {}, "", "anat_ingress"
            )

        if self.cfg.surface_analysis["freesurfer"]["ingress_reconall"]:  # type: ignore[attr-defined]
            self.ingress_freesurfer()


class StratPool(_Pool):
    """A pool of ResourcePools keyed by strategy."""

    def __init__(
        self,
        rpool: Optional[dict] = None,
        *,
        name: str = "",
    ) -> None:
        """Initialize a StratPool."""
        super().__init__(name=name)
        if not rpool:
            self.rpool = STRAT_DICT({})
        else:
            self.rpool = STRAT_DICT(rpool)

    def append_name(self, name):
        self.name.append(name)

    @overload
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[PIPE_IDX] = None,
        report_fetched: Literal[False] = False,
        optional: bool = False,
    ) -> Optional[Resource]: ...
    @overload
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[PIPE_IDX],
        report_fetched: Literal[True],
        optional: bool = False,
    ) -> tuple[Optional[Resource], Optional[str]]: ...
    @overload
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[PIPE_IDX] = None,
        *,
        report_fetched: Literal[True],
        optional: bool = False,
    ) -> tuple[Optional[Resource], Optional[str]]: ...
    @overload
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[PIPE_IDX] = None,
        report_fetched: bool = False,
        optional: bool = False,
    ) -> Optional[Resource] | tuple[Optional[Resource], Optional[str]]: ...
    def get(
        self,
        resource: list[str] | str,
        pipe_idx: Optional[PIPE_IDX] = None,
        report_fetched: bool = False,
        optional: bool = False,
    ):
        """Return a Resource."""
        return super().get(resource, pipe_idx, report_fetched, optional)

    @overload
    def get_data(
        self, resource: list[str] | str, report_fetched: Literal[True]
    ) -> tuple[ResourceData, str]: ...
    @overload
    def get_data(
        self, resource: list[str] | str, report_fetched: Literal[False] = False
    ) -> ResourceData: ...
    def get_data(self, resource, report_fetched):
        """Get ResourceData from a StratPool."""
        _resource = self.get(resource, report_fetched=report_fetched)
        if report_fetched:
            assert isinstance(_resource, tuple)
            connect, fetched = _resource
            assert isinstance(connect, Resource) and isinstance(fetched, str)
            return connect.data, fetched
        assert isinstance(_resource, Resource)
        return _resource.data