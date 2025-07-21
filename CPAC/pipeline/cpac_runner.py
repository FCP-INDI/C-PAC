# Copyright (C) 2022-2024  C-PAC Developers

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
"""Run C-PAC."""

from multiprocessing import Process
import os
from pathlib import Path
from time import strftime
import warnings

from voluptuous.error import Invalid
import yaml

from CPAC.longitudinal.wf.anat import anat_longitudinal_wf
from CPAC.pipeline.utils import get_shell
from CPAC.utils.configuration import check_pname, Configuration, set_subject
from CPAC.utils.configuration.yaml_template import upgrade_pipeline_to_1_8
from CPAC.utils.ga import track_run
from CPAC.utils.monitoring import failed_to_start, init_loggers, log_nodes_cb, WFLOGGER


def run_condor_jobs(c, config_file, subject_list_file, p_name):
    """Run condor jobs."""
    # Import packages
    import subprocess
    from time import strftime

    try:
        sublist = yaml.safe_load(open(os.path.realpath(subject_list_file), "r"))
    except:
        msg = "Subject list is not in proper YAML format. Please check your file"
        raise Exception(msg)

    cluster_files_dir = os.path.join(os.getcwd(), "cluster_files")
    subject_bash_file = os.path.join(
        cluster_files_dir, "submit_%s.condor" % str(strftime("%Y_%m_%d_%H_%M_%S"))
    )
    f = open(subject_bash_file, "w")

    print("Executable = /usr/bin/python", file=f)
    print("Universe = vanilla", file=f)
    print("transfer_executable = False", file=f)
    print("getenv = True", file=f)
    print(
        "log = %s"
        % os.path.join(
            cluster_files_dir, "c-pac_%s.log" % str(strftime("%Y_%m_%d_%H_%M_%S"))
        ),
        file=f,
    )

    for sidx in range(1, len(sublist) + 1):
        print(
            "error = %s"
            % os.path.join(
                cluster_files_dir,
                "c-pac_%s.%s.err" % (str(strftime("%Y_%m_%d_%H_%M_%S")), str(sidx)),
            ),
            file=f,
        )
        print(
            "output = %s"
            % os.path.join(
                cluster_files_dir,
                "c-pac_%s.%s.out" % (str(strftime("%Y_%m_%d_%H_%M_%S")), str(sidx)),
            ),
            file=f,
        )

        print(
            "arguments = \"-c 'import CPAC; CPAC.pipeline.cpac_pipeline.run( ''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'')'\""
            % (
                str(config_file),
                subject_list_file,
                str(sidx),
                c.maskSpecificationFile,
                c.roiSpecificationFile,
                c.templateSpecificationFile,
                p_name,
            ),
            file=f,
        )
        print("queue", file=f)

    f.close()

    # commands.getoutput('chmod +x %s' % subject_bash_file )
    WFLOGGER.info(subprocess.getoutput("condor_submit %s ", subject_bash_file))


# Create and run script for CPAC to run on cluster
def run_cpac_on_cluster(config_file, subject_list_file, cluster_files_dir):
    """Build a batch job submission script and submit to the scheduler."""
    # Import packages
    import getpass
    import re
    import subprocess
    from time import strftime

    from indi_schedulers import cluster_templates

    # Load in pipeline config
    try:
        pipeline_dict = yaml.safe_load(open(os.path.realpath(config_file), "r"))
        pipeline_config = Configuration(pipeline_dict)
    except:
        msg = "Pipeline config is not in proper YAML format. Please check your file"
        raise Exception(msg)
    # Load in the subject list
    try:
        sublist = yaml.safe_load(open(os.path.realpath(subject_list_file), "r"))
    except:
        msg = "Subject list is not in proper YAML format. Please check your file"
        raise Exception(msg)

    # Init variables
    timestamp = str(strftime("%Y_%m_%d_%H_%M_%S"))
    job_scheduler = pipeline_config.pipeline_setup["system_config"]["on_grid"][
        "resource_manager"
    ].lower()

    # For SLURM time limit constraints only, hh:mm:ss
    hrs_limit = 8 * len(sublist)
    time_limit = "%d:00:00" % hrs_limit

    # Batch file variables
    user_account = getpass.getuser()
    num_subs = len(sublist)

    # Run CPAC via python -c command
    python_cpac_str = (
        'python -c "from CPAC.pipeline.cpac_pipeline import run; '
        "run('%(config_file)s', '%(subject_list_file)s', "
        "%(env_arr_idx)s, '%(pipeline_name)s', "
        "plugin='MultiProc', plugin_args=%(plugin_args)s)\""
    )

    # Init plugin arguments
    plugin_args = {
        "n_procs": pipeline_config.pipeline_setup["system_config"][
            "max_cores_per_participant"
        ],
        "memory_gb": pipeline_config.pipeline_setup["system_config"][
            "maximum_memory_per_participant"
        ],
        "raise_insufficient": pipeline_config.pipeline_setup["system_config"][
            "raise_insufficient"
        ],
        "status_callback": log_nodes_cb,
    }

    # Set up run command dictionary
    run_cmd_dict = {
        "config_file": config_file,
        "subject_list_file": subject_list_file,
        "pipeline_name": pipeline_config.pipeline_setup["pipeline_name"],
        "plugin_args": plugin_args,
    }

    # Set up config dictionary
    config_dict = {
        "timestamp": timestamp,
        "shell": get_shell(),
        "job_name": "CPAC_" + pipeline_config.pipeline_setup["pipeline_name"],
        "num_tasks": num_subs,
        "queue": pipeline_config.pipeline_setup["system_config"]["on_grid"]["SGE"][
            "queue"
        ],
        "par_env": pipeline_config.pipeline_setup["system_config"]["on_grid"]["SGE"][
            "parallel_environment"
        ],
        "cores_per_task": pipeline_config.pipeline_setup["system_config"][
            "max_cores_per_participant"
        ],
        "user": user_account,
        "work_dir": cluster_files_dir,
        "time_limit": time_limit,
    }

    # Get string template for job scheduler
    if job_scheduler == "pbs":
        env_arr_idx = "$PBS_ARRAYID"
        batch_file_contents = cluster_templates.pbs_template
        confirm_str = r"(?<=Your job-array )\d+"
        exec_cmd = "qsub"
    elif job_scheduler == "sge":
        env_arr_idx = "$SGE_TASK_ID"
        batch_file_contents = cluster_templates.sge_template
        confirm_str = r"(?<=Your job-array )\d+"
        exec_cmd = "qsub"
    elif job_scheduler == "slurm":
        env_arr_idx = "$SLURM_ARRAY_TASK_ID"
        batch_file_contents = cluster_templates.slurm_template
        confirm_str = r"(?<=Submitted batch job )\d+"
        exec_cmd = "sbatch"

    # Populate rest of dictionary
    config_dict["env_arr_idx"] = env_arr_idx
    run_cmd_dict["env_arr_idx"] = env_arr_idx
    config_dict["run_cmd"] = python_cpac_str % run_cmd_dict

    # Populate string from config dict values
    batch_file_contents = batch_file_contents % config_dict
    # Write file
    batch_filepath = os.path.join(
        cluster_files_dir, "cpac_submit_%s.%s" % (timestamp, job_scheduler)
    )
    with open(batch_filepath, "w") as f:
        f.write(batch_file_contents)

    # Get output response from job submission
    out = subprocess.getoutput("%s %s" % (exec_cmd, batch_filepath))

    # Check for successful qsub submission
    if re.search(confirm_str, out) is None:
        err_msg = "Error submitting C-PAC pipeline run to %s queue" % job_scheduler
        raise Exception(err_msg)

    # Get pid and send to pid file
    pid = re.search(confirm_str, out).group(0)
    pid_file = os.path.join(cluster_files_dir, "pid.txt")
    with open(pid_file, "w") as f:
        f.write(pid)


def run_T1w_longitudinal(sublist, cfg: Configuration, dry_run: bool = False):
    subject_id_dict = {}

    for sub in sublist:
        if sub["subject_id"] in subject_id_dict:
            subject_id_dict[sub["subject_id"]].append(sub)
        else:
            subject_id_dict[sub["subject_id"]] = [sub]

    # subject_id_dict has the subject_id as keys and a list of
    # sessions for each participant as value
    for subject_id, sub_list in subject_id_dict.items():
        if len(sub_list) > 1:
            log_dir: str
            _, _, log_dir = set_subject(sub_list[0], cfg)
            log_dir = str(Path(log_dir).parent / f"{subject_id}_longitudinal")
            init_loggers(subject_id, cfg, log_dir, mock=True)
            anat_longitudinal_wf(subject_id, sub_list, cfg, dry_run=dry_run)
        elif len(sub_list) == 1:
            warnings.warn(
                "\n\nThere is only one anatomical session "
                "for sub-%s. Longitudinal preprocessing "
                "will be skipped for this subject."
                "\n\n" % subject_id
            )


def run(
    subject_list_file,
    config_file=None,
    p_name=None,
    plugin=None,
    plugin_args=None,
    tracking=True,
    num_subs_at_once=None,
    debug=False,
    test_config=False,
) -> int:
    """Run C-PAC subjects via job queue.

    Returns
    -------
    int
        exit code
    """
    exitcode = 0

    # Import packages
    import os
    import time

    from CPAC.pipeline.cpac_pipeline import run_workflow

    WFLOGGER.info("Run called with config file %s", config_file)

    if plugin_args is None:
        plugin_args = {"status_callback": log_nodes_cb}

    if not config_file:
        import pkg_resources as p

        config_file = p.resource_filename(
            "CPAC", os.path.join("resources", "configs", "pipeline_config_template.yml")
        )

    # Init variables
    sublist = None
    if ".yaml" in subject_list_file or ".yml" in subject_list_file:
        subject_list_file = os.path.realpath(subject_list_file)
    else:
        from CPAC.utils.bids_utils import (
            bids_gen_cpac_sublist,
            collect_bids_files_configs,
        )

        (file_paths, config) = collect_bids_files_configs(subject_list_file, None)
        sublist = bids_gen_cpac_sublist(subject_list_file, file_paths, config, None)
        if not sublist:
            WFLOGGER.error("Did not find data in %s", subject_list_file)
            return 1

    # take date+time stamp for run identification purposes
    unique_pipeline_id = strftime("%Y%m%d%H%M%S")
    pipeline_start_stamp = strftime("%Y-%m-%d_%H:%M:%S")

    # Load in pipeline config file
    config_file = os.path.realpath(config_file)
    try:
        if not os.path.exists(config_file):
            raise IOError
        else:
            try:
                c = Configuration(yaml.safe_load(open(config_file, "r")))
            except Invalid:
                try:
                    upgrade_pipeline_to_1_8(config_file)
                    c = Configuration(yaml.safe_load(open(config_file, "r")))
                except Exception as e:
                    msg = (
                        "C-PAC could not upgrade pipeline configuration file "
                        f"{config_file} to v1.8 syntax"
                    )
                    raise RuntimeError(msg) from e
            except Exception as e:
                raise e
    except IOError as e:
        msg = f"config file {config_file} doesn't exist"
        raise FileNotFoundError(msg) from e
    except yaml.parser.ParserError as e:
        error_detail = '"%s" at line %d' % (e.problem, e.problem_mark.line)
        msg = (
            f"Error parsing config file: {config_file}\n\n"
            "Error details:\n"
            f"    {error_detail}"
            "\n\n"
        )
        raise Exception(msg)
    except Exception as e:
        msg = (
            f"Error parsing config file: {config_file}\n\n"
            "Error details:\n"
            f"    {e}"
            "\n\n"
        )
        raise Exception(msg)

    c.pipeline_setup["log_directory"]["path"] = os.path.abspath(
        c.pipeline_setup["log_directory"]["path"]
    )
    c.pipeline_setup["working_directory"]["path"] = os.path.abspath(
        c.pipeline_setup["working_directory"]["path"]
    )
    if "s3://" not in c.pipeline_setup["output_directory"]["path"]:
        c.pipeline_setup["output_directory"]["path"] = os.path.abspath(
            c.pipeline_setup["output_directory"]["path"]
        )
    if debug:
        c.pipeline_setup["output_directory"]["path"]["write_debugging_outputs"] = "[1]"

    if num_subs_at_once:
        if not str(num_subs_at_once).isdigit():
            msg = "[!] Value entered for --num_cores not a digit."
            raise Exception(msg)
        c.pipeline_setup["system_config"]["num_participants_at_once"] = int(
            num_subs_at_once
        )

    # Do some validation
    if not c.pipeline_setup["working_directory"]["path"]:
        msg = "Working directory not specified"
        raise Exception(msg)

    if len(c.pipeline_setup["working_directory"]["path"]) > 70:
        warnings.warn(
            "We recommend that the working directory full path "
            "should have less then 70 characters. "
            "Long paths might not work in your operating system."
        )
        warnings.warn(
            "Current working directory: "
            f"{c.pipeline_setup['working_directory']['path']}"
        )

    # Get the pipeline name
    p_name = check_pname(p_name, c)

    # Load in subject list
    try:
        if not sublist:
            sublist = yaml.safe_load(open(subject_list_file, "r"))
    except:
        msg = "Subject list is not in proper YAML format. Please check your file"
        raise FileNotFoundError(msg)

    # Populate subject scan map
    sub_scan_map = {}
    try:
        for sub in sublist:
            if sub["unique_id"]:
                s = sub["subject_id"] + "_" + sub["unique_id"]
            else:
                s = sub["subject_id"]
            scan_ids = ["scan_anat"]

            if "func" in sub:
                for id in sub["func"]:
                    scan_ids.append("scan_" + str(id))

            if "rest" in sub:
                for id in sub["rest"]:
                    scan_ids.append("scan_" + str(id))

            sub_scan_map[s] = scan_ids
    except Exception as e:
        msg = (
            "\n\nERROR: Subject list file not in proper format - check if you loaded"
            " the correct file?\nError name: cpac_runner_0001\n\n"
        )
        raise ValueError(msg) from e

    pipeline_timing_info = []
    pipeline_timing_info.append(unique_pipeline_id)
    pipeline_timing_info.append(pipeline_start_stamp)
    pipeline_timing_info.append(len(sublist))

    if tracking:
        try:
            track_run(
                level="participant" if not test_config else "test",
                participants=len(sublist),
            )
        except:
            WFLOGGER.error("Usage tracking failed for this run.")

    # If we're running on cluster, execute job scheduler
    if c.pipeline_setup["system_config"]["on_grid"]["run"]:
        # Create cluster log dir
        cluster_files_dir = os.path.join(
            c.pipeline_setup["log_directory"]["path"], "cluster_files"
        )
        if not os.path.exists(cluster_files_dir):
            os.makedirs(cluster_files_dir)

        # Check if its a condor job, and run that
        if (
            "condor"
            in c.pipeline_setup["system_config"]["on_grid"]["resource_manager"].lower()
        ):
            run_condor_jobs(c, config_file, subject_list_file, p_name)
        # All other schedulers are supported
        else:
            run_cpac_on_cluster(config_file, subject_list_file, cluster_files_dir)

    # Run on one computer
    else:
        # Create working dir
        if not os.path.exists(c.pipeline_setup["working_directory"]["path"]):
            try:
                os.makedirs(c.pipeline_setup["working_directory"]["path"])
            except:
                err = (
                    "\n\n[!] CPAC says: Could not create the working "
                    "directory: %s\n\nMake sure you have permissions "
                    "to write to this directory.\n\n"
                    % c.pipeline_setup["working_directory"]["path"]
                )
                raise Exception(err)
        """
        if not os.path.exists(c.pipeline_setup['log_directory']['path']):
            try:
                os.makedirs(c.pipeline_setup['log_directory']['path'])
            except:
                err = "\n\n[!] CPAC says: Could not create the log " \
                      "directory: %s\n\nMake sure you have permissions " \
                      "to write to this directory.\n\n" % c.pipeline_setup['log_directory']['path']
                raise Exception(err)
        """

        # BEGIN LONGITUDINAL TEMPLATE PIPELINE
        if c["longitudinal_template_generation", "run"]:
            run_T1w_longitudinal(sublist, c, dry_run=test_config)
            # TODO functional longitudinal pipeline

        # If it only allows one, run it linearly
        if c.pipeline_setup["system_config"]["num_participants_at_once"] == 1:
            for sub in sublist:
                try:
                    exitcode = run_workflow(
                        sub,
                        c,
                        True,
                        pipeline_timing_info,
                        p_name,
                        plugin,
                        plugin_args,
                        test_config,
                    )
                except Exception as exception:  # pylint: disable=broad-except
                    exitcode = 1
                    failed_to_start(set_subject(sub, c)[2], exception)
            return exitcode

        # Init job queue
        job_queue = []

        # Allocate processes
        processes = [
            Process(
                target=run_workflow,
                args=(
                    sub,
                    c,
                    True,
                    pipeline_timing_info,
                    p_name,
                    plugin,
                    plugin_args,
                    test_config,
                ),
            )
            for sub in sublist
        ]
        working_dir = os.path.join(
            c["pipeline_setup", "working_directory", "path"], p_name
        )
        # Create pipeline-specific working dir if not exists
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        # Set PID context to pipeline-specific file
        with open(os.path.join(working_dir, "pid.txt"), "w", encoding="utf-8") as pid:
            # If we're allocating more processes than are subjects, run
            # them all
            if (
                len(sublist)
                <= c.pipeline_setup["system_config"]["num_participants_at_once"]
            ):
                for i, _p in enumerate(processes):
                    try:
                        _p.start()
                        print(_p.pid, file=pid)
                    # pylint: disable=broad-except
                    except Exception as exception:
                        exitcode = 1
                        failed_to_start(set_subject(sublist[i], c)[2], exception)
            # Otherwise manage resources to run processes incrementally
            else:
                idx = 0
                while idx < len(sublist):
                    # If the job queue is empty and we haven't started indexing
                    if len(job_queue) == 0 and idx == 0:
                        # Init subject process index
                        idc = idx
                        # Launch processes (one for each subject)
                        for _p in processes[
                            idc : idc
                            + c.pipeline_setup["system_config"][
                                "num_participants_at_once"
                            ]
                        ]:
                            try:
                                _p.start()
                                print(_p.pid, file=pid)
                                job_queue.append(_p)
                                idx += 1
                            # pylint: disable=broad-except
                            except Exception as exception:
                                exitcode = 1
                                failed_to_start(
                                    set_subject(sublist[idx], c)[2], exception
                                )
                    # Otherwise, jobs are running - check them
                    else:
                        # Check every job in the queue's status
                        for job in job_queue:
                            # If the job is not alive
                            if not job.is_alive():
                                WFLOGGER.warning("found dead job %s", job)
                                # Find job and delete it from queue
                                loc = job_queue.index(job)
                                del job_queue[loc]
                                # ...and start the next available
                                # process (subject)
                                try:
                                    processes[idx].start()
                                    # Append this to job queue and
                                    # increment index
                                    # pylint: disable=modified-iterating-list
                                    job_queue.append(processes[idx])
                                    idx += 1
                                # pylint: disable=broad-except
                                except Exception as exception:
                                    exitcode = 1
                                    failed_to_start(
                                        set_subject(sublist[idx], c)[2], exception
                                    )
                        # Add sleep so while loop isn't consuming 100% of CPU
                        time.sleep(2)
            # set exitcode to 1 if any exception
            if hasattr(pid, "exitcode"):
                exitcode = exitcode or pid.exitcode
            # Close PID txt file to indicate finish
            pid.close()
    return exitcode
