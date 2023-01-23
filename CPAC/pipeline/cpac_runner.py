"""Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with C-PAC. If not, see <https://www.gnu.org/licenses/>."""
import os
import sys
import warnings
from multiprocessing import Process
from time import strftime
import yaml
from voluptuous.error import Invalid
from CPAC.utils.configuration import check_pname, Configuration, set_subject
from CPAC.utils.ga import track_run
from CPAC.utils.monitoring import failed_to_start, log_nodes_cb
from CPAC.longitudinal_pipeline.longitudinal_workflow import \
    anat_longitudinal_wf
from CPAC.utils.configuration.yaml_template import upgrade_pipeline_to_1_8


# Run condor jobs
def run_condor_jobs(c, config_file, subject_list_file, p_name):

    # Import packages
    import subprocess
    from time import strftime

    try:
        sublist = yaml.safe_load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        raise Exception("Subject list is not in proper YAML format. Please check your file")

    cluster_files_dir = os.path.join(os.getcwd(), 'cluster_files')
    subject_bash_file = os.path.join(cluster_files_dir, 'submit_%s.condor' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    f = open(subject_bash_file, 'w')

    print("Executable = /usr/bin/python", file=f)
    print("Universe = vanilla", file=f)
    print("transfer_executable = False", file=f)
    print("getenv = True", file=f)
    print("log = %s" % os.path.join(cluster_files_dir, 'c-pac_%s.log' % str(strftime("%Y_%m_%d_%H_%M_%S"))), file=f)

    for sidx in range(1,len(sublist)+1):
        print("error = %s" % os.path.join(cluster_files_dir, 'c-pac_%s.%s.err' % (str(strftime("%Y_%m_%d_%H_%M_%S")), str(sidx))), file=f)
        print("output = %s" % os.path.join(cluster_files_dir, 'c-pac_%s.%s.out' % (str(strftime("%Y_%m_%d_%H_%M_%S")), str(sidx))), file=f)

        print("arguments = \"-c 'import CPAC; CPAC.pipeline.cpac_pipeline.run( ''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'')\'\"" % (str(config_file), subject_list_file, str(sidx), c.maskSpecificationFile, c.roiSpecificationFile, c.templateSpecificationFile, p_name), file=f)
        print("queue", file=f)

    f.close()

    #commands.getoutput('chmod +x %s' % subject_bash_file )
    print(subprocess.getoutput("condor_submit %s " % (subject_bash_file)))


# Create and run script for CPAC to run on cluster
def run_cpac_on_cluster(config_file, subject_list_file,
                        cluster_files_dir):
    '''
    Function to build a SLURM batch job submission script and
    submit it to the scheduler via 'sbatch'
    '''

    # Import packages
    import subprocess
    import getpass
    import re
    from time import strftime
    from indi_schedulers import cluster_templates

    # Load in pipeline config
    try:
        pipeline_dict = yaml.safe_load(open(os.path.realpath(config_file), 'r'))
        pipeline_config = Configuration(pipeline_dict)
    except:
        raise Exception('Pipeline config is not in proper YAML format. '\
                        'Please check your file')
    # Load in the subject list
    try:
        sublist = yaml.safe_load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        raise Exception('Subject list is not in proper YAML format. '\
                        'Please check your file')

    # Init variables
    timestamp = str(strftime("%Y_%m_%d_%H_%M_%S"))
    job_scheduler = pipeline_config.pipeline_setup['system_config']['on_grid']['resource_manager'].lower()

    # For SLURM time limit constraints only, hh:mm:ss
    hrs_limit = 8 * len(sublist)
    time_limit = '%d:00:00' % hrs_limit

    # Batch file variables
    shell = subprocess.getoutput('echo $SHELL')
    user_account = getpass.getuser()
    num_subs = len(sublist)

    # Run CPAC via python -c command
    python_cpac_str = 'python -c "from CPAC.pipeline.cpac_pipeline import run; '\
                      'run(\'%(config_file)s\', \'%(subject_list_file)s\', '\
                      '%(env_arr_idx)s, \'%(pipeline_name)s\', '\
                      'plugin=\'MultiProc\', plugin_args=%(plugin_args)s)"'

    # Init plugin arguments
    plugin_args = {
        'n_procs': pipeline_config.pipeline_setup['system_config'][
            'max_cores_per_participant'],
        'memory_gb': pipeline_config.pipeline_setup['system_config'][
            'maximum_memory_per_participant'],
        'raise_insufficient': pipeline_config.pipeline_setup['system_config'][
            'raise_insufficient'],
        'status_callback': log_nodes_cb}

    # Set up run command dictionary
    run_cmd_dict = {'config_file': config_file,
                    'subject_list_file': subject_list_file,
                    'pipeline_name': pipeline_config.pipeline_setup[
                        'pipeline_name'],
                    'plugin_args': plugin_args}

    # Set up config dictionary
    config_dict = {'timestamp': timestamp,
                   'shell': shell,
                   'job_name': 'CPAC_' + pipeline_config.pipeline_setup[
                       'pipeline_name'],
                   'num_tasks': num_subs,
                   'queue': pipeline_config.pipeline_setup['system_config'][
                       'on_grid']['SGE']['queue'],
                   'par_env': pipeline_config.pipeline_setup['system_config'][
                       'on_grid']['SGE']['parallel_environment'],
                   'cores_per_task': pipeline_config.pipeline_setup[
                       'system_config']['max_cores_per_participant'],
                   'user': user_account,
                   'work_dir': cluster_files_dir,
                   'time_limit': time_limit}

    # Get string template for job scheduler
    if job_scheduler == 'pbs':
        env_arr_idx = '$PBS_ARRAYID'
        batch_file_contents = cluster_templates.pbs_template
        confirm_str = '(?<=Your job-array )\d+'
        exec_cmd = 'qsub'
    elif job_scheduler == 'sge':
        env_arr_idx = '$SGE_TASK_ID'
        batch_file_contents = cluster_templates.sge_template
        confirm_str = '(?<=Your job-array )\d+'
        exec_cmd = 'qsub'
    elif job_scheduler == 'slurm':
        env_arr_idx = '$SLURM_ARRAY_TASK_ID'
        batch_file_contents = cluster_templates.slurm_template
        confirm_str = '(?<=Submitted batch job )\d+'
        exec_cmd = 'sbatch'

    # Populate rest of dictionary
    config_dict['env_arr_idx'] = env_arr_idx
    run_cmd_dict['env_arr_idx'] = env_arr_idx
    config_dict['run_cmd'] = python_cpac_str % run_cmd_dict

    # Populate string from config dict values
    batch_file_contents = batch_file_contents % config_dict
    # Write file
    batch_filepath = os.path.join(cluster_files_dir, 'cpac_submit_%s.%s' \
                                  % (timestamp, job_scheduler))
    with open(batch_filepath, 'w') as f:
        f.write(batch_file_contents)

    # Get output response from job submission
    out = subprocess.getoutput('%s %s' % (exec_cmd, batch_filepath))

    # Check for successful qsub submission
    if re.search(confirm_str, out) == None:
        err_msg = 'Error submitting C-PAC pipeline run to %s queue' \
                  % job_scheduler
        raise Exception(err_msg)

    # Get pid and send to pid file
    pid = re.search(confirm_str, out).group(0)
    pid_file = os.path.join(cluster_files_dir, 'pid.txt')
    with open(pid_file, 'w') as f:
        f.write(pid)


def run_T1w_longitudinal(sublist, cfg):
    subject_id_dict = {}

    for sub in sublist:
        if sub['subject_id'] in subject_id_dict:
            subject_id_dict[sub['subject_id']].append(sub)
        else:
            subject_id_dict[sub['subject_id']] = [sub]

    # subject_id_dict has the subject_id as keys and a list of
    # sessions for each participant as value
    valid_longitudinal_data = False
    for subject_id, sub_list in subject_id_dict.items():
        if len(sub_list) > 1:
            valid_longitudinal_data = True
            anat_longitudinal_wf(subject_id, sub_list, cfg)
        elif len(sub_list) == 1:
            warnings.warn("\n\nThere is only one anatomical session "
                          "for sub-%s. Longitudinal preprocessing "
                          "will be skipped for this subject."
                          "\n\n" % subject_id)


# Run C-PAC subjects via job queue
def run(subject_list_file, config_file=None, p_name=None, plugin=None,
        plugin_args=None, tracking=True, num_subs_at_once=None, debug=False,
        test_config=False) -> int:
    """
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

    print('Run called with config file {0}'.format(config_file))

    if plugin_args is None:
        plugin_args = {'status_callback': log_nodes_cb}

    if not config_file:
        import pkg_resources as p
        config_file = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "pipeline_config_template.yml"))

    # Init variables
    sublist = None
    if '.yaml' in subject_list_file or '.yml' in subject_list_file:
        subject_list_file = os.path.realpath(subject_list_file)
    else:
        from CPAC.utils.bids_utils import collect_bids_files_configs, \
            bids_gen_cpac_sublist
        (file_paths, config) = collect_bids_files_configs(subject_list_file,
                                                          None)
        sublist = bids_gen_cpac_sublist(subject_list_file, file_paths,
                                        config, None)
        if not sublist:
            print(f"Did not find data in {subject_list_file}")
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
                c = Configuration(yaml.safe_load(open(config_file, 'r')))
            except Invalid:
                try:
                    upgrade_pipeline_to_1_8(config_file)
                    c = Configuration(yaml.safe_load(open(config_file, 'r')))
                except Exception as e:
                    print(
                        'C-PAC could not upgrade pipeline configuration file '
                        f'{config_file} to v1.8 syntax',
                        file=sys.stderr
                    )
                    raise e
            except Exception as e:
                raise e
    except IOError:
        print("config file %s doesn't exist" % config_file)
        raise
    except yaml.parser.ParserError as e:
        error_detail = "\"%s\" at line %d" % (
            e.problem,
            e.problem_mark.line
        )
        raise Exception(
            "Error parsing config file: {0}\n\n"
            "Error details:\n"
            "    {1}"
            "\n\n".format(config_file, error_detail)
        )
    except Exception as e:
        raise Exception(
            "Error parsing config file: {0}\n\n"
            "Error details:\n"
            "    {1}"
            "\n\n".format(config_file, e)
        )

    c.pipeline_setup['log_directory']['path'] = os.path.abspath(c.pipeline_setup['log_directory']['path'])
    c.pipeline_setup['working_directory']['path'] = os.path.abspath(c.pipeline_setup['working_directory']['path'])
    if 's3://' not in c.pipeline_setup['output_directory']['path']:
        c.pipeline_setup['output_directory']['path'] = os.path.abspath(c.pipeline_setup['output_directory']['path'])
    if debug:
        c.pipeline_setup['output_directory']['path']['write_debugging_outputs'] = "[1]"

    if num_subs_at_once:
        if not str(num_subs_at_once).isdigit():
            raise Exception('[!] Value entered for --num_cores not a digit.')
        c.pipeline_setup['system_config']['num_participants_at_once'] = int(num_subs_at_once)

    # Do some validation
    if not c.pipeline_setup['working_directory']['path']:
        raise Exception('Working directory not specified')

    if len(c.pipeline_setup['working_directory']['path']) > 70:
        warnings.warn("We recommend that the working directory full path "
                      "should have less then 70 characters. "
                      "Long paths might not work in your operating system.")
        warnings.warn("Current working directory: "
                      f"{c.pipeline_setup['working_directory']['path']}")

    # Get the pipeline name
    p_name = check_pname(p_name, c)

    # Load in subject list
    try:
        if not sublist:
            sublist = yaml.safe_load(open(subject_list_file, 'r'))
    except:
        print("Subject list is not in proper YAML format. Please check " \
              "your file")
        raise Exception

    # Populate subject scan map
    sub_scan_map = {}
    try:
        for sub in sublist:
            if sub['unique_id']:
                s = sub['subject_id'] + "_" + sub["unique_id"]
            else:
                s = sub['subject_id']
            scan_ids = ['scan_anat']

            if 'func' in sub:
                for id in sub['func']:
                    scan_ids.append('scan_'+ str(id))

            if 'rest' in sub:
                for id in sub['rest']:
                    scan_ids.append('scan_'+ str(id))

            sub_scan_map[s] = scan_ids
    except:
        print("\n\n" + "ERROR: Subject list file not in proper format - " \
              "check if you loaded the correct file?" + "\n" + \
              "Error name: cpac_runner_0001" + "\n\n")
        raise Exception

    pipeline_timing_info = []
    pipeline_timing_info.append(unique_pipeline_id)
    pipeline_timing_info.append(pipeline_start_stamp)
    pipeline_timing_info.append(len(sublist))

    if tracking:
        try:
            track_run(
                level='participant' if not test_config else 'test',
                participants=len(sublist)
            )
        except:
            print("Usage tracking failed for this run.")

    # If we're running on cluster, execute job scheduler
    if c.pipeline_setup['system_config']['on_grid']['run']:

        # Create cluster log dir
        cluster_files_dir = os.path.join(c.pipeline_setup['log_directory']['path'], 'cluster_files')
        if not os.path.exists(cluster_files_dir):
            os.makedirs(cluster_files_dir)

        # Check if its a condor job, and run that
        if 'condor' in c.pipeline_setup['system_config']['on_grid']['resource_manager'].lower():
            run_condor_jobs(c, config_file, subject_list_file, p_name)
        # All other schedulers are supported
        else:
            run_cpac_on_cluster(config_file, subject_list_file, cluster_files_dir)

    # Run on one computer
    else:
        # Create working dir
        if not os.path.exists(c.pipeline_setup['working_directory']['path']):
            try:
                os.makedirs(c.pipeline_setup['working_directory']['path'])
            except:
                err = "\n\n[!] CPAC says: Could not create the working " \
                      "directory: %s\n\nMake sure you have permissions " \
                      "to write to this directory.\n\n" % c.pipeline_setup['working_directory']['path']
                raise Exception(err)
        '''
        if not os.path.exists(c.pipeline_setup['log_directory']['path']):
            try:
                os.makedirs(c.pipeline_setup['log_directory']['path'])
            except:
                err = "\n\n[!] CPAC says: Could not create the log " \
                      "directory: %s\n\nMake sure you have permissions " \
                      "to write to this directory.\n\n" % c.pipeline_setup['log_directory']['path']
                raise Exception(err)
        '''

        # BEGIN LONGITUDINAL TEMPLATE PIPELINE
        if hasattr(c, 'longitudinal_template_generation') and \
                        c.longitudinal_template_generation['run']:

            run_T1w_longitudinal(sublist, c)
            # TODO functional longitudinal pipeline

        '''
            if valid_longitudinal_data:
                rsc_file_list = []
                for dirpath, dirnames, filenames in os.walk(c.pipeline_setup[
                                                                'output_directory']['path']):
                    for f in filenames:
                        # TODO is there a better way to check output folder name?
                        if f != '.DS_Store' and 'T1w_longitudinal_pipeline' in dirpath:
                            rsc_file_list.append(os.path.join(dirpath, f))

                subject_specific_dict = {subj: [] for subj in subject_id_dict.keys()}
                session_specific_dict = {os.path.join(session['subject_id'], session['unique_id']): [] for session in sublist}
                for rsc_path in rsc_file_list:
                    key = [s for s in session_specific_dict.keys() if s in rsc_path]
                    if key:
                        session_specific_dict[key[0]].append(rsc_path)
                    else:
                        subj = [s for s in subject_specific_dict.keys() if s in rsc_path]
                        if subj:
                            subject_specific_dict[subj[0]].append(rsc_path)
                
                # update individual-specific outputs: 
                # anatomical_brain, anatomical_brain_mask and anatomical_reorient
                for key in session_specific_dict.keys():
                    for f in session_specific_dict[key]:
                        sub, ses = key.split('/')
                        ses_list = [subj for subj in sublist if sub in subj['subject_id'] and ses in subj['unique_id']]
                        if len(ses_list) > 1:
                            raise Exception("There are several files containing " + f)
                        if len(ses_list) == 1:
                            ses = ses_list[0]
                            subj_id = ses['subject_id']
                            tmp = f.split(c.pipeline_setup['output_directory']['path'])[-1]
                            keys = tmp.split(os.sep)
                            if keys[0] == '':
                                keys = keys[1:]
                            if len(keys) > 1:
                                if ses.get('resource_pool') is None:
                                    ses['resource_pool'] = {
                                        keys[0].split(c.pipeline_setup['pipeline_name'] + '_')[-1]: {
                                            keys[-2]: f
                                        }
                                    }
                                else:
                                    strat_key = keys[0].split(c.pipeline_setup['pipeline_name'] + '_')[-1]
                                    if ses['resource_pool'].get(strat_key) is None:
                                        ses['resource_pool'].update({
                                            strat_key: {
                                                keys[-2]: f
                                            }
                                        })
                                    else:
                                        ses['resource_pool'][strat_key].update({
                                                keys[-2]: f
                                            })

                for key in subject_specific_dict:
                    for f in subject_specific_dict[key]:
                        ses_list = [subj for subj in sublist if key in subj['anat']]
                        for ses in ses_list:
                            tmp = f.split(c.pipeline_setup['output_directory']['path'])[-1]
                            keys = tmp.split(os.sep)
                            if keys[0] == '':
                                keys = keys[1:]
                            if len(keys) > 1:
                                if ses.get('resource_pool') is None:
                                    ses['resource_pool'] = {
                                        keys[0].split(c.pipeline_setup['pipeline_name'] + '_')[-1]: {
                                            keys[-2]: f
                                        }
                                    }
                                else:
                                    strat_key = keys[0].split(c.pipeline_setup['pipeline_name'] + '_')[-1]
                                    if ses['resource_pool'].get(strat_key) is None:
                                        ses['resource_pool'].update({
                                            strat_key: {
                                                keys[-2]: f
                                            }
                                        })
                                    else:
                                        if keys[-2] == 'anatomical_brain' or keys[-2] == 'anatomical_brain_mask' or keys[-2] == 'anatomical_skull_leaf':
                                            pass
                                        elif 'apply_warp_anat_longitudinal_to_standard' in keys[-2] or 'fsl_apply_xfm_longitudinal' in keys[-2]:
                                            # TODO update!!!
                                            # it assumes session id == last key (ordered by session count instead of session id) + 1
                                            # might cause problem if session id is not continuous
                                            def replace_index(target1, target2, file_path):
                                                index1 = file_path.index(target1)+len(target1)
                                                index2 = file_path.index(target2)+len(target2)
                                                file_str_list = list(file_path)
                                                file_str_list[index1] = "*"
                                                file_str_list[index2] = "*"
                                                file_path_updated = "".join(file_str_list)
                                                file_list = glob.glob(file_path_updated)
                                                file_list.sort()
                                                return file_list
                                            if ses['unique_id'] == str(int(keys[-2][-1])+1):
                                                if keys[-3] == 'seg_probability_maps':
                                                    f_list = replace_index('seg_probability_maps_', 'segment_prob_', f)
                                                    ses['resource_pool'][strat_key].update({
                                                        keys[-3]: f_list
                                                    })
                                                elif keys[-3] == 'seg_partial_volume_files':
                                                    f_list = replace_index('seg_partial_volume_files_', 'segment_pve_', f)
                                                    ses['resource_pool'][strat_key].update({
                                                        keys[-3]: f_list
                                                    })
                                                else:
                                                    ses['resource_pool'][strat_key].update({
                                                        keys[-3]: f # keys[-3]: 'anatomical_to_standard'
                                                    })
                                        elif keys[-2] != 'warp_list':
                                            ses['resource_pool'][strat_key].update({
                                                    keys[-2]: f
                                                })
                                        elif keys[-2] == 'warp_list':
                                            if 'ses-'+ses['unique_id'] in tmp:
                                                ses['resource_pool'][strat_key].update({
                                                    keys[-2]: f
                                                })
                for key in subject_specific_dict:
                    ses_list = [subj for subj in sublist if key in subj['anat']]
                    for ses in ses_list:
                        for reg_strat in strat_list:
                            try: 
                                ss_strat_list = list(ses['resource_pool'])
                                for strat_key in ss_strat_list:
                                    try:
                                        ses['resource_pool'][strat_key].update({
                                            'registration_method': reg_strat['registration_method']
                                        })
                                    except KeyError:
                                        pass
                            except KeyError:
                                pass

                yaml.dump(sublist, open(os.path.join(c.pipeline_setup['working_directory']['path'],'data_config_longitudinal.yml'), 'w'), default_flow_style=False)
                print('\n\n' + 'Longitudinal pipeline completed.' + '\n\n')

                # skip main preprocessing
                if (
                    not c.anatomical_preproc['run'] and
                    not c.functional_preproc['run']
                ):
                    sys.exit()
        '''
        # END LONGITUDINAL TEMPLATE PIPELINE

        # If it only allows one, run it linearly
        if c.pipeline_setup['system_config']['num_participants_at_once'] == 1:
            for sub in sublist:
                try:
                    run_workflow(sub, c, True, pipeline_timing_info,
                                 p_name, plugin, plugin_args, test_config)
                except Exception as exception:  # pylint: disable=broad-except
                    exitcode = 1
                    failed_to_start(set_subject(sub, c)[2], exception)
            return exitcode

        # Init job queue
        job_queue = []

        # Allocate processes
        processes = [Process(target=run_workflow,
                             args=(sub, c, True, pipeline_timing_info, p_name,
                                   plugin, plugin_args, test_config)) for
                     sub in sublist]
        working_dir = os.path.join(c['pipeline_setup', 'working_directory',
                                     'path'], p_name)
        # Create pipeline-specific working dir if not exists
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        # Set PID context to pipeline-specific file
        with open(os.path.join(working_dir, 'pid.txt'), 'w', encoding='utf-8'
                  ) as pid:
            # If we're allocating more processes than are subjects, run
            # them all
            if len(sublist) <= c.pipeline_setup['system_config'][
                    'num_participants_at_once']:
                for i, _p in enumerate(processes):
                    try:
                        _p.start()
                        print(_p.pid, file=pid)
                    # pylint: disable=broad-except
                    except Exception as exception:
                        exitcode = 1
                        failed_to_start(set_subject(sublist[i], c)[2],
                                        exception)
            # Otherwise manage resources to run processes incrementally
            else:
                idx = 0
                while idx < len(sublist):
                    # If the job queue is empty and we haven't started indexing
                    if len(job_queue) == 0 and idx == 0:
                        # Init subject process index
                        idc = idx
                        # Launch processes (one for each subject)
                        for _p in processes[idc: idc + c.pipeline_setup[
                                'system_config']['num_participants_at_once']]:
                            try:
                                _p.start()
                                print(_p.pid, file=pid)
                                job_queue.append(_p)
                                idx += 1
                            # pylint: disable=broad-except
                            except Exception as exception:
                                exitcode = 1
                                failed_to_start(set_subject(sublist[idx],
                                                            c)[2], exception)
                    # Otherwise, jobs are running - check them
                    else:
                        # Check every job in the queue's status
                        for job in job_queue:
                            # If the job is not alive
                            if not job.is_alive():
                                # Find job and delete it from queue
                                print('found dead job ', job)
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
                                    failed_to_start(set_subject(sublist[idx],
                                                                c)[2],
                                                    exception)
                        # Add sleep so while loop isn't consuming 100% of CPU
                        time.sleep(2)
            # set exitcode to 1 if any exception
            if hasattr(pid, 'exitcode'):
                exitcode = exitcode or pid.exitcode
            # Close PID txt file to indicate finish
            pid.close()
    return exitcode
