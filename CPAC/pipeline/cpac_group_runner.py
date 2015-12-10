import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
from time import strftime

from multiprocessing import Process

import re
import os
import sys
import glob
import time
import csv

from nipype import logging

from CPAC.utils import Configuration


def split_folders(path):
    folders = []
    
    while 1:
        path, folder = os.path.split(path)
        if folder != "":
            folders.append(folder)
        else:
            if path != "":
                folders.append(path)
            break

    folders.reverse()
    #print folders
    return folders



def run_sge_jobs(c, config_file, resource, subject_infos):


    import commands
    import pickle
    from time import strftime

    temp_files_dir = os.path.join(os.getcwd(), 'cluster_temp_files')
    resource_file = os.path.join(temp_files_dir, 'resource.obj')
    f = open(resource_file, 'w')
    pickle.dump(resource, f)
    f.close()

    subject_infos_file = os.path.join(temp_files_dir, 'subject_infos.obj')
    f = open(subject_infos_file, 'w')
    pickle.dump(subject_infos, f)
    f.close()

    shell = commands.getoutput('echo $SHELL')

    subject_bash_file = ''
    if c.runBASC:
        subject_bash_file = os.path.join(temp_files_dir, 'submit_BASC_%s.sge' % str(strftime("%Y_%m_%d_%H_%M_%S")))

    if c.runCWAS:
        subject_bash_file = os.path.join(temp_files_dir, 'submit_CWAS_%s.sge' % str(strftime("%Y_%m_%d_%H_%M_%S")))

    if c.runGroupAnalysis:
        subject_bash_file = os.path.join(temp_files_dir, 'submit_GroupAnalysis_%s.sge' % str(strftime("%Y_%m_%d_%H_%M_%S")))

    f = open(subject_bash_file, 'w')
    print >>f, '#! %s' % shell
    print >>f, '#$ -cwd'
    print >>f, '#$ -S %s' % shell
    print >>f, '#$ -V'
    print >>f, '#$ -q %s' % c.queue
    print >>f, '#$ -pe %s %d' % (c.parallelEnvironment, c.numCoresPerSubject)
    print >>f, '#$ -e %s' % os.path.join(temp_files_dir, 'c-pac_%s.err' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    print >>f, '#$ -o %s' % os.path.join(temp_files_dir, 'c-pac_%s.out' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    print >>f, 'source ~/.bashrc'

    if c.runBASC:

        print >>f, "python -c \"import CPAC; CPAC.pipeline.cpac_basc_pipeline.run(\\\"%s\\\" , \\\"%s\\\") \" " % (str(config_file), subject_infos_file)

    elif c.runCWAS:

        print >>f, "python -c \"import CPAC; CPAC.pipeline.cpac_cwas_pipeline.run(\\\"%s\\\" , \\\"%s\\\") \" " % (str(config_file), subject_infos_file)

    elif c.runGroupAnalysis:

        print >>f, "python -c \"import CPAC; CPAC.pipeline.cpac_group_analysis_pipeline.run(\\\"%s\\\" , \\\"%s\\\", \\\"%s\\\") \" " % (str(config_file), subject_infos_file, resource_file)


    f.close()

    commands.getoutput('chmod +x %s' % subject_bash_file )
#    print commands.getoutput('qsub  %s ' % (subject_bash_file))
    print 'qsub  %s ' % (subject_bash_file)



def run_pbs_jobs(c, config_file, resource, subject_infos):

    import commands
    import pickle
    from time import strftime

    temp_files_dir = os.path.join(os.getcwd(), 'cluster_temp_files')
    resource_file = os.path.join(temp_files_dir, 'resource.obj')
    f = open(resource_file, 'w')
    pickle.dump(resource, f)
    f.close()

    subject_infos_file = os.path.join(temp_files_dir, 'subject_infos.obj')
    f = open(subject_infos_file, 'w')
    pickle.dump(subject_infos, f)
    f.close()

    subject_bash_file = ''
    shell = commands.getoutput('echo $SHELL')
    if c.runBASC:
        subject_bash_file = os.path.join(temp_files_dir, 'submit_BASC_%s.sge' % str(strftime("%Y_%m_%d_%H_%M_%S")))

    if c.runCWAS:
        subject_bash_file = os.path.join(temp_files_dir, 'submit_CWAS_%s.sge' % str(strftime("%Y_%m_%d_%H_%M_%S")))

    if c.runGroupAnalysis:
        subject_bash_file = os.path.join(temp_files_dir, 'submit_GroupAnalysis_%s.sge' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    f = open(subject_bash_file, 'w')
    print >>f, '#! %s' % shell
    print >>f, '#PBS -S %s' % shell
    print >>f, '#PBS -V'
    print >>f, '#PBS -q %s' % c.queue
    print >>f, '#PBS -l nodes=1:ppn=%d' % c.numCoresPerSubject
    print >>f, '#PBS -e %s' % os.path.join(temp_files_dir, 'c-pac_%s.err' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    print >>f, '#PBS -o %s' % os.path.join(temp_files_dir, 'c-pac_%s.out' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    print >>f, 'source ~/.bashrc'

    if c.runBASC:

        print >>f, "python -c \"import CPAC; CPAC.pipeline.cpac_basc_pipeline.run(\\\"%s\\\" , \\\"%s\\\") \" " % (str(config_file), subject_infos_file)

    elif c.runCWAS:

        print >>f, "python -c \"import CPAC; CPAC.pipeline.cpac_cwas_pipeline.run(\\\"%s\\\" , \\\"%s\\\") \" " % (str(config_file), subject_infos_file)

    elif c.runGroupAnalysis:

        print >>f, "python -c \"import CPAC; CPAC.pipeline.cpac_group_analysis_pipeline.run(\\\"%s\\\" , \\\"%s\\\", \\\"%s\\\") \" " % (str(config_file), subject_infos_file, resource_file)


    f.close()

    commands.getoutput('chmod +x %s' % subject_bash_file )
    print commands.getoutput('qsub  %s ' % (subject_bash_file))



def load_group_config(group_config_file):

    import yaml

    try:
        ga_config_path = os.path.realpath(group_config_file)

        with open(ga_config_path,"r") as f:
            ga_config = Configuration(yaml.load(f))
    except:
        raise Exception("\n\nError in reading %s configuration file\n\n" \
                        % group_config_file)

    if len(ga_config.derivative_list) == 0:
        print '[!] CPAC says: You do not have any derivatives selected ' \
                'to run for group-level analysis. Return to your group-' \
                'analysis configuration file and select at least one.'
        print 'Group analysis configuration file: %s\n\n' \
                % group_config_file
        raise Exception

    return ga_config



def load_group_subject_list(ga_config):

    import pandas as pd

    # get the group subject list
    if not ga_config.subject_list.endswith(".csv"):
        err = "\n\n[!] CPAC says: The group-level analysis participant " \
                "list should be a CSV file (.csv).\n\n"
        raise Exception(err)

    with open(ga_config.subject_list,"r") as f:
        ga_sublist = pd.read_csv(f)


    if "participant" not in ga_sublist.columns:
        err = "\n\n[!] CPAC says: Your group-level analysis subject "\
                "list CSV is missing a 'participant' column.\n\n"
        raise Exception(err)

    if ga_config.repeated_measures == True:
        if "session" not in ga_sublist.columns and \
            "series" not in ga_sublist.columns:
            err = "\n\n[!] CPAC says: You have selected to run " \
                    "repeated measures in your group-level analysis, " \
                    "but you do not have either a 'session' or " \
                    "'series' column in your group analysis subject " \
                    "list CSV.\n\n"
            raise Exception(err)
    else:
        if "session" in ga_sublist.columns or \
            "series" in ga_sublist.columns:
            err = "\n\n[!] CPAC says: Your group-level analysis " \
                    "subject list CSV is formatted for repeated " \
                    "measures, but you do not have repeated measures " \
                    "selected in your group model configuration.\n\n"
            raise Exception(err)


    return ga_sublist



def wildcards_into_filepath(ga_sublist, subject_path, matched_subs, output):

    # 1. take in the "subject_path", which is the filepath of the current
    #    output file being iterated over in run()
    # 2. replace the participant IDs and potentially the session and series
    #    IDs appropriately with wildcards (*)

    ''' WHEN THIS LOOP RUNS THROUGH WITHOUT EVER HITTING SOMETHING IN THE '''
    ''' SUBJECT LIST, THAT MEANS IT WAS GIVEN A SUBJECT PATH FOR A SUBJECT '''
    ''' THAT WASN'T IN THE SUBLIST- HOWEVER, WE NEED SOMETHING THAT WILL '''
    ''' TAKE THE SUBS IN THE SUBLIST THAT DID NOT FIND A DERIVATIVE; AND '''
    ''' TO DO THIS FOR EACH DERIVATIVE '''

    session_id = None
    series_id = None

    if ("session" in ga_sublist.columns) and \
        ("series" in ga_sublist.columns):

        for subject, session, series in zip(ga_sublist.participant, \
            ga_sublist.session, ga_sublist.series):

            subject = str(subject)
            session = str(session)
            series = str(series)

            if (subject in subject_path) and \
                (session in subject_path) and \
                (series in subject_path):

                key = subject_path.replace(subject, "*")
                key = key.replace(session, "*")
                key = key.replace(series, "*")
                subject_id = subject
                session_id = session
                series_id = series

                if output not in matched_subs.keys():
                    matched_subs[output] = []

                matched_subs[output].append((subject, session, series))
                
                break

        else:
            # the subject path is for an output file of a subject that
            # isn't in the group subject list
            return None, None, None, None, matched_subs

    elif ("session" in ga_sublist.columns) and \
        ("series" not in ga_sublist.columns):

        for subject, session in zip(ga_sublist.participant, \
            ga_sublist.session):

            subject = str(subject)
            session = str(session)

            if (subject in subject_path) and \
                (session in subject_path):

                key = subject_path.replace(subject, "*")
                key = key.replace(session, "*")
                subject_id = subject
                session_id = session

                if output not in matched_subs.keys():
                    matched_subs[output] = []

                matched_subs[output].append((str(subject), session))

                break

        else:
            return None, None, None, None, matched_subs

    elif ("series" in ga_sublist.columns) and \
        ("session" not in ga_sublist.columns):

        for subject, series in zip(ga_sublist.participant, \
            ga_sublist.series):

            subject = str(subject)
            series = str(series)

            if (subject in subject_path) and \
                (series in subject_path):

                key = subject_path.replace(subject, "*")
                key = key.replace(series, "*")
                subject_id = subject
                series_id = series

                if output not in matched_subs.keys():
                    matched_subs[output] = []

                if (subject, series) not in matched_subs[output]:
                    matched_subs[output].append((str(subject), series))
                else:
                    # this fires if there are multiple sessions, but the user
                    # did not denote multiple sessions under a 'sessions'
                    # column in the group analysis subject list
                    err = "\n\n[!] CPAC says: It appears there are multiple "\
                          "sessions in your data, but you do not have any " \
                          "sessions denoted in your group analysis subject " \
                          "list.\n\n"
                    raise Exception(err)

                break

        else:
            return None, None, None, None, matched_subs

    elif "participant" in ga_sublist.columns:
        # each group of subjects from each session go into their own
        # separate model, instead of combining all sessions into one
        #     i.e. session_1 will have its own set of group
        #          analysis outputs, session_2 will have another..
        #          they will be separated by session just like
        #          they would be separated by scan
                    
        for subject in ga_sublist.participant:
            if str(subject) in subject_path:

                key = subject_path.replace(str(subject), '*')
                subject_id = str(subject)

                if output not in matched_subs.keys():
                    matched_subs[output] = []

                matched_subs[output].append((str(subject)))
                
                break
        
        else:     
            # the subject path is for an output file of a subject that
            # isn't in the group subject list     
            return None, None, None, None, matched_subs

    else:

        err = "\n\n[!] CPAC says: Your group-level analysis subject list " \
              "is missing a column labeled 'Participant'.\n\n"
        raise Exception(err)


    return key, subject_id, session_id, series_id, matched_subs



def run(config_file, output_path_file):
    
    # Runs group analysis

    import yaml

    # Load the config file into 'c'
    with open(os.path.realpath(config_file),"r") as f:
        c = Configuration(yaml.load(f))

    if len(c.modelConfigs) == 0:
        print '[!] CPAC says: You do not have any models selected ' \
              'to run for group-level analysis. Return to your pipeline ' \
              'configuration file and create or select at least one.\n\n'
        raise Exception

    subject_paths = []
    ind_outputs = []
    ga_configs = {}
    ga_sublists = {}
    matched_subs = {}

    # load the group model configs and the group participant lists into
    # dictionaries, with the model names as keys
    for group_config_file in c.modelConfigs:

        ga_config = load_group_config(group_config_file)

        # gather which outputs the user has selected to run group analysis for
        for output_type in ga_config.derivative_list:
            if output_type not in ind_outputs:
                ind_outputs.append(output_type)

        if ga_config.model_name in ga_configs.keys():
            err = "\n\n[!] CPAC says: You have more than one group-level " \
                  "analysis model configuration with the same name!\n\n" \
                  "Duplicate model name: %s\n\n" % ga_config.model_name
            raise Exception(err)

        ga_configs[ga_config.model_name] = ga_config

    for ga_config in ga_configs.values():

        ga_sublist = load_group_subject_list(ga_config)

        ga_sublists[ga_config.model_name] = ga_sublist

    
    # collect all of the output paths from individual-level analysis   
    for root, folders, files in os.walk(output_path_file):
    
        split_output_dir_path = output_path_file.split("/")
    
        for filename in files:
        
            if filename.endswith("nii.gz"):
    
                fullpath = os.path.join(root, filename)
            
                split_fullpath = fullpath.split("/")
                
                deriv_folder_name = \
                    split_fullpath[len(split_output_dir_path)+1]
            
                for output_name in ind_outputs:
            
                    if output_name == deriv_folder_name:
        
                        subject_paths.append(fullpath)


    if len(subject_paths) == 0:
        print '[!] CPAC says: No individual-level analysis outputs were ' \
              'found given the path file you provided.\n\nPipeline Output ' \
              'Directory provided: ', output_path_file, '\n\nEither make ' \
              'sure your Output Directory path is correct, or that ' \
              'individual-level analysis completed successfully.\n\n'
        raise Exception


    # 'subject_paths' is a list of every output file from every subject in
    # the output folder of the run, for the outputs selected for grp analysis

    # converts the subject_paths list into a set to enforce no duplicates
    set_subject_paths = set(subject_paths)

    # converts the set back into a list
    subject_paths = list(set_subject_paths)


    base_path = c.outputDirectory


    from collections import defaultdict
    analysis_map_gp = defaultdict(list)


    print "Parsing through output paths. This may take a little while " \
          "depending on how many subjects, group analysis models, or " \
          "selected derivatives you have..\n"

    count = 0

    for subject_path in subject_paths:

        # get classifying information for the current subject_path (which is
        # an output file path from individual analysis)

        split_output_dir_path = output_path_file.split("/")
        split_fullpath = subject_path.split("/")

        pipeline_id = split_fullpath[len(split_output_dir_path)-1]
        subject_unique_id = split_fullpath[len(split_output_dir_path)]
        resource_id = split_fullpath[len(split_output_dir_path)+1]
        scan_id = split_fullpath[len(split_output_dir_path)+2]

        
        # add auxiliary stuff to resource_id if applicable
        
        if ("_mask_" in subject_path) and (("sca_roi" in subject_path) or \
            ("sca_tempreg" in subject_path)):
            
            for dirname in split_fullpath:
                if "_mask_" in dirname:
                    maskname = dirname
                    
            filename = split_fullpath[-1]
            
            if ".nii.gz" in filename:
                filename = filename.replace(".nii.gz","")
            elif ".nii" in filename:
                filename = filename.replace(".nii","")
            
            resource_name = resource_id + "_%s_%s" % (maskname, filename)

            
        elif ("_spatial_map_" in subject_path) and \
            ("dr_tempreg" in subject_path):
            
            for dirname in split_fullpath:
                if "_spatial_map_" in dirname:
                    mapname = dirname
                    
            filename = split_fullpath[-1]
            
            if ".nii.gz" in filename:
                filename = filename.replace(".nii.gz","")
            elif ".nii" in filename:
                filename = filename.replace(".nii","")
            
            resource_name = resource_id + "_%s_%s" % (mapname, filename)
            
            
        elif ("_mask_" in subject_path) and ("centrality" in subject_path):
            
            for dirname in split_fullpath:
                if "_mask_" in dirname:
                    maskname = dirname
                    
            filename = split_fullpath[-1]
            
            if ".nii.gz" in filename:
                filename = filename.replace(".nii.gz","")
            elif ".nii" in filename:
                filename = filename.replace(".nii","")
            
            resource_name = resource_id + "_%s_%s" % (maskname, filename)
            
            
        else:
        
            resource_name = resource_id


        for model in ga_configs:

            if resource_id in list(ga_configs[model].derivative_list):

                if model not in matched_subs.keys():
                    matched_subs[model] = {}

                key, subject_id, session_id, series_id, matched_subs[model] =\
                    wildcards_into_filepath(ga_sublists[model], subject_path,\
                        matched_subs[model], resource_id)

                # toss it in the pile!
                if key != None:
                    analysis_map_gp[(resource_name, group_config_file, model, key)].append((pipeline_id, subject_id, session_id, series_id, scan_id, subject_path))



        count += 1

        if count == int(len(subject_paths)*0.7):
            print "Almost finished parsing output paths.."

        # with this loop, 'analysis_map_gp' is a dictionary with a key for
        # each individual output file - and each entry is a list of tuples,
        # one tuple for each subject in the subject list, containing
        # 'subject_path', which is a full path to that output file for that
        # one particular subject


    # let's find missing subjects
    sublist_tuples = []
    missing_subs = {}

    # create a list of tuples based on the group analysis participant list
    for model in ga_sublists.keys():
        for sub_row in ga_sublists[model].values:
            sub_row = [str(i) for i in sub_row]
            if len(sub_row) > 1:
                sublist_tuples.append(tuple(sub_row))
            else:
                sublist_tuples.append(sub_row[0])


    sublist_tuples = list(set(sublist_tuples))


    # compare what's been matched to what's in the participant list
    for group_model in matched_subs.keys():
        for output_type in matched_subs[group_model].keys():
            for sub_tuple in sublist_tuples:
                if sub_tuple not in matched_subs[group_model][output_type]:
            
                    if group_model not in missing_subs.keys():
                        missing_subs[group_model] = {}

                    if output_type not in missing_subs[group_model].keys():
                        missing_subs[group_model][output_type] = []

                    missing_subs[group_model][output_type].append(sub_tuple)

    if len(missing_subs) > 0:

        miss_msg = "\n\n[!] CPAC warns: Output files missing for " \
                   "participants included in group-level analysis:\n"

        for group_model in missing_subs.keys():
            for output_type in missing_subs[group_model].keys():
                miss_msg = miss_msg + "\n" + output_type + "\n"
                for sub_tuple in missing_subs[group_model][output_type]:
                    miss_msg = miss_msg + str(sub_tuple) + "\n"

        print miss_msg


    # let's create new subject lists
    new_sublists = {}

    for group_model in matched_subs.keys():

        if group_model not in new_sublists.keys():
            new_sublists[group_model] = {}

        for output_type in matched_subs[group_model].keys():

            if output_type not in new_sublists[group_model].keys():
                new_sublists[group_model][output_type] = []

            new_list = []

            tuple_list = matched_subs[group_model][output_type]

            new_list.append(ga_sublists[group_model].columns)

            tuple_list = list(set(tuple_list))

            for sub_tuple in tuple_list:

                if type(sub_tuple) is tuple:
                    new_list.append(list(sub_tuple))
                elif type(sub_tuple) is str:
                    new_list.append(sub_tuple)

            new_sublists[group_model][output_type] = new_list


    # check in case there are no output files

    if len(analysis_map_gp) == 0:
        err = "\n\n[!] CPAC says: No output files from individual-level " \
              "analysis were found for the subjects or sessions/series " \
              "in the group analysis subject list(s)!\n\n"
        raise Exception(err)



    print "Finished parsing through output paths!\n"


    '''
    for resource, group_model, glob_key in analysis_map_gp.keys():
        if resource == 'functional_mni':


            if 1 in c.runBASC:

                if not c.runOnGrid:
                    from CPAC.pipeline.cpac_basc_pipeline import prep_basc_workflow
                    prep_basc_workflow(c, analysis_map_gp[(resource, group_model, glob_key)])
                else:
                    if 'sge' in c.resourceManager.lower():
                        run_sge_jobs(c, config_file, resource, analysis_map_gp[(resource, group_model, glob_key)])

                    elif 'pbs' in c.resourceManager.lower():
                        run_pbs_jobs(c, config_file, resource, analysis_map_gp[(resource, group_model, glob_key)])


            if 1 in c.runCWAS:

                if not c.runOnGrid:

                    from CPAC.pipeline.cpac_cwas_pipeline import prep_cwas_workflow
                    prep_cwas_workflow(c, analysis_map[(resource, group_model, glob_key)])

                else:
                    if 'sge' in c.resourceManager.lower():
                        run_sge_jobs(c, config_file, resource, analysis_map[(resource, group_model, glob_key)])

                    elif 'pbs' in c.resourceManager.lower():
                        run_pbs_jobs(c, config_file, resource, analysis_map[(resource, group_model, glob_key)])
    '''

    procss = []
    

    for resource, group_model_config, group_model, glob_key in analysis_map_gp.keys():

        # 'resource' is each type of output
        # 'glob_key' is a path to each and every individual output file,
        # except with the subject ID replaced with a wildcard (*)
                      
        #get all the motion parameters across subjects
        from CPAC.utils import extract_parameters
        scrub_threshold = extract_parameters.run(c.outputDirectory, c.runScrubbing)

        if not c.runOnGrid:

            from CPAC.pipeline.cpac_ga_model_generator import prep_group_analysis_workflow
            procss.append(Process(target=prep_group_analysis_workflow, args=(c, group_model_config, resource, new_sublists[group_model][resource], analysis_map_gp[(resource, group_model_config, group_model, glob_key)], scrub_threshold)))
            
        else:
        
            print "\n\n[!] CPAC says: Group-level analysis has not yet " \
                  "been implemented to handle runs on a cluster or grid.\n\n"\
                  "Please turn off 'Run CPAC On A Cluster/Grid' in order " \
                  "to continue with group-level analysis. This will submit " \
                  "the job to only one node, however.\n\nWe will update " \
                  "users on when this feature will be available through " \
                  "release note announcements.\n\n"

    
          
    pid = open(os.path.join(c.outputDirectory, 'pid_group.txt'), 'w')
                        
    jobQueue = []
    if len(procss) <= c.numGPAModelsAtOnce:
        """
        Stream all the subjects as sublist is
        less than or equal to the number of 
        subjects that need to run
        """
        for p in procss:
            p.start()
            print >>pid,p.pid
                
    else:
        """
        Stream the subject workflows for preprocessing.
        At Any time in the pipeline c.numSubjectsAtOnce
        will run, unless the number remaining is less than
        the value of the parameter stated above
        """
        idx = 0
        while(idx < len(procss)):
                
            if len(jobQueue) == 0 and idx == 0:
                
                idc = idx
                    
                for p in procss[idc: idc + c.numGPAModelsAtOnce]:
                
                    p.start()
                    print >>pid,p.pid
                    jobQueue.append(p)
                    idx += 1
                
            else:
                
                for job in jobQueue:
                
                    if not job.is_alive():
                        print 'found dead job ', job
                        loc = jobQueue.index(job)
                        del jobQueue[loc]
                        procss[idx].start()
                
                        jobQueue.append(procss[idx])
                        idx += 1

                
    pid.close()
    
    

