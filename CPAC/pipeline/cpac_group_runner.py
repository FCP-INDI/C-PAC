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



def run(config_file, subject_list_file, output_path_file):
    
    # Runs group analysis

    import yaml

    # Load the config file into 'c'
    c = Configuration(yaml.load(open(os.path.realpath(config_file), 'r')))


    # load the subject list (in the main GUI window, not the group analysis
    # one), and parse the yaml so that the subIDs and session IDs can be
    # accessed for below
    try:
        sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        print "Subject list is not in proper YAML format. Please check your file"
        raise Exception


    subject_paths = []


    # 'output_path_file' is the wildcard-filled path to the 'Derivative Path
    # File' provided in the dialog box when group analysis is first run
    #for file in glob.glob(os.path.abspath(output_path_file)):
    #    path_list = open(file, 'r').readlines()
    #    subject_paths.extend([s.rstrip('\r\n') for s in path_list])
        
           
    ind_outputs = ['alff_to_standard_zstd', 'alff_to_standard_smooth_zstd', 'falff_to_standard_zstd', 'falff_to_standard_smooth_zstd', 'reho_to_standard_zstd', 'reho_to_standard_smooth_zstd', 'sca_roi_files_to_standard_fisher_zstd', 'sca_roi_files_to_standard_smooth_fisher_zstd', 'sca_seed_to_standard_fisher_zstd', 'sca_seed_to_standard_smooth_fisher_zstd', 'sca_tempreg_maps_zstat_files_smooth', 'vmhc_fisher_zstd', 'vmhc_fisher_zstd_zstat_map', 'centrality_outputs_zstd', 'centrality_outputs_smoothed_zstd', 'dr_tempreg_maps_files_to_standard', 'dr_tempreg_maps_files_to_standard_smooth', 'dr_tempreg_maps_zstat_files_to_standard', 'dr_tempreg_maps_zstat_files_to_standard_smooth', 'alff_to_standard', 'alff_to_standard_smooth', 'falff_to_standard', 'falff_to_standard_smooth', 'reho_to_standard', 'reho_to_standard_smooth', 'sca_roi_files_to_standard', 'sca_roi_files_to_standard_smooth', 'sca_seed_to_standard', 'sca_seed_to_standard_smooth', 'sca_tempreg_maps_files', 'sca_tempreg_maps_files_smooth', 'sca_tempreg_maps_zstat_files', 'sca_tempreg_maps_zstat_files_smooth', 'vmhc_raw_score', 'centrality_outputs', 'centrality_outputs_smoothed', 'dr_tempreg_maps_files_to_standard', 'dr_tempreg_maps_files_to_standard_smooth', 'dr_tempreg_maps_zstat_files_to_standard', 'dr_tempreg_maps_zstat_files_to_standard_smooth']
            
            
    
    # collect all of the output paths
    
    for root, folders, files in os.walk(output_path_file):
    
        split_output_dir_path = output_path_file.split("/")
    
        for filename in files:
        
            if filename.endswith("nii.gz"):
    
                fullpath = os.path.join(root, filename)
            
                split_fullpath = fullpath.split("/")
                
                #subID = split_fullpath[len(split_output_dir_path)]
                deriv_folder_name = split_fullpath[len(split_output_dir_path)+1]
            
                #second_half_filepath = fullpath.split(subID)
            
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


    if len(c.modelConfigs) == 0:
        print '[!] CPAC says: You do not have any models selected ' \
              'to run for group-level analysis. Return to your pipeline ' \
              'configuration file and create or select at least one.\n\n'
        raise Exception



    # 'subject_paths' is a list of every output from every subject included
    # in the output folder of the run

    # converts the subject_paths list into a set to enforce no duplicates
    set_subject_paths = set(subject_paths)

    # converts the set back into a list
    subject_paths = list(set_subject_paths)


    #base_path = os.path.dirname(os.path.commonprefix(subject_paths))
    base_path = c.outputDirectory


    from collections import defaultdict
    analysis_map = defaultdict(list)
    analysis_map_gp = defaultdict(list)


    print "Parsing through output paths. This may take a little while " \
          "depending on how many subjects, group analysis models, or " \
          "selected derivatives you have..\n"

    count = 0

    for subject_path in subject_paths:

        # each 'subject_path' is a full filepath to one of the output files

        # Remove the base bath offset
        #rs_path = subject_path.replace(base_path, "", 1)
        #rs_path = rs_path.lstrip('/')

        # rs_path is now the path to the output file, except everything before
        # the pipeline folder (named with the pipeline ID) is stripped from
        # the path

        #folders = split_folders(rs_path)
 
        #pipeline_id = folders[0]
        #subject_unique_id = folders[1]
        #resource_id = folders[2]
        #scan_id = folders[3]


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


        # get list of all unique IDs (session IDs)
        # loop through them and check subject_path for existence of any of the
        # session IDs
        # if it exists, load it into unique_id

        # init subject_id to None
        subject_id = None
        for sub in sublist:
            if sub['subject_id'] in subject_unique_id:
                subject_id = sub['subject_id']

        # If subject_id never gets set for this specific subject, move on to next subject
        if not subject_id:
            continue

        # 'resource_id' is each type of output
        # 'key' is a path to each and every individual output file,
        # except with the subject ID replaced with a wildcard (*)

        # loop here to replace the one below it:
        #     go through model configs, make a list of all ders included
        #     enumerate list of selected derivatives and the models they are in
        #     like: (resource_id, group_model, key)
        for group_config_file in c.modelConfigs:

            try:
                ga_config = Configuration(yaml.load(open(os.path.realpath(group_config_file), 'r')))
            except:
                raise Exception("\n\nError in reading %s configuration file\n\n" % group_config_file)

            if len(ga_config.derivative_list) == 0:
                print '[!] CPAC says: You do not have any derivatives selected ' \
                      'to run for group-level analysis. Return to your group-analysis ' \
                      'configuration file and select at least one.'
                print 'Group analysis configuration file: %s\n\n' % group_config_file
                raise Exception


            if resource_id in ga_config.derivative_list:

                # include all of the scans and sessions in one model if True
                if ga_config.repeated_measures == True:
                    key = subject_path.replace(subject_unique_id, '*')
                    key = key.replace(scan_id, '*')
                else:
                    # each group of subjects from each session go into their own
                    # separate model, instead of combining all sessions into one
                    try:
                        key = subject_path.replace(subject_id, '*')
                    except:
                        # this fires if 'subject_id' was never given a value basically
                        print '\n\n[!] CPAC says: Either the derivative path file ' \
                              'you provided does not contain the output directory ' \
                              'given in the pipeline configuration file.\n'
                        print 'Derivative path file: ', output_path_file, '\n'
                        print 'Output directory: ', c.outputDirectory, '\n'
                        print '- OR -\n'
                        print 'Your subject list does not contain all of the ' \
                              'subjects you wish to run group-level analysis on.\n'
                        print 'Please correct this and try again.\n\n\n'
                        raise Exception


                analysis_map[(resource_name, group_config_file, key)].append((pipeline_id, subject_id, scan_id, subject_path))

                analysis_map_gp[(resource_name, group_config_file, key)].append((pipeline_id, subject_id, scan_id, subject_path))

        count += 1

        if count == int(len(subject_paths)*0.7):
            print "Almost finished parsing output paths.."     

        # with this loop, 'analysis_map_gp' is a dictionary with a key for
        # each individual output file - and each entry is a list of tuples,
        # one tuple for each subject in the subject list, containing
        # 'subject_path', which is a full path to that output file for that
        # one particular subject


    print "Finished parsing through output paths!\n"



    for resource, group_model, glob_key in analysis_map.keys():
        if resource == 'functional_mni':


            if 1 in c.runBASC:

                if not c.runOnGrid:
                    from CPAC.pipeline.cpac_basc_pipeline import prep_basc_workflow
                    prep_basc_workflow(c, analysis_map[(resource, group_model, glob_key)])
                else:
                    if 'sge' in c.resourceManager.lower():
                        run_sge_jobs(c, config_file, resource, analysis_map[(resource, group_model, glob_key)])

                    elif 'pbs' in c.resourceManager.lower():
                        run_pbs_jobs(c, config_file, resource, analysis_map[(resource, group_model, glob_key)])


            if 1 in c.runCWAS:

                if not c.runOnGrid:

                    from CPAC.pipeline.cpac_cwas_pipeline import prep_cwas_workflow
                    prep_cwas_workflow(c, analysis_map[(resource, group_model, glob_key)])

                else:
                    if 'sge' in c.resourceManager.lower():
                        run_sge_jobs(c, config_file, resource, analysis_map[(resource, group_model, glob_key)])

                    elif 'pbs' in c.resourceManager.lower():
                        run_pbs_jobs(c, config_file, resource, analysis_map[(resource, group_model, glob_key)])



    procss = []
    

    for resource, group_model, glob_key in analysis_map_gp.keys():

        # 'resource' is each type of output
        # 'glob_key' is a path to each and every individual output file,
        # except with the subject ID replaced with a wildcard (*)
                      
        #get all the motion parameters across subjects

        print "Pulling motion parameters for all subjects..\n"

        from CPAC.utils import extract_parameters
        scrub_threshold = extract_parameters.run(c.outputDirectory, c.runScrubbing)

        if not c.runOnGrid:
                    
            print "Starting group analysis pipeline setup..\n"

            from CPAC.pipeline.cpac_ga_model_generator import prep_group_analysis_workflow
            procss.append(Process(target=prep_group_analysis_workflow, args=(c, group_model, resource, analysis_map_gp[(resource, group_model, glob_key)], scrub_threshold)))
            
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
    
    

