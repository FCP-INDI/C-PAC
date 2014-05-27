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



def run(config_file, output_path_file):
    
    # Runs group analysis

    import yaml

    outDump = open('group_runner_outDump.txt', 'a')

    # Load the config file into 'c'
    c = Configuration(yaml.load(open(os.path.realpath(config_file), 'r')))

    subject_paths = []


    # 'output_path_file' is the wildcard-filled path to the 'Derivative Path
    # File' provided in the dialog box when group analysis is first run

    for file in glob.glob(os.path.abspath(output_path_file)):
        path_list = open(file, 'r').readlines()
        subject_paths.extend([s.rstrip('\r\n') for s in path_list])

#
    print >>outDump, 'path_list: ', path_list, '\n\n'
    print >>outDump, 'subject_paths: ', subject_paths, '\n\n'

    # 'subject_paths' is a list of every output from every subject included
    # in the output folder of the run

    # converts the subject_paths list into a set to enforce no duplicates
    set_subject_paths = set(subject_paths)

    # converts the set back into a list
    subject_paths = list(set_subject_paths)


    #base_path = os.path.dirname(os.path.commonprefix(subject_paths))
    base_path = c.outputDirectory

#
    print >>outDump, 'subject_paths, post list(): ', subject_paths, '\n\n'
    outDump.close()

    from collections import defaultdict
    analysis_map = defaultdict(list)
    analysis_map_gp = defaultdict(list)


    for subject_path in subject_paths:

        #Remove the base bath offset
        rs_path = subject_path.replace(base_path, "", 1)

        rs_path = rs_path.lstrip('/')

        folders = split_folders(rs_path)
        
        pipeline_id = folders[0]
        subject_id = folders[1]
        resource_id = folders[2]
        scan_id = folders[3]

        key = subject_path.replace(subject_id, '*')

        # include all of the scans in one model if True
        if c.repeatedMeasures == True:
            key = key.replace(scan_id, '*')


        # 'resource_id' is each type of output
        # 'key' is a path to each and every individual output file,
        # except with the subject ID replaced with a wildcard (*)

        if resource_id in c.derivativeList:

            analysis_map[(resource_id, key)].append((pipeline_id, subject_id, scan_id, subject_path))

            analysis_map_gp[(resource_id, key)].append((pipeline_id, subject_id, scan_id, subject_path))


        # with this loop, 'analysis_map_gp' is a dictionary with a key for
        # each individual output file - and each entry is a list of tuples,
        # one tuple for each subject in the subject list, containing
        # 'subject_path', which is a full path to that output file for that
        # one particular subject



    for resource, glob_key in analysis_map.keys():
        if resource == 'functional_mni':


            if 1 in c.runBASC:

                if not c.runOnGrid:
                    from CPAC.pipeline.cpac_basc_pipeline import prep_basc_workflow
                    prep_basc_workflow(c, analysis_map[(resource, glob_key)])
                else:
                    if 'sge' in c.resourceManager.lower():
                        run_sge_jobs(c, config_file, resource, analysis_map[(resource, glob_key)])

                    elif 'pbs' in c.resourceManager.lower():
                        run_pbs_jobs(c, config_file, resource, analysis_map[(resource, glob_key)])


            if 1 in c.runCWAS:

                if not c.runOnGrid:

                    from CPAC.pipeline.cpac_cwas_pipeline import prep_cwas_workflow
                    prep_cwas_workflow(c, analysis_map[(resource, glob_key)])

                else:
                    if 'sge' in c.resourceManager.lower():
                        run_sge_jobs(c, config_file, resource, analysis_map[(resource, glob_key)])

                    elif 'pbs' in c.resourceManager.lower():
                        run_pbs_jobs(c, config_file, resource, analysis_map[(resource, glob_key)])



    procss = []
        
    for resource, glob_key in analysis_map_gp.keys():

        # 'resource' is each type of output
        # 'glob_key' is a path to each and every individual output file,
        # except with the subject ID replaced with a wildcard (*)
        
        if resource in c.derivativeList:
            
            if 1 in c.runGroupAnalysis:
            
                #get all the motion parameters across subjects
                try:

                    from CPAC.utils import extract_parameters
                    extract_parameters.run(c.outputDirectory)

                except Exception:
                    print "Extract parameters script did not run correctly"

                if not c.runOnGrid:
                    
                    from CPAC.pipeline.cpac_group_analysis_pipeline import prep_group_analysis_workflow
                    procss.append(Process(target=prep_group_analysis_workflow, args=(c, resource, analysis_map_gp[(resource, glob_key)])))
          
    
    
          
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
    
    

