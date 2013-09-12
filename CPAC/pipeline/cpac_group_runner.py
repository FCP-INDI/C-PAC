import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

import re
import os
import sys
import glob

import time

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

    import re
    import os
    import glob
    import yaml

    # Load the config file into 'c'
    c = Configuration(yaml.load(open(os.path.realpath(config_file), 'r')))

    #diag = open(os.path.join('/home/data/Projects/CPAC_Regression_Test/2013-08-19-20_v0-3-1/fsl-model/2013-09-03', 'group_runner_diagnostic.txt'), 'wt')

    #print >>diag, "Config file: ", c
    #print >>diag, ""
    #print >>diag, "Output path file: ", output_path_file
    #print >>diag, ""

    subject_paths = []

    for file in glob.glob(os.path.abspath(output_path_file)):
        path_list = open(file, 'r').readlines()
        subject_paths.extend([s.rstrip('\r\n') for s in path_list])

    #print >>diag, "Subject paths list size: "
    #print >>diag, len(subject_paths)
    #print >>diag, ""

    #print >>diag, "First subject path: "
    #print >>diag, subject_paths[0]
    #print >>diag, ""


    set_subject_paths = set(subject_paths)
    subject_paths = list(set_subject_paths)
    #base_path = os.path.dirname(os.path.commonprefix(subject_paths))
    base_path = c.outputDirectory

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


        #if scan_id == '_scan_rest_1_rest':

        key = subject_path.replace(subject_id, '*')
        analysis_map[(resource_id, key)].append((pipeline_id, subject_id, scan_id, subject_path))

        # separate map for group analysis
        #if c.mixedScanAnalysis == True:
        #    key = key.replace(scan_id, '*')

        analysis_map_gp[(resource_id, key)].append((pipeline_id, subject_id, scan_id, subject_path))


    #print >>diag, ""
    #print >>diag, "Analysis_map_gp dictionary size: "
    #print >>diag, len(analysis_map_gp)
    #print >>diag, ""


    #print >>diag, "Derivative list: "
    #print >>diag, c.derivativeList
    #print >>diag, ""


    timing = open(os.path.join(c.outputDirectory, 'group_analysis_timing.txt'), 'wt')
    #timing = open(os.path.join('/home/data/Projects/CPAC_Regression_Test/2013-08-19-20_v0-3-1/fsl-model/2013-09-03', 'group_analysis_timing.txt'), 'wt')

    sca_roi_runs = 0
    sca_roi_time = 0
    sca_seed_runs = 0
    sca_seed_time = 0
    sca_tempreg_runs = 0
    sca_tempreg_time = 0
    dr_tempreg_runs = 0
    dr_tempreg_time = 0
    vmhc_z_runs = 0
    vmhc_z_time = 0
    alff_Z_runs = 0
    alff_Z_time = 0
    falff_Z_runs = 0
    falff_Z_time = 0
    reho_Z_runs = 0
    reho_Z_time = 0
    centrality_outputs_runs = 0
    centrality_outputs_time = 0

    # Start timing here
    gpa_start_time = time.time()


    for resource, glob_key in analysis_map.keys():
        if resource == 'functional_mni':

            wf_start_time = time.time()

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

            print >>timing, "Group analysis workflow completed for resource: ", resource
            print >>timing, "Elapsed run time (minutes): ", ((time.time() - wf_start_time)/60)
            print >>timing, ""




    for resource, glob_key in analysis_map_gp.keys():

        if resource in c.derivativeList:

            wf_start_time = time.time()

            #print >>diag, "Resource: "
            #print >>diag, resource
            #print >>diag, ""

            #print >>diag, "glob key: "
            #print >>diag, glob_key
            #print >>diag, ""

            #print >>diag, "Analysis map gp entry: "
            #print >>diag, analysis_map_gp[(resource,glob_key)]
            #print >>diag, ""

            if 1 in c.runGroupAnalysis:
              
                #get all the motion parameters across subjects

                try:

                    from CPAC.utils import extract_parameters
                    extract_parameters.run(c.outputDirectory)

                except Exception:

                    print "Extract parameters script did not run correctly"


                if not c.runOnGrid:

                    from CPAC.pipeline.cpac_group_analysis_pipeline import prep_group_analysis_workflow

                    #print c, "   ", resource, "   ", analysis_map_gp[(resource, glob_key)], "   ", glob_key
                    prep_group_analysis_workflow(c, resource, analysis_map_gp[(resource, glob_key)])



                else:

                    if 'sge' in c.resourceManager.lower():
                        
                        run_sge_jobs(c, config_file, resource, analysis_map_gp[(resource, glob_key)])
                       
                    elif 'pbs' in c.resourceManager.lower():
                     
                        run_pbs_jobs(c, config_file, resource, analysis_map_gp[(resource, glob_key)])


            print >>timing, "Group analysis workflow completed for resource: ", resource
            print >>timing, "Elapsed run time (minutes): ", ((time.time() - wf_start_time)/60)
            print >>timing, ""

            # This can be implemented more sleekly using a dictionary, have to do this at some point
            if resource == 'sca_roi_Z_to_standard_smooth':
                sca_roi_runs += 1
                sca_roi_time = sca_roi_time + ((time.time() - wf_start_time)/60)
            elif resource == 'sca_seed_Z_to_standard_smooth':
                sca_seed_runs += 1
                sca_seed_time = sca_seed_time + ((time.time() - wf_start_time)/60)
            elif resource == 'sca_tempreg_maps_z_files_smooth':
                sca_tempreg_runs += 1
                sca_tempreg_time = sca_tempreg_time + ((time.time() - wf_start_time)/60)
            elif resource == 'dr_tempreg_maps_z_files_smooth':
                dr_tempreg_runs += 1
                dr_tempreg_time = dr_tempreg_time + ((time.time() - wf_start_time)/60)
            elif resource == 'vmhc_z_score_stat_map':
                vmhc_z_runs += 1
                vmhc_z_time = vmhc_z_time + ((time.time() - wf_start_time)/60)
            elif resource == 'alff_Z_to_standard_smooth':
                alff_Z_runs += 1
                alff_Z_time = alff_Z_time + ((time.time() - wf_start_time)/60)
            elif resource == 'falff_Z_to_standard_smooth':
                falff_Z_runs += 1
                falff_Z_time = falff_Z_time + ((time.time() - wf_start_time)/60)
            elif resource == 'reho_Z_to_standard_smooth':
                reho_Z_runs += 1
                reho_Z_time = reho_Z_time + ((time.time() - wf_start_time)/60)
            elif resource == 'centrality_outputs_smoothed':
                centrality_outputs_runs += 1
                centrality_outputs_time = centrality_outputs_time + ((time.time() - wf_start_time)/60)
            
    print >>timing, "Entire group analysis run complete."
    print >>timing, "Elapsed run time (minutes): ", ((time.time() - gpa_start_time)/60)
    print >>timing, ""

    print >>timing, "sca_roi_Z_to_standard_smooth"
    print >>timing, "Number of runs: ", sca_roi_runs
    print >>timing, "Total run time (minutes): ", sca_roi_time
    print >>timing, ""

    print >>timing, "sca_seed_Z_to_standard_smooth"
    print >>timing, "Number of runs: ", sca_seed_runs
    print >>timing, "Total run time (minutes): ", sca_seed_time
    print >>timing, ""

    print >>timing, "sca_tempreg_maps_z_files_smooth"
    print >>timing, "Number of runs: ", sca_tempreg_runs
    print >>timing, "Total run time (minutes): ", sca_tempreg_time
    print >>timing, ""

    print >>timing, "dr_tempreg_maps_z_files_smooth"
    print >>timing, "Number of runs: ", dr_tempreg_runs
    print >>timing, "Total run time (minutes): ", dr_tempreg_time
    print >>timing, ""

    print >>timing, "vmhc_z_score_stat_map"
    print >>timing, "Number of runs: ", vmhc_z_runs
    print >>timing, "Total run time (minutes): ", vmhc_z_time
    print >>timing, ""

    print >>timing, "alff_Z_to_standard_smooth"
    print >>timing, "Number of runs: ", alff_Z_runs
    print >>timing, "Total run time (minutes): ", alff_Z_time
    print >>timing, ""

    print >>timing, "falff_Z_to_standard_smooth"
    print >>timing, "Number of runs: ", falff_Z_runs
    print >>timing, "Total run time (minutes): ", falff_Z_time
    print >>timing, ""

    print >>timing, "reho_Z_to_standard_smooth"
    print >>timing, "Number of runs: ", reho_Z_runs
    print >>timing, "Total run time (minutes): ", reho_Z_time
    print >>timing, ""

    print >>timing, "centrality_outputs_smoothed"
    print >>timing, "Number of runs: ", centrality_outputs_runs
    print >>timing, "Total run time (minutes): ", centrality_outputs_time
    print >>timing, ""



    timing.close()
    #diag.close()



