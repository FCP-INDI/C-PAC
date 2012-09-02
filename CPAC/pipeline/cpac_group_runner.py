import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

import re
import os
import sys
import glob

from CPAC.utils.datasource import create_gpa_dataflow

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

        print >>f, "python CPAC.pipeline.cpac_basc_pipeline.py -c ", str(config_file), " -i ", subject_infos_file

    elif c.runCWAS:

        print >>f, "python CPAC.pipeline.cpac_cwas_pipeline.py -c ", str(config_file), " -i ", subject_infos_file

    elif c.runGroupAnalysis:

        print >>f, "python CPAC.pipeline.cpac_cwas_pipeline.py -c ", str(config_file), " -i ", subject_infos_file, " -r ", resource_file


    f.close()

    commands.getoutput('chmod +x %s' % subject_bash_file )
    print commands.getoutput('qsub  %s ' % (subject_bash_file))


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

        print >>f, "python CPAC.pipeline.cpac_basc_pipeline.py -c ", str(config_file), " -i ", subject_infos_file

    elif c.runCWAS:

        print >>f, "python CPAC.pipeline.cpac_cwas_pipeline.py -c ", str(config_file), " -i ", subject_infos_file

    elif c.runGroupAnalysis:

        print >>f, "python CPAC.pipeline.cpac_cwas_pipeline.py -c ", str(config_file), " -i ", subject_infos_file, " -r ", resource_file


    f.close()

    commands.getoutput('chmod +x %s' % subject_bash_file )
    print commands.getoutput('qsub  %s ' % (subject_bash_file))



if __name__ == "__main__":
    import argparse
    import re

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--config',
                        dest='config',
                        required=True,
                        help='location of config file'
                        )

    parser.add_argument('-s', '--subjects',
                        dest='subjects',
                        required=True,
                        help='location of subjects file'
                        )

    args = parser.parse_args()
    path, fname = os.path.split(os.path.realpath(args.config))
    sys.path.append(path)
    c = __import__(fname.split('.')[0])

    subject_paths = open(args.subjects, 'r').readlines()
    subject_paths = [s.rstrip('\r\n') for s in subject_paths]

    #base_path = os.path.dirname(os.path.commonprefix(subject_paths))
    base_path = c.sinkDirectory

    from collections import defaultdict
    analysis_map = defaultdict(list)
    analysis_map_gp = defaultdict(list)

    for subject_path in subject_paths:
        #Remove the base bath offset
        rs_path = subject_path.replace(base_path, "", 1)
        folders = split_folders(rs_path)

        pipeline_id = folders[1]
        subject_id = folders[2]
        resource_id = folders[3]
        scan_id = folders[4]

        key = subject_path.replace(subject_id, '*')
        analysis_map[(resource_id, key)].append((pipeline_id, subject_id, scan_id, subject_path))

        # separate map for group analysis
        if c.mixed_scan_analysis == True:
            key = key.replace(scan_id, '*')

        analysis_map_gp[(resource_id, key)].append((pipeline_id, subject_id, scan_id, subject_path))

    for resource, glob_key in analysis_map.keys():
        if resource == 'functional_mni':
            if 1 in c.runBASC and not c.runOnGrid:
                from CPAC.cpac_basc_pipeline import prep_basc_workflow
                prep_basc_workflow(c, analysis_map[(resource, glob_key)])

            if 1 in c.runCWAS and not c.runOnGrid:
                from CPAC.cpac_cwas_pipeline import prep_cwas_workflow
                prep_cwas_workflow(c, analysis_map[(resource, glob_key)])


    for resource, glob_key in analysis_map_gp.keys():
        if resource in c.derivative_list:

            if 1 in c.runGroupAnalysis and not c.runOnGrid:
                from CPAC.cpac_group_analysis_pipeline import prep_group_analysis_workflow
                prep_group_analysis_workflow(c, resource, analysis_map_gp[(resource, glob_key)])
        #if re.match(r"^falff", resource):


