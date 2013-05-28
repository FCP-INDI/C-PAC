from multiprocessing import Pool
from multiprocessing import Process
import sys
import os
import argparse
from CPAC.pipeline import cpac_pipeline
from CPAC.utils.utils import create_seeds_
from CPAC.utils import Configuration
import yaml

def get_vectors(strat):

    paths = []
    def dfs(val_list, path):

        if val_list == []:
            paths.append(path)

        else:
            vals = val_list.pop()
            for val in vals:
                if path == '':
                    dfs(list(val_list), str(val))
                else:
                    dfs(list(val_list), str(val) + '#' + path)


    val_list = []

    for key in sorted(strat.keys()):
        val_list.append(strat[key])

    dfs(val_list, '')

    print paths
    return paths


def make_enteries(paths, path_iterables):

    enteries = []
    idx = 1
    for path in sorted(paths):

        sub_enteries = []
        values = path.split('#')

        indx = 0
        for value in values:

            if '[' or '(' in value:

                value = value.strip('[]')
                value = value.strip('()')

            if ',' in value:
                import re
                value = re.sub(r',', '.', value)
                value = re.sub(r' ', '', value)
            sub_enteries.append(path_iterables[indx] + '_' + value)
            indx += 1


        enteries.append(sub_enteries)

    return enteries



def build_strategies(configuration):

    import collections

    path_iterables = ['_gm_threshold', '_wm_threshold', '_csf_threshold', '_threshold', '_compcor', '_target_angle_deg']
    non_strategy_iterables = ['_fwhm', '_hp', '_lp', '_bandpass_freqs']

    proper_names = {'_threshold':'Scrubbing Threshold = ', '_csf_threshold':'Cerebral Spinal Fluid Threshold = ',
                    '_gm_threshold':'Gray Matter Threshold = ',
                    'nc':'Compcor: Number Of Components = ', '_compcor':'Nuisance Signal Corrections = ',
                    '_target_angle_deg':'Median Angle Correction: Traget Angle in Degree = ', '_wm_threshold':'White Matter Threshold = '}


    config_iterables = {'_gm_threshold': eval('configuration.grayMatterThreshold'), '_wm_threshold': eval('configuration.whiteMatterThreshold'), '_csf_threshold': eval('configuration.cerebralSpinalFluidThreshold'), '_threshold': eval('configuration.scrubbingThreshold'), '_compcor': eval('configuration.Corrections'), '_target_angle_deg': eval('configuration.targetAngleDeg')}


    ### This is really dirty code and ordering of corrections in 
    ### in output directory is dependant on the nuisance workflow
    ### when the workflow is changed , change this section as well
    corrections_order = ['pc1', 'linear', 'wm', 'global', 'motion', 'quadratic', 'gm', 'compcor', 'csf']


    corrections_dict_list = config_iterables['_compcor']


    print corrections_dict_list

    main_all_options = []
    for corrections_dict in corrections_dict_list:
        string = ""
        for correction in corrections_order:

            string += correction + str(corrections_dict[correction]) + '.'

        string = string[0:len(string) -1]



        cmpcor_components = eval('configuration.nComponents')


        all_options = []
        for comp in cmpcor_components:

            all_options.append('ncomponents_%d' %comp + '_selector_' + string)

        main_all_options.append(str(str(all_options).strip('[]')).strip('\'\''))


    config_iterables['_compcor'] = main_all_options


    ############

    paths = get_vectors(config_iterables)

    strategy_enteries = make_enteries(paths, sorted(path_iterables))

    return strategy_enteries




def run_sge_jobs(c, config_file, strategies_file, subject_list_file):


    import commands
    import pickle
    from time import strftime

    try:
        sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        raise Exception ("Subject list is not in proper YAML format. Please check your file")

    shell = commands.getoutput('echo $SHELL')

    temp_files_dir = os.path.join(os.getcwd(), 'cluster_temp_files')
    subject_bash_file = os.path.join(temp_files_dir, 'submit_%s.sge' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    f = open(subject_bash_file, 'w')
    print >>f, '#! %s' % shell
    print >>f, '#$ -cwd'
    print >>f, '#$ -S %s' % shell
    print >>f, '#$ -V'
    print >>f, '#$ -t 1-%d' % len(sublist)
    print >>f, '#$ -q %s' % c.queue
    print >>f, '#$ -pe %s %d' % (c.parallelEnvironment, c.numCoresPerSubject)
    print >>f, '#$ -e %s' % os.path.join(temp_files_dir, 'c-pac_%s.err' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    print >>f, '#$ -o %s' % os.path.join(temp_files_dir, 'c-pac_%s.out' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    print >>f, 'source ~/.bashrc'

#    print >>f, "python CPAC.pipeline.cpac_pipeline.py -c ", str(config_file), " -s ", subject_list_file, " -indx $SGE_TASK_ID  -strategies ", strategies_file
    print >>f, "python -c \"import CPAC; CPAC.pipeline.cpac_pipeline.run(\\\"%s\\\" , \\\"%s\\\", \\\"$SGE_TASK_ID\\\" , \\\"%s\\\", \\\"%s\\\" , \\\"%s\\\", \\\"%s\\\") \" " % (str(config_file), \
        subject_list_file, strategies_file, c.maskSpecificationFile, c.roiSpecificationFile, c.templateSpecificationFile)

    f.close()

    commands.getoutput('chmod +x %s' % subject_bash_file )
    print commands.getoutput('qsub  %s ' % (subject_bash_file))


def run_condor_jobs(c, config_file, strategies_file, subject_list_file):


    import commands
    import pickle
    from time import strftime

    try:
        sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        raise Exception ("Subject list is not in proper YAML format. Please check your file")

    temp_files_dir = os.path.join(os.getcwd(), 'cluster_temp_files')
    subject_bash_file = os.path.join(temp_files_dir, 'submit_%s.condor' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    f = open(subject_bash_file, 'w')

    print >>f, "Executable = /usr/bin/python"
    print >>f, "Universe = vanilla"
    print >>f, "transfer_executable = False"
    print >>f, "getenv = True"
    print >>f, "log = %s" % os.path.join(temp_files_dir, 'c-pac_%s.log' % str(strftime("%Y_%m_%d_%H_%M_%S")))

    sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    for sidx in range(1,len(sublist)+1):
        print >>f, "error = %s" % os.path.join(temp_files_dir, 'c-pac_%s.%s.err' % (str(strftime("%Y_%m_%d_%H_%M_%S")), str(sidx)))
        print >>f, "output = %s" % os.path.join(temp_files_dir, 'c-pac_%s.%s.out' % (str(strftime("%Y_%m_%d_%H_%M_%S")), str(sidx)))

        print >>f, "arguments = \"-c 'import CPAC; CPAC.pipeline.cpac_pipeline.run( ''%s'',''%s'',''%s'',''%s'', ''%s'',''%s'',''%s'')\'\"" % (str(config_file), subject_list_file, str(sidx), strategies_file, c.maskSpecificationFile, c.roiSpecificationFile, c.templateSpecificationFile)
        print >>f, "queue"

    f.close()

    #commands.getoutput('chmod +x %s' % subject_bash_file )
    print commands.getoutput("condor_submit %s " % (subject_bash_file))





def run_pbs_jobs(c, config_file, strategies_file, subject_list_file):



    import commands
    import pickle
    from time import strftime


    try:
        sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        raise Exception ("Subject list is not in proper YAML format. Please check your file")
    
    temp_files_dir = os.path.join(os.getcwd(), 'cluster_temp_files')
    shell = commands.getoutput('echo $SHELL')
    subject_bash_file = os.path.join(temp_files_dir, 'submit_%s.pbs' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    f = open(subject_bash_file, 'w')
    print >>f, '#! %s' % shell
    print >>f, '#PBS -S %s' % shell
    print >>f, '#PBS -V'
    print >>f, '#PBS -t 1-%d' % len(sublist)
    print >>f, '#PBS -q %s' % c.queue
    print >>f, '#PBS -l nodes=1:ppn=%d' % c.numCoresPerSubject
    print >>f, '#PBS -e %s' % os.path.join(temp_files_dir, 'c-pac_%s.err' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    print >>f, '#PBS -o %s' % os.path.join(temp_files_dir, 'c-pac_%s.out' % str(strftime("%Y_%m_%d_%H_%M_%S")))
    print >>f, 'source ~/.bashrc'

    print >>f, "python -c \"import CPAC; CPAC.pipeline.cpac_pipeline.run(\\\"%s\\\",\\\"%s\\\",\\\"${PBS_ARRAYID}\\\",\\\"%s\\\", \\\"%s\\\" , \\\"%s\\\", \\\"%s\\\") \" " % (str(config_file), \
        subject_list_file, strategies_file, c.maskSpecificationFile, c.roiSpecificationFile, c.templateSpecificationFile)
#    print >>f, "python -c \"import CPAC; CPAC.pipeline.cpac_pipeline.py -c %s -s %s -indx ${PBS_ARRAYID} -strategies %s \" " %(str(config_file), subject_list_file, strategies_file)
    #print >>f, "python CPAC.pipeline.cpac_pipeline.py -c ", str(config_file), "-s ", subject_list_file, " -indx ${PBS_ARRAYID} -strategies ", strategies_file
    f.close()

    commands.getoutput('chmod +x %s' % subject_bash_file )
    print commands.getoutput('qsub  %s ' % (subject_bash_file))


def append_seeds_to_file(working_dir, seed_list, seed_file):

    existing_seeds = []
    filtered_list = []
    try:

        if os.path.isfile(seed_file):
            existing_seeds += [line.rstrip('\r\n') for line in open(seed_file, 'r').readlines() if not (line.startswith('#') and line == '\n')]

            for seed in seed_list:
                if not seed in existing_seeds:
                    filtered_list.append(seed)

            if not len(filtered_list) == 0:
                f = open(seed_file, 'a')
                for seed in filtered_list:
                    f.write("%s\n" % seed)
                f.close()

            return seed_file

        else:
            raise


    except:
        #make tempfile and add seeds to it
        import tempfile

        try:
            if not os.path.exists(working_dir):
                os.makedirs(working_dir)

        except Exception, e:

            print 'error encountered : ', e
            raise

        some_number, f_name = tempfile.mkstemp(suffix='.txt', prefix='temp_roi_seeds', dir=working_dir, text=True)

        f_handle = open(f_name, 'w')

        for seed in seed_list:
            f_handle.write('%s\n' % seed)

        f_handle.close()
        return f_name



def run(config_file, subject_list_file):
    
    try:
    
        if not os.path.exists(config_file):
            raise IOError("config file %s doesn't exist" %config_file)
        else:
            c = Configuration(yaml.load(open(os.path.realpath(config_file), 'r')))
    
    except Exception:
        print "Error reading config file - %s"%config_file 

    try:
        sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        raise Exception ("Subject list is not in proper YAML format. Please check your file")

    strategies = sorted(build_strategies(c))

    seeds_created = []
    if not (c.seedSpecificationFile is None):

        try:
            if os.path.exists(c.seedSpecificationFile):
                seeds_created = create_seeds_(c.seedOutputLocation, c.seedSpecificationFile, c.FSLDIR)
                print 'seeds created %s -> ' % seeds_created
        except:
            raise IOError('Problem in seedSpecificationFile')

    if 1 in c.runVoxelTimeseries:

        if 2 in c.useSeedInAnalysis:

            c.maskSpecificationFile = append_seeds_to_file(c.workingDirectory, seeds_created, c.maskSpecificationFile)

    if 1 in c.runROITimeseries:

        if 1 in c.useSeedInAnalysis:

            c.roiSpecificationFile = append_seeds_to_file(c.workingDirectory, seeds_created, c.roiSpecificationFile)

    if 1 in c.runNetworkCentrality:

        if 3 in c.useSeedInAnalysis:

            c.templateSpecificationFile = append_seeds_to_file(c.workingDirectory, seeds_created, c.templateSpecificationFile)



    if not c.runOnGrid:

        from CPAC.pipeline.cpac_pipeline import prep_workflow
        procss = [Process(target=prep_workflow, args=(sub, c, strategies)) for sub in sublist]

        jobQueue = []
        if len(sublist) <= c.numSubjectsAtOnce:
            """
            Stream all the subjects as sublist is
            less than or equal to the number of 
            subjects that need to run
            """
            for p in procss:
                p.start()


        else:

            """
            Stream the subject worlflows for preprocessing.
            At Any time in the pipeline c.numSubjectsAtOnce
            will run, unless the number remaining is less than
            the value of the parameter stated above
            """
            idx = 0
            while(idx < len(sublist)):

                if len(jobQueue) == 0 and idx == 0:

                    idc = idx
                    for p in procss[idc: idc + c.numSubjectsAtOnce]:

                        p.start()
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



    else:

        import commands
        import pickle
        from time import strftime

        temp_files_dir = os.path.join(os.getcwd(), 'cluster_temp_files')
        print commands.getoutput("mkdir -p %s" % temp_files_dir)


        strategies_file = os.path.join(temp_files_dir, 'strategies.obj')
        f = open(strategies_file, 'w')
        pickle.dump(strategies, f)
        f.close()




        if 'sge' in c.resourceManager.lower():

            run_sge_jobs(c, config_file, strategies_file, subject_list_file)


        elif 'pbs' in c.resourceManager.lower():

            run_pbs_jobs(c, config_file, strategies_file, subject_list_file)

        elif 'condor' in c.resourceManager.lower():

            run_condor_jobs(c, config_file, strategies_file, subject_list_file)
