import argparse
import e_afni
import sys
import os
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from guppy import hpy
import cProfile
import psutil
from sca import create_sca




def get_workflow(c):

    """
    Setup standard file paths and initialize the sca workflow object

    Parameters
    ----------

    c : an object of the config python file which is standard to C-PAC
        Contains all the variables defined in the Config python file

    Returns
    -------
    sca : A workflow
        Seed Based Correlation Analysis Workflow


    """
    standard = os.path.join(c.FSLDIR,
            'data/standard/MNI152_T1_%s.nii.gz' % (c.resolution_for_func))

    sca = None
    sca = create_sca(c.correlationSpace)
    sca.inputs.fwhm_input.fwhm = c.fwhm
    sca.get_node('fwhm_input').iterables = ('fwhm',
                                                c.fwhm)
    sca.inputs.inputspec.standard = standard
    return sca


def prep_workflow(c):

    """
    Calls the get_workflow function to intialize the sca workflow and then it runs it

    Parameters
    ----------

    c : an object of the config python file which is standard to C-PAC
        Contains all the variables defined in the Config python file

    Returns
    -------

    None
    
    Notes
    -----

    All the inputs to the  sca workflow are hard coded , ideally the inputs should be from the benchmark data.

    """

    sca = get_workflow(c)
    sca.base_dir = c.workingDirectory
    sca.crash_dir = c.crashLogDirectory
    sca.config['execution'] = {'hash_method': 'timestamp'}

    sca.inputs.seed_list_input.seed_list = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/settings/seeds/rest_Dickstein_accumbens.nii.gz')
    sca.inputs.inputspec.premat = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/reg/example_func2highres.mat')
    sca.inputs.inputspec.postmat = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/reg/highres2example_func.mat')
    sca.inputs.inputspec.rest_res_filt = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/rest_res_filt.nii.gz')
    sca.inputs.inputspec.fieldcoeff_file = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/anat/reg/highres2standard_warp.nii.gz')
    sca.inputs.inputspec.rest_mask2standard = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/rest_mask2standard.nii.gz')
    sca.inputs.inputspec.ref = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/example_func.nii.gz')
    sca.run()


def main():

    """
    The main purpose of this function is to call the functions that setup and run sca workflow
    This function on its own does the profiling of the workflow

    Parameters
    ----------

    None

    Returns
    -------

    None

    Notes
    -----

    sca_profiler reports Memory Usage, the CPU usage and the IO usage of the sca workflow


    """

    parser = argparse.ArgumentParser(description="example: \
                        run sca_profiler.py -c config.py")
    parser.add_argument('-c', '--config',
                        dest='config',
                        required=True,
                        help='location of config file'
                        )
    args = parser.parse_args()
    path, fname = os.path.split(os.path.realpath(args.config))
    sys.path.append(path)
    c = __import__(fname.split('.')[0])

    p = psutil.Process(os.getpid())
    p.get_cpu_percent(interval=0)
    cpu_before = p.get_cpu_times()
    prep_workflow(c)

    print '\n\nSCA Memory CPU and IO Stats'
    print     '---------------------------'
    print 'Memory Usage: ', p.get_memory_info()
    print 'CPU Times before worlflow is run: ', cpu_before, ' & after workflow is run: ', p.get_cpu_times()
    print 'CPU Usage: ', p.get_cpu_percent(interval=1)
    print 'IO Usage: ', p.get_io_counters()
    print     '---------------------------\n\n'

if __name__ == "__main__":

    main()
