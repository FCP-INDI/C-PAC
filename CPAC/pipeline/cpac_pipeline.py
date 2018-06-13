# CPAC/pipeline/cpac_pipeline.py

'''
This module prepares and executes the main C-PAC workflow
'''

# Import packages
import os
import time
from time import strftime
import zlib
import linecache
import csv
import pickle
import pandas as pd
import pkg_resources as p

# Nipype packages
import nipype
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
from nipype.interfaces.afni import preprocess
from nipype.pipeline.engine.utils import format_dot
import nipype.interfaces.ants as ants
import nipype.interfaces.c3 as c3
from nipype import config
from nipype import logging

# INDI-Tools
from indi_aws import aws_utils, fetch_creds

# CPAC packages
import CPAC
from CPAC import network_centrality
from CPAC.network_centrality.utils import merge_lists
from CPAC.anat_preproc.anat_preproc import create_anat_preproc
from CPAC.EPI_DistCorr.EPI_DistCorr import create_EPI_DistCorr
from CPAC.func_preproc.func_preproc import create_func_preproc, \
    create_wf_edit_func
from CPAC.seg_preproc.seg_preproc import create_seg_preproc

from CPAC.registration import create_nonlinear_register, \
    create_register_func_to_anat, \
    create_bbregister_func_to_anat, \
    create_wf_calculate_ants_warp, \
    create_wf_apply_ants_warp, \
    create_wf_c3d_fsl_to_itk, \
    create_wf_collect_transforms
from CPAC.nuisance import create_nuisance, bandpass_voxels

from CPAC.median_angle import create_median_angle_correction
from CPAC.generate_motion_statistics import motion_power_statistics
from CPAC.generate_motion_statistics import fristons_twenty_four
from CPAC.scrubbing import create_scrubbing_preproc
from CPAC.timeseries import create_surface_registration, get_roi_timeseries, \
    get_voxel_timeseries, get_vertices_timeseries, \
    get_spatial_map_timeseries
from CPAC.network_centrality import create_resting_state_graphs, \
    get_cent_zscore
from CPAC.utils.datasource import *
from CPAC.utils import Configuration, create_all_qc

# TODO - QA pages - re-introduce these
from CPAC.qc.qc import create_montage, create_montage_gm_wm_csf
from CPAC.qc.utils import register_pallete, make_edge, drop_percent_, \
    gen_histogram, gen_plot_png, gen_motion_plt, \
    gen_std_dev, gen_func_anat_xfm, gen_snr, \
    generateQCPages, cal_snr_val

from CPAC.utils.utils import extract_one_d, set_gauss, \
    process_outputs, get_scan_params, \
    get_tr, extract_txt, create_log, \
    create_log_template, extract_output_mean, \
    create_output_mean_csv, get_zscore, \
    get_fisher_zscore, dbg_file_lineno, add_afni_prefix
from CPAC.vmhc.vmhc import create_vmhc
from CPAC.reho.reho import create_reho
from CPAC.alff.alff import create_alff
from CPAC.sca.sca import create_sca, create_temporal_reg

# Init variables
logger = logging.getLogger('workflow')


class strategy:
    def __init__(self):
        self.resource_pool = {}
        self.leaf_node = None
        self.leaf_out_file = None
        self.name = []

    def append_name(self, name):
        self.name.append(name)

    def get_name(self):
        return self.name

    def set_leaf_properties(self, node, out_file):
        self.leaf_node = node
        self.leaf_out_file = out_file

    def get_leaf_properties(self):
        return self.leaf_node, self.leaf_out_file

    def get_resource_pool(self):
        return self.resource_pool

    def get_node_from_resource_pool(self, resource_key):
        try:
            return self.resource_pool[resource_key]
        except:
            logger.info('no node for output: ')
            logger.info(resource_key)
            raise

    def update_resource_pool(self, resources):
        for key, value in resources.items():
            if key in self.resource_pool:
                logger.info('Warning key %s already exists in resource' \
                            ' pool, replacing with %s ' % (key, value))
            self.resource_pool[key] = value


def create_new_fork(strat):
    """Create a new strategy fork in the pipeline by producing a copy of the
    current strategy (with resource pool)."""

    tmp = strategy()
    tmp.resource_pool = dict(strat.resource_pool)
    tmp.leaf_node = strat.leaf_node
    tmp.leaf_out_file = str(strat.leaf_out_file)
    tmp.name = list(strat.name)
    strat = tmp

    return strat


# Create and prepare C-PAC pipeline workflow
def prep_workflow(sub_dict, c, strategies, run, pipeline_timing_info=None,
                  p_name=None, plugin='MultiProc', plugin_args=None):
    '''
    Function to prepare and, optionally, run the C-PAC workflow

    Parameters
    ----------
    sub_dict : dictionary
        subject dictionary with anatomical and functional image paths
    c : Configuration object
        CPAC pipelin configuration dictionary object
    strategies : obj
        strategies object describing what strategies to run the pipeline
        through
    run : boolean
        flag to indicate whether to run the prepared workflow
    pipeline_timing_info : list (optional); default=None
        list of pipeline info for reporting timing information
    p_name : string (optional); default=None
        name of pipeline
    plugin : string (optional); defaule='MultiProc'
        nipype plugin to utilize when the workflow is ran
    plugin_args : dictionary (optional); default=None
        plugin-specific arguments for the workflow plugin

    Returns
    -------
    workflow : nipype workflow
        the prepared nipype workflow object containing the parameters
        specified in the config
    '''

    # Import packages
    from CPAC.utils.utils import check_config_resources, check_system_deps

    # Settle some things about the resource pool keys and the output directory
    keys_csv = p.resource_filename('CPAC', 'resources/cpac_outputs.csv')

    try:
        keys = pd.read_csv(keys_csv)
    except Exception as e:
        err = "\n[!] Could not access or read the cpac_outputs.csv " \
              "resource file:\n{0}\n\nError details {1}\n".format(keys_csv, e)
        raise Exception(err)

    # TODO: this way of pulling from the dataframe produces a warning
    # TODO: also may want to put these into a function of some sort

    # outputs marked as optional in the matrix file, but we want them to be
    # written out no matter what for a specific reason
    override_optional = list(keys[keys['Override optional'] == 'yes']['Resource'])

    # extra outputs that we don't write to the output directory, unless the
    # user selects to do so
    debugging_outputs = list(keys[keys['Optional outputs: Debugging outputs'] == 'yes']['Resource'])

    # outputs to write out if the user selects to write all the functional
    # resources and files CPAC generates
    extra_functional_outputs = list(keys[keys['Optional outputs: Extra functionals'] == 'yes']['Resource'])

    # outputs to send into smoothing, if smoothing is enabled, and
    # outputs to write out if the user selects to write non-smoothed outputs
    # "_mult" is for items requiring mapnodes
    outputs_native_nonsmooth = list(keys[keys['Optional outputs: Native space'] == 'yes'][keys['Optional outputs: Non-smoothed'] == 'yes'][keys['Multiple outputs'] != 'yes']['Resource'])
    outputs_native_nonsmooth_mult = list(keys[keys['Optional outputs: Native space'] == 'yes'][keys['Optional outputs: Non-smoothed'] == 'yes'][keys['Multiple outputs'] == 'yes']['Resource'])
    outputs_template_nonsmooth = list(keys[keys['Space'] == 'template'][keys['Optional outputs: Non-smoothed'] == 'yes'][keys['Multiple outputs'] != 'yes']['Resource'])
    outputs_template_nonsmooth_mult = list(keys[keys['Space'] == 'template'][keys['Optional outputs: Non-smoothed'] == 'yes'][keys['Multiple outputs'] == 'yes']['Resource'])

    # don't write these, unless the user selects to write native-space outputs
    outputs_native_smooth = list(keys[keys['Space'] != 'template'][keys['Derivative'] == 'yes'][keys['Optional outputs: Non-smoothed'] != 'yes']['Resource'])

    # ever used??? contains template-space, smoothed, both raw and z-scored
    outputs_template_smooth = list(keys[keys['Space'] == 'template'][keys['Derivative'] == 'yes'][keys['Optional outputs: Non-smoothed'] != 'yes']['Resource'])

    # outputs to send into z-scoring, if z-scoring is enabled, and
    # outputs to write out if user selects to write non-z-scored outputs
    # "_mult" is for items requiring mapnodes
    outputs_template_raw = list(keys[keys['Space'] == 'template'][keys['Multiple outputs'] != 'yes'][keys['Optional outputs: Raw scores'] == 'yes']['Resource'])
    outputs_template_raw_mult = list(keys[keys['Space'] == 'template'][keys['Multiple outputs'] == 'yes'][keys['Optional outputs: Raw scores'] == 'yes']['Resource'])

    # outputs to send into the average calculation nodes
    # "_mult" is for items requiring mapnodes
    outputs_average = list(keys[keys['Calculate averages'] == 'yes'][keys['Multiple outputs'] != 'yes']['Resource'])
    outputs_average_mult = list(keys[keys['Calculate averages'] == 'yes'][keys['Multiple outputs'] == 'yes']['Resource'])

    # Start timing here
    pipeline_start_time = time.time()
    # at end of workflow, take timestamp again, take time elapsed and check
    # tempfile add time to time data structure inside tempfile, and increment
    # number of subjects

    cores_msg = 'VERSION: CPAC %s' % CPAC.__version__

    # Check pipeline config resources
    sub_mem_gb, num_cores_per_sub, num_ants_cores = \
        check_config_resources(c)

    if plugin_args:
        plugin_args['memory_gb'] = sub_mem_gb
        plugin_args['n_procs'] = num_cores_per_sub
    else:
        plugin_args = {'memory_gb': sub_mem_gb, 'n_procs': num_cores_per_sub}

    # perhaps in future allow user to set threads maximum
    # this is for centrality mostly    
    # import mkl
    numThreads = '1'

    os.environ['OMP_NUM_THREADS'] = '1'  # str(num_cores_per_sub)
    os.environ['MKL_NUM_THREADS'] = '1'  # str(num_cores_per_sub)
    os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = str(num_ants_cores)

    # calculate maximum potential use of cores according to current pipeline
    # configuration
    max_core_usage = int(c.maxCoresPerParticipant) * \
                     int(c.numParticipantsAtOnce)

    cores_msg = cores_msg + "\n\nSetting maximum number of cores per " \
                            "participant to %s\n" % c.maxCoresPerParticipant

    cores_msg = cores_msg + 'Setting number of participants at once to %s\n' \
                            % c.numParticipantsAtOnce

    cores_msg = cores_msg + 'Setting OMP_NUM_THREADS to %s\n' % numThreads
    cores_msg = cores_msg + 'Setting MKL_NUM_THREADS to %s\n' % numThreads

    if 'ANTS' in c.regOption:
        cores_msg = cores_msg + 'Setting ANTS/ITK thread usage to %d\n\n' \
                                % c.num_ants_threads

    cores_msg = cores_msg + 'Maximum potential number of cores that might ' \
                            'be used during this run: %d\n\n' % max_core_usage

    logger.info(cores_msg)

    qc_montage_id_a = {}
    qc_montage_id_s = {}
    qc_plot_id = {}
    qc_hist_id = {}
    if sub_dict['unique_id']:
        subject_id = sub_dict['subject_id'] + "_" + sub_dict['unique_id']
    else:
        subject_id = sub_dict['subject_id']

    log_dir = os.path.join(c.logDirectory, subject_id)

    if not os.path.exists(log_dir):
        os.makedirs(os.path.join(log_dir))

    # temp
    already_skullstripped = c.already_skullstripped[0]
    if already_skullstripped == 2:
        already_skullstripped = 0
    elif already_skullstripped == 3:
        already_skullstripped = 1

    subject_info = {}
    subject_info['subject_id'] = subject_id
    subject_info['start_time'] = pipeline_start_time
    subject_info['strategies'] = strategies

    '''
    input filepaths check and tool setup check
    '''

    # this section checks all of the file paths provided in the pipeline
    # config yaml file and ensures the files exist and are accessible

    pipeline_config_map = c.return_config_elements()

    wrong_filepath_list = []

    for config_item in pipeline_config_map:

        label = config_item[0]
        val = config_item[1]

        # ensures it is only checking file paths
        if isinstance(val, str) and '/' in val:
            if ('.txt' in val) or ('.nii' in val) or ('.nii.gz' in val) \
                    or ('.mat' in val) or ('.cnf' in val) or ('.sch' in val):
                if not os.path.isfile(val):
                    wrong_filepath_list.append((label, val))

    if len(wrong_filepath_list) > 0:
        print '\n\n'
        print 'Whoops! - Filepaths provided do not exist:\n'

        for file_tuple in wrong_filepath_list:
            print file_tuple[0], ' - ', file_tuple[1]

        print '\nPlease double-check your pipeline configuration file.\n\n'

        # VERY TEMPORARY
        if (len(wrong_filepath_list) == 1) and \
                (wrong_filepath_list[0][0] == "dilated_symmetric_brain_mask"):
            pass
        else:
            raise Exception

    # Check system dependencies
    if 'ANTS' in c.regOption:
        check_ants = True
    else:
        check_ants = False
    check_system_deps(check_ants)

    '''
    workflow preliminary setup
    '''

    wfname = 'resting_preproc_' + str(subject_id)
    workflow = pe.Workflow(name=wfname)
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {'hash_method': 'timestamp',
                                    'crashdump_dir': os.path.abspath(
                                        c.crashLogDirectory)}

    try:
        if c.run_logging == True:
            config.update_config(
                {'logging': {'log_directory': log_dir, 'log_to_file': True}})
        else:
            config.update_config(
                {'logging': {'log_to_file': False}})
    except AttributeError:
        config.update_config(
            {'logging': {'log_directory': log_dir, 'log_to_file': True}})

    logging.update_logging(config)

    if c.reGenerateOutputs is True:
        import commands
        cmd = "find %s -name \'*sink*\' -exec rm -rf {} \\;" % os.path.join(
            c.workingDirectory, wfname)
        logger.info(cmd)
        commands.getoutput(cmd)
        cmd = "find %s -name \'*link*\' -exec rm -rf {} \\;" % os.path.join(
            c.workingDirectory, wfname)
        logger.info(cmd)
        commands.getoutput(cmd)
        cmd = "find %s -name \'*log*\' -exec rm -rf {} \\;" % os.path.join(
            c.workingDirectory, wfname)
        logger.info(cmd)
        commands.getoutput(cmd)

    def create_log_node(wflow, output, indx, scan_id=None):
        # call logging workflow

        if wflow:
            log_wf = create_log(wf_name='log_%s' % wflow.name)
            log_wf.inputs.inputspec.workflow = wflow.name
            log_wf.inputs.inputspec.index = indx
            log_wf.inputs.inputspec.log_dir = log_dir
            workflow.connect(wflow, output, log_wf, 'inputspec.inputs')
        else:
            log_wf = create_log(wf_name='log_done_%s' % scan_id,
                                scan_id=scan_id)
            log_wf.base_dir = log_dir
            log_wf.inputs.inputspec.workflow = 'DONE'
            log_wf.inputs.inputspec.index = indx
            log_wf.inputs.inputspec.log_dir = log_dir
            log_wf.inputs.inputspec.inputs = log_dir
            return log_wf

    def logStandardError(sectionName, errLine, errNum, errInfo=None):
        logger.info("\n\nERROR: {0} - {1}\n\nError name: cpac_pipeline"
                    "_{2}\n\n".format(sectionName, errLine, errNum))
        if errInfo:
            logger.info("Error details: {0}\n\n".format(errInfo))

    def logConnectionError(workflow_name, numStrat, resourcePool, errNum,
                           errInfo=None):
        logger.info(
            "\n\n" + 'ERROR: Invalid Connection: %s: %s, resource_pool: %s' \
            % (workflow_name, numStrat,
               resourcePool) + "\n\n" + "Error name: cpac_pipeline_%s" % (
            errNum) + \
            "\n\n" + "This is a pipeline creation error - the workflows "
                     "have not started yet." + "\n\n")
        if errInfo:
            logger.info(str(errInfo))

    def logStandardWarning(sectionName, warnLine):
        logger.info(
            "\n\n" + 'WARNING: %s - %s' % (sectionName, warnLine) + "\n\n")

    def getNodeList(strategy):
        nodes = []
        for node in strategy.name:
            nodes.append(node[:-2])
        return nodes

    strat_list = []

    workflow_bit_id = {}
    workflow_counter = 0

    """""""""""""""""""""""""""""""""""""""""""""""""""
     PREPROCESSING
    """""""""""""""""""""""""""""""""""""""""""""""""""

    '''
    Initialize Anatomical Input Data Flow
    '''

    num_strat = 0

    strat_initial = strategy()

    # bundle_node = pe.Node(util.Function(input_names=))

    # Extract credentials path if it exists
    try:
        creds_path = sub_dict['creds_path']
        if creds_path and 'none' not in creds_path.lower():
            if os.path.exists(creds_path):
                input_creds_path = os.path.abspath(creds_path)
            else:
                err_msg = 'Credentials path: "%s" for subject "%s" was not ' \
                          'found. Check this path and try again.' % (
                          creds_path, subject_id)
                raise Exception(err_msg)
        else:
            input_creds_path = None
    except KeyError:
        input_creds_path = None

    flow = create_anat_datasource()
    flow.inputs.inputnode.subject = subject_id
    flow.inputs.inputnode.anat = sub_dict['anat']
    flow.inputs.inputnode.creds_path = input_creds_path
    flow.inputs.inputnode.dl_dir = c.workingDirectory

    anat_flow = flow.clone('anat_gather_%d' % num_strat)

    strat_initial.set_leaf_properties(anat_flow, 'outputspec.anat')

    num_strat += 1

    strat_list.append(strat_initial)

    '''
    Inserting Anatomical Preprocessing workflow
    '''
    new_strat_list = []
    num_strat = 0

    workflow_bit_id['anat_preproc'] = workflow_counter

    for strat in strat_list:
        # create a new node, Remember to change its name!
        # anat_preproc = create_anat_preproc(use_AFNI, already_skullstripped).clone(
        #    'anat_preproc_%d' % num_strat)
      if already_skullstripped == 0:
        if "AFNI" in c.skullstrip_option:
            anat_preproc = create_anat_preproc(True,already_skullstripped,wf_name = 'anat_preproc_%d' % num_strat)
            
            try:
               node, out_file = strat.get_leaf_properties()
               workflow.connect(node, out_file, anat_preproc, 'inputspec.anat')
               anat_preproc.inputs.AFNI_options.shrink_factor = c.shrink_factor
               anat_preproc.inputs.AFNI_options.var_shrink_fac = c.var_shrink_fac
               anat_preproc.inputs.AFNI_options.shrink_factor_bot_lim = c.shrink_factor_bot_lim
               anat_preproc.inputs.AFNI_options.avoid_vent = c.avoid_vent
               anat_preproc.inputs.AFNI_options.niter = c.n_iterations
               anat_preproc.inputs.AFNI_options.pushout = c.pushout
               anat_preproc.inputs.AFNI_options.touchup = c.touchup
               anat_preproc.inputs.AFNI_options.fill_hole = c.fill_hole
               anat_preproc.inputs.AFNI_options.avoid_eyes = c.avoid_eyes
               anat_preproc.inputs.AFNI_options.use_edge = c.use_edge
               anat_preproc.inputs.AFNI_options.exp_frac = c.exp_frac
               anat_preproc.inputs.AFNI_options.smooth_final = c.smooth_final
               anat_preproc.inputs.AFNI_options.push_to_edge = c.push_to_edge
               anat_preproc.inputs.AFNI_options.use_skull = c.use_skull
               anat_preproc.inputs.AFNI_options.perc_init = c.perc_init
               anat_preproc.inputs.AFNI_options.max_inter_iter = c.max_inter_iter
               anat_preproc.inputs.AFNI_options.blur_fwhm = c.blur_fwhm
               anat_preproc.inputs.AFNI_options.fac = c.fac
            
            
               anat_preproc.get_node('AFNI_options.shrink_factor').iterables = ('shrink_factor',c.shrink_factor)
               anat_preproc.get_node('AFNI_options.var_shrink_fac').iterables = ('var_shrink_fac',c.var_shrink_fac)
               anat_preproc.get_node('AFNI_options.shrink_factor_bot_lim').iterables = ('shrink_factor_bot_lim',c.shrink_factor_bottom_lim)
               anat_preproc.get_node('AFNI_options.avoid_vent').iterables = ('avoid_vent',c.avoid_vent)
               anat_preproc.get_node('AFNI_options.niter').iterables = ('niter',c.n_iterations)
               anat_preproc.get_node('AFNI_options.pushout').iterables = ('pushout',c.pushout)
               anat_preproc.get_node('AFNI_options.touchup').iterables = ('touchup',c.touchup)
               anat_preproc.get_node('AFNI_options.fill_hole').iterables = ('fill_hole',c.fill_hole)
               anat_preproc.get_node('AFNI_options.avoid_eyes').iterables = ('avoid_eyes',c.avoid_eyes)
               anat_preproc.get_node('AFNI_options.use_edge').iterables = ('use_edge',c.use_edge)
               anat_preproc.get_node('AFNI_options.exp_frac').iterables = ('exp_frac',c.exp_frac)
               anat_preproc.get_node('AFNI_options.smooth_final').iterables = ('smooth_final',c.smooth_final)
               anat_preproc.get_node('AFNI_options.push_to_edge').iterables = ('push_to_edge',c.push_to_edge)
               anat_preproc.get_node('AFNI_options.use_skull').iterables = ('use_skull',c.use_skull)
               anat_preproc.get_node('AFNI_options.PercInit').iterables = ('perc_init',c.perc_init)
               anat_preproc.get_node('AFNI_options.max_inter_iter').iterables = ('max_inter_iter',c.max_inter_iter)
               anat_preproc.get_node('AFNI_options.blur_FWHM').iterables = ('blur_FWHM',c.blurFWHM)
               anat_preproc.get_node('AFNI_options.fac').iterables = ('fac',c.fac)
            except:
                logConnectionError('Anatomical Preprocessing No valid Previous for strat',
                num_strat, strat.get_resource_pool(), '0001')
                continue
             
           
            
       
            if "BET" in c.skullstrip_option:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = strat.leaf_node
                tmp.leaf_out_file = (strat.leaf_out_file)
                tmp.name=list(strat.name)
                strat=tmp
                new_strat_list.append(strat)
            
            strat.append_name(anat_preproc.name)    
            strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
            # add stuff to resource pool if we need it

            strat.update_resource_pool(
            {'anatomical_brain': (anat_preproc, 'outputspec.brain')})
            strat.update_resource_pool(
            {'anatomical_reorient': (anat_preproc, 'outputspec.reorient')})
         # write to log
            create_log_node(anat_preproc, 'outputspec.brain', num_strat)  
            num_strat += 1 
    
    strat_list += new_strat_list
    
    new_strat_list = []
    
    for strat in strat_list:

        nodes = getNodeList(strat)
        if ("BET" in c.skullstrip_option) and ('anat_preproc' not in nodes):
               anat_preproc = create_anat_preproc(False,already_skullstripped,wf_name = 'anat_preproc_%d' % num_strat)
               
               try:
                   node, out_file = strat.get_leaf_properties()
                   workflow.connect(node, out_file, anat_preproc, 'inputspec.anat')
                   
                   anat_preproc.inputs.BET_options.center = c.center
                  # anat_preproc.inputs.BET_options.center_y = c.center_y
                   #anat_preproc.inputs.BET_options.center_z = c.center_z
                   anat_preproc.inputs.BET_options.frac = c.frac
                   anat_preproc.inputs.BET_options.mask_boolean = c.mask_boolean
                   anat_preproc.inputs.BET_options.mesh_boolean = c.mesh_boolean
                   anat_preproc.inputs.BET_options.outline = c.outline
                   anat_preproc.inputs.BET_options.padding = c.padding
                   anat_preproc.inputs.BET_options.radius = c.radius
                   anat_preproc.inputs.BET_options.reduce_bias = c.reduce_bias
                   anat_preproc.inputs.BET_options.remove_eyes = c.remove_eyes
                   anat_preproc.inputs.BET_options.robust = c.robust
                   anat_preproc.inputs.BET_options.skull = c.skull
                   anat_preproc.inputs.BET_options.surfaces = c.surfaces
                   anat_preproc.inputs.BET_options.threshold = c.threshold
                   anat_preproc.inputs.BET_options.vertical_gradient = c.vertical_gradient
                   
                   anat_preproc.get_node('BET_options.center').iterables = ('center',c.center)
                   #anat_preproc.get_node('BET_options.center_y').iterables = ('center_y',c.center_y)
                   #anat_preproc.get_node('BET_options.center_z').iterables = ('center_z',c.center_z)
                   #anat_preproc.get_node('AFNI_options.shrink_factor').iterables = ('shrink_factor',c.shrink_factor)
                   anat_preproc.get_node('BET_options.mask_boolean').iterables = ('mask_boolean',c.mask_boolean)
                   anat_preproc.get_node('BET_options.mesh_boolean').iterables = ('mesh_boolean',c.mesh_boolean)
                   anat_preproc.get_node('BET_options.outline').iterables = ('outline',c.outline)
                   anat_preproc.get_node('BET_options.padding').iterables = ('padding',c.padding)
                   anat_preproc.get_node('BET_options.radius').iterables = ('radius',c.radius)
                   anat_preproc.get_node('BET_options.reduce_bias').iterables = ('reduce_bias',c.reduce_bias)
                   anat_preproc.get_node('BET_options.remove_eyes').iterables = ('remove_eyes',c.remove_eyes)
                   anat_preproc.get_node('BET_options.robust').iterables = ('robust',c.robust)
                   anat_preproc.get_node('BET_options.skull').iterables = ('skull',c.skull)
                   anat_preproc.get_node('BET_options.surfaces').iterables = ('surfaces',c.surfaces)
                   anat_preproc.get_node('BET_options.threshold').iterables = ('threshold',c.threshold)
                   anat_preproc.get_node('BET_options.vertical_gradient').iterables = ('vertical_gradient',c.vertical_gradient) 

               except:
                   logConnectionError('Anatomical Preprocessing No valid Previous for strat',num_strat, strat.get_resource_pool(), '0001')
                   raise
 
               strat.append_name(anat_preproc.name)    
               strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
               # add stuff to resource pool if we need it

               strat.update_resource_pool({'anatomical_brain': (anat_preproc, 'outputspec.brain')})
               strat.update_resource_pool({'anatomical_reorient': (anat_preproc, 'outputspec.reorient')})
               # write to log
               create_log_node(anat_preproc, 'outputspec.brain', num_strat)  
               num_strat += 1 
    
    strat_list += new_strat_list

    '''
    Set Up FWHM iterable
    '''

    inputnode_fwhm = None
    if c.fwhm != None:
        inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                                 name='fwhm_input')
        inputnode_fwhm.iterables = ("fwhm", c.fwhm)

    '''
    T1 -> Template, Non-linear registration (FNIRT or ANTS)
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1

    # either run FSL anatomical-to-MNI registration, or...

    workflow_bit_id['anat_mni_register'] = workflow_counter
    for strat in strat_list:
          if 'FSL' in c.regOption:

            # this is to prevent the user from running FNIRT if they are
            # providing already-skullstripped inputs. this is because
            # FNIRT requires an input with the skull still on
            if already_skullstripped == 1:
                err_msg = '\n\n[!] CPAC says: FNIRT (for anatomical ' \
                          'registration) will not work properly if you ' \
                          'are providing inputs that have already been ' \
                          'skull-stripped.\n\nEither switch to using ' \
                          'ANTS for registration or provide input ' \
                          'images that have not been already ' \
                          'skull-stripped.\n\n'

                logger.info(err_msg)
                raise Exception

            fnirt_reg_anat_mni = create_nonlinear_register(
                'anat_mni_fnirt_register_%d' % num_strat)

            try:
                node, out_file = strat.get_node_from_resource_pool(
                    'anatomical_brain')
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.input_brain')

                node, out_file = strat.get_node_from_resource_pool(
                    'anatomical_reorient')
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.input_skull')

                # pass the reference files
                fnirt_reg_anat_mni.inputs.inputspec.reference_brain = c.template_brain_only_for_anat
                fnirt_reg_anat_mni.inputs.inputspec.reference_skull = c.template_skull_for_anat
                fnirt_reg_anat_mni.inputs.inputspec.ref_mask = c.ref_mask

                # assign the FSL FNIRT config file specified in pipeline
                # config.yml
                fnirt_reg_anat_mni.inputs.inputspec.fnirt_config = c.fnirtConfig


            except:
                logConnectionError('Anatomical Registration (FSL)', num_strat,
                                   strat.get_resource_pool(), '0002')
                raise

            if 'ANTS' in c.regOption:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(fnirt_reg_anat_mni.name)
            strat.set_leaf_properties(fnirt_reg_anat_mni,
                                      'outputspec.output_brain')

            strat.update_resource_pool({'anatomical_to_mni_linear_xfm': (fnirt_reg_anat_mni, 'outputspec.linear_xfm'),
                                        'anatomical_to_mni_nonlinear_xfm': (fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm'),
                                        'mni_to_anatomical_linear_xfm': (fnirt_reg_anat_mni, 'outputspec.invlinear_xfm'),
                                        'anatomical_to_standard': (fnirt_reg_anat_mni, 'outputspec.output_brain')})

            create_log_node(fnirt_reg_anat_mni, 'outputspec.output_brain',
                            num_strat)

            num_strat += 1

    strat_list += new_strat_list

    new_strat_list = []

    for strat in strat_list:

        nodes = getNodeList(strat)

        # or run ANTS anatomical-to-MNI registration instead
        if ('ANTS' in c.regOption) and \
                ('anat_mni_fnirt_register' not in nodes):

            ants_reg_anat_mni = \
                create_wf_calculate_ants_warp(
                    'anat_mni_ants_register_%d' % num_strat,
                    c.regWithSkull[0],
                    num_threads=num_ants_cores)

            try:
                # calculating the transform with the skullstripped is
                # reported to be better, but it requires very high
                # quality skullstripping. If skullstripping is imprecise
                # registration with skull is preferred
                if 1 in c.regWithSkull:

                    if already_skullstripped == 1:
                        err_msg = '\n\n[!] CPAC says: You selected ' \
                                  'to run anatomical registration with ' \
                                  'the skull, but you also selected to ' \
                                  'use already-skullstripped images as ' \
                                  'your inputs. This can be changed ' \
                                  'in your pipeline configuration ' \
                                  'editor.\n\n'

                        logger.info(err_msg)
                        raise Exception

                    # get the skull-stripped anatomical from resource pool
                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_brain')

                    # pass the anatomical to the workflow
                    workflow.connect(node, out_file,
                                     ants_reg_anat_mni,
                                     'inputspec.anatomical_brain')

                    # pass the reference file
                    ants_reg_anat_mni.inputs.inputspec.reference_brain = \
                        c.template_brain_only_for_anat

                    # get the reorient skull-on anatomical from resource pool
                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_reorient')

                    # pass the anatomical to the workflow
                    workflow.connect(node, out_file,
                                     ants_reg_anat_mni,
                                     'inputspec.anatomical_skull')

                    # pass the reference file
                    ants_reg_anat_mni.inputs.inputspec.reference_skull = \
                        c.template_skull_for_anat

                else:
                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_brain')

                    workflow.connect(node, out_file, ants_reg_anat_mni,
                                     'inputspec.anatomical_brain')

                    # pass the reference file
                    ants_reg_anat_mni.inputs.inputspec. \
                        reference_brain = c.template_brain_only_for_anat

                ants_reg_anat_mni.inputs.inputspec.dimension = 3
                ants_reg_anat_mni.inputs.inputspec. \
                    use_histogram_matching = True
                ants_reg_anat_mni.inputs.inputspec. \
                    winsorize_lower_quantile = 0.01
                ants_reg_anat_mni.inputs.inputspec. \
                    winsorize_upper_quantile = 0.99
                ants_reg_anat_mni.inputs.inputspec. \
                    metric = ['MI', 'MI', 'CC']
                ants_reg_anat_mni.inputs.inputspec.metric_weight = [1, 1, 1]
                ants_reg_anat_mni.inputs.inputspec. \
                    radius_or_number_of_bins = [32, 32, 4]
                ants_reg_anat_mni.inputs.inputspec. \
                    sampling_strategy = ['Regular', 'Regular', None]
                ants_reg_anat_mni.inputs.inputspec. \
                    sampling_percentage = [0.25, 0.25, None]
                ants_reg_anat_mni.inputs.inputspec. \
                    number_of_iterations = [[1000, 500, 250, 100],
                                            [1000, 500, 250, 100],
                                            [100, 100, 70, 20]]
                ants_reg_anat_mni.inputs.inputspec. \
                    convergence_threshold = [1e-8, 1e-8, 1e-9]
                ants_reg_anat_mni.inputs.inputspec. \
                    convergence_window_size = [10, 10, 15]
                ants_reg_anat_mni.inputs.inputspec. \
                    transforms = ['Rigid', 'Affine', 'SyN']
                ants_reg_anat_mni.inputs.inputspec. \
                    transform_parameters = [[0.1], [0.1], [0.1, 3, 0]]
                ants_reg_anat_mni.inputs.inputspec. \
                    shrink_factors = [[8, 4, 2, 1], [8, 4, 2, 1],
                                      [6, 4, 2, 1]]
                ants_reg_anat_mni.inputs.inputspec. \
                    smoothing_sigmas = [[3, 2, 1, 0], [3, 2, 1, 0],
                                        [3, 2, 1, 0]]

            except:
                logConnectionError('Anatomical Registration (ANTS)',
                                   num_strat, strat.get_resource_pool(),
                                   '0003')
                raise

            strat.append_name(ants_reg_anat_mni.name)
            strat.set_leaf_properties(ants_reg_anat_mni,
                                      'outputspec.normalized_output_brain')

            strat.update_resource_pool({'ants_initial_xfm': (ants_reg_anat_mni, 'outputspec.ants_initial_xfm'),
                                        'ants_rigid_xfm': (ants_reg_anat_mni, 'outputspec.ants_rigid_xfm'),
                                        'ants_affine_xfm': (ants_reg_anat_mni, 'outputspec.ants_affine_xfm'),
                                        'anatomical_to_mni_nonlinear_xfm': (ants_reg_anat_mni, 'outputspec.warp_field'),
                                        'mni_to_anatomical_nonlinear_xfm': (ants_reg_anat_mni, 'outputspec.inverse_warp_field'),
                                        'anat_to_mni_ants_composite_xfm': (ants_reg_anat_mni, 'outputspec.composite_transform'),
                                        'anatomical_to_standard': (ants_reg_anat_mni, 'outputspec.normalized_output_brain')})

            create_log_node(ants_reg_anat_mni,
                            'outputspec.normalized_output_brain', num_strat)

            num_strat += 1

    strat_list += new_strat_list

    '''
    [SYMMETRIC] T1 -> Symmetric Template, Non-linear registration (FNIRT/ANTS)
    '''
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1

    # either run FSL anatomical-to-MNI registration, or...

    if 1 in c.runVMHC:

        if not os.path.exists(c.template_symmetric_brain_only):
            logger.info("\n\n" + (
            "ERROR: Missing file - %s" % c.template_symmetric_brain_only) + "\n\n" + \
                        "Error name: cpac_pipeline_0017" + "\n\n")
            raise Exception

        if not os.path.exists(c.template_symmetric_skull):
            logger.info("\n\n" + (
            "ERROR: Missing file - %s" % c.template_symmetric_skull) + "\n\n" + \
                        "Error name: cpac_pipeline_0018" + "\n\n")
            raise Exception

        workflow_bit_id['anat_mni_symmetric_register'] = workflow_counter
        for strat in strat_list:

            nodes = getNodeList(strat)

            if 'FSL' in c.regOption and \
                    ('anat_mni_ants_register' not in nodes):

                # this is to prevent the user from running FNIRT if they are
                # providing already-skullstripped inputs. this is because
                # FNIRT requires an input with the skull still on
                if already_skullstripped == 1:
                    err_msg = '\n\n[!] CPAC says: FNIRT (for anatomical ' \
                              'registration) will not work properly if you ' \
                              'are providing inputs that have already been ' \
                              'skull-stripped.\n\nEither switch to using ' \
                              'ANTS for registration or provide input ' \
                              'images that have not been already ' \
                              'skull-stripped.\n\n'

                    logger.info(err_msg)
                    raise Exception

                fnirt_reg_anat_symm_mni = create_nonlinear_register(
                    'anat_symmetric_mni_fnirt_register_%d' % num_strat)

                try:
                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_brain')
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.input_brain')

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_reorient')
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.input_skull')

                    # pass the reference files
                    fnirt_reg_anat_symm_mni.inputs.inputspec.reference_brain = c.template_symmetric_brain_only
                    fnirt_reg_anat_symm_mni.inputs.inputspec.reference_skull = c.template_symmetric_skull
                    fnirt_reg_anat_symm_mni.inputs.inputspec.ref_mask = c.dilated_symmetric_brain_mask

                    # assign the FSL FNIRT config file specified in pipeline
                    # config.yml
                    fnirt_reg_anat_symm_mni.inputs.inputspec.fnirt_config = c.configFileTwomm

                    # node, out_file = strat.get_node_from_resource_pool('mni_normalized_anatomical')
                    # workflow.connect(node, out_file,
                    #                 fnirt_reg_anat_symm_mni, 'inputspec.wait')

                except:
                    logConnectionError(
                        'Symmetric Anatomical Registration (FSL)', num_strat,
                        strat.get_resource_pool(), '0002')
                    raise

                strat.append_name(fnirt_reg_anat_symm_mni.name)
                strat.set_leaf_properties(fnirt_reg_anat_symm_mni,
                                          'outputspec.output_brain')

                strat.update_resource_pool({'anatomical_to_symmetric_mni_linear_xfm': (fnirt_reg_anat_symm_mni, 'outputspec.linear_xfm'),
                                            'anatomical_to_symmetric_mni_nonlinear_xfm': (fnirt_reg_anat_symm_mni, 'outputspec.nonlinear_xfm'),
                                            'symmetric_mni_to_anatomical_linear_xfm': (fnirt_reg_anat_symm_mni, 'outputspec.invlinear_xfm'),
                                            'symmetric_anatomical_to_standard': (fnirt_reg_anat_symm_mni, 'outputspec.output_brain')})

                create_log_node(fnirt_reg_anat_symm_mni,
                                'outputspec.output_brain', num_strat)

                num_strat += 1

        strat_list += new_strat_list

        new_strat_list = []

        for strat in strat_list:

            nodes = getNodeList(strat)

            # or run ANTS anatomical-to-MNI registration instead
            if ('ANTS' in c.regOption) and \
                    ('anat_mni_fnirt_register' not in nodes) and \
                    ('anat_symmetric_mni_fnirt_register' not in nodes):

                ants_reg_anat_symm_mni = \
                    create_wf_calculate_ants_warp(
                        'anat_symmetric_mni_ants_register_%d' % num_strat,
                        c.regWithSkull[0],
                        num_threads=num_ants_cores)

                try:
                    # calculating the transform with the skullstripped is
                    # reported to be better, but it requires very high
                    # quality skullstripping. If skullstripping is imprecise
                    # registration with skull is preferred
                    if 1 in c.regWithSkull:

                        if already_skullstripped == 1:
                            err_msg = '\n\n[!] CPAC says: You selected ' \
                                      'to run anatomical registration with ' \
                                      'the skull, but you also selected to ' \
                                      'use already-skullstripped images as ' \
                                      'your inputs. This can be changed ' \
                                      'in your pipeline configuration ' \
                                      'editor.\n\n'

                            logger.info(err_msg)
                            raise Exception

                        # get the skullstripped anatomical from resource pool
                        node, out_file = strat.get_node_from_resource_pool(
                            'anatomical_brain')

                        # pass the anatomical to the workflow
                        workflow.connect(node, out_file,
                                         ants_reg_anat_symm_mni,
                                         'inputspec.anatomical_brain')

                        # pass the reference file
                        ants_reg_anat_symm_mni.inputs.inputspec.reference_brain = \
                            c.template_symmetric_brain_only

                        # get the reorient skull-on anatomical from resource
                        # pool
                        node, out_file = strat.get_node_from_resource_pool(
                            'anatomical_reorient')

                        # pass the anatomical to the workflow
                        workflow.connect(node, out_file,
                                         ants_reg_anat_symm_mni,
                                         'inputspec.anatomical_skull')

                        # pass the reference file
                        ants_reg_anat_symm_mni.inputs.inputspec.reference_skull = \
                            c.template_symmetric_skull

                    else:
                        # get the skullstripped anatomical from resource pool
                        node, out_file = strat.get_node_from_resource_pool(
                            'anatomical_brain')

                        workflow.connect(node, out_file,
                                         ants_reg_anat_symm_mni,
                                         'inputspec.anatomical_brain')

                        # pass the reference file
                        ants_reg_anat_symm_mni.inputs.inputspec. \
                            reference_brain = c.template_symmetric_brain_only

                    ants_reg_anat_symm_mni.inputs.inputspec.dimension = 3
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        use_histogram_matching = True
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        winsorize_lower_quantile = 0.01
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        winsorize_upper_quantile = 0.99
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        metric = ['MI', 'MI', 'CC']
                    ants_reg_anat_symm_mni.inputs.inputspec.metric_weight = [
                        1, 1, 1]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        radius_or_number_of_bins = [32, 32, 4]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        sampling_strategy = ['Regular', 'Regular', None]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        sampling_percentage = [0.25, 0.25, None]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        number_of_iterations = [[1000, 500, 250, 100], \
                                                [1000, 500, 250, 100],
                                                [100, 100, 70, 20]]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        convergence_threshold = [1e-8, 1e-8, 1e-9]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        convergence_window_size = [10, 10, 15]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        transforms = ['Rigid', 'Affine', 'SyN']
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        transform_parameters = [[0.1], [0.1], [0.1, 3, 0]]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        shrink_factors = [[8, 4, 2, 1], [8, 4, 2, 1],
                                          [6, 4, 2, 1]]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        smoothing_sigmas = [[3, 2, 1, 0], [3, 2, 1, 0],
                                            [3, 2, 1, 0]]

                except:
                    logConnectionError(
                        'Symmetric Anatomical Registration (ANTS)', num_strat,
                        strat.get_resource_pool(), '0003')
                    raise

                strat.append_name(ants_reg_anat_symm_mni.name)
                strat.set_leaf_properties(ants_reg_anat_symm_mni,
                                          'outputspec.normalized_output_brain')

                strat.update_resource_pool({'ants_symmetric_initial_xfm': (ants_reg_anat_symm_mni, 'outputspec.ants_initial_xfm'),
                                            'ants_symmetric_rigid_xfm': (ants_reg_anat_symm_mni, 'outputspec.ants_rigid_xfm'),
                                            'ants_symmetric_affine_xfm': (ants_reg_anat_symm_mni, 'outputspec.ants_affine_xfm'),
                                            'anatomical_to_symmetric_mni_nonlinear_xfm': (ants_reg_anat_symm_mni, 'outputspec.warp_field'),
                                            'symmetric_mni_to_anatomical_nonlinear_xfm': (ants_reg_anat_symm_mni, 'outputspec.inverse_warp_field'),
                                            'anat_to_symmetric_mni_ants_composite_xfm': (ants_reg_anat_symm_mni, 'outputspec.composite_transform'),
                                            'symmetric_anatomical_to_standard': (ants_reg_anat_symm_mni, 'outputspec.normalized_output_brain')})

                create_log_node(ants_reg_anat_symm_mni,
                                'outputspec.normalized_output_brain',
                                num_strat)

                num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting Segmentation Preprocessing
    Workflow
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runSegmentationPreprocessing:
        workflow_bit_id['seg_preproc'] = workflow_counter
        for strat in strat_list:

            nodes = getNodeList(strat)

            if 'anat_mni_fnirt_register' in nodes:
                seg_preproc = create_seg_preproc(False,
                                                 'seg_preproc_%d' % num_strat)
            elif 'anat_mni_ants_register' in nodes:
                seg_preproc = create_seg_preproc(True,
                                                 'seg_preproc_%d' % num_strat)

            try:
                node, out_file = strat.get_node_from_resource_pool(
                    'anatomical_brain')
                workflow.connect(node, out_file,
                                 seg_preproc, 'inputspec.brain')

                if 'anat_mni_fnirt_register' in nodes:
                    node, out_file = strat.get_node_from_resource_pool(
                        'mni_to_anatomical_linear_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc,
                                     'inputspec.standard2highres_mat')
                elif 'anat_mni_ants_register' in nodes:
                    node, out_file = strat.get_node_from_resource_pool(
                        'ants_initial_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc,
                                     'inputspec.standard2highres_init')
                    node, out_file = strat.get_node_from_resource_pool(
                        'ants_rigid_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc,
                                     'inputspec.standard2highres_rig')

                    node, out_file = strat.get_node_from_resource_pool(
                        'ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc,
                                     'inputspec.standard2highres_mat')

                seg_preproc.inputs.inputspec.PRIOR_CSF = c.PRIORS_CSF
                seg_preproc.inputs.inputspec.PRIOR_GRAY = c.PRIORS_GRAY
                seg_preproc.inputs.inputspec.PRIOR_WHITE = c.PRIORS_WHITE


            except:
                logConnectionError('Segmentation Preprocessing', num_strat,
                                   strat.get_resource_pool(), '0004')
                raise

            if 0 in c.runSegmentationPreprocessing:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(seg_preproc.name)
            strat.update_resource_pool(
                {'anatomical_gm_mask': (seg_preproc, 'outputspec.gm_mask'),
                 'anatomical_csf_mask': (seg_preproc, 'outputspec.csf_mask'),
                 'anatomical_wm_mask': (seg_preproc, 'outputspec.wm_mask'),
                 'seg_probability_maps': (
                 seg_preproc, 'outputspec.probability_maps'),
                 'seg_mixeltype': (seg_preproc, 'outputspec.mixeltype'),
                 'seg_partial_volume_map': (
                 seg_preproc, 'outputspec.partial_volume_map'),
                 'seg_partial_volume_files': (
                 seg_preproc, 'outputspec.partial_volume_files')})

            create_log_node(seg_preproc, 'outputspec.partial_volume_map',
                            num_strat)
            num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting Functional Data workflow
    '''

    num_strat = 0

    for strat in strat_list:
        # create a new node, Remember to change its name!
        # Flow = create_func_datasource(sub_dict['rest'])
        # Flow.inputs.inputnode.subject = subject_id

        # keep this in so that older participant lists that still have the
        # "rest" flag will still work
        try:
            # func_paths_dict is a dictionary of paths to the functional scans
            func_paths_dict = sub_dict['func']
        except KeyError:
            func_paths_dict = sub_dict['rest']

        try:
            funcFlow = create_func_datasource(func_paths_dict,
                                              'func_gather_%d' % num_strat)
            funcFlow.inputs.inputnode.subject = subject_id
            funcFlow.inputs.inputnode.creds_path = input_creds_path
            funcFlow.inputs.inputnode.dl_dir = c.workingDirectory
            funcFlow.get_node('inputnode').iterables = \
                ("scan", func_paths_dict.keys())
        except Exception as xxx:
            logger.info("Error create_func_datasource failed. "
                        "(%s:%d)" % dbg_file_lineno())
            raise

        """
        Add in nodes to get parameters from configuration file
        """

        scan_imports = ['import os', 'import json', 'import warnings',
                        'from CPAC.utils import check']

        try:
            # a node which checks if scan_parameters are present for each scan
            scan_params = \
                pe.Node(util.Function(input_names=['data_config_scan_params',
                                                   'subject_id',
                                                   'scan',
                                                   'pipeconfig_tr',
                                                   'pipeconfig_tpattern',
                                                   'pipeconfig_start_indx',
                                                   'pipeconfig_stop_indx'],
                                      output_names=['tr',
                                                    'tpattern',
                                                    'ref_slice',
                                                    'start_indx',
                                                    'stop_indx'],
                                      function=get_scan_params,
                                      imports=scan_imports),
                        name='scan_params_%d' % num_strat)
        except Exception as xxx:
            logger.info("Error creating scan_params node. (%s:%d)"
                        % dbg_file_lineno())
            raise

        if "Selected Functional Volume" in c.func_reg_input:

            try:
                get_func_volume = pe.Node(interface=preprocess.Calc(),
                                          name='get_func_volume_%d' % num_strat)

                get_func_volume.inputs.expr = 'a'
                get_func_volume.inputs.single_idx = c.func_reg_input_volume
                get_func_volume.inputs.outputtype = 'NIFTI_GZ'

            except Exception as xxx:
                logger.info("Error creating get_func_volume node." + \
                            " (%s:%d)" % dbg_file_lineno())
                raise

            try:
                workflow.connect(funcFlow, 'outputspec.rest',
                                 get_func_volume, 'in_file_a')

            except Exception as xxx:
                logger.info("Error connecting get_func_volume node." + \
                            " (%s:%d)" % dbg_file_lineno())
                raise

        # wire in the scan parameter workflow
        try:
            workflow.connect(funcFlow, 'outputspec.scan_params',
                             scan_params, 'data_config_scan_params')
        except Exception as xxx:
            logger.info("Error connecting scan_params 'data_config_"
                        "scan_params' input. (%s:%d)" % dbg_file_lineno())
            raise

        try:
            workflow.connect(funcFlow, 'outputspec.subject',
                             scan_params, 'subject_id')
        except Exception as xxx:
            logger.info("Error connecting scan_params 'subject_id' input." + \
                        " (%s:%d)" % dbg_file_lineno())
            raise

        try:
            workflow.connect(funcFlow, 'outputspec.scan',
                             scan_params, 'scan')
        except Exception as xxx:
            logger.info("Error connecting scan_params 'scan' input." + \
                        " (%s:%d)" % dbg_file_lineno())
            raise

        # connect in constants
        scan_params.inputs.pipeconfig_tr = c.TR
        scan_params.inputs.pipeconfig_tpattern = c.slice_timing_pattern
        scan_params.inputs.pipeconfig_start_indx = c.startIdx
        scan_params.inputs.pipeconfig_stop_indx = c.stopIdx

        # node to convert TR between seconds and milliseconds
        try:
            convert_tr = pe.Node(util.Function(input_names=['tr'],
                                               output_names=['tr'],
                                               function=get_tr),
                                 name='convert_tr_%d' % num_strat)
        except Exception as xxx:
            logger.info("Error creating convert_tr node." + \
                        " (%s:%d)" % dbg_file_lineno())
            raise

        try:
            workflow.connect(scan_params, 'tr',
                             convert_tr, 'tr')
        except Exception as xxx:
            logger.info("Error connecting convert_tr 'tr' input." + \
                        " (%s:%d)" % dbg_file_lineno())
            raise

        strat.set_leaf_properties(funcFlow, 'outputspec.rest')
        strat.update_resource_pool(
            {'raw_functional': (funcFlow, 'outputspec.rest')})

        if 1 in c.runEPI_DistCorr:
            try:
                strat.update_resource_pool(
                    {"fmap_phase_diff": (funcFlow, 'outputspec.phase_diff'),
                     "fmap_magnitude": (funcFlow, 'outputspec.magnitude')})
            except:
                err = "\n\n[!] You have selected to run field map " \
                      "distortion correction, but at least one of your " \
                      "scans listed in your data configuration file is " \
                      "missing either a field map phase difference file " \
                      "or a field map magnitude file, or both.\n\n"
                raise Exception(err)

        if "Selected Functional Volume" in c.func_reg_input:
            strat.update_resource_pool(
                {'selected_func_volume': (get_func_volume, 'out_file')})

        num_strat += 1

    """
    Truncate scan length based on configuration information
    """

    num_strat = 0

    for strat in strat_list:
        try:
            trunc_wf = create_wf_edit_func(
                wf_name="edit_func_%d" % (num_strat))
        except Exception as xxx:
            logger.info("Error create_wf_edit_func failed." + \
                        " (%s:%d)" % (dbg_file_lineno()))
            raise

        # find the output data on the leaf node
        try:
            node, out_file = strat.get_leaf_properties()
        except Exception as xxx:
            logger.info("Error  get_leaf_properties failed." + \
                        " (%s:%d)" % dbg_file_lineno())
            raise

        # connect the functional data from the leaf node into the wf
        try:
            workflow.connect(node, out_file, trunc_wf, 'inputspec.func')
        except Exception as xxx:
            logger.info("Error connecting input 'func' to trunc_wf." + \
                        " (%s:%d)" % dbg_file_lineno())
            print xxx
            raise

        # connect the other input parameters
        try:
            workflow.connect(scan_params, 'start_indx',
                             trunc_wf, 'inputspec.start_idx')
        except Exception as xxx:
            logger.info("Error connecting input 'start_indx' to trunc_wf." + \
                        " (%s:%d)" % dbg_file_lineno())
            raise

        try:
            workflow.connect(scan_params, 'stop_indx',
                             trunc_wf, 'inputspec.stop_idx')
        except Exception as xxx:
            logger.info("Error connecting input 'stop_idx' to trunc_wf." + \
                        " (%s:%d)" % dbg_file_lineno())
            raise

        # replace the leaf node with the output from the recently added
        # workflow
        strat.set_leaf_properties(trunc_wf, 'outputspec.edited_func')
        num_strat = num_strat + 1

    """
    EPI Field-Map based Distortion Correction
    """
 
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
   
    if 1 in c.runEPI_DistCorr:

        # if not fmap_phasediff:
        #     err = "\n\n[!] Field-map distortion correction is enabled, but " \
        #           "there is no field map phase difference file set in the " \
        #           "data configuration YAML file for participant {0}." \
        #           "\n\n".format(sub_dict['subject_id'])
        #     raise Exception(err)
        #
        # if not fmap_mag:
        #     err = "\n\n[!] Field-map distortion correction is enabled, but " \
        #           "there is no field map magnitude file set in the " \
        #           "data configuration YAML file for participant {0}." \
        #           "\n\n".format(sub_dict['subject_id'])
        #     raise Exception(err)

        workflow_bit_id['epi_distcorr'] = workflow_counter
    
        for strat in strat_list:

            if 'BET' in c.fmap_distcorr_skullstrip:
               epi_distcorr = create_EPI_DistCorr(use_BET = True, wf_name='epi_distcorr_%d' % (num_strat))
            else:
               epi_distcorr = create_EPI_DistCorr(use_BET = False, wf_name='epi_distcorr_%d' % (num_strat))

            epi_distcorr.inputs.bet_frac_input.bet_frac = c.fmap_distcorr_frac
            epi_distcorr.inputs.deltaTE_input.deltaTE = c.fmap_distcorr_deltaTE
            epi_distcorr.inputs.dwellT_input.dwellT = c.fmap_distcorr_dwell_time
            epi_distcorr.inputs.dwell_asym_ratio_input.dwell_asym_ratio = c.fmap_distcorr_dwell_asym_ratio

            epi_distcorr.get_node('bet_frac_input').iterables = ('bet_frac',c.fmap_distcorr_frac)
            epi_distcorr.get_node('deltaTE_input').iterables = ('deltaTE',
                                                   c.fmap_distcorr_deltaTE)
            epi_distcorr.get_node('dwellT_input').iterables = ('dwellT',
                                                   c.fmap_distcorr_dwell_time)
            epi_distcorr.get_node('dwell_asym_ratio_input').iterables = ('dwell_asym_ratio',c.fmap_distcorr_dwell_asym_ratio)

            try:
                node,out_file = strat.get_leaf_properties()
                workflow.connect(node,out_file,epi_distcorr,'inputspec.func_file')

                node,out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                workflow.connect(node,out_file,epi_distcorr,'inputspec.anat_file')

                node, out_file = strat.get_node_from_resource_pool('fmap_phase_diff')
                workflow.connect(node, out_file, epi_distcorr, 'inputspec.fmap_pha')

                node,out_file = strat.get_node_from_resource_pool('fmap_magnitude')
                workflow.connect(node,out_file,epi_distcorr, 'inputspec.fmap_mag')

            except:
                logConnectionError('EPI_DistCorr Workflow', num_strat,strat.get_resource_pool(), '0004')

            if 0 in c.runEPI_DistCorr:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file=str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(epi_distcorr.name)

            strat.update_resource_pool({'despiked_fieldmap':(epi_distcorr,'outputspec.fmap_despiked')})
            strat.update_resource_pool({'fieldmap_mask':(epi_distcorr,'outputspec.fieldmapmask')})
            strat.update_resource_pool({'prepared_fieldmap_map':(epi_distcorr,'outputspec.fieldmap')})
           
            num_strat += 1

    strat_list += new_strat_list

    """
    Inserting slice timing correction
    Workflow
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.slice_timing_correction:

        for strat in strat_list:

            # create TShift AFNI node
            try:
                func_slice_timing_correction = pe.Node(
                    interface=preprocess.TShift(),
                    name='func_slice_timing_correction_%d' % (num_strat))
                func_slice_timing_correction.inputs.outputtype = 'NIFTI_GZ'
            except Exception as xxx:
                logger.info("Error connecting input 'stop_idx' to trunc_wf. "
                            "(%s:%d)" % dbg_file_lineno())
                raise

            # find the output data on the leaf node
            try:
                node, out_file = strat.get_leaf_properties()
            except Exception as xxx:
                logger.info("Error  get_leaf_properties failed. "
                            "(%s:%d)" % dbg_file_lineno())
                raise

            # connect the output of the leaf node as the in_file
            try:
                workflow.connect(node, out_file,
                                 func_slice_timing_correction, 'in_file')
            except Exception as xxx:
                logger.info(
                    "Error connecting input 'infile' to func_slice_timing_"
                    "correction afni node. (%s:%d)" % dbg_file_lineno())
                raise

            logger.info("connected input to slc")
            # we might prefer to use the TR stored in the NIFTI header
            # if not, use the value in the scan_params node
            try:
                if c.TR:
                    if isinstance(c.TR, str):
                        if "None" in c.TR or "none" in c.TR:
                            pass
                        else:
                            workflow.connect(scan_params, 'tr',
                                             func_slice_timing_correction, 'tr')
                    else:
                        workflow.connect(scan_params, 'tr',
                                         func_slice_timing_correction, 'tr')
            except Exception as xxx:
                logger.info(
                    "Error connecting input 'tr' to func_slice_timing_"
                    "correction afni node. (%s:%d)" % dbg_file_lineno())
                print xxx
                raise
            logger.info("connected TR")

            # we might prefer to use the slice timing information stored in
            # the NIFTI header if not, use the value in the scan_params node
            logger.info("slice timing pattern %s" % c.slice_timing_pattern)
            try:
                if not "Use NIFTI Header" in c.slice_timing_pattern:
                    try:
                        logger.info("connecting slice timing pattern %s" %
                                    c.slice_timing_pattern)

                        # add the @ prefix to the tpattern file going into
                        # AFNI 3dTshift - needed this so the tpattern file
                        # output from get_scan_params would be tied downstream
                        # via a connection (to avoid poofing)
                        add_prefix = pe.Node(util.Function(input_names=['tpattern'],
                                                           output_names=['afni_prefix'],
                                                           function=add_afni_prefix),
                                             name='func_slice_timing_correction_add_afni_prefix_%d' % num_strat)
                        workflow.connect(scan_params, 'tpattern',
                                         add_prefix, 'tpattern')
                        workflow.connect(add_prefix, 'afni_prefix',
                                         func_slice_timing_correction,
                                         'tpattern')

                        logger.info("connected slice timing pattern %s" %
                                    c.slice_timing_patter)
                    except Exception as xxx:
                        logger.info(
                            "Error connecting input 'acquisition' to "
                            "func_slice_timing_correction afni node. "
                            "(%s:%d)" % dbg_file_lineno())
                        print xxx
                        raise
                    logger.info("connected slice timing pattern %s" %
                                c.slice_timing_pattern)
            except Exception as xxx:
                logger.info(
                    "Error connecting input 'acquisition' to "
                    "func_slice_timing_correction afni node. "
                    "(%s:%d)" % dbg_file_lineno())
                print xxx
                raise

            # add the name of the node to the strat name
            strat.append_name(func_slice_timing_correction.name)

            # set the leaf node
            strat.set_leaf_properties(func_slice_timing_correction,
                                      'out_file')

            # add the outputs to the resource pool
            strat.update_resource_pool({'slice_time_corrected':
                                            (func_slice_timing_correction, 'out_file')})
            num_strat += 1

    # add new strats (if forked)
    strat_list += new_strat_list

    logger.info(" finished connecting slice timing pattern")

    """
    Inserting Functional Image Preprocessing
    Workflow
    """

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1

    workflow_bit_id['func_preproc'] = workflow_counter

    for strat in strat_list:

        if '3dAutoMask' in c.functionalMasking:

            try:
                func_preproc = create_func_preproc(use_bet=False,
                                                   wf_name='func_preproc_automask_%d' % num_strat)
            except Exception as xxx:
                logger.info("Error allocating func_preproc." + \
                            " (%s:%d)" % dbg_file_lineno())
                raise

            try:
                node, out_file = strat.get_leaf_properties()
                try:
                    workflow.connect(node, out_file, func_preproc,
                                     'inputspec.func')
                except Exception as xxx:
                    logger.info(
                        "Error connecting leafnode to func, func_preproc." + \
                        " (%s:%d)" % (dbg_file_lineno()))
                    print xxx
                    raise
                logger.info("infile rest connected")
            except Exception as xxx:
                logConnectionError('Functional Preprocessing', num_strat,
                                   strat.get_resource_pool(), '0005_automask')
                num_strat += 1
                raise

            if 'BET' in c.functionalMasking:
                # we are forking so create a new node
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(func_preproc.name)

            strat.set_leaf_properties(func_preproc, 'outputspec.preprocessed')

            # add stuff to resource pool if we need it
            strat.update_resource_pool(
                {'mean_functional': (func_preproc, 'outputspec.example_func')})
            strat.update_resource_pool(
                {'functional_preprocessed_mask': (func_preproc, 'outputspec.preprocessed_mask')})
            strat.update_resource_pool(
                {'movement_parameters': (func_preproc, 'outputspec.movement_parameters')})
            strat.update_resource_pool(
                {'max_displacement': (func_preproc, 'outputspec.max_displacement')})
            strat.update_resource_pool(
                {'functional_preprocessed': (func_preproc, 'outputspec.preprocessed')})
            strat.update_resource_pool(
                {'functional_brain_mask': (func_preproc, 'outputspec.mask')})
            strat.update_resource_pool(
                {'motion_correct': (func_preproc, 'outputspec.motion_correct')})
            strat.update_resource_pool(
                {'coordinate_transformation': (func_preproc, 'outputspec.oned_matrix_save')})

            create_log_node(func_preproc, 'outputspec.preprocessed',
                            num_strat)
            num_strat += 1

    strat_list += new_strat_list

    new_strat_list = []

    for strat in strat_list:

        nodes = getNodeList(strat)

        if ('BET' in c.functionalMasking) and ('func_preproc_automask' not in nodes):

            func_preproc = create_func_preproc(use_bet=True,
                                               wf_name='func_preproc_bet_%d' % num_strat)
            node = None
            out_file = None
            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file, func_preproc,
                                 'inputspec.func')

            except Exception as xxx:
                logConnectionError('Functional Preprocessing', num_strat,
                                   strat.get_resource_pool(), '0005_bet')
                num_strat += 1
                raise

            strat.append_name(func_preproc.name)

            strat.set_leaf_properties(func_preproc, 'outputspec.preprocessed')

            # TODO redundant with above resource pool additions?

            # add stuff to resource pool if we need it
            strat.update_resource_pool({'mean_functional': (
            func_preproc, 'outputspec.example_func')})
            strat.update_resource_pool({'functional_preprocessed_mask': (
            func_preproc, 'outputspec.preprocessed_mask')})
            strat.update_resource_pool({'movement_parameters': (
            func_preproc, 'outputspec.movement_parameters')})
            strat.update_resource_pool({'max_displacement': (
            func_preproc, 'outputspec.max_displacement')})
            strat.update_resource_pool(
                {'functional_preprocessed': (func_preproc, 'outputspec.preprocessed')})
            strat.update_resource_pool(
                {'functional_brain_mask': (func_preproc, 'outputspec.mask')})
            strat.update_resource_pool({'motion_correct': (
            func_preproc, 'outputspec.motion_correct')})
            strat.update_resource_pool({'coordinate_transformation': (
            func_preproc, 'outputspec.oned_matrix_save')})

            create_log_node(func_preproc, 'outputspec.preprocessed',
                            num_strat)
            num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting Friston's 24 parameter Workflow
    In case this workflow runs , it overwrites the movement_parameters file
    So the file contains 24 parameters for motion and that gets wired to all the workflows
    that depend on. The effect should be seen when regressing out nuisance signals and motion
    is used as one of the regressors
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runFristonModel:
        workflow_bit_id['fristons_parameter_model'] = workflow_counter

        for strat in strat_list:

            fristons_model = fristons_twenty_four(
                wf_name='fristons_parameter_model_%d' % num_strat)

            try:

                node, out_file = strat.get_node_from_resource_pool(
                    'movement_parameters')
                workflow.connect(node, out_file,
                                 fristons_model, 'inputspec.movement_file')

            except Exception as xxx:
                logConnectionError('Friston\'s Parameter Model', num_strat,
                                   strat.get_resource_pool(), '0006')
                raise

            if 0 in c.runFristonModel:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)
            strat.append_name(fristons_model.name)

            strat.update_resource_pool({'movement_parameters': (fristons_model, 'outputspec.movement_file')})

            create_log_node(fristons_model, 'outputspec.movement_file',
                            num_strat)

            num_strat += 1

    strat_list += new_strat_list

    '''
    Func -> T1 Registration (Initial Linear reg)
    '''

    # Depending on configuration, either passes output matrix to
    # Func -> Template ApplyWarp, or feeds into linear reg of BBReg operation
    # (if BBReg is enabled)

    new_strat_list = []
    num_strat = 0
    workflow_counter += 1

    if 1 in c.runRegisterFuncToAnat:

        workflow_bit_id['func_to_anat'] = workflow_counter

        for strat in strat_list:

            nodes = getNodeList(strat)

            # if field map-based distortion correction is on, but BBR is off,
            # send in the distortion correction files here
            # TODO: is this robust to the possibility of forking both
            # TODO: distortion correction and BBR at the same time?
            # TODO: (note if you are forking with BBR on/off, at this point
            # TODO:  there is still only one strat, so you would have to fork
            # TODO:  here instead to have a func->anat with fieldmap and
            # TODO:  without, and send the without-fieldmap to the BBR fork)
            dist_corr = False
            if 'epi_distcorr' in nodes and 1 not in c.runBBReg:
                dist_corr = True
                # TODO: for now, disabling dist corr when BBR is disabled
                err = "\n\n[!] Field map distortion correction is enabled, " \
                      "but Boundary-Based Registration is off- BBR is " \
                      "required for distortion correction.\n\n"
                raise Exception(err)

            func_to_anat = create_register_func_to_anat(dist_corr,
                                                        'func_to_anat_FLIRT'
                                                        '_%d' % num_strat)

            # Input registration parameters
            func_to_anat.inputs.inputspec.interp = 'trilinear'

            try:
                def pick_wm(seg_prob_list):
                    seg_prob_list.sort()
                    return seg_prob_list[-1]

                if 'Mean Functional' in c.func_reg_input:
                    # Input functional image (mean functional)
                    node, out_file = strat.get_node_from_resource_pool(
                        'mean_functional')
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.func')

                elif 'Selected Functional Volume' in c.func_reg_input:
                    # Input functional image (specific volume)
                    node, out_file = strat.get_node_from_resource_pool(
                        'selected_func_volume')
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.func')

                # Input skull-stripped anatomical (anat.nii.gz)
                node, out_file = strat.get_node_from_resource_pool(
                    'anatomical_brain')
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.anat')

                if dist_corr:
                    # apply field map distortion correction outputs to
                    # the func->anat registration

                    func_to_anat.inputs.echospacing_input.echospacing = c.fmap_distcorr_dwell_time[0]
                    func_to_anat.inputs.pedir_input.pedir = c.fmap_distcorr_pedir

                    node, out_file = strat.get_node_from_resource_pool("despiked_fieldmap")
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.fieldmap')

                    node, out_file = strat.get_node_from_resource_pool("fieldmap_mask")
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.fieldmapmask')

            except:
                logConnectionError(
                    'Register Functional to Anatomical (pre BBReg)',
                    num_strat, strat.get_resource_pool(), '0007')
                raise

            if 0 in c.runRegisterFuncToAnat:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(func_to_anat.name)
            # strat.set_leaf_properties(func_mni_warp, 'out_file')

            strat.update_resource_pool({'mean_functional_in_anat': (func_to_anat, 'outputspec.anat_func_nobbreg'),
                                        'functional_to_anat_linear_xfm': (func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg')})

            # Outputs:
            # functional_to_anat_linear_xfm = func-t1.mat, linear, sent to 'premat' of post-FNIRT applywarp,
            #                                 or to the input of the post-ANTS c3d_affine_tool

            # create_log_node(func_to_anat, 'outputspec.mni_func', num_strat)
            num_strat += 1

    strat_list += new_strat_list

    '''
    Func -> T1 Registration (BBREG)
    '''

    # Outputs 'functional_to_anat_linear_xfm', a matrix file of the
    # functional-to-anatomical registration warp to be applied LATER in
    # func_mni_warp, which accepts it as input 'premat'

    new_strat_list = []
    num_strat = 0
    workflow_counter += 1

    if (1 in c.runRegisterFuncToAnat) and (1 in c.runBBReg):

        workflow_bit_id['func_to_anat_bbreg'] = workflow_counter

        for strat in strat_list:

            nodes = getNodeList(strat)

            # this is needed here in case tissue segmentation is set on/off
            # and you have bbreg enabled- this will ensure bbreg will run for
            # the strat that has segmentation but will not run (thus avoiding
            # a crash) on the strat without segmentation
            if 'seg_preproc' in nodes:

                dist_corr = False
                if 'epi_distcorr' in nodes:
                    dist_corr = True

                func_to_anat_bbreg = create_bbregister_func_to_anat(dist_corr,
                                                                    'func_to_anat_bbreg_%d' % num_strat)

                # Input registration parameters
                func_to_anat_bbreg.inputs.inputspec.bbr_schedule = c.boundaryBasedRegistrationSchedule

                try:
                    def pick_wm(seg_prob_list):
                        seg_prob_list.sort()
                        return seg_prob_list[-1]

                    if 'Mean Functional' in c.func_reg_input:
                        # Input functional image (mean functional)
                        node, out_file = strat.get_node_from_resource_pool(
                            'mean_functional')
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg, 'inputspec.func')

                    elif 'Selected Functional Volume' in c.func_reg_input:
                        # Input functional image (specific volume)
                        node, out_file = strat.get_node_from_resource_pool(
                            'selected_func_volume')
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg, 'inputspec.func')

                    # Input anatomical whole-head image (reoriented)
                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_reorient')
                    workflow.connect(node, out_file,
                                     func_to_anat_bbreg,
                                     'inputspec.anat_skull')

                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     func_to_anat_bbreg,
                                     'inputspec.linear_reg_matrix')

                    # Input segmentation probability maps for white matter
                    # segmentation
                    node, out_file = strat.get_node_from_resource_pool(
                        'seg_probability_maps')
                    workflow.connect(node, (out_file, pick_wm),
                                     func_to_anat_bbreg,
                                     'inputspec.anat_wm_segmentation')

                    if dist_corr:
                        # apply field map distortion correction outputs to
                        # the func->anat registration

                        func_to_anat_bbreg.inputs.echospacing_input.echospacing = c.fmap_distcorr_dwell_time[0]
                        func_to_anat_bbreg.inputs.pedir_input.pedir = c.fmap_distcorr_pedir

                        node, out_file = strat.get_node_from_resource_pool(
                            "despiked_fieldmap")
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg,
                                         'inputspec.fieldmap')

                        node, out_file = strat.get_node_from_resource_pool(
                            "fieldmap_mask")
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg,
                                         'inputspec.fieldmapmask')

                except:
                    logConnectionError(
                        'Register Functional to Anatomical (BBReg)',
                        num_strat, strat.get_resource_pool(), '0008')
                    raise

                if 0 in c.runBBReg:
                    tmp = strategy()
                    tmp.resource_pool = dict(strat.resource_pool)
                    tmp.leaf_node = (strat.leaf_node)
                    tmp.out_file = str(strat.leaf_out_file)

                    # This line is needed here for some reason, otherwise the
                    # connection between func_preproc and nuisance will
                    # break - even though this workflow has nothing to do with
                    # it - but excluding this line below removes the leaf node
                    # from the new forked strat
                    tmp.leaf_out_file = str(strat.leaf_out_file)

                    tmp.name = list(strat.name)
                    strat = tmp
                    new_strat_list.append(strat)

                strat.append_name(func_to_anat_bbreg.name)

                strat.update_resource_pool({'mean_functional_in_anat': (func_to_anat_bbreg, 'outputspec.anat_func'),
                                            'functional_to_anat_linear_xfm': (func_to_anat_bbreg, 'outputspec.func_to_anat_linear_xfm')})

                # Outputs:
                # functional_to_anat_linear_xfm = func-t1.mat, linear, sent to 'premat' of post-FNIRT applywarp,
                #                                 or to the input of the post-ANTS c3d_affine_tool

                # create_log_node(func_to_anat, 'outputspec.mni_func', num_strat)
                num_strat += 1

            else:
                # anatomical segmentation is not being run in this particular
                # strategy/fork - we don't want this to stop workflow building
                # unless there is only one strategy
                if len(strat_list) > 1:
                    pass
                else:
                    err = "\n\n[!] Boundary-based registration (BBR) for " \
                          "functional-to-anatomical registration is " \
                          "enabled, but anatomical segmentation is not. " \
                          "BBR requires the outputs of segmentation. " \
                          "Please modify your pipeline configuration and " \
                          "run again.\n\n"
                    raise Exception(err)

    strat_list += new_strat_list

    '''
    Inserting Generate Motion Statistics Workflow
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1

    workflow_bit_id['gen_motion_stats'] = workflow_counter
    for strat in strat_list:

        gen_motion_stats = motion_power_statistics(c.fdCalc[0],
                                                   'gen_motion_stats_%d'
                                                   % num_strat)
        gen_motion_stats.inputs.scrubbing_input.threshold = c.spikeThreshold
        gen_motion_stats.inputs.scrubbing_input.remove_frames_before = c.numRemovePrecedingFrames
        gen_motion_stats.inputs.scrubbing_input.remove_frames_after = c.numRemoveSubsequentFrames
        gen_motion_stats.get_node('scrubbing_input').iterables = ('threshold',
                                                                  c.spikeThreshold)

        try:
            # **special case where the workflow is not getting outputs from
            # resource pool but is connected to functional datasource
            workflow.connect(funcFlow, 'outputspec.subject',
                             gen_motion_stats, 'inputspec.subject_id')

            workflow.connect(funcFlow, 'outputspec.scan',
                             gen_motion_stats, 'inputspec.scan_id')

            node, out_file = strat.get_node_from_resource_pool(
                'motion_correct')
            workflow.connect(node, out_file,
                             gen_motion_stats, 'inputspec.motion_correct')

            node, out_file = strat.get_node_from_resource_pool(
                'movement_parameters')
            workflow.connect(node, out_file,
                             gen_motion_stats,
                             'inputspec.movement_parameters')

            node, out_file = strat.get_node_from_resource_pool(
                'max_displacement')
            workflow.connect(node, out_file,
                             gen_motion_stats, 'inputspec.max_displacement')

            node, out_file = strat.get_node_from_resource_pool(
                'functional_brain_mask')
            workflow.connect(node, out_file,
                             gen_motion_stats, 'inputspec.mask')

            node, out_file = strat.get_node_from_resource_pool(
                'coordinate_transformation')
            workflow.connect(node, out_file,
                             gen_motion_stats, 'inputspec.oned_matrix_save')

        except:
            logConnectionError('Generate Motion Statistics', num_strat,
                               strat.get_resource_pool(), '0009')
            raise

        strat.append_name(gen_motion_stats.name)

        strat.update_resource_pool({'frame_wise_displacement_power': (
                                        gen_motion_stats, 'outputspec.FDP_1D'),
                                    'frame_wise_displacement_jenkinson': (
                                        gen_motion_stats, 'outputspec.FDJ_1D'),
                                    'power_params': (gen_motion_stats,
                                                     'outputspec.power_params'),
                                    'motion_params': (gen_motion_stats,
                                                      'outputspec.motion_params')})

        if "De-Spiking" in c.runMotionSpike and 1 in c.runNuisance:
            strat.update_resource_pool({'despiking_frames_excluded': (
                                            gen_motion_stats, 'outputspec.frames_ex_1D'),
                                        'despiking_frames_included': (
                                            gen_motion_stats, 'outputspec.frames_in_1D')})

            '''
            if "Scrubbing" in c.runMotionSpike:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)
            '''

        if "Scrubbing" in c.runMotionSpike and 1 in c.runNuisance:
            strat.update_resource_pool({'scrubbing_frames_excluded': (
                                            gen_motion_stats, 'outputspec.frames_ex_1D'),
                                        'scrubbing_frames_included': (
                                            gen_motion_stats, 'outputspec.frames_in_1D')})

            '''
            if "De-Spiking" in c.runMotionSpike:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)
            '''

        create_log_node(gen_motion_stats, 'outputspec.motion_params',
                        num_strat)
        num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting Nuisance Workflow
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1

    if 1 in c.runNuisance:

        workflow_bit_id['nuisance'] = workflow_counter

        for strat in strat_list:

            nodes = getNodeList(strat)

            # this is needed here in case tissue segmentation is set on/off
            # and you have nuisance enabled- this will ensure nuisance will
            # run for the strat that has segmentation but will not run (thus
            # avoiding a crash) on the strat without segmentation
            if 'seg_preproc' in nodes:

                subwf_name = "nuisance"
                if "De-Spiking" in c.runMotionSpike:
                    subwf_name = "nuisance_with_despiking"

                if 'anat_mni_fnirt_register' in nodes:
                    nuisance = create_nuisance(False,
                                               '{0}_{1}'.format(subwf_name,
                                                                num_strat))
                else:
                    nuisance = create_nuisance(True,
                                               '{0}_{1}'.format(subwf_name,
                                                                num_strat))

                nuisance.get_node('residuals').iterables = (
                    [('selector', c.Regressors),
                     ('compcor_ncomponents', c.nComponents)])

                nuisance.inputs.inputspec.lat_ventricles_mask = c.lateral_ventricles_mask

                try:
                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.subject')

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_gm_mask')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.gm_mask')

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_wm_mask')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.wm_mask')

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_csf_mask')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.csf_mask')

                    node, out_file = strat.get_node_from_resource_pool(
                        'movement_parameters')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.motion_components')

                    if "De-Spiking" in c.runMotionSpike:
                        node, out_file = strat.get_node_from_resource_pool(
                            'despiking_frames_excluded')
                        workflow.connect(node, out_file,
                                         nuisance, 'inputspec.frames_ex')
                    else:
                        nuisance.inputs.inputspec.frames_ex = None

                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     nuisance,
                                     'inputspec.func_to_anat_linear_xfm')

                    if 'anat_mni_fnirt_register' in nodes:
                        node, out_file = strat.get_node_from_resource_pool(
                            'mni_to_anatomical_linear_xfm')
                        workflow.connect(node, out_file,
                                         nuisance,
                                         'inputspec.mni_to_anat_linear_xfm')
                    else:
                        # pass the ants_affine_xfm to the input for the
                        # INVERSE transform, but ants_affine_xfm gets inverted
                        # within the workflow

                        node, out_file = strat.get_node_from_resource_pool(
                            'ants_initial_xfm')
                        workflow.connect(node, out_file,
                                         nuisance,
                                         'inputspec.anat_to_mni_initial_xfm')

                        node, out_file = strat.get_node_from_resource_pool(
                            'ants_rigid_xfm')
                        workflow.connect(node, out_file,
                                         nuisance,
                                         'inputspec.anat_to_mni_rigid_xfm')

                        node, out_file = strat.get_node_from_resource_pool(
                            'ants_affine_xfm')
                        workflow.connect(node, out_file,
                                         nuisance,
                                         'inputspec.anat_to_mni_affine_xfm')

                except Exception as e:
                    logConnectionError('Nuisance', num_strat,
                                       strat.get_resource_pool(), '0010', e)
                    raise

                if 0 in c.runNuisance:
                    new_strat_list.append(create_new_fork(strat))

                if 1 in c.runNuisance and "De-Spiking" in c.runMotionSpike and \
                        "Scrubbing" in c.runMotionSpike:
                    # create a new fork that will run nuisance like above but
                    # without the de-spiking
                    new_strat_list.append(create_new_fork(strat))

                if 1 in c.runNuisance and "De-Spiking" in c.runMotionSpike and \
                        "None" in c.runMotionSpike:
                    # create a new fork that will run nuisance like above but
                    # without the de-spiking
                    new_strat_list.append(create_new_fork(strat))

                strat.append_name(nuisance.name)

                strat.set_leaf_properties(nuisance, 'outputspec.subject')

                strat.update_resource_pool({'functional_nuisance_residuals': (nuisance, 'outputspec.subject')})
                strat.update_resource_pool({'functional_nuisance_regressors': (nuisance, 'outputspec.regressors')})

                create_log_node(nuisance, 'outputspec.subject', num_strat)

                num_strat += 1

    strat_list += new_strat_list

    # set a flag in case we're doing nuisance on/off
    non_nuisance_strat = False

    for strat in strat_list:

        nodes = getNodeList(strat)

        if 0 in c.runNuisance and \
                ("nuisance" not in nodes and "nuisance_with_despiking" not in nodes):
            if not non_nuisance_strat:
                # save one of the strats so that it won't have any nuisance
                # at all - this only fires if nuisance is on/off
                non_nuisance_strat = True
                continue

        if 1 in c.runNuisance and "De-Spiking" in c.runMotionSpike and \
                "nuisance_with_despiking" not in nodes and \
                ("Scrubbing" in c.runMotionSpike or "None" in c.runMotionSpike):
            # run nuisance in the new fork (if created), without de-spiking,
            # so that we can have nuisance and then scrubbing, or a nuisance
            # strat without de-spiking if doing de-spiking on/off
            #     this only runs if we have ["De-Spiking", "Scrubbing"] or
            #     ["De-Spiking", "Scrubbing", "Off"] in c.runMotionSpike

            # this is needed here in case tissue segmentation is set on/off
            # and you have nuisance enabled- this will ensure nuisance will
            # run for the strat that has segmentation but will not run (thus
            # avoiding a crash) on the strat without segmentation
            if 'seg_preproc' in nodes:

                if 'anat_mni_fnirt_register' in nodes:
                    nuisance = create_nuisance(False,
                                               'nuisance_no_despiking_%d' % num_strat)
                else:
                    nuisance = create_nuisance(True,
                                               'nuisance_no_despiking_%d' % num_strat)

                nuisance.get_node('residuals').iterables = (
                    [('selector', c.Regressors),
                     ('compcor_ncomponents', c.nComponents)])

                nuisance.inputs.inputspec.lat_ventricles_mask = c.lateral_ventricles_mask

                try:
                    # enforcing no de-spiking here!
                    # TODO: when condensing these sub-wf builders, pass
                    # TODO: something so that the check in the nuisance strat
                    # TODO: above can be modified for this version down here
                    nuisance.inputs.inputspec.frames_ex = None

                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.subject')

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_gm_mask')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.gm_mask')

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_wm_mask')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.wm_mask')

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_csf_mask')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.csf_mask')

                    node, out_file = strat.get_node_from_resource_pool(
                        'movement_parameters')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.motion_components')

                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     nuisance,
                                     'inputspec.func_to_anat_linear_xfm')

                    if 'anat_mni_fnirt_register' in nodes:
                        node, out_file = strat.get_node_from_resource_pool(
                            'mni_to_anatomical_linear_xfm')
                        workflow.connect(node, out_file,
                                         nuisance,
                                         'inputspec.mni_to_anat_linear_xfm')
                    else:
                        # pass the ants_affine_xfm to the input for the
                        # INVERSE transform, but ants_affine_xfm gets inverted
                        # within the workflow

                        node, out_file = strat.get_node_from_resource_pool(
                            'ants_initial_xfm')
                        workflow.connect(node, out_file,
                                         nuisance,
                                         'inputspec.anat_to_mni_initial_xfm')

                        node, out_file = strat.get_node_from_resource_pool(
                            'ants_rigid_xfm')
                        workflow.connect(node, out_file,
                                         nuisance,
                                         'inputspec.anat_to_mni_rigid_xfm')

                        node, out_file = strat.get_node_from_resource_pool(
                            'ants_affine_xfm')
                        workflow.connect(node, out_file,
                                         nuisance,
                                         'inputspec.anat_to_mni_affine_xfm')

                except Exception as e:
                    logConnectionError('Nuisance', num_strat,
                                       strat.get_resource_pool(), '0010b', e)
                    raise

                strat.append_name(nuisance.name)

                strat.set_leaf_properties(nuisance, 'outputspec.subject')

                strat.update_resource_pool(
                    {'functional_nuisance_residuals': (nuisance, 'outputspec.subject')})
                strat.update_resource_pool(
                    {'functional_nuisance_regressors': (nuisance, 'outputspec.regressors')})

                create_log_node(nuisance, 'outputspec.subject', num_strat)

                num_strat += 1

    '''
    Inserting Median Angle Correction Workflow
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runMedianAngleCorrection:
        workflow_bit_id['median_angle_corr'] = workflow_counter
        for strat in strat_list:
            median_angle_corr = create_median_angle_correction(
                'median_angle_corr_%d' % num_strat)

            median_angle_corr.get_node('median_angle_correct').iterables = (
            'target_angle_deg', c.targetAngleDeg)
            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 median_angle_corr, 'inputspec.subject')
            except:
                logConnectionError('Median Angle Correction', num_strat,
                                   strat.get_resource_pool(), '0011')
                raise

            if 0 in c.runMedianAngleCorrection:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(median_angle_corr.name)

            strat.set_leaf_properties(median_angle_corr, 'outputspec.subject')

            strat.update_resource_pool({'functional_median_angle_corrected': (median_angle_corr, 'outputspec.subject')})

            create_log_node(median_angle_corr, 'outputspec.subject',
                            num_strat)

            num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting ALFF/fALFF
    Workflow
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runALFF:
        for strat in strat_list:
            alff = create_alff('alff_falff_%d' % num_strat)
            alff.inputs.hp_input.hp = c.highPassFreqALFF
            alff.inputs.lp_input.lp = c.lowPassFreqALFF
            alff.get_node('hp_input').iterables = ('hp',
                                                   c.highPassFreqALFF)
            alff.get_node('lp_input').iterables = ('lp',
                                                   c.lowPassFreqALFF)

            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 alff, 'inputspec.rest_res')
                node, out_file = strat.get_node_from_resource_pool(
                    'functional_brain_mask')
                workflow.connect(node, out_file,
                                 alff, 'inputspec.rest_mask')

            except:
                logConnectionError('ALFF', num_strat,
                                   strat.get_resource_pool(), '0012')
                raise

            strat.append_name(alff.name)

            strat.update_resource_pool({'alff': (alff, 'outputspec.alff_img')})
            strat.update_resource_pool({'falff': (alff, 'outputspec.falff_img')})

            create_log_node(alff, 'outputspec.falff_img', num_strat)

            num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting Frequency Filtering Node
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runFrequencyFiltering:
        workflow_bit_id['frequency_filter'] = workflow_counter
        filter_imports = ['import os', 'import nibabel as nb',
                          'import numpy as np',
                          'from scipy.fftpack import fft, ifft']
        for strat in strat_list:
            frequency_filter = pe.Node(
                util.Function(input_names=['realigned_file',
                                           'bandpass_freqs',
                                           'sample_period'],
                              output_names=['bandpassed_file'],
                              function=bandpass_voxels,
                              imports=filter_imports),
                name='frequency_filter_%d' % num_strat)

            frequency_filter.iterables = ('bandpass_freqs', c.nuisanceBandpassFreq)
            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 frequency_filter, 'realigned_file')

            except:
                logConnectionError('Frequency Filtering', num_strat,
                                   strat.get_resource_pool(), '0013')
                raise

            if 0 in c.runFrequencyFiltering:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(frequency_filter.name)

            strat.set_leaf_properties(frequency_filter, 'bandpassed_file')

            strat.update_resource_pool({'functional_freq_filtered': (frequency_filter, 'bandpassed_file')})

            create_log_node(frequency_filter, 'bandpassed_file', num_strat)

            num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting Scrubbing Workflow
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1

    if "Scrubbing" in c.runMotionSpike and 1 in c.runNuisance:

        workflow_bit_id['scrubbing'] = workflow_counter

        # set a flag in case we're doing nuisance on/off
        non_nuisance_strat = False

        for strat in strat_list:

            nodes = getNodeList(strat)

            if 0 in c.runNuisance and \
                    "nuisance" not in nodes and \
                        "nuisance_with_despiking" not in nodes and \
                            "nuisance_no_despiking" not in nodes:
                if not non_nuisance_strat:
                    # save one of the strats so that it won't have any
                    # nuisance at all - this only fires if nuisance is on/off
                    non_nuisance_strat = True
                    continue

            # skip if this strat had de-spiking (mutually exclusive)
            if "nuisance_with_despiking" in nodes:
                continue

            if 'gen_motion_stats' in nodes:
                scrubbing = \
                    create_scrubbing_preproc('scrubbing_%d' % num_strat)

                try:
                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file,
                                     scrubbing, 'inputspec.preprocessed')

                    node, out_file = strat.get_node_from_resource_pool(
                        'scrubbing_frames_included')
                    workflow.connect(node, out_file,
                                     scrubbing, 'inputspec.frames_in_1D')

                    node, out_file = strat.get_node_from_resource_pool(
                        'movement_parameters')
                    workflow.connect(node, out_file,
                                     scrubbing,
                                     'inputspec.movement_parameters')

                except:
                    logConnectionError('Scrubbing Workflow', num_strat,
                                       strat.get_resource_pool(), '0014')
                    raise

                if "None" in c.runMotionSpike:
                    new_strat_list.append(create_new_fork(strat))

                strat.append_name(scrubbing.name)

                strat.set_leaf_properties(scrubbing,
                                          'outputspec.preprocessed')

                strat.update_resource_pool({'scrubbing_movement_parameters': (
                                                scrubbing, 'outputspec.scrubbed_movement_parameters'),
                                            'scrubbed_preprocessed': (
                                                scrubbing, 'outputspec.preprocessed')})

                create_log_node(scrubbing, 'outputspec.preprocessed',
                                num_strat)

                num_strat += 1

    strat_list += new_strat_list

    '''
    Func -> Template, uses antsApplyTransforms (ANTS) or ApplyWarp (FSL) to
    apply the warp; also includes mean functional warp
    '''

    new_strat_list = []
    num_strat = 0
    if 1 in c.runRegisterFuncToMNI:

        for strat in strat_list:

            nodes = getNodeList(strat)

            # Run FSL ApplyWarp
            if 'anat_mni_fnirt_register' in nodes:

                func_mni_warp = pe.Node(interface=fsl.ApplyWarp(),
                                        name='func_mni_fsl_warp_%d' % num_strat)
                func_mni_warp.inputs.ref_file = c.template_brain_only_for_func

                functional_brain_mask_to_standard = pe.Node(
                    interface=fsl.ApplyWarp(),
                    name='func_mni_fsl_warp_mask_%d' % num_strat)
                functional_brain_mask_to_standard.inputs.interp = 'nn'
                functional_brain_mask_to_standard.inputs.ref_file = c.template_skull_for_func

                mean_functional_warp = pe.Node(interface=fsl.ApplyWarp(),
                                               name='mean_func_fsl_warp_%d' % num_strat)
                mean_functional_warp.inputs.ref_file = c.template_brain_only_for_func

                motion_correct_warp = pe.Node(interface=fsl.ApplyWarp(),
                                              name="motion_correct_fsl_warp_%d" % num_strat)
                motion_correct_warp.inputs.ref_file = c.template_brain_only_for_func

                try:

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     func_mni_warp, 'field_file')

                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     func_mni_warp, 'premat')

                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file,
                                     func_mni_warp, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     functional_brain_mask_to_standard,
                                     'field_file')
                    workflow.connect(node, out_file,
                                     mean_functional_warp, 'field_file')
                    workflow.connect(node, out_file,
                                     motion_correct_warp, 'field_file')

                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     functional_brain_mask_to_standard,
                                     'premat')
                    workflow.connect(node, out_file,
                                     mean_functional_warp, 'premat')
                    workflow.connect(node, out_file,
                                     motion_correct_warp, 'premat')

                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_brain_mask')
                    workflow.connect(node, out_file,
                                     functional_brain_mask_to_standard,
                                     'in_file')

                    node, out_file = strat.get_node_from_resource_pool(
                        'mean_functional')
                    workflow.connect(node, out_file, mean_functional_warp,
                                     'in_file')

                    node, out_file = strat.get_node_from_resource_pool(
                        'motion_correct')
                    workflow.connect(node, out_file, motion_correct_warp,
                                     'in_file')

                except:
                    logConnectionError(
                        'Functional Timeseries Registration to MNI space (FSL)',
                        num_strat, strat.get_resource_pool(), '0015')
                    raise

                strat.update_resource_pool(
                    {'functional_to_standard': (func_mni_warp, 'out_file'),
                     'functional_brain_mask_to_standard': (
                     functional_brain_mask_to_standard, 'out_file'),
                     'mean_functional_to_standard': (
                     mean_functional_warp, 'out_file'),
                     'motion_correct_to_standard': (
                     motion_correct_warp, 'out_file')})

                strat.append_name(func_mni_warp.name)
                create_log_node(func_mni_warp, 'out_file', num_strat)

                num_strat += 1

        strat_list += new_strat_list

        for strat in strat_list:

            nodes = getNodeList(strat)

            if 'ANTS' in c.regOption and \
                    'anat_mni_fnirt_register' not in nodes:

                # ANTS warp application

                def fsl_to_itk_conversion(source_file, reference, func_name):

                    # converts FSL-format .mat affine xfm into ANTS-format
                    # .txt; .mat affine comes from Func->Anat registration
                    fsl_to_itk_func_mni = create_wf_c3d_fsl_to_itk(name='fsl_to_itk_%s_%d' % (func_name, num_strat))

                    try:

                        # convert the .mat from linear Func->Anat to
                        # ANTS format
                        node, out_file = strat.get_node_from_resource_pool( \
                            'functional_to_anat_linear_xfm')
                        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                                         'inputspec.affine_file')

                        node, out_file = strat.get_node_from_resource_pool( \
                            reference)
                        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                                         'inputspec.reference_file')

                        node, out_file = strat.get_node_from_resource_pool( \
                            source_file)
                        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                                         'inputspec.source_file')

                    except:
                        logConnectionError('Functional Timeseries ' \
                                           'Registration to MNI space (ANTS)',
                                           num_strat, \
                                           strat.get_resource_pool(), '0016')
                        raise

                    strat.update_resource_pool({'itk_func_anat_affine_%s' % (func_name): (fsl_to_itk_func_mni, 'outputspec.itk_transform')})

                    strat.append_name(fsl_to_itk_func_mni.name)
                    create_log_node(fsl_to_itk_func_mni, 'outputspec.itk_transform', num_strat)

                def collect_transforms_func_mni(func_name):

                    # collects series of warps to be applied
                    collect_transforms_func_mni = \
                        create_wf_collect_transforms(name='collect_transforms_%s_%d' % (func_name, num_strat))

                    try:

                        # Field file from anatomical nonlinear registration
                        node, out_file = strat.get_node_from_resource_pool( \
                            'anatomical_to_mni_nonlinear_xfm')
                        workflow.connect(node, out_file,
                                         collect_transforms_func_mni,
                                         'inputspec.warp_file')

                        # initial transformation from anatomical registration
                        node, out_file = strat.get_node_from_resource_pool( \
                            'ants_initial_xfm')
                        workflow.connect(node, out_file,
                                         collect_transforms_func_mni,
                                         'inputspec.linear_initial')

                        # affine transformation from anatomical registration
                        node, out_file = strat.get_node_from_resource_pool( \
                            'ants_affine_xfm')
                        workflow.connect(node, out_file,
                                         collect_transforms_func_mni,
                                         'inputspec.linear_affine')

                        # rigid transformation from anatomical registration
                        node, out_file = strat.get_node_from_resource_pool( \
                            'ants_rigid_xfm')
                        workflow.connect(node, out_file,
                                         collect_transforms_func_mni,
                                         'inputspec.linear_rigid')

                        # Premat from Func->Anat linear reg and bbreg
                        # (if bbreg is enabled)
                        node, out_file = strat.get_node_from_resource_pool( \
                            'itk_func_anat_affine_%s' % func_name)
                        workflow.connect(node, out_file,
                                         collect_transforms_func_mni,
                                         'inputspec.fsl_to_itk_affine')

                    except:
                        logConnectionError('Functional Timeseries ' \
                                           'Registration to MNI space (ANTS)',
                                           num_strat, \
                                           strat.get_resource_pool(), '0016')
                        raise

                    strat.update_resource_pool({'itk_collected_warps_%s' % \
                                                (func_name): (
                    collect_transforms_func_mni, \
                    'outputspec.transformation_series')})

                    strat.append_name(collect_transforms_func_mni.name)
                    create_log_node(collect_transforms_func_mni, \
                                    'outputspec.transformation_series',
                                    num_strat)

                def ants_apply_warps_func_mni(input_node, input_outfile,
                                              ref_node, ref_outfile, standard,
                                              func_name, interp,
                                              input_image_type):

                    # converts FSL-format .mat affine xfm into ANTS-format
                    # .txt; .mat affine comes from Func->Anat registration
                    fsl_to_itk_func_mni = create_wf_c3d_fsl_to_itk(name= 'fsl_to_itk_%s_%d' % (func_name, num_strat))

                    # collects series of warps to be applied
                    collect_transforms_func_mni = \
                        create_wf_collect_transforms(name= 'collect_transforms_%s_%d' % (func_name, num_strat))

                    # apply ants warps
                    apply_ants_warp_func_mni = \
                        create_wf_apply_ants_warp(name='apply_ants_warp_%s_%d' % (func_name, num_strat),
                                                  ants_threads=int(num_ants_cores))

                    apply_ants_warp_func_mni.inputs.inputspec.reference_image = standard
                    apply_ants_warp_func_mni.inputs.inputspec.dimension = 3
                    apply_ants_warp_func_mni.inputs.inputspec.interpolation = interp
                    # input_image_type:
                    # (0 or 1 or 2 or 3)
                    # Option specifying the input image type of scalar
                    # (default), vector, tensor, or time series.
                    apply_ants_warp_func_mni.inputs.inputspec. \
                        input_image_type = input_image_type

                    try:

                        # convert the .mat from linear Func->Anat to
                        # ANTS format
                        node, out_file = strat.get_node_from_resource_pool( \
                            'functional_to_anat_linear_xfm')
                        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                                         'inputspec.affine_file')

                        node, out_file = strat.get_node_from_resource_pool(
                            "anatomical_brain")
                        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                                         'inputspec.reference_file')

                        workflow.connect(ref_node, ref_outfile,
                                         fsl_to_itk_func_mni,
                                         'inputspec.source_file')

                        # Field file from anatomical nonlinear registration
                        node, out_file = strat.get_node_from_resource_pool( \
                            'anatomical_to_mni_nonlinear_xfm')
                        workflow.connect(node, out_file,
                                         collect_transforms_func_mni,
                                         'inputspec.warp_file')

                        # initial transformation from anatomical registration
                        node, out_file = strat.get_node_from_resource_pool( \
                            'ants_initial_xfm')
                        workflow.connect(node, out_file,
                                         collect_transforms_func_mni,
                                         'inputspec.linear_initial')

                        # affine transformation from anatomical registration
                        node, out_file = strat.get_node_from_resource_pool( \
                            'ants_affine_xfm')
                        workflow.connect(node, out_file,
                                         collect_transforms_func_mni,
                                         'inputspec.linear_affine')

                        # rigid transformation from anatomical registration
                        node, out_file = strat.get_node_from_resource_pool( \
                            'ants_rigid_xfm')
                        workflow.connect(node, out_file,
                                         collect_transforms_func_mni,
                                         'inputspec.linear_rigid')

                        # Premat from Func->Anat linear reg and bbreg
                        # (if bbreg is enabled)
                        workflow.connect(fsl_to_itk_func_mni,
                                         'outputspec.itk_transform',
                                         collect_transforms_func_mni,
                                         'inputspec.fsl_to_itk_affine')

                        # this <node, out_file> pulls in directly because
                        # it pulls in the leaf in some instances
                        workflow.connect(input_node, input_outfile,
                                         apply_ants_warp_func_mni,
                                         'inputspec.input_image')

                        workflow.connect(collect_transforms_func_mni,
                                         'outputspec.transformation_series',
                                         apply_ants_warp_func_mni,
                                         'inputspec.transforms')


                    except:
                        logConnectionError('Functional Timeseries ' \
                                           'Registration to MNI space (ANTS)',
                                           num_strat, \
                                           strat.get_resource_pool(), '0016')
                        raise

                    strat.update_resource_pool({func_name: \
                                                    (apply_ants_warp_func_mni, \
                                                     'outputspec.output_image')})

                    strat.append_name(apply_ants_warp_func_mni.name)
                    create_log_node(apply_ants_warp_func_mni, \
                                    'outputspec.output_image', num_strat)

                # 4D FUNCTIONAL apply warp
                node, out_file = strat.get_leaf_properties()
                node2, out_file2 = \
                    strat.get_node_from_resource_pool("mean_functional")
                ants_apply_warps_func_mni(node, out_file,
                                          node2, out_file2,
                                          c.template_brain_only_for_func,
                                          "functional_to_standard",
                                          "Linear", 3)

                # 4D FUNCTIONAL MOTION-CORRECTED apply warp
                node, out_file = \
                    strat.get_node_from_resource_pool('motion_correct')
                node2, out_file2 = \
                    strat.get_node_from_resource_pool("mean_functional")
                ants_apply_warps_func_mni(node, out_file,
                                          node2, out_file2,
                                          c.template_brain_only_for_func,
                                          "motion_correct_to_standard",
                                          "Linear", 3)

                # FUNCTIONAL BRAIN MASK (binary, no timeseries) apply warp
                node, out_file = \
                    strat.get_node_from_resource_pool("functional_brain_mask")
                ants_apply_warps_func_mni(node, out_file,
                                          node, out_file,
                                          c.template_brain_only_for_func,
                                          "functional_brain_mask_to_standard",
                                          "NearestNeighbor", 0)

                # FUNCTIONAL MEAN (no timeseries) apply warp
                node, out_file = \
                    strat.get_node_from_resource_pool("mean_functional")
                ants_apply_warps_func_mni(node, out_file,
                                          node, out_file,
                                          c.template_brain_only_for_func,
                                          "mean_functional_to_standard",
                                          "Linear", 0)

                num_strat += 1

    strat_list += new_strat_list

    """""""""""""""""""""""""""""""""""""""""""""""""""
     OUTPUTS
    """""""""""""""""""""""""""""""""""""""""""""""""""

    '''
    Inserting VMHC
    Workflow
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runVMHC:

        for strat in strat_list:

            nodes = getNodeList(strat)

            if 'func_mni_fsl_warp' in nodes:
                vmhc = create_vmhc(False, 'vmhc_%d' % num_strat)
            else:
                vmhc = create_vmhc(True, 'vmhc_%d' % num_strat,
                                   int(num_ants_cores))

            vmhc.inputs.inputspec.standard_for_func = c.template_skull_for_func
            vmhc.inputs.fwhm_input.fwhm = c.fwhm
            vmhc.get_node('fwhm_input').iterables = ('fwhm', c.fwhm)

            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.rest_res')

                node, out_file = strat.get_node_from_resource_pool(
                    'functional_to_anat_linear_xfm')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.example_func2highres_mat')

                node, out_file = strat.get_node_from_resource_pool(
                    'functional_brain_mask')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.rest_mask')

                node, out_file = strat.get_node_from_resource_pool(
                    'mean_functional')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.mean_functional')

                node, out_file = strat.get_node_from_resource_pool(
                    'anatomical_brain')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.brain')

                if ('ANTS' in c.regOption) and \
                        ('anat_mni_fnirt_register' not in nodes) and \
                        ('anat_symmetric_mni_fnirt_register' not in nodes):

                    node, out_file = strat.get_node_from_resource_pool(
                        'ants_symmetric_initial_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.ants_symm_initial_xfm')

                    node, out_file = strat.get_node_from_resource_pool(
                        'ants_symmetric_rigid_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.ants_symm_rigid_xfm')

                    node, out_file = strat.get_node_from_resource_pool(
                        'ants_symmetric_affine_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.ants_symm_affine_xfm')

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_to_symmetric_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.ants_symm_warp_field')

                else:

                    node, out_file = strat.get_node_from_resource_pool(
                        'anatomical_to_symmetric_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.fnirt_nonlinear_warp')


            except:
                logConnectionError('VMHC', num_strat,
                                   strat.get_resource_pool(), '0019')
                raise

            strat.update_resource_pool(
                {'vmhc_raw_score': (vmhc, 'outputspec.VMHC_FWHM_img')})
            strat.update_resource_pool(
                {'vmhc_fisher_zstd': (vmhc, 'outputspec.VMHC_Z_FWHM_img')})
            strat.update_resource_pool({'vmhc_fisher_zstd_zstat_map': (
            vmhc, 'outputspec.VMHC_Z_stat_FWHM_img')})
            strat.append_name(vmhc.name)

            create_log_node(vmhc, 'outputspec.VMHC_FWHM_img', num_strat)

            num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting REHO
    Workflow
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runReHo:
        for strat in strat_list:

            preproc = create_reho()
            cluster_size = c.clusterSize
            # Check the cluster size is supported
            if not (cluster_size == 27 or \
                                cluster_size == 19 or \
                                cluster_size == 7):
                err_msg = 'Cluster size specified: %d, is not supported. ' \
                          'Change to 7, 19, or 27 and try again' % cluster_size
                raise Exception(err_msg)
            else:
                preproc.inputs.inputspec.cluster_size = cluster_size
                reho = preproc.clone('reho_%d' % num_strat)

            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 reho, 'inputspec.rest_res_filt')

                node, out_file = strat.get_node_from_resource_pool(
                    'functional_brain_mask')
                workflow.connect(node, out_file,
                                 reho, 'inputspec.rest_mask')
            except:
                logConnectionError('ReHo', num_strat,
                                   strat.get_resource_pool(), '0020')
                raise

            strat.update_resource_pool(
                {'reho': (reho, 'outputspec.raw_reho_map')})
            strat.append_name(reho.name)

            create_log_node(reho, 'outputspec.raw_reho_map', num_strat)

            num_strat += 1
    strat_list += new_strat_list

    '''
    Timeseries and SCA config selections processing
    '''

    ts_analysis_dict = {}
    sca_analysis_dict = {}

    # c.tsa_roi_paths and c.sca_roi_paths come in a format as such:
    # a list containing a dictionary
    # [{'/path/to/rois1.nii.gz': 'Avg, MultReg',
    # '/path/to/rois2.nii.gz': 'Avg, MultReg',
    # '/path/to/rois3.nii.gz': 'Avg, MultReg',
    # '/path/to/rois4.nii.gz': 'DualReg'}]

    if 1 in c.runROITimeseries:

        if c.tsa_roi_paths:
            tsa_roi_dict = c.tsa_roi_paths[0]
        else:
            err = "\n\n[!] CPAC says: Time Series Extraction is " \
                  "set to run, but no ROI NIFTI file paths were provided!" \
                  "\n\n"
            raise Exception(err)

        # flip the dictionary
        for roi_path in tsa_roi_dict.keys():

            for analysis_type in tsa_roi_dict[roi_path].split(","):

                analysis_type = analysis_type.replace(" ", "")

            if analysis_type not in ts_analysis_dict.keys():
                    ts_analysis_dict[analysis_type] = []

            ts_analysis_dict[analysis_type].append(roi_path)

    if 1 in c.runSCA:

        if c.sca_roi_paths:
            sca_roi_dict = c.sca_roi_paths[0]
        else:
            err = "\n\n[!] CPAC says: Seed-based Correlation Analysis is " \
                  "set to run, but no ROI NIFTI file paths were provided!" \
                  "\n\n"
            raise Exception(err)

        # flip the dictionary
        for roi_path in sca_roi_dict.keys():

            for analysis_type in sca_roi_dict[roi_path].split(","):

                analysis_type = analysis_type.replace(" ", "")

                if analysis_type not in sca_analysis_dict.keys():
                    sca_analysis_dict[analysis_type] = []

                sca_analysis_dict[analysis_type].append(roi_path)

    '''
    Spatial Regression Based Time Series
    '''

    new_strat_list = []
    num_strat = 0

    # if 1 in c.runSpatialRegression:

    if ("SpatialReg" in ts_analysis_dict.keys()) or \
            ("DualReg" in sca_analysis_dict.keys()):

        for strat in strat_list:

            if "SpatialReg" in ts_analysis_dict.keys():
                resample_spatial_map_to_native_space = pe.Node(
                    interface=fsl.FLIRT(),
                    name='resample_spatial_map_to_native_space_%d' % num_strat)
                resample_spatial_map_to_native_space.inputs.interp = 'nearestneighbour'
                resample_spatial_map_to_native_space.inputs.apply_xfm = True
                resample_spatial_map_to_native_space.inputs.in_matrix_file = c.identityMatrix

                spatial_map_dataflow = create_spatial_map_dataflow(
                    ts_analysis_dict["SpatialReg"],
                    'spatial_map_dataflow_%d' % num_strat)

                spatial_map_timeseries = get_spatial_map_timeseries(
                    'spatial_map_timeseries_%d' % num_strat)
                spatial_map_timeseries.inputs.inputspec.demean = True  # c.spatialDemean

            if "DualReg" in sca_analysis_dict.keys():
                resample_spatial_map_to_native_space_for_dr = pe.Node(
                    interface=fsl.FLIRT(),
                    name='resample_spatial_map_to_native_space_for_DR_%d' % num_strat)
                resample_spatial_map_to_native_space_for_dr.inputs.interp = 'nearestneighbour'
                resample_spatial_map_to_native_space_for_dr.inputs.apply_xfm = True
                resample_spatial_map_to_native_space_for_dr.inputs.in_matrix_file = c.identityMatrix

                spatial_map_dataflow_for_dr = create_spatial_map_dataflow(
                    sca_analysis_dict["DualReg"],
                    'spatial_map_dataflow_for_DR_%d' % num_strat)

                spatial_map_timeseries_for_dr = get_spatial_map_timeseries(
                    'spatial_map_timeseries_for_DR_%d' % num_strat)
                spatial_map_timeseries_for_dr.inputs.inputspec.demean = True  # c.spatialDemean

            try:

                if "SpatialReg" in ts_analysis_dict.keys():
                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_standard')
                    node2, out_file2 = strat.get_node_from_resource_pool(
                        'functional_brain_mask_to_standard')

                    # resample the input functional file and functional mask
                    # to spatial map
                    workflow.connect(node, out_file,
                                     resample_spatial_map_to_native_space,
                                     'reference')
                    workflow.connect(spatial_map_dataflow,
                                     'select_spatial_map.out_file',
                                     resample_spatial_map_to_native_space,
                                     'in_file')

                    # connect it to the spatial_map_timeseries
                    workflow.connect(resample_spatial_map_to_native_space,
                                     'out_file',
                                     spatial_map_timeseries,
                                     'inputspec.spatial_map')
                    workflow.connect(node2, out_file2,
                                     spatial_map_timeseries,
                                     'inputspec.subject_mask')
                    workflow.connect(node, out_file,
                                     spatial_map_timeseries,
                                     'inputspec.subject_rest')

                if "DualReg" in sca_analysis_dict.keys():
                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_standard')
                    node2, out_file2 = strat.get_node_from_resource_pool(
                        'functional_brain_mask_to_standard')

                    # resample the input functional file and functional mask
                    # to spatial map
                    workflow.connect(node, out_file,
                                     resample_spatial_map_to_native_space_for_dr,
                                     'reference')
                    workflow.connect(spatial_map_dataflow_for_dr,
                                     'select_spatial_map.out_file',
                                     resample_spatial_map_to_native_space_for_dr,
                                     'in_file')

                    # connect it to the spatial_map_timeseries
                    workflow.connect(
                        resample_spatial_map_to_native_space_for_dr,
                        'out_file',
                        spatial_map_timeseries_for_dr,
                        'inputspec.spatial_map')
                    workflow.connect(node2, out_file2,
                                     spatial_map_timeseries_for_dr,
                                     'inputspec.subject_mask')
                    workflow.connect(node, out_file,
                                     spatial_map_timeseries_for_dr,
                                     'inputspec.subject_rest')

            except Exception as e:
                logConnectionError('Spatial map timeseries extraction',
                                   num_strat, strat.get_resource_pool(),
                                   '0029')
                print "\nError message: %s\n\n" % e
                raise Exception

            if "SpatialReg" in ts_analysis_dict.keys():
                strat.append_name(spatial_map_timeseries.name)
                strat.update_resource_pool({'spatial_map_timeseries': (spatial_map_timeseries, 'outputspec.subject_timeseries')})
                create_log_node(spatial_map_timeseries, 'outputspec.subject_timeseries', num_strat)

            if "DualReg" in sca_analysis_dict.keys():
                strat.append_name(spatial_map_timeseries_for_dr.name)
                strat.update_resource_pool({'spatial_map_timeseries_for_DR': (spatial_map_timeseries_for_dr, 'outputspec.subject_timeseries')})
                create_log_node(spatial_map_timeseries_for_dr, 'outputspec.subject_timeseries', num_strat)

            if ("SpatialReg" in ts_analysis_dict.keys()) or \
                    ("DualReg" in sca_analysis_dict.keys()):
                num_strat += 1

    strat_list += new_strat_list

    '''
    ROI Based Time Series
    '''

    new_strat_list = []
    num_strat = 0

    if ("Avg" in ts_analysis_dict.keys()) or \
            ("Avg" in sca_analysis_dict.keys()) or \
            ("MultReg" in sca_analysis_dict.keys()):

        for strat in strat_list:

            if "Avg" in ts_analysis_dict.keys():
                resample_functional_to_roi = pe.Node(interface=fsl.FLIRT(),
                                                     name='resample_functional_to_roi_%d' % num_strat)
                resample_functional_to_roi.inputs.interp = 'trilinear'
                resample_functional_to_roi.inputs.apply_xfm = True
                resample_functional_to_roi.inputs.in_matrix_file = c.identityMatrix

                roi_dataflow = create_roi_mask_dataflow(
                    ts_analysis_dict["Avg"], 'roi_dataflow_%d' % num_strat)

                roi_timeseries = get_roi_timeseries(
                    'roi_timeseries_%d' % num_strat)

            if "Avg" in sca_analysis_dict.keys():
                # same workflow, except to run TSE and send it to the resource
                # pool so that it will not get sent to SCA
                resample_functional_to_roi_for_sca = pe.Node(
                    interface=fsl.FLIRT(),
                    name='resample_functional_to_roi_for_sca_%d' % num_strat)
                resample_functional_to_roi_for_sca.inputs.interp = 'trilinear'
                resample_functional_to_roi_for_sca.inputs.apply_xfm = True
                resample_functional_to_roi_for_sca.inputs.in_matrix_file = c.identityMatrix

                roi_dataflow_for_sca = create_roi_mask_dataflow(
                    sca_analysis_dict["Avg"],
                    'roi_dataflow_for_sca_%d' % num_strat)

                roi_timeseries_for_sca = get_roi_timeseries(
                    'roi_timeseries_for_sca_%d' % num_strat)

            if "MultReg" in sca_analysis_dict.keys():
                # same workflow, except to run TSE and send it to the resource
                # pool so that it will not get sent to SCA
                resample_functional_to_roi_for_multreg = pe.Node(
                    interface=fsl.FLIRT(),
                    name='resample_functional_to_roi_for_mult_reg_%d' % num_strat)
                resample_functional_to_roi_for_multreg.inputs.interp = 'trilinear'
                resample_functional_to_roi_for_multreg.inputs.apply_xfm = True
                resample_functional_to_roi_for_multreg.inputs.in_matrix_file = c.identityMatrix

                roi_dataflow_for_multreg = create_roi_mask_dataflow(
                    sca_analysis_dict["MultReg"],
                    'roi_dataflow_for_mult_reg_%d' % num_strat)

                roi_timeseries_for_multreg = get_roi_timeseries(
                    'roi_timeseries_for_mult_reg_%d' % num_strat)

            try:
                if "Avg" in ts_analysis_dict.keys():
                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_standard')

                    # resample the input functional file to roi
                    workflow.connect(node, out_file,
                                     resample_functional_to_roi, 'in_file')
                    workflow.connect(roi_dataflow, 'outputspec.out_file',
                                     resample_functional_to_roi, 'reference')

                    # connect it to the roi_timeseries
                    workflow.connect(roi_dataflow, 'outputspec.out_file',
                                     roi_timeseries, 'input_roi.roi')
                    workflow.connect(resample_functional_to_roi, 'out_file',
                                     roi_timeseries, 'inputspec.rest')

                if "Avg" in sca_analysis_dict.keys():
                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_standard')

                    # resample the input functional file to roi
                    workflow.connect(node, out_file,
                                     resample_functional_to_roi_for_sca,
                                     'in_file')
                    workflow.connect(roi_dataflow_for_sca,
                                     'outputspec.out_file',
                                     resample_functional_to_roi_for_sca,
                                     'reference')

                    # connect it to the roi_timeseries
                    workflow.connect(roi_dataflow_for_sca,
                                     'outputspec.out_file',
                                     roi_timeseries_for_sca, 'input_roi.roi')
                    workflow.connect(resample_functional_to_roi_for_sca,
                                     'out_file',
                                     roi_timeseries_for_sca, 'inputspec.rest')

                if "MultReg" in sca_analysis_dict.keys():
                    node, out_file = strat.get_node_from_resource_pool(
                        'functional_to_standard')

                    # resample the input functional file to roi
                    workflow.connect(node, out_file,
                                     resample_functional_to_roi_for_multreg,
                                     'in_file')
                    workflow.connect(roi_dataflow_for_multreg,
                                     'outputspec.out_file',
                                     resample_functional_to_roi_for_multreg,
                                     'reference')

                    # connect it to the roi_timeseries
                    workflow.connect(roi_dataflow_for_multreg,
                                     'outputspec.out_file',
                                     roi_timeseries_for_multreg,
                                     'input_roi.roi')
                    workflow.connect(resample_functional_to_roi_for_multreg,
                                     'out_file',
                                     roi_timeseries_for_multreg,
                                     'inputspec.rest')

            except:
                logConnectionError('ROI Timeseries analysis', num_strat,
                                   strat.get_resource_pool(), '0031')
                raise

            if "Avg" in ts_analysis_dict.keys():
                strat.append_name(roi_timeseries.name)
                strat.update_resource_pool({'roi_timeseries': (roi_timeseries, 'outputspec.roi_outputs')})
                create_log_node(roi_timeseries, 'outputspec.roi_outputs',
                                num_strat)

            if "Avg" in sca_analysis_dict.keys():
                strat.append_name(roi_timeseries_for_sca.name)
                strat.update_resource_pool({'roi_timeseries_for_SCA': (roi_timeseries_for_sca, 'outputspec.roi_outputs')})
                create_log_node(roi_timeseries_for_sca,
                                'outputspec.roi_outputs', num_strat)

            if "MultReg" in sca_analysis_dict.keys():
                strat.append_name(roi_timeseries_for_multreg.name)
                strat.update_resource_pool({'roi_timeseries_for_SCA_multreg': (roi_timeseries_for_multreg, 'outputspec.roi_outputs')})
                create_log_node(roi_timeseries_for_multreg,
                                'outputspec.roi_outputs', num_strat)

            if ("Avg" in ts_analysis_dict.keys()) or \
                    ("Avg" in sca_analysis_dict.keys()) or \
                    ("MultReg" in sca_analysis_dict.keys()):
                num_strat += 1

    strat_list += new_strat_list

    '''
    Voxel Based Time Series
    '''

    new_strat_list = []
    num_strat = 0

    # if 1 in c.runVoxelTimeseries:

    if "Voxel" in ts_analysis_dict.keys():

        for strat in strat_list:

            resample_functional_to_mask = pe.Node(interface=fsl.FLIRT(),
                                                  name='resample_functional_to_mask_%d' % num_strat)
            resample_functional_to_mask.inputs.interp = 'trilinear'
            resample_functional_to_mask.inputs.apply_xfm = True
            resample_functional_to_mask.inputs.in_matrix_file = c.identityMatrix

            mask_dataflow = create_roi_mask_dataflow(
                ts_analysis_dict["Voxel"], 'mask_dataflow_%d' % num_strat)

            voxel_timeseries = get_voxel_timeseries(
                'voxel_timeseries_%d' % num_strat)
            voxel_timeseries.inputs.inputspec.output_type = c.roiTSOutputs

            try:

                node, out_file = strat.get_node_from_resource_pool(
                    'functional_to_standard')

                # resample the input functional file to mask
                workflow.connect(node, out_file,
                                 resample_functional_to_mask, 'in_file')
                workflow.connect(mask_dataflow, 'outputspec.out_file',
                                 resample_functional_to_mask, 'reference')

                # connect it to the voxel_timeseries
                workflow.connect(mask_dataflow, 'outputspec.out_file',
                                 voxel_timeseries, 'input_mask.mask')
                workflow.connect(resample_functional_to_mask, 'out_file',
                                 voxel_timeseries, 'inputspec.rest')


            except:
                logConnectionError('Voxel timeseries analysis', num_strat,
                                   strat.get_resource_pool(), '0030')
                raise

            strat.append_name(voxel_timeseries.name)
            strat.update_resource_pool({'voxel_timeseries': (voxel_timeseries, 'outputspec.mask_outputs')})
            create_log_node(voxel_timeseries, 'outputspec.mask_outputs',
                            num_strat)
            num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting SCA
    Workflow for ROI INPUT
    '''

    new_strat_list = []
    num_strat = 0

    # if 1 in c.runSCA and (1 in c.runROITimeseries):

    if "Avg" in sca_analysis_dict.keys():

        for strat in strat_list:

            sca_roi = create_sca('sca_roi_%d' % num_strat)

            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 sca_roi, 'inputspec.functional_file')

                node, out_file = strat.get_node_from_resource_pool(
                    'roi_timeseries_for_SCA')
                workflow.connect(node, (out_file, extract_one_d),
                                 sca_roi, 'inputspec.timeseries_one_d')
            except:
                logConnectionError('SCA ROI', num_strat,
                                   strat.get_resource_pool(), '0032')
                raise

            strat.update_resource_pool({'sca_roi_files': (sca_roi, 'outputspec.correlation_files')})

            create_log_node(sca_roi, 'outputspec.correlation_stack',
                            num_strat)

            strat.append_name(sca_roi.name)
            num_strat += 1

    strat_list += new_strat_list

    '''
    (Dual Regression) Temporal Regression for Dual Regression
    '''

    new_strat_list = []
    num_strat = 0

    # if 1 in c.runDualReg and (1 in c.runSpatialRegression):

    if "DualReg" in sca_analysis_dict.keys():

        for strat in strat_list:

            dr_temp_reg = create_temporal_reg(
                'temporal_dual_regression_%d' % num_strat)
            dr_temp_reg.inputs.inputspec.normalize = c.mrsNorm
            dr_temp_reg.inputs.inputspec.demean = True  # c.mrsDemean

            try:
                node, out_file = strat.get_node_from_resource_pool(
                    'spatial_map_timeseries_for_DR')

                node2, out_file2 = strat.get_leaf_properties()
                node3, out_file3 = strat.get_node_from_resource_pool(
                    'functional_brain_mask')

                workflow.connect(node2, out_file2,
                                 dr_temp_reg, 'inputspec.subject_rest')

                workflow.connect(node, out_file,
                                 dr_temp_reg, 'inputspec.subject_timeseries')

                workflow.connect(node3, out_file3,
                                 dr_temp_reg, 'inputspec.subject_mask')

            except:
                logConnectionError(
                    'Temporal multiple regression for dual regression',
                    num_strat, strat.get_resource_pool(), '0033')
                raise

            strat.update_resource_pool({'dr_tempreg_maps_files': (dr_temp_reg, 'outputspec.temp_reg_map_files')})
            strat.update_resource_pool({'dr_tempreg_maps_zstat_files': (dr_temp_reg, 'outputspec.temp_reg_map_z_files')})

            strat.append_name(dr_temp_reg.name)

            create_log_node(dr_temp_reg, 'outputspec.temp_reg_map', num_strat)

            num_strat += 1

    '''
    elif 1 in c.runDualReg and (0 in c.runSpatialRegression):
        logger.info("\n\n" + "WARNING: Dual Regression - Spatial regression was turned off for at least one of the strategies.")
        logger.info("Spatial regression is required for dual regression." + "\n\n")
    '''

    strat_list += new_strat_list

    '''
    (Multiple Regression) Temporal Regression for SCA
    '''

    new_strat_list = []
    num_strat = 0

    # if 1 in c.runMultRegSCA and (1 in c.runROITimeseries):

    if "MultReg" in sca_analysis_dict.keys():

        for strat in strat_list:

            sc_temp_reg = create_temporal_reg(
                'temporal_regression_sca_%d' % num_strat, which='RT')
            sc_temp_reg.inputs.inputspec.normalize = c.mrsNorm
            sc_temp_reg.inputs.inputspec.demean = True  # c.mrsDemean

            try:
                node, out_file = strat.get_node_from_resource_pool(
                    'functional_to_standard')
                node2, out_file2 = strat.get_node_from_resource_pool(
                    'roi_timeseries_for_SCA_multreg')
                node3, out_file3 = strat.get_node_from_resource_pool(
                    'functional_brain_mask_to_standard')

                workflow.connect(node, out_file,
                                 sc_temp_reg, 'inputspec.subject_rest')

                workflow.connect(node2, (out_file2, extract_one_d),
                                 sc_temp_reg, 'inputspec.subject_timeseries')

                workflow.connect(node3, out_file3,
                                 sc_temp_reg, 'inputspec.subject_mask')

            except:
                logConnectionError(
                    'Temporal multiple regression for seed based connectivity',
                    num_strat, strat.get_resource_pool(), '0037')
                raise

            strat.update_resource_pool({'sca_tempreg_maps_files': (sc_temp_reg, 'outputspec.temp_reg_map_files')})
            strat.update_resource_pool({'sca_tempreg_maps_zstat_files': (sc_temp_reg, 'outputspec.temp_reg_map_z_files')})

            create_log_node(sc_temp_reg, 'outputspec.temp_reg_map', num_strat)

            strat.append_name(sc_temp_reg.name)
            num_strat += 1

    strat_list += new_strat_list

    '''
    Inserting Network centrality
    '''

    new_strat_list = []
    num_strat = 0

    # If we're running centrality
    if 1 in c.runNetworkCentrality:

        # validate the mask file path
        if not c.templateSpecificationFile.endswith(".nii") and \
                not c.templateSpecificationFile.endswith(".nii.gz"):
            err = "\n\n[!] CPAC says: The Network Centrality mask " \
                  "specification file must be a NIFTI file (ending in .nii " \
                  "or .nii.gz).\nFile path you provided: %s\n\n" \
                  % c.templateSpecificationFile

            raise Exception(err)

        # Check for the existence of AFNI 3dDegreeCentrality/LFCD binaries
        import subprocess
        try:
            ret_code = subprocess.check_call(['which', '3dDegreeCentrality'],
                                             stdout=open(os.devnull, 'wb'))
            if ret_code == 0:
                afni_centrality_found = True
                logger.info('Using AFNI centrality function')
        except subprocess.CalledProcessError as exc:
            afni_centrality_found = False
            logger.info('Using C-PAC centrality function')
        try:
            ret_code = subprocess.check_call(['which', '3dLFCD'],
                                             stdout=open(os.devnull, 'wb'))
            if ret_code == 0:
                afni_lfcd_found = True
                logger.info('Using AFNI LFCD function')
        except subprocess.CalledProcessError as exc:
            afni_lfcd_found = False
            logger.info('Using C-PAC LFCD function')

        # For each desired strategy
        for strat in strat_list:

            # Resample the functional mni to the centrality mask resolution
            resample_functional_to_template = pe.Node(interface=fsl.FLIRT(),
                                                      name='resample_functional_to_template_%d' % num_strat)
            resample_functional_to_template.inputs.interp = 'trilinear'
            resample_functional_to_template.inputs.in_matrix_file = c.identityMatrix
            resample_functional_to_template.inputs.apply_xfm = True

            # Get nipype  node and out file of the func mni img
            node, out_file = strat.get_node_from_resource_pool(
                'functional_to_standard')

            # Resample the input functional file to template(roi/mask)
            workflow.connect(node, out_file,
                             resample_functional_to_template, 'in_file')

            resample_functional_to_template.inputs.reference = \
                c.templateSpecificationFile

            # Init merge node for appending method output lists to one another
            merge_node = pe.Node(util.Function(input_names=['deg_list',
                                                            'eig_list',
                                                            'lfcd_list'],
                                               output_names=['merged_list'],
                                               function=merge_lists),
                                 name='merge_node_%d' % num_strat)

            # Function to connect the CPAC centrality python workflow
            # into pipeline
            def connectCentralityWorkflow(methodOption,
                                          thresholdOption,
                                          threshold,
                                          weightOptions,
                                          mList):

                # Create centrality workflow
                network_centrality = \
                    create_resting_state_graphs(
                        wf_name='network_centrality_%d-%s' \
                                % (num_strat, methodOption),
                        allocated_memory=c.memoryAllocatedForDegreeCentrality)

                # Connect resampled (to template/mask resolution)
                # functional_mni to inputspec
                workflow.connect(resample_functional_to_template, 'out_file',
                                 network_centrality, 'inputspec.in_file')
                # Subject mask/parcellation image
                network_centrality.inputs.inputspec.template = \
                    c.templateSpecificationFile
                # Give which method we're doing
                network_centrality.inputs.inputspec.method_option = \
                    methodOption
                # Type of threshold
                network_centrality.inputs.inputspec.threshold_option = \
                    thresholdOption
                # Connect threshold value (float)
                network_centrality.inputs.inputspec.threshold = threshold

                # Merge output with others via merge_node connection
                workflow.connect(network_centrality,
                                 'outputspec.centrality_outputs',
                                 merge_node,
                                 mList)
                # Append this as a strategy
                strat.append_name(network_centrality.name)
                # Create log node for strategy
                create_log_node(network_centrality,
                                'outputspec.centrality_outputs',
                                num_strat)

            # Function to connect the afni 3dDegreeCentrality workflow
            # into pipeline
            def connect_afni_centrality_wf(method_option, threshold_option,
                                           threshold):

                # Import packages
                from CPAC.network_centrality.afni_network_centrality \
                    import create_afni_centrality_wf
                import CPAC.network_centrality.utils as cent_utils

                # Init variables
                # Set method_options variables
                if method_option == 'degree':
                    out_list = 'deg_list'
                elif method_option == 'eigenvector':
                    out_list = 'eig_list'
                elif method_option == 'lfcd':
                    out_list = 'lfcd_list'

                # Init workflow name and resource limits
                wf_name = 'afni_centrality_%d_%s' % (num_strat, method_option)
                num_threads = c.maxCoresPerParticipant
                memory = c.memoryAllocatedForDegreeCentrality

                # Format method and threshold options properly and check for
                # errors
                method_option, threshold_option = \
                    cent_utils.check_centrality_params(method_option,
                                                       threshold_option,
                                                       threshold)

                # Change sparsity thresholding to % to work with afni
                if threshold_option == 'sparsity':
                    threshold = threshold * 100

                # Init the workflow
                afni_centrality_wf = \
                    create_afni_centrality_wf(wf_name, method_option,
                                              threshold_option,
                                              threshold, num_threads, memory)

                # Connect pipeline resources to workflow
                # Dataset
                workflow.connect(resample_functional_to_template, 'out_file',
                                 afni_centrality_wf, 'inputspec.in_file')
                # Mask
                afni_centrality_wf.inputs.inputspec.template = \
                    c.templateSpecificationFile

                # Connect outputs to merge node
                workflow.connect(afni_centrality_wf,
                                 'outputspec.outfile_list',
                                 merge_node,
                                 out_list)

            # Degree/eigen check
            if afni_centrality_found:
                if c.degWeightOptions.count(True) > 0:
                    connect_afni_centrality_wf('degree',
                                               c.degCorrelationThresholdOption,
                                               c.degCorrelationThreshold)
                if c.eigWeightOptions.count(True) > 0:
                    connect_afni_centrality_wf('eigenvector',
                                               c.eigCorrelationThresholdOption,
                                               c.eigCorrelationThreshold)
            # Otherwise run the CPAC python workflow
            else:
                # If we're calculating degree centrality
                if c.degWeightOptions.count(True) > 0:
                    connectCentralityWorkflow('degree',
                                              c.degCorrelationThresholdOption,
                                              c.degCorrelationThreshold,
                                              c.degWeightOptions,
                                              'deg_list')
                # If we're calculating eigenvector centrality
                if c.eigWeightOptions.count(True) > 0:
                    connectCentralityWorkflow('eigenvector',
                                              c.eigCorrelationThresholdOption,
                                              c.eigCorrelationThreshold,
                                              c.eigWeightOptions,
                                              'eig_list')
            # LFCD check
            if afni_lfcd_found:
                # If we're calculating lFCD
                if c.lfcdWeightOptions.count(True) > 0:
                    connect_afni_centrality_wf('lfcd',
                                               c.lfcdCorrelationThresholdOption,
                                               c.lfcdCorrelationThreshold)
            # Otherwise run the CPAC python workflow
            else:
                # If we're calculating lFCD
                if c.lfcdWeightOptions.count(True) > 0:
                    connectCentralityWorkflow('lfcd',
                                              c.lfcdCorrelationThresholdOption,
                                              c.lfcdCorrelationThreshold,
                                              c.lfcdWeightOptions,
                                              'lfcd_list')

            # Update resource pool with centrality outputs
            try:
                strat.update_resource_pool(
                    {'centrality_outputs': (merge_node, 'merged_list')})
            except:
                logConnectionError('Network Centrality', num_strat,
                                   strat.get_resource_pool(), '0050')
                raise

            if 0 in c.runNetworkCentrality:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            num_strat += 1

    strat_list += new_strat_list

    num_strat = 0

    """""""""""""""""""""""""""""""""""""""""""""""""""
     Apply warps, Z-scoring, Smoothing, Averages
    """""""""""""""""""""""""""""""""""""""""""""""""""

    ''''''
    ''' APPLY WARPS '''

    def output_to_standard(output_name, strat, num_strat, pipeline_config_obj,
                           map_node=False, input_image_type=0):

        nodes = getNodeList(strat)

        if 'apply_ants_warp_functional_to_standard' in nodes:

            # ANTS WARP APPLICATION

            # convert the func-to-anat linear warp from FSL FLIRT to
            # ITK (ANTS) format
            fsl_to_itk_convert = create_wf_c3d_fsl_to_itk(input_image_type,
                                                          map_node,
                                                          name='{0}_fsl_to_itk_{1}'.format(output_name, num_strat))

            # collect the list of warps into a single stack to feed into the
            # ANTS warp apply tool
            collect_transforms = create_wf_collect_transforms(map_node,
                                                              name='{0}_collect_transforms_{1}'.format(output_name, num_strat))

            # ANTS apply warp
            apply_ants_warp = create_wf_apply_ants_warp(map_node,
                                                        name='{0}_to_standard_{1}'.format(output_name, num_strat),
                                                        ants_threads=int(pipeline_config_obj.num_ants_threads))

            apply_ants_warp.inputs.inputspec.dimension = 3
            apply_ants_warp.inputs.inputspec.interpolation = 'Linear'
            apply_ants_warp.inputs.inputspec.reference_image = \
                pipeline_config_obj.template_brain_only_for_func

            apply_ants_warp.inputs.inputspec.input_image_type = \
                input_image_type

            try:
                # affine from FLIRT func->anat linear registration
                node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                workflow.connect(node, out_file, fsl_to_itk_convert,
                                 'inputspec.affine_file')

                # reference used in FLIRT func->anat linear registration
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file, fsl_to_itk_convert,
                                 'inputspec.reference_file')

                # output file to be converted
                node, out_file = \
                    strat.get_node_from_resource_pool(output_name)
                workflow.connect(node, out_file, fsl_to_itk_convert,
                                 'inputspec.source_file')

                # nonlinear warp from anatomical->template ANTS registration
                node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                workflow.connect(node, out_file, collect_transforms,
                                 'inputspec.warp_file')

                # linear initial from anatomical->template ANTS registration
                node, out_file = strat.get_node_from_resource_pool('ants_initial_xfm')
                workflow.connect(node, out_file, collect_transforms,
                                 'inputspec.linear_initial')

                # linear affine from anatomical->template ANTS registration
                node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                workflow.connect(node, out_file, collect_transforms,
                                 'inputspec.linear_affine')

                # rigid affine from anatomical->template ANTS registration
                node, out_file = strat.get_node_from_resource_pool('ants_rigid_xfm')
                workflow.connect(node, out_file, collect_transforms,
                                 'inputspec.linear_rigid')

                # converted FLIRT func->anat affine, now in ITK (ANTS) format
                workflow.connect(fsl_to_itk_convert,
                                 'outputspec.itk_transform',
                                 collect_transforms,
                                 'inputspec.fsl_to_itk_affine')

                # output file to be converted
                node, out_file = strat.get_node_from_resource_pool(output_name)
                workflow.connect(node, out_file, apply_ants_warp,
                                 'inputspec.input_image')

                # collection of warps to be applied to the output file
                workflow.connect(collect_transforms,
                                 'outputspec.transformation_series',
                                 apply_ants_warp,
                                 'inputspec.transforms')

            except:
                logConnectionError('{0} to MNI (ANTS)'.format(output_name),
                                   num_strat, strat.get_resource_pool(),
                                   '0022')
                raise

            strat.update_resource_pool({'{0}_to_standard'.format(output_name): (apply_ants_warp, 'outputspec.output_image')})
            strat.append_name(apply_ants_warp.name)

            num_strat += 1

        else:
            # FSL WARP APPLICATION
            if map_node:
                apply_fsl_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                                            name='{0}_to_standard_{1}'.format(output_name, num_strat),
                                            iterfield=['in_file'])
            else:
                apply_fsl_warp = pe.Node(interface=fsl.ApplyWarp(),
                                         name='{0}_to_standard_{1}'.format(output_name,
                                                                           num_strat))

            apply_fsl_warp.inputs.ref_file = \
                pipeline_config_obj.template_skull_for_func

            try:
                # output file to be warped
                node, out_file = strat.get_node_from_resource_pool(output_name)
                workflow.connect(node, out_file, apply_fsl_warp, 'in_file')

                # linear affine from func->anat linear FLIRT registration
                node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                workflow.connect(node, out_file, apply_fsl_warp, 'premat')

                # nonlinear warp from anatomical->template FNIRT registration
                node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                workflow.connect(node, out_file, apply_fsl_warp, 'field_file')

            except:
                logConnectionError('{0} to MNI (FSL)'.format(output_name),
                                   num_strat, strat.get_resource_pool(),
                                   '0021')
                raise Exception

            strat.update_resource_pool({'{0}_to_standard'.format(output_name): (apply_fsl_warp, 'out_file')})
            strat.append_name(apply_fsl_warp.name)

        return strat

    ''''''
    ''' Z-SCORING '''

    def z_score_standardize(output_name, mask_name, strat, num_strat,
                            map_node=False):

        # call the z-scoring sub-workflow builder
        z_score_std = get_zscore(output_name, map_node,
                                 'z_score_std_%s_%d' % (output_name, num_strat))

        try:
            node, out_file = strat.get_node_from_resource_pool(output_name)
            workflow.connect(node, out_file, z_score_std,
                             'inputspec.input_file')

            # get the mask
            if "/" in mask_name:
                # mask_name is a direct file path and not the name of a
                # resource pool key
                z_score_std.inputs.inputspec.mask_file = mask_name
            else:
                node, out_file = strat.get_node_from_resource_pool(mask_name)
                workflow.connect(node, out_file,
                                 z_score_std, 'inputspec.mask_file')

        except Exception as e:
            logConnectionError('%s z-score standardize' % output_name,
                               num_strat, strat.get_resource_pool(),
                               '0127', e)
            raise

        strat.append_name(z_score_std.name)
        strat.update_resource_pool({'{0}_zstd'.format(output_name): (z_score_std, 'outputspec.z_score_img')})

        return strat

    def fisher_z_score_standardize(output_name, timeseries_oned_file, strat,
                                   num_strat, map_node=False):

        # call the fisher r-to-z sub-workflow builder
        fisher_z_score_std = get_fisher_zscore(output_name, map_node,
                                               'fisher_z_score_std_%s_%d' \
                                               % (output_name, num_strat))

        try:
            node, out_file = strat. \
                get_node_from_resource_pool(output_name)

            workflow.connect(node, out_file, fisher_z_score_std,
                             'inputspec.correlation_file')

            node, out_file = strat. \
                get_node_from_resource_pool(timeseries_oned_file)
            workflow.connect(node, out_file, fisher_z_score_std,
                             'inputspec.timeseries_one_d')

        except:
            logConnectionError('%s fisher z-score standardize' % output_name,
                               num_strat, strat.get_resource_pool(), '0128')
            raise

        strat.append_name(fisher_z_score_std.name)
        strat.update_resource_pool({'{0}_fisher_zstd'.format(output_name): (fisher_z_score_std, 'outputspec.fisher_z_score_img')})

        return strat

    ''''''
    ''' SMOOTHING '''

    def output_smooth(output_name, mask_name, strat, num_strat,
                      map_node=False):

        if map_node:
            output_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                       name='{0}_smooth_{1}'.format(output_name,
                                                                    num_strat),
                                       iterfield=['in_file'])
        else:
            output_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                                    name='{0}_smooth_{1}'.format(output_name,
                                                                    num_strat))

        try:
            # get the resource to be smoothed
            node, out_file = strat.get_node_from_resource_pool(output_name)

            workflow.connect(node, out_file, output_smooth, 'in_file')

            # get the parameters for fwhm
            workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                             output_smooth, 'op_string')

            # get the mask
            if "/" in mask_name:
                # mask_name is a direct file path and not the name of a
                # resource pool key
                output_smooth.inputs.operand_files = mask_name
            else:
                node, out_file = strat.get_node_from_resource_pool(mask_name)
                workflow.connect(node, out_file,
                                 output_smooth, 'operand_files')

        except:
            logConnectionError('{0} smooth'.format(output_name), num_strat,
                               strat.get_resource_pool(), '0027')
            raise

        strat.append_name(output_smooth.name)
        strat.update_resource_pool({'{0}_smooth'.format(output_name): (output_smooth, 'out_file')})

        return strat

    ''''''
    ''' AVERAGING '''

    def calc_avg(output_name, strat, num_strat, map_node=False):
        """Calculate the average of an output using AFNI 3dmaskave."""

        extract_imports = ['import os']

        if map_node:
            calc_average = pe.MapNode(interface=preprocess.Maskave(),
                                      name='{0}_mean_{1}'.format(output_name,
                                                                 num_strat),
                                      iterfield=['in_file'])

            mean_to_csv = pe.MapNode(util.Function(input_names=['in_file',
                                                                'output_name'],
                                                   output_names=['output_mean'],
                                                   function=extract_output_mean,
                                                   imports=extract_imports),
                                     name='{0}_mean_to_txt_{1}'.format(output_name,
                                                                       num_strat),
                                     iterfield=['in_file'])
        else:
            calc_average = pe.Node(interface=preprocess.Maskave(),
                                   name='{0}_mean_{1}'.format(output_name,
                                                              num_strat))

            mean_to_csv = pe.Node(util.Function(input_names=['in_file',
                                                             'output_name'],
                                                output_names=['output_mean'],
                                                function=extract_output_mean,
                                                imports=extract_imports),
                                  name='{0}_mean_to_txt_{1}'.format(output_name,
                                                                    num_strat))

        mean_to_csv.inputs.output_name = output_name

        try:
            node, out_file = strat.get_node_from_resource_pool(output_name)
            workflow.connect(node, out_file, calc_average, 'in_file')
            workflow.connect(calc_average, 'out_file', mean_to_csv, 'in_file')
        except:
            logConnectionError('{0} calc average'.format(output_name),
                               num_strat, strat.get_resource_pool(), '0128')
            raise

        strat.append_name(calc_average.name)
        strat.update_resource_pool({'output_means.@{0}_average'.format(output_name): (mean_to_csv, 'output_mean')})

        return strat

    '''
    Loop through the resource pool and connect the nodes for:
        - applying warps to standard
        - z-score standardization
        - smoothing
        - calculating output averages
    '''

    num_strat = 0
    for strat in strat_list:

        if 1 in c.runRegisterFuncToMNI:
            rp = strat.get_resource_pool()
            for key in sorted(rp.keys()):
                # connect nodes to apply warps to template
                if key in outputs_native_nonsmooth:
                    # smoothing happens at the end, so only the non-smooth
                    # named output labels for the native-space outputs
                    strat = output_to_standard(key, strat, num_strat, c)
                elif key in outputs_native_nonsmooth_mult:
                    strat = output_to_standard(key, strat, num_strat, c,
                                               map_node=True)

        if "Before" in c.smoothing_order:
            # run smoothing before Z-scoring
            if 1 in c.run_smoothing:
                rp = strat.get_resource_pool()
                for key in sorted(rp.keys()):
                    # connect nodes for smoothing
                    if "centrality" in key:
                        # centrality needs its own mask
                        strat = output_smooth(key,
                                              c.templateSpecificationFile,
                                              strat, num_strat, map_node=True)
                    elif key in outputs_native_nonsmooth:
                        # native space
                        strat = output_smooth(key, "functional_brain_mask",
                                              strat, num_strat)
                    elif key in outputs_native_nonsmooth_mult:
                        # native space with multiple files (map nodes)
                        strat = output_smooth(key, "functional_brain_mask",
                                              strat, num_strat, map_node=True)
                    elif key in outputs_template_nonsmooth:
                        # template space
                        strat = output_smooth(key,
                                              "functional_brain_mask_to_standard",
                                              strat, num_strat)
                    elif key in outputs_template_nonsmooth_mult:
                        # template space with multiple files (map nodes)
                        strat = output_smooth(key,
                                              "functional_brain_mask_to_standard",
                                              strat, num_strat, map_node=True)

            if 1 in c.runZScoring:
                rp = strat.get_resource_pool()
                for key in sorted(rp.keys()):
                    # connect nodes for z-score standardization
                    if "sca_roi_files_to_standard" in key:
                        # correlation files need the r-to-z
                        strat = fisher_z_score_standardize(key,
                                                           "roi_timeseries_for_SCA",
                                                           strat, num_strat,
                                                           map_node=True)
                    elif "centrality" in key:
                        # specific mask
                        strat = z_score_standardize(key,
                                                    c.templateSpecificationFile,
                                                    strat, num_strat,
                                                    map_node=True)
                    elif key in outputs_template_raw:
                        # raw score, in template space
                        strat = z_score_standardize(key,
                                                    "functional_brain_mask_to_standard",
                                                    strat, num_strat)
                    elif key in outputs_template_raw_mult:
                        # same as above but multiple files so mapnode required
                        strat = z_score_standardize(key,
                                                    "functional_brain_mask_to_standard",
                                                    strat, num_strat,
                                                    map_node=True)

        elif "After" in c.smoothing_order:
            # run smoothing after Z-scoring
            if 1 in c.runZScoring:
                rp = strat.get_resource_pool()
                for key in sorted(rp.keys()):
                    # connect nodes for z-score standardization
                    if "sca_roi_files_to_standard" in key:
                        # correlation files need the r-to-z
                        strat = fisher_z_score_standardize(key,
                                                           "roi_timeseries_for_SCA",
                                                           strat, num_strat,
                                                           map_node=True)
                    elif "centrality" in key:
                        # specific mask
                        strat = z_score_standardize(key,
                                                    c.templateSpecificationFile,
                                                    strat, num_strat,
                                                    map_node=True)
                    elif key in outputs_template_raw:
                        # raw score, in template space
                        strat = z_score_standardize(key,
                                                    "functional_brain_mask_to_standard",
                                                    strat, num_strat)
                    elif key in outputs_template_raw_mult:
                        # same as above but multiple files so mapnode required
                        strat = z_score_standardize(key,
                                                    "functional_brain_mask_to_standard",
                                                    strat, num_strat,
                                                    map_node=True)

            if 1 in c.run_smoothing:
                rp = strat.get_resource_pool()
                for key in sorted(rp.keys()):
                    # connect nodes for smoothing
                    if "centrality" in key:
                        # centrality needs its own mask
                        strat = output_smooth(key,
                                              c.templateSpecificationFile,
                                              strat, num_strat, map_node=True)
                    elif key in outputs_native_nonsmooth:
                        # native space
                        strat = output_smooth(key, "functional_brain_mask",
                                              strat, num_strat)
                    elif key in outputs_native_nonsmooth_mult:
                        # native space with multiple files (map nodes)
                        strat = output_smooth(key, "functional_brain_mask",
                                              strat, num_strat, map_node=True)
                    elif key in outputs_template_nonsmooth:
                        # template space
                        strat = output_smooth(key,
                                              "functional_brain_mask_to_standard",
                                              strat, num_strat)
                    elif key in outputs_template_nonsmooth_mult:
                        # template space with multiple files (map nodes)
                        strat = output_smooth(key,
                                              "functional_brain_mask_to_standard",
                                              strat, num_strat, map_node=True)

        rp = strat.get_resource_pool()
        for key in sorted(rp.keys()):
            # connect nodes to calculate averages
            if key in outputs_average:
                # the outputs we need the averages for
                strat = calc_avg(key, strat, num_strat)
            elif key in outputs_average_mult:
                # those outputs, but the ones with multiple files (map nodes)
                strat = calc_avg(key, strat, num_strat, map_node=True)

        num_strat += 1

    logger.info('\n\n' + 'Pipeline building completed.' + '\n\n')

    """""""""""""""""""""""""""""""""""""""""""""""""""
     QUALITY CONTROL
    """""""""""""""""""""""""""""""""""""""""""""""""""
   
    if 1 in c.generateQualityControlImages:

        preproc, out_file = strat.get_node_from_resource_pool('functional_preprocessed')
        brain_mask, mask_file = strat.get_node_from_resource_pool('functional_brain_mask')
        func_to_anat_xfm, xfm_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
        anat_ref, ref_file = strat.get_node_from_resource_pool('anatomical_brain')
        mfa, mfa_file = strat.get_node_from_resource_pool('mean_functional_in_anat')

        # register color palettes
        register_pallete(os.path.realpath(
                os.path.join(CPAC.__path__[0], 'qc', 'red.py')), 'red')
        register_pallete(os.path.realpath(
                os.path.join(CPAC.__path__[0], 'qc', 'green.py')), 'green')
        register_pallete(os.path.realpath(
                os.path.join(CPAC.__path__[0], 'qc', 'blue.py')), 'blue')
        register_pallete(os.path.realpath(
                os.path.join(CPAC.__path__[0], 'qc', 'red_to_blue.py')), 'red_to_blue')
        register_pallete(os.path.realpath(
                os.path.join(CPAC.__path__[0], 'qc', 'cyan_to_yellow.py')), 'cyan_to_yellow')
                
        hist = pe.Node(util.Function(input_names=['measure_file','measure'],output_names = ['hist_path'],function = gen_histogram),name = 'histogram')

        for strat in strat_list:

            nodes = getNodeList(strat)

            # make SNR plot
            try:
                std_dev_imports = ['import os', 'import subprocess']
                std_dev = pe.Node(util.Function(input_names=['mask_',
                                                             'func_'],
                                                output_names=['new_fname'],
                                                function=gen_std_dev,
                                                imports=std_dev_imports),
                                  name='std_dev_%d' % num_strat)

                workflow.connect(preproc, out_file,
                                 std_dev, 'func_')

                workflow.connect(brain_mask, mask_file,
                                 std_dev, 'mask_')

                std_dev_anat_imports = ['import os',
                                        'import subprocess']
                std_dev_anat = pe.Node(util.Function(input_names=['func_',
                                                                  'ref_',
                                                                  'xfm_',
                                                                  'interp_'],
                                                     output_names=['new_fname'],
                                                     function=gen_func_anat_xfm,
                                                     imports=std_dev_anat_imports),
                                       name='std_dev_anat_%d' % num_strat)

                workflow.connect(std_dev, 'new_fname',
                                 std_dev_anat, 'func_')

                workflow.connect(func_to_anat_xfm, xfm_file,
                                 std_dev_anat, 'xfm_')

                workflow.connect(anat_ref, ref_file,
                                 std_dev_anat, 'ref_')

                snr_imports = ['import os', 'import subprocess']
                snr = pe.Node(util.Function(input_names=['std_dev',
                                                         'mean_func_anat'],
                                            output_names=['new_fname'],
                                            function=gen_snr,
                                            imports=snr_imports),
                              name='snr_%d' % num_strat)

                workflow.connect(std_dev_anat, 'new_fname', snr, 'std_dev')
                workflow.connect(mfa, mfa_file, snr, 'mean_func_anat')

                snr_val_imports = ['import os',
                                   'import nibabel as nb',
                                   'import numpy.ma as ma']
                snr_val = pe.Node(util.Function(input_names=['measure_file'],
                                                output_names=['snr_storefl'],
                                                function=cal_snr_val,
                                                imports=snr_val_imports),
                                  name='snr_val%d' % num_strat)
                std_dev_anat.inputs.interp_ = 'trilinear'

                workflow.connect(snr, 'new_fname',
                                 snr_val, 'measure_file')

                hist_ = hist.clone('hist_snr_%d' % num_strat)
                hist_.inputs.measure = 'snr'

                workflow.connect(snr, 'new_fname',
                                 hist_, 'measure_file')

                drop_percent = pe.Node(
                    util.Function(input_names=['measure_file',
                                               'percent_'],
                                  output_names=['modified_measure_file'],
                                  function=drop_percent_),
                    name='dp_snr_%d' % num_strat)
                drop_percent.inputs.percent_ = 99

                workflow.connect(snr, 'new_fname',
                                 drop_percent, 'measure_file')

                montage_snr = create_montage('montage_snr_%d' % num_strat,
                                             'red_to_blue', 'snr')

                workflow.connect(drop_percent, 'modified_measure_file',
                                 montage_snr, 'inputspec.overlay')
                                                        
                workflow.connect(anat_ref, ref_file,
                                 montage_snr, 'inputspec.underlay')

                strat.update_resource_pool({'qc___snr_a': (montage_snr, 'outputspec.axial_png'),
                                            'qc___snr_s': (montage_snr, 'outputspec.sagittal_png'),
                                            'qc___snr_hist': (hist_, 'hist_path'),
                                            'qc___snr_val': (snr_val, 'snr_storefl')})
                if not 3 in qc_montage_id_a:
                    qc_montage_id_a[3] = 'snr_a'
                    qc_montage_id_s[3] = 'snr_s'
                    qc_hist_id[3] = 'snr_hist'

            except:
                logStandardError('QC', 'unable to get resources for SNR plot', '0051')
                raise

            # make motion parameters plot
            try:
                mov_param, out_file = strat.get_node_from_resource_pool('movement_parameters')

                mov_plot_imports = ['import os', 'import math',
                                    'import numpy as np',
                                    'from matplotlib import pyplot as plt']
                mov_plot = pe.Node(util.Function(input_names=['motion_parameters'],
                                                 output_names=['translation_plot',
                                                               'rotation_plot'],
                                                 function=gen_motion_plt,
                                                 imports=mov_plot_imports),
                                   name='motion_plt_%d' % num_strat)

                workflow.connect(mov_param, out_file, mov_plot, 'motion_parameters')

                strat.update_resource_pool({'qc___movement_trans_plot': (mov_plot, 'translation_plot'),'qc___movement_rot_plot': (mov_plot, 'rotation_plot')})
            
                if not 6 in qc_plot_id:
                    qc_plot_id[6] = 'movement_trans_plot'
            
                if not 7 in qc_plot_id:
                    qc_plot_id[7] = 'movement_rot_plot'

            except:
                logStandardError('QC', 'unable to get resources for Motion Parameters plot', '0052')
                raise

            # make FD plot and volumes removed
            if 'gen_motion_stats' in nodes:
                if 1 in c.runNuisance:
                    try:
                        if c.fdCalc == 'Power':
                            fd, out_file = strat.get_node_from_resource_pool('frame_wise_displacement_power')
                        else:
                            fd, out_file = strat.get_node_from_resource_pool('frame_wise_displacement_jenkinson')

                        fd_plot_imports = ['import os', 'import numpy as np',
                                           'import matplotlib',
                                           'from matplotlib import pyplot']

                        fd_plot = \
                            pe.Node(util.Function(input_names=['arr',
                                                               'measure',
                                                               'ex_vol'],
                                                  output_names=['hist_path'],
                                                  function=gen_plot_png,
                                                  imports=fd_plot_imports),
                                                  name='fd_plot_%d' % num_strat)
                        fd_plot.inputs.measure = 'FD'

                        workflow.connect(fd, out_file,fd_plot, 'arr')

                        if "De-Spiking" in c.runMotionSpike:
                            excluded, out_file_ex = strat.get_node_from_resource_pool('despiking_frames_excluded')
                            workflow.connect(excluded, out_file_ex,
                                             fd_plot, 'ex_vol')
                        elif "Scrubbing" in c.runMotionSpike:
                            excluded, out_file_ex = strat.get_node_from_resource_pool('scrubbing_frames_excluded')
                            workflow.connect(excluded, out_file_ex,
                                             fd_plot, 'ex_vol')

                        strat.update_resource_pool({'qc___fd_plot': (fd_plot, 'hist_path')})

                        if not 8 in qc_plot_id:
                            qc_plot_id[8] = 'fd_plot'

                    except:
                        logStandardError('QC', 'unable to get resources for FD plot', '0053')
                        raise
                            
            # make QC montages for Skull Stripping Visualization
            try:
                anat_underlay, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                skull, out_file_s = strat.get_node_from_resource_pool('anatomical_reorient')

                skull_edge = make_edge(wf_name='skull_edge_%d' % num_strat)

                workflow.connect(skull, out_file_s,
                                 skull_edge, 'inputspec.file_')

                montage_skull = create_montage('montage_skull_%d' % num_strat,
                                               'red', 'skull_vis')

                workflow.connect(skull_edge, 'outputspec.new_fname',
                                 montage_skull, 'inputspec.overlay')
                workflow.connect(anat_underlay, out_file,
                                 montage_skull, 'inputspec.underlay')

                strat.update_resource_pool({'qc___skullstrip_vis_a': (montage_skull, 'outputspec.axial_png'),'qc___skullstrip_vis_s': (montage_skull, 'outputspec.sagittal_png')})
                #skull_edge = #pe.Node(util.Function(input_names=['file_'],output_names=['new_fname#'],function=make_edge),name='skull_edge_%d' % num_strat)

                #workflow.connect(skull, out_file_s,skull_edge, 'file_')

                #workflow.connect(anat_underlay, out_file,montage_skull,'inputspec.underlay')

                # workflow.connect(skull_edge, 'new_fname',montage_skull,'inputspec.overlay')

                # strat.update_resource_pool({'qc___skullstrip_vis_a': (montage_skull, #'outputspec.axial_png'),'qc___skullstrip_vis_s': (montage_skull, #'outputspec.sagittal_png')})

                if not 1 in qc_montage_id_a:
                    qc_montage_id_a[1] = 'skullstrip_vis_a'
                    qc_montage_id_s[1] = 'skullstrip_vis_s'
                
            except:
                    logStandardError('QC', 'Cannot generate QC montages for Skull Stripping: Resources Not Found', '0054')
                    raise

            # make QC montages for mni normalized anatomical image
            try:
                mni_anat_underlay, out_file = strat.get_node_from_resource_pool('mean_functional_in_anat')

                montage_mni_anat = create_montage('montage_mni_anat_%d' % num_strat,
                                                  'red', 'mni_anat')

                workflow.connect(mni_anat_underlay, out_file,
                                 montage_mni_anat, 'inputspec.underlay')

                montage_mni_anat.inputs.inputspec.overlay = p.resource_filename('CPAC','resources/templates/MNI152_Edge_AllTissues.nii.gz')

                strat.update_resource_pool({'qc___mni_normalized_anatomical_a': (montage_mni_anat, 'outputspec.axial_png'),
                                            'qc___mni_normalized_anatomical_s': (montage_mni_anat, 'outputspec.sagittal_png')})

                if not 6 in qc_montage_id_a:
                        qc_montage_id_a[6] = 'mni_normalized_anatomical_a'
                        qc_montage_id_s[6] = 'mni_normalized_anatomical_s'

            except:
                logStandardError('QC', 'Cannot generate QC montages for MNI normalized anatomical: Resources Not Found', '0054')
                raise

            # make QC montages for CSF WM GM
            if 'seg_preproc' in nodes:

                try:
                    anat_underlay, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    csf_overlay, out_file_csf = strat.get_node_from_resource_pool('anatomical_csf_mask')
                    wm_overlay, out_file_wm = strat.get_node_from_resource_pool('anatomical_wm_mask')
                    gm_overlay, out_file_gm = strat.get_node_from_resource_pool('anatomical_gm_mask')

                    montage_csf_gm_wm = create_montage_gm_wm_csf('montage_csf_gm_wm_%d' % num_strat,
                                        'montage_csf_gm_wm')

                    workflow.connect(anat_underlay, out_file,
                                     montage_csf_gm_wm, 'inputspec.underlay')

                    workflow.connect(csf_overlay, out_file_csf,
                                     montage_csf_gm_wm, 'inputspec.overlay_csf')

                    workflow.connect(wm_overlay, out_file_wm,
                                     montage_csf_gm_wm, 'inputspec.overlay_wm')

                    workflow.connect(gm_overlay, out_file_gm,
                                     montage_csf_gm_wm, 'inputspec.overlay_gm')

                    strat.update_resource_pool({'qc___csf_gm_wm_a': (montage_csf_gm_wm, 'outputspec.axial_png'),
                                                'qc___csf_gm_wm_s': (montage_csf_gm_wm, 'outputspec.sagittal_png')})

                    if not 2 in qc_montage_id_a:
                            qc_montage_id_a[2] = 'csf_gm_wm_a'
                            qc_montage_id_s[2] = 'csf_gm_wm_s'

                except:
                    logStandardError('QC', 'Cannot generate QC montages for WM GM CSF masks: Resources Not Found', '0055')
                    raise

            # make QC montage for Mean Functional in T1 with T1 edge
            try:
                anat, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                m_f_a, out_file_mfa = strat.get_node_from_resource_pool('mean_functional_in_anat')

                anat_edge = make_edge(wf_name='anat_edge_%d' % num_strat)

                workflow.connect(anat, out_file, anat_edge, 'inputspec.file_')

                montage_anat = create_montage('montage_anat_%d' % num_strat,
                                              'red', 't1_edge_on_mean_func_in_t1')

                workflow.connect(anat_edge, 'outputspec.new_fname',
                                 montage_anat, 'inputspec.overlay')

                workflow.connect(m_f_a, out_file_mfa,
                                 montage_anat, 'inputspec.underlay')

                # workflow.connect(anat_edge, 'new_fname', montage_anat,
                #                  'inputspec.overlay')
                                    
                strat.update_resource_pool({'qc___mean_func_with_t1_edge_a': (montage_anat, 'outputspec.axial_png'),'qc___mean_func_with_t1_edge_s': (montage_anat, 'outputspec.sagittal_png')})
                                    
                #anat_edge = pe.Node(util.Function(input_names=['file_'],
                #                                                   output_names=['new_fname'],
                #                                                  function=make_edge),
                #                                     name='anat_edge_%d' % num_strat)

                #                workflow.connect(anat, out_file,
                #                                 anat_edge, 'file_')

                #
                #                workflow.connect(m_f_a, out_file_mfa,
                #                                 montage_anat, 'inputspec.underlay')
                #
                #                workflow.connect(anat_edge, 'new_fname',
                #                                 montage_anat, 'inputspec.overlay')
                #
                #                strat.update_resource_pool({'qc___mean_func_with_t1_edge_a': (montage_anat, #'outputspec.axial_png'),
                #                                        'qc___mean_func_with_t1_edge_s': (montage_anat, #'outputspec.sagittal_png')})

                if not 4 in qc_montage_id_a:
                        qc_montage_id_a[4] = 'mean_func_with_t1_edge_a'
                        qc_montage_id_s[4] = 'mean_func_with_t1_edge_s'

            except:
                logStandardError('QC', 'Cannot generate QC montages for Mean Functional in T1 with T1 edge: Resources Not Found', '0056')
                raise

            # make QC montage for Mean Functional in MNI with MNI edge
            try:
                m_f_i, out_file = strat.get_node_from_resource_pool('mean_functional_to_standard')

                montage_mfi = create_montage('montage_mfi_%d' % num_strat,
                                    'red', 'MNI_edge_on_mean_func_mni')   ###

                #   MNI_edge = pe.Node(util.Function(input_names=['file_'],
                #                                      output_names=['new_fname'],
                #                                      function=make_edge),
                #                        name='MNI_edge_%d' % num_strat)
                #   #MNI_edge.inputs.file_ = c.template_brain_only_for_func
                #  workflow.connect(MNI_edge, 'new_fname',
                #                   montage_mfi, 'inputspec.overlay')

                workflow.connect(m_f_i, out_file,
                                 montage_mfi, 'inputspec.underlay')

                montage_mfi.inputs.inputspec.overlay = p.resource_filename('CPAC','resources/templates/MNI152_Edge_AllTissues.nii.gz')


                strat.update_resource_pool({'qc___mean_func_with_mni_edge_a': (montage_mfi, 'outputspec.axial_png'),
                                            'qc___mean_func_with_mni_edge_s': (montage_mfi, 'outputspec.sagittal_png')})

                if not 5 in qc_montage_id_a:
                        qc_montage_id_a[5] = 'mean_func_with_mni_edge_a'
                        qc_montage_id_s[5] = 'mean_func_with_mni_edge_s'


            except:
                logStandardError('QC', 'Cannot generate QC montages for Mean Functional in MNI with MNI edge: Resources Not Found', '0057')
                raise

            # QA pages function
            def QA_montages(measure, idx):
                try:
                    overlay, out_file = strat.get_node_from_resource_pool(measure)

                    drop_percent = pe.MapNode(util.Function(input_names=['measure_file',
                                                                         'percent_'],
                                                            output_names=['modified_measure_file'],
                                                            function=drop_percent_),
                                              name='dp_%s_%d' % (measure, num_strat),
                                              iterfield=['measure_file'])
                    drop_percent.inputs.percent_ = 99.999

                    workflow.connect(overlay, out_file,
                                     drop_percent, 'measure_file')

                    montage = create_montage('montage_%s_%d' % (measure, num_strat),'cyan_to_yellow', measure)
                    montage.inputs.inputspec.underlay = c.template_brain_only_for_func

                    workflow.connect(drop_percent, 'modified_measure_file',
                                     montage, 'inputspec.overlay')

                    histogram = hist.clone('hist_%s_%d' % (measure, num_strat))
                    histogram.inputs.measure = measure

                    workflow.connect(overlay, out_file,
                                     histogram, 'measure_file')

                    strat.update_resource_pool({'qc___%s_a' % measure: (montage, 'outputspec.axial_png'),
                                                'qc___%s_s' % measure: (montage, 'outputspec.sagittal_png'),
                                                'qc___%s_hist' % measure: (histogram, 'hist_path')})

                    if not idx in qc_montage_id_a:
                        qc_montage_id_a[idx] = '%s_a' % measure
                        qc_montage_id_s[idx] = '%s_s' % measure
                        qc_hist_id[idx] = '%s_hist' % measure

                except Exception as e:
                    print "[!] Connection of QA montages workflow for %s " \
                          "has failed.\n" % measure
                    print "Error: %s" % e
                    pass

            # ALFF and f/ALFF QA montages
            if 1 in c.runALFF:
                if 1 in c.runRegisterFuncToMNI:
                    QA_montages('alff_to_standard', 7)
                    QA_montages('falff_to_standard', 8)

                    if c.fwhm != None:
                        QA_montages('alff_to_standard_smooth', 9)
                        QA_montages('falff_to_standard_smooth', 10)

                    if 1 in c.runZScoring:

                        if c.fwhm != None:
                            QA_montages('alff_to_standard_zstd_smooth', 11)
                            QA_montages('falff_to_standard_zstd_smooth', 12)

                        else:
                            QA_montages('alff_to_standard_zstd', 13)
                            QA_montages('falff_to_standard_zstd', 14)

            # ReHo QA montages
            if 1 in c.runReHo:
                if 1 in c.runRegisterFuncToMNI:
                    QA_montages('reho_to_standard', 15)

                    if c.fwhm != None:
                        QA_montages('reho_to_standard_smooth', 16)

                    if 1 in c.runZScoring:

                        if c.fwhm != None:
                            QA_montages('reho_to_standard_zstd_smooth', 17)

                        else:
                            QA_montages('reho_to_standard_fisher_zstd', 18)

            # SCA ROI QA montages
            if (1 in c.runSCA) and (1 in c.runROITimeseries):
                if 1 in c.runRegisterFuncToMNI:
                    QA_montages('sca_roi_to_standard', 19)

                    if c.fwhm != None:
                        QA_montages('sca_roi_to_standard_smooth', 20)

                    if 1 in c.runZScoring:

                        if c.fwhm != None:
                            QA_montages('sca_roi_to_standard_zstd_fisher_smooth', 22)

                        else:
                            QA_montages('sca_roi_to_standard_fisher_zstd', 21)

            # SCA Seed QA montages
            if (1 in c.runSCA) and ("Voxel" in ts_analysis_dict.keys()): #(1 in c.runVoxelTimeseries):

                if 1 in c.runRegisterFuncToMNI:
                    QA_montages('sca_seed_to_standard', 23)

                    if c.fwhm != None:
                        QA_montages('sca_seed_to_standard_smooth', 24)

                    if 1 in c.runZScoring:

                        if c.fwhm != None:
                            QA_montages('sca_seed_to_standard_zstd_fisher_smooth', 26)

                        else:
                            QA_montages('sca_seed_to_standard_fisher_zstd', 25)
                            
            # SCA Multiple Regression
            if "MultReg" in sca_analysis_dict.keys(): #(1 in c.runMultRegSCA) and (1 in c.runROITimeseries):

                if 1 in c.runRegisterFuncToMNI:
                   QA_montages('sca_tempreg_maps_files', 27)
                   QA_montages('sca_tempreg_maps_zstat_files', 28)

                   if c.fwhm != None:
                        QA_montages('sca_tempreg_maps_files_smooth', 29)
                        QA_montages('sca_tempreg_maps_zstat_files_smooth', 30)

            # Dual Regression QA montages
            if ("DualReg" in sca_analysis_dict.keys()) and ("SpatialReg" in ts_analysis_dict.keys()):

                QA_montages('dr_tempreg_maps_files', 31)
                QA_montages('dr_tempreg_maps_zstat_files', 32)

                if 1 in c.runRegisterFuncToMNI:
                    QA_montages('dr_tempreg_maps_files_to_standard', 33)
                    QA_montages('dr_tempreg_maps_zstat_files_to_standard', 34)

                    if c.fwhm != None:
                        QA_montages('dr_tempreg_maps_files_to_standard_smooth', 35)
                        QA_montages('dr_tempreg_maps_zstat_files_to_standard_smooth', 36)

            # VMHC QA montages
            if 1 in c.runVMHC:

                QA_montages('vmhc_raw_score', 37)
                QA_montages('vmhc_fisher_zstd', 38)
                QA_montages('vmhc_fisher_zstd_zstat_map', 39)

            # Network Centrality QA montages
            if 1 in c.runNetworkCentrality:

                QA_montages('centrality_outputs', 40)
                QA_montages('centrality_outputs_zstd', 41)

                if c.fwhm != None:
                    QA_montages('centrality_outputs_smoothed', 42)
                    QA_montages('centrality_outputs_smoothed_zstd', 43)

            num_strat += 1

    logger.info('\n\n' + 'Pipeline building completed.' + '\n\n')

    # end of workflow

    # Run the pipeline only if the user signifies.
    # otherwise, only construct the pipeline (above)
    if run == 1:
        try:
            workflow.write_graph(graph2use='orig')
        except:
            pass

        # this section creates names for the different branched strategies.
        # it identifies where the pipeline has forked and then appends the
        # name of the forked nodes to the branch name in the output directory
        renamedStrats = []
        forkPoints = []
        forkPointsDict = {}

        def is_number(s):
            # function which returns boolean checking if a character
            # is a number or not
            try:
                float(s)
                return True
            except ValueError:
                return False

        for strat in strat_list:

            # load list of nodes in this one particular
            # strat into the list "nodeList"
            nodeList = strat.name
            renamedNodesList = []

            # strip the _n (n being the strat number) from
            # each node name and return to a list
            for node in nodeList:

                renamedNode = node
                lastNodeChar = node[len(node) - 1]

                while lastNodeChar == '_' or lastNodeChar == '-' or is_number(
                        lastNodeChar):
                    # make 'renamedNode' the node name with the last character
                    # stripped off, continue this until the _# at the end
                    # of it is gone - does it this way instead of just cutting
                    # off the last two characters in case of a large amount of
                    # strats which can reach double digits
                    renamedNode = renamedNode[:-1]
                    lastNodeChar = renamedNode[len(renamedNode) - 1]

                renamedNodesList.append(renamedNode)

            renamedStrats.append(renamedNodesList)

        # here, renamedStrats is a list containing each strat (forks)
        for strat in renamedStrats:

            tmpForkPoint = []

            # here, 'strat' is a list of node names within one of the forks
            for nodeName in strat:

                # compare each strat against the first one in the strat list,
                # and if any node names in the new strat are not present in
                # the 'original' one, then append to a list of 'fork points'
                for renamedStratNodes in renamedStrats:

                    if nodeName not in renamedStratNodes and \
                                    nodeName not in tmpForkPoint:
                        tmpForkPoint.append(nodeName)

            forkPoints.append(tmpForkPoint)

        # forkPoints is a list of lists, each list containing node names of
        # nodes run in that strat/fork that are unique to that strat/fork

        forkNames = []

        # here 'forkPoint' is an individual strat with its unique nodes
        for forkPoint in forkPoints:

            forkName = ''
            forklabel = ''
            for fork in forkPoint:
                forklabel = ''
                if 'ants' in fork:
                    forklabel = 'ANTS'
                if 'fnirt' in fork:
                    forklabel = 'FNIRT'
                if 'automask' in fork:
                    forklabel = '3dAutoMask(func)'
                if 'bet' in fork:
                    forklabel = 'BET(func)'
                if 'epi_distcorr' in fork:
                    forklabel = 'dist_corr'
                if 'bbreg' in fork:
                    forklabel = 'bbreg'
                if 'frequency' in fork:
                    forklabel = 'freq-filter'
                if 'nuisance_with_despiking' in fork:
                    forklabel = 'nuisance_with_despiking'
                elif 'nuisance_no_despiking' in fork:
                    forklabel = 'nuisance_no_despiking'
                elif 'nuisance' in fork:
                    forklabel = 'nuisance'
                if 'median' in fork:
                    forklabel = 'median'
                if 'friston' in fork:
                    forklabel = 'friston'
                if 'motion_stats' in fork:
                    forklabel = 'motion'
                if 'scrubbing' in fork:
                    forklabel = 'scrub'
                if 'slice' in fork:
                    forklabel = 'slice'
                if 'AFNI' in fork:
                    forklabel = 'afni(anat)'
                if 'BET' in fork:
                    forklabel = 'BET(anat)'
                    #if 'epi_distcorr' in fork:
                    #forklabel = 'epi_distcorr'
                if forklabel not in forkName:
                    forkName = forkName + '__' + forklabel

            forkNames.append(forkName)

        # match each strat_list with fork point list
        # this is for the datasink
        for x in range(len(strat_list)):
            forkPointsDict[strat_list[x]] = forkNames[x]

        '''
        Datasink
        '''
        import networkx as nx
        num_strat = 0
        sink_idx = 0
        pip_ids = []

        wf_names = []
        scan_ids = ['scan_anat']

        try:
            for scanID in sub_dict['func']:
                scan_ids.append('scan_' + str(scanID))
        except KeyError:
            for scanID in sub_dict['rest']:
                scan_ids.append('scan_' + str(scanID))

        pipes = []
        origStrat = 0

        for strat in strat_list:
            rp = strat.get_resource_pool()

            # build helper dictionary to assist with a clean strategy label
            # for symlinks
            strategy_tag_helper_symlinks = {}

            if any('scrubbing' in name for name in strat.get_name()):
                strategy_tag_helper_symlinks['_threshold'] = 1
            else:
                strategy_tag_helper_symlinks['_threshold'] = 0

            if any('seg_preproc' in name for name in strat.get_name()):
                strategy_tag_helper_symlinks['_csf_threshold'] = 1
                strategy_tag_helper_symlinks['_wm_threshold'] = 1
                strategy_tag_helper_symlinks['_gm_threshold'] = 1
            else:
                strategy_tag_helper_symlinks['_csf_threshold'] = 0
                strategy_tag_helper_symlinks['_wm_threshold'] = 0
                strategy_tag_helper_symlinks['_gm_threshold'] = 0

            if any('median_angle_corr' in name for name in strat.get_name()):
                strategy_tag_helper_symlinks['_target_angle_deg'] = 1
            else:
                strategy_tag_helper_symlinks['_target_angle_deg'] = 0

            if any('nuisance' in name for name in strat.get_name()):
                strategy_tag_helper_symlinks['nuisance'] = 1
            else:
                strategy_tag_helper_symlinks['nuisance'] = 0

            strat_tag = ""
            hash_val = 0

            for name in strat.get_name():
                import re
                extra_string = re.search('_\d+', name).group(0)

                if extra_string:
                    name = name.split(extra_string)[0]

                if workflow_bit_id.get(name) is not None:
                    strat_tag += name + '_'
                    print name, ' --- ', 2 ** workflow_bit_id[name]
                    hash_val += 2 ** workflow_bit_id[name]

            if p_name is None or p_name == 'None':
                if forkPointsDict[strat]:
                    pipeline_id = c.pipelineName + forkPointsDict[strat]
                else:
                    pipeline_id = c.pipelineName
                    # if running multiple pipelines with gui, need to change
                    # this in future
                    p_name = None
            else:
                if forkPointsDict[strat]:
                    pipeline_id = c.pipelineName + forkPointsDict[strat]
                else:
                    pipeline_id = p_name
                    # if running multiple pipelines with gui, need to change
                    # this in future
                    p_name = None

            pip_ids.append(pipeline_id)
            wf_names.append(strat.get_name())

            # Extract credentials path for output if it exists
            s3_str = 's3://'
            try:
                # Get path to creds file
                creds_path = None
                if c.awsOutputBucketCredentials:
                    creds_path = str(c.awsOutputBucketCredentials)
                    creds_path = os.path.abspath(creds_path)
                    
                # Test for s3 write access
                s3_write_access = \
                    aws_utils.test_bucket_access(creds_path,
                                                 c.outputDirectory)
                if not s3_write_access:
                    raise Exception('Not able to write to bucket!')
            except Exception as exc:
                if c.outputDirectory.lower().startswith(s3_str):
                    err_msg = 'There was an error processing credentials or ' \
                              'accessing the S3 bucket. Check and try again.\n' \
                              'Error: %s' % exc
                    raise Exception(err_msg)
            try:
                encrypt_data = bool(c.s3Encryption[0])
            except Exception as exc:
                encrypt_data = False

            # TODO: remove this once forking for despiking/scrubbing is
            # TODO: modified at the gen motion params level
            # ensure X_frames_included/excluded only gets sent to output dir
            # for appropriate strats
            nodes = getNodeList(strat)
            if "nuisance_with_despiking" not in nodes:
                if "despiking_frames_included" in rp.keys():
                    del rp["despiking_frames_included"]
                if "despiking_frames_excluded" in rp.keys():
                    del rp["despiking_frames_excluded"]
            if "scrubbing" not in nodes:
                if "scrubbing_frames_included" in rp.keys():
                    del rp["scrubbing_frames_included"]
                if "scrubbing_frames_excluded" in rp.keys():
                    del rp["scrubbing_frames_excluded"]

            for key in sorted(rp.keys()):

                if key not in override_optional:

                    if 1 not in c.write_func_outputs:
                        if key in extra_functional_outputs:
                            continue

                    if 1 not in c.write_debugging_outputs:
                        if key in debugging_outputs:
                            continue

                    if 0 not in c.runRegisterFuncToMNI:
                        if key in outputs_native_nonsmooth or \
                            key in outputs_native_nonsmooth_mult or \
                                key in outputs_native_smooth:
                            continue

                    if 0 not in c.runZScoring:
                        # write out only the z-scored outputs
                        if key in outputs_template_raw or \
                                key in outputs_template_raw_mult:
                            continue

                    if 0 not in c.run_smoothing:
                        # write out only the smoothed outputs
                        if key in outputs_native_nonsmooth or \
                                key in outputs_template_nonsmooth or \
                                    key in outputs_native_nonsmooth_mult or \
                                        key in outputs_template_nonsmooth_mult:
                            continue

                ds = pe.Node(nio.DataSink(), name='sinker_%d' % sink_idx)

                # Write QC outputs to log directory
                if 'qc' in key.lower():
                    ds.inputs.base_directory = c.logDirectory
                else:
                    ds.inputs.base_directory = c.outputDirectory
                    # For each pipeline ID, generate the QC pages
                    #  for pip_id in pip_ids:
                        # Define pipeline-level logging for QC
                        #    pipeline_out_base = os.path.join(c.logDirectory, 'pipeline_%s' % pip_id)
                        #qc_output_folder = os.path.join(pipeline_out_base, subject_id, 'qc_files_here')
                        #For each subject, create a QC index.html page
                        #make_QC_html_pages(qc_output_folder)

                ds.inputs.creds_path = creds_path
                ds.inputs.encrypt_bucket_keys = encrypt_data
                ds.inputs.container = os.path.join('pipeline_%s' % pipeline_id, subject_id)
                ds.inputs.regexp_substitutions = [(r"/_sca_roi(.)*[/]", '/'),
                                                  (r"/_smooth_centrality_(\d)+[/]", '/'),
                                                  (r"/_z_score(\d)+[/]", "/"),
                                                  (r"/_dr_tempreg_maps_zstat_files_smooth_(\d)+[/]", "/"),
                                                  (r"/_sca_tempreg_maps_zstat_files_smooth_(\d)+[/]", "/"),
                                                  (r"/qc___", '/qc/')]

                node, out_file = rp[key]
                workflow.connect(node, out_file, ds, key)

                link_node = pe.Node(interface=util.Function(
                                        input_names=['in_file', 'strategies',
                                                     'subject_id', 'pipeline_id', 'helper',
                                                     'create_sym_links'],
                                        output_names=[],
                                        function=process_outputs),
                                    name='process_outputs_%d' % sink_idx)

                link_node.inputs.strategies = strategies
                link_node.inputs.subject_id = subject_id
                link_node.inputs.pipeline_id = 'pipeline_%s' % pipeline_id
                link_node.inputs.helper = dict(strategy_tag_helper_symlinks)

                if 1 in c.runSymbolicLinks:
                    link_node.inputs.create_sym_links = True
                else:
                    link_node.inputs.create_sym_links = False

                workflow.connect(ds, 'out_file', link_node, 'in_file')

                sink_idx += 1
                logger.info('sink index: %s' % sink_idx)

            d_name = os.path.join(c.logDirectory, ds.inputs.container)

            if not os.path.exists(d_name):
                os.makedirs(d_name)

            try:
                G = nx.DiGraph()
                strat_name = strat.get_name()
                G.add_edges_from([(strat_name[s], strat_name[s + 1]) for s in
                                  range(len(strat_name) - 1)])
                dotfilename = os.path.join(d_name, 'strategy.dot')
                nx.drawing.nx_pydot.write_dot(G, dotfilename)
                format_dot(dotfilename, 'png')
            except:
                logStandardWarning('Datasink',
                                   'Cannot Create the strategy and pipeline '
                                   'graph, dot or/and pygraphviz is not '
                                   'installed')
                pass

            logger.info('%s*' % d_name)
            num_strat += 1

            pipes.append(pipeline_id)

        # creates the HTML files used to represent the logging-based status
        create_log_template(pip_ids, wf_names, scan_ids, subject_id, log_dir)

        logger.info("\n\nStrategy forks: {0}\n\n".format(str(set(pipes))))

        pipeline_start_date = strftime("%Y-%m-%d")
        pipeline_start_datetime = strftime("%Y-%m-%d %H:%M:%S")
        pipeline_starttime_string = pipeline_start_datetime.replace(' ', '_')
        pipeline_starttime_string = pipeline_starttime_string.replace(':',
                                                                      '-')

        strat_no = 0

        subject_info['resource_pool'] = []

        for strat in strat_list:
            strat_label = 'strat_%d' % strat_no
            subject_info[strat_label] = strat.get_name()
            subject_info['resource_pool'].append(strat.get_resource_pool())
            strat_no += 1

        subject_info['status'] = 'Running'

        '''
        subject_info_pickle = open(os.getcwd() + '/subject_info.p', 'wb')

        pickle.dump(subject_info, subject_info_pickle)

        subject_info_pickle.close()
        '''

        # TODO:set memory and num_threads of critical nodes if running
        # MultiProcPlugin

        # Create callback logger
        import logging as cb_logging
        cb_log_filename = os.path.join(log_dir,
                                       'callback_%s.log' % sub_dict[
                                           'subject_id'])

        try:
            if not os.path.exists(os.path.dirname(cb_log_filename)):
                os.makedirs(os.path.dirname(cb_log_filename))

        except IOError:
            pass

        # Add handler to callback log file
        cb_logger = cb_logging.getLogger('callback')
        cb_logger.setLevel(cb_logging.DEBUG)
        handler = cb_logging.FileHandler(cb_log_filename)
        cb_logger.addHandler(handler)

        # Add status callback function that writes in callback log
        if nipype.__version__ not in ('0.13.1', '0.14.0'):
            err_msg = "This version of nipype may not be compatible with " \
                      "CPAC v%s, please install version 0.14.0\n" \
                       % (CPAC.__version__)
            logger.error(err_msg)
        else:
            if nipype.__version__ == '0.13.1':
                from nipype.pipeline.plugins.callback_log import log_nodes_cb
                plugin_args['status_callback'] = log_nodes_cb
            else:
                from nipype.utils.profiler import log_nodes_cb
                plugin_args['status_callback'] = log_nodes_cb

        # Actually run the pipeline now, for the current subject
        workflow.run(plugin=plugin, plugin_args=plugin_args)

        # Dump subject info pickle file to subject log dir
        subject_info['status'] = 'Completed'
        subject_info_pickle = open(
            os.path.join(log_dir, 'subject_info_%s.p' % subject_id), 'wb')
        pickle.dump(subject_info, subject_info_pickle)
        subject_info_pickle.close()

        '''
        # Actually run the pipeline now
        try:

            workflow.run(plugin='MultiProc', plugin_args={'n_procs': c.numCoresPerSubject})

        except:

            crashString = "\n\n" + "ERROR: CPAC run stopped prematurely with an error - see above.\n" + ("pipeline configuration- %s \n" % c.pipelineName) + \
            ("subject workflow- %s \n\n" % wfname) + ("Elapsed run time before crash (minutes): %s \n\n" % ((time.time() - pipeline_start_time)/60)) + \
            ("Timing information saved in %s/cpac_timing_%s_%s.txt \n" % (c.outputDirectory, c.pipelineName, pipeline_starttime_string)) + \
            ("System time of start:      %s \n" % pipeline_start_datetime) + ("System time of crash: %s" % strftime("%Y-%m-%d %H:%M:%S")) + "\n\n"

            logger.info(crashString)

            print >>timing, "ERROR: CPAC run stopped prematurely with an error."
            print >>timing, "Pipeline configuration: %s" % c.pipelineName
            print >>timing, "Subject workflow: %s" % wfname
            print >>timing, "\n" + "Elapsed run time before crash (minutes): ", ((time.time() - pipeline_start_time)/60)
            print >>timing, "System time of crash: ", strftime("%Y-%m-%d %H:%M:%S")
            print >>timing, "\n\n"

            timing.close()

            raise Exception
        '''

        '''
        try:

            workflow.run(plugin='MultiProc', plugin_args={'n_procs': c.numCoresPerSubject})

        except Exception as e:

            print "Error: CPAC Pipeline has failed."
            print ""
            print e
            print type(e)
            ###raise Exception
        '''

        # subject_dir = os.path.join(c.outputDirectory, 'pipeline_' + pipeline_id, subject_id)
        # create_output_mean_csv(subject_dir)

        for count, scanID in enumerate(pip_ids):
            for scan in scan_ids:
                create_log_node(None, None, count, scan).run()
        if 1 in c.generateQualityControlImages:
            for pip_id in pip_ids:
                try:
                    pipeline_base = os.path.join(c.outputDirectory, 'pipeline_%s' % pip_id)
                    qc_output_folder = os.path.join(pipeline_base, subject_id, 'qc_files_here')
                    generateQCPages(qc_output_folder,qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id)
            #create_all_qc.run(pipeline_base)
                except Exception as e:
                    print "Error: The QC function page generation is not running"
                    print ""
                    print e
                    print type(e)
                    raise Exception
                    


            # Generate the QC pages -- this function isn't even running, because there is noparameter for qc_montage_id_a/qc_montage_id_s/qc_plot_id,qc_hist_id
                #two methods can be done here:
                #i) group all the qc_montage_ids in the resource pool or
                #add a loop for all the files in the qc output folder, generate the html pages using the same functions, but with different parameters
                # generateQCPages(qc_output_folder, qc_montage_id_a,
                #  qc_montage_id_s, qc_plot_id, qc_hist_id)
                # Automatically generate QC index page
                #create_all_qc.run(pipeline_out_base)
        

        # pipeline timing code starts here

        # have this check in case the user runs cpac_runner from terminal and
        # the timing parameter list is not supplied as usual by the GUI
        if pipeline_timing_info != None:

            # pipeline_timing_info list:
            #  [0] - unique pipeline ID
            #  [1] - pipeline start time stamp (first click of 'run' from GUI)
            #  [2] - number of subjects in subject list
            unique_pipeline_id = pipeline_timing_info[0]
            pipeline_start_stamp = pipeline_timing_info[1]
            num_subjects = pipeline_timing_info[2]

            # elapsed time data list:
            #  [0] - elapsed time in minutes
            elapsed_time_data = []

            elapsed_time_data.append(
                int(((time.time() - pipeline_start_time) / 60)))

            # elapsedTimeBin list:
            #  [0] - cumulative elapsed time (minutes) across all subjects
            #  [1] - number of times the elapsed time has been appended
            #        (effectively a measure of how many subjects have run)

            # needs to happen:
            # write more doc for all this
            # warning in .csv that some runs may be partial
            # code to delete .tmp file

            timing_temp_file_path = os.path.join(c.logDirectory,
                                                 '%s_pipeline_timing.tmp' % unique_pipeline_id)

            if not os.path.isfile(timing_temp_file_path):
                elapsedTimeBin = []
                elapsedTimeBin.append(0)
                elapsedTimeBin.append(0)

                with open(timing_temp_file_path, 'wb') as handle:
                    pickle.dump(elapsedTimeBin, handle)

            with open(timing_temp_file_path, 'rb') as handle:
                elapsedTimeBin = pickle.loads(handle.read())

            elapsedTimeBin[0] = elapsedTimeBin[0] + elapsed_time_data[0]
            elapsedTimeBin[1] = elapsedTimeBin[1] + 1

            with open(timing_temp_file_path, 'wb') as handle:
                pickle.dump(elapsedTimeBin, handle)

            # this happens once the last subject has finished running!
            if elapsedTimeBin[1] == num_subjects:

                pipelineTimeDict = {}
                pipelineTimeDict['Pipeline'] = c.pipelineName
                pipelineTimeDict[
                    'Cores_Per_Subject'] = c.maxCoresPerParticipant
                pipelineTimeDict[
                    'Simultaneous_Subjects'] = c.numParticipantsAtOnce
                pipelineTimeDict['Number_of_Subjects'] = num_subjects
                pipelineTimeDict['Start_Time'] = pipeline_start_stamp
                pipelineTimeDict['End_Time'] = strftime("%Y-%m-%d_%H:%M:%S")
                pipelineTimeDict['Elapsed_Time_(minutes)'] = elapsedTimeBin[0]
                pipelineTimeDict['Status'] = 'Complete'

                gpaTimeFields = ['Pipeline', 'Cores_Per_Subject',
                                 'Simultaneous_Subjects',
                                 'Number_of_Subjects', 'Start_Time',
                                 'End_Time', 'Elapsed_Time_(minutes)',
                                 'Status']
                timeHeader = dict((n, n) for n in gpaTimeFields)

                timeCSV = open(os.path.join(c.logDirectory,
                                            'cpac_individual_timing_%s.csv' % c.pipelineName),
                               'a')
                readTimeCSV = open(os.path.join(c.logDirectory,
                                                'cpac_individual_timing_%s.csv' % c.pipelineName),
                                   'rb')
                timeWriter = csv.DictWriter(timeCSV, fieldnames=gpaTimeFields)
                timeReader = csv.DictReader(readTimeCSV)

                headerExists = False
                for line in timeReader:
                    if 'Start_Time' in line:
                        headerExists = True

                if headerExists == False:
                    timeWriter.writerow(timeHeader)

                timeWriter.writerow(pipelineTimeDict)
                timeCSV.close()
                readTimeCSV.close()

                # remove the temp timing file now that it is no longer needed
                os.remove(timing_temp_file_path)

        # Upload logs to s3 if s3_str in output directory
        if c.outputDirectory.lower().startswith(s3_str):
            try:
                # Store logs in s3 output director/logs/...
                s3_log_dir = c.outputDirectory + '/logs/' + \
                             os.path.basename(log_dir)
                bucket_name = c.outputDirectory.split('/')[2]
                bucket = fetch_creds.return_bucket(creds_path, bucket_name)
                # Collect local log files
                local_log_files = []
                for root, dirs, files in os.walk(log_dir):
                    local_log_files.extend([os.path.join(root, fil) \
                                            for fil in files])
                    # Form destination keys
                s3_log_files = [loc.replace(log_dir, s3_log_dir) \
                                for loc in local_log_files]
                # Upload logs
                aws_utils.s3_upload(bucket, (local_log_files, s3_log_files),
                                    encrypt=encrypt_data)
                # Delete local log files
                for log_f in local_log_files:
                    os.remove(log_f)
            except Exception as exc:
                err_msg = 'Unable to upload CPAC log files in: %s.\nError: %s' \
                          % (log_dir, exc)
                logger.error(err_msg)

        # Remove working directory when done
        sub_w_path = os.path.join(c.workingDirectory, wfname)
        if c.removeWorkingDir:
            try:
                if os.path.exists(sub_w_path):
                    import shutil
                    logger.info("removing dir -> %s" % sub_w_path)
                    shutil.rmtree(sub_w_path)
            except:
                logStandardWarning('Datasink', (
                'Couldn\'t remove subjects %s working directory' % wfname))
                pass

        endString = (
                    "End of subject workflow %s \n\n" % wfname) + "CPAC run complete:\n" + (
                    "pipeline configuration- %s \n" % c.pipelineName) + \
                    ("subject workflow- %s \n\n" % wfname) + (
                    "Elapsed run time (minutes): %s \n\n" % (
                    (time.time() - pipeline_start_time) / 60)) + \
                    (
                    "Timing information saved in %s/cpac_individual_timing_%s.csv \n" % (
                    c.logDirectory, c.pipelineName)) + \
                    (
                    "System time of start:      %s \n" % pipeline_start_datetime) + (
                    "System time of completion: %s" % strftime(
                        "%Y-%m-%d %H:%M:%S"))

        logger.info(endString)

    return workflow


# Run the prep_workflow function with specific arguments
def run(config, subject_list_file, indx, strategies, p_name=None, \
        plugin=None, plugin_args=None):
    '''
    Function to build and execute the complete workflow

    Parameters
    ----------
    config: string
        filepath to a C-PAC config file
    subject_list_file : string
        filepath to a C-PAC subject list file
    indx : integer
        index of the subject in the subject list to run
    strategies : string
        filepath to a C-PAC strategies file
    maskSpecificationFile : string
        filepath to the mask-specification file
    roiSpecificationFile : string
        filepath to the roi-specification file
    templateSpecificationFile : string
        filepath to the template-specification file
    p_name : string (optional)
        name of the pipeline configuration
    plugin : string (optional)
        name of the plugin  used to schedule nodes
    plugin_args : dict (optional)
        arguments of plugin
    '''

    # Import packages
    import commands
    commands.getoutput('source ~/.bashrc')
    import yaml

    # Import configuration file
    c = Configuration(yaml.load(open(os.path.realpath(config), 'r')))

    # Try and load in the subject list
    try:
        sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        raise Exception(
            "Subject list is not in proper YAML format. Please check your file")

    # Grab the subject of interest
    sub_dict = sublist[int(indx) - 1]
    sub_id = sub_dict['subject_id']

    try:
        # Build and run the pipeline
        prep_workflow(sub_dict, c, pickle.load(open(strategies, 'r')), 1,
                      p_name, plugin=plugin, plugin_args=plugin_args)
    except Exception as e:
        print 'Could not complete cpac run for subject: %s!' % sub_id
        print 'Error: %s' % e
