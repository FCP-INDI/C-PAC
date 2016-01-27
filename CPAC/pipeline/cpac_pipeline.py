# CPAC/pipeline/cpac_pipeline.py
# 

'''
This module prepares and executes the main C-PAC workflow
'''

# Import packages
import os
import time
from time import strftime
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
from nipype.interfaces.afni import preprocess
from   nipype.pipeline.utils import format_dot
import nipype.interfaces.ants as ants
import nipype.interfaces.c3 as c3
from nipype import config
from nipype import logging
from CPAC import network_centrality
from CPAC.network_centrality.utils import merge_lists
logger = logging.getLogger('workflow')
import pkg_resources as p
import CPAC
from CPAC.anat_preproc.anat_preproc import create_anat_preproc
from CPAC.func_preproc.func_preproc import create_func_preproc, create_wf_edit_func
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
### no create_log_template here, move in CPAC/utils/utils.py
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
                             get_fisher_zscore, dbg_file_lineno
from CPAC.vmhc.vmhc import create_vmhc
from CPAC.reho.reho import create_reho
from CPAC.alff.alff import create_alff
from CPAC.sca.sca import create_sca, create_temporal_reg
import zlib
import linecache
import csv
import pickle


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

    

def prep_workflow(sub_dict, c, strategies, run, pipeline_timing_info=None, p_name=None):


    """""""""""""""""""""""""""""""""""""""""""""""""""
     SETUP
    """""""""""""""""""""""""""""""""""""""""""""""""""

    '''
    preliminaries
    '''

    # Start timing here
    pipeline_start_time = time.time()
    # at end of workflow, take timestamp again, take time elapsed and check
    # tempfile add time to time data structure inside tempfile, and increment
    # number of subjects


    cores_msg = 'VERSION: CPAC %s' % CPAC.__version__


    # perhaps in future allow user to set threads maximum
    # this is for centrality mostly    
    # import mkl
    numThreads = '1'

    os.environ['OMP_NUM_THREADS'] = numThreads
    os.environ['MKL_NUM_THREADS'] = numThreads
    os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = str(c.num_ants_threads)

    # calculate maximum potential use of cores according to current pipeline
    # configuration
    if 'ANTS' in c.regOption:
        max_core_usage = int(c.numCoresPerSubject) * \
                             int(c.numSubjectsAtOnce) * int(numThreads) * \
                             int(c.num_ants_threads)
    else:
        max_core_usage = int(c.numCoresPerSubject) * \
                             int(c.numSubjectsAtOnce) * int(numThreads)

    cores_msg = cores_msg + '\n\nSetting number of cores per subject to %s\n' \
                            % c.numCoresPerSubject

    cores_msg = cores_msg + 'Setting number of subjects at once to %s\n' \
                            % c.numSubjectsAtOnce

    cores_msg = cores_msg + 'Setting OMP_NUM_THREADS to %s\n' % numThreads
    cores_msg = cores_msg + 'Setting MKL_NUM_THREADS to %s\n' % numThreads

    if 'ANTS' in c.regOption:
        cores_msg = cores_msg + 'Setting ANTS/ITK thread usage to %d\n\n' \
                    % c.num_ants_threads

    cores_msg = cores_msg + 'Maximum potential number of cores that might ' \
                'be used during this run: %d\n\n' % max_core_usage

    cores_msg = cores_msg + 'If that\'s more cores than you have, better ' \
                'fix that quick! Hint: This can be changed via the settings '\
                '\'Number of Cores Per Subject\', \'Number of Subjects ' \
                'to Run Simultaneously\', and \'Number of Cores for ' \
                'Anatomical Registration (ANTS only)\' in the pipeline ' \
                'configuration editor under the tab \'Computer Settings' \
                '\'.\n\n'

    logger.info(cores_msg)


    qc_montage_id_a = {}
    qc_montage_id_s = {}
    qc_plot_id = {}
    qc_hist_id = {}
    if sub_dict['unique_id']:
        subject_id = sub_dict['subject_id'] + "_" + sub_dict['unique_id']
    else:
        subject_id = sub_dict['subject_id']
        
    log_dir = os.path.join(c.outputDirectory, 'logs', subject_id)

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
        if (len(wrong_filepath_list) == 1) and (wrong_filepath_list[0][0] == "dilated_symmetric_brain_mask"):
            pass
        else:
            raise Exception
            
            
    # this checks to make sure the user has appropriately installed and
    # configured necessary tools (i.e. AFNI, FSL, ANTS..)
    
    missing_install = []
    
    if os.system("3dcalc >/dev/null") == 32512:
        missing_install.append("AFNI")
          
    if os.system("fslmaths >/dev/null") == 32512:
        missing_install.append("FSL")
    
    if "ANTS" in c.regOption:
    
        if os.system("c3d_affine_tool >/dev/null") == 32512:
            missing_install.append("C3D")
    
        if os.system("antsRegistration >/dev/null") == 32512:
            missing_install.append("ANTS")
            
            
    if len(missing_install) > 0:
   
        missing_string = ""
        
        for string in missing_install:
            missing_string = missing_string + string + "\n"
   
        err = "\n\n[!] CPAC says: It appears the following software " \
              "packages are not installed or configured properly:\n\n%s\n" \
              "Consult the CPAC Installation Guide for instructions.\n\n" \
              % missing_string
        raise Exception(err)
    
                
    

    '''
    workflow preliminary setup
    '''

    wfname = 'resting_preproc_' + str(subject_id)
    workflow = pe.Workflow(name=wfname)
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {'hash_method': 'timestamp', 'crashdump_dir': os.path.abspath(c.crashLogDirectory)}
    config.update_config({'logging': {'log_directory': log_dir, 'log_to_file': True}})
    logging.update_logging(config)


    if c.reGenerateOutputs is True:

        import commands
        cmd = "find %s -name \'*sink*\' -exec rm -rf {} \\;" % os.path.join(c.workingDirectory, wfname)
        logger.info(cmd)
        commands.getoutput(cmd)
        cmd = "find %s -name \'*link*\' -exec rm -rf {} \\;" % os.path.join(c.workingDirectory, wfname)
        logger.info(cmd)
        commands.getoutput(cmd)
        cmd = "find %s -name \'*log*\' -exec rm -rf {} \\;" % os.path.join(c.workingDirectory, wfname)
        logger.info(cmd)
        commands.getoutput(cmd)


    def create_log_node(wflow, output, indx, scan_id = None):
        #call logging workflow

        if wflow: 
            log_wf = create_log(wf_name = 'log_%s' %wflow.name)
            log_wf.inputs.inputspec.workflow = wflow.name
            log_wf.inputs.inputspec.index = indx
            log_wf.inputs.inputspec.log_dir = log_dir
            workflow.connect(wflow, output, log_wf, 'inputspec.inputs')
        else:
            log_wf = create_log(wf_name = 'log_done_%s'%scan_id, scan_id= scan_id)
            log_wf.inputs.inputspec.workflow = 'DONE'
            log_wf.inputs.inputspec.index = indx
            log_wf.inputs.inputspec.log_dir = log_dir
            log_wf.inputs.inputspec.inputs = log_dir
            return log_wf
        
    def logStandardError(sectionName, errLine, errNum):
        
        logger.info("\n\n" + 'ERROR: %s - %s' % (sectionName, errLine) + "\n\n" + \
                    "Error name: cpac_pipeline_%s" % (errNum) + "\n\n")
        
    def logConnectionError(workflow_name, numStrat, resourcePool, errNum):
        
        logger.info("\n\n" + 'ERROR: Invalid Connection: %s: %s, resource_pool: %s' \
                    % (workflow_name, numStrat, resourcePool) + "\n\n" + "Error name: cpac_pipeline_%s" % (errNum) + \
                    "\n\n" + "This is a pipeline creation error - the workflows have not started yet." + "\n\n")
        
    def logStandardWarning(sectionName, warnLine):
        
        logger.info("\n\n" + 'WARNING: %s - %s' % (sectionName, warnLine) + "\n\n")
        
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

    strat_initial = None

    for gather_anat in c.runAnatomicalDataGathering:
        strat_initial = strategy()

        if gather_anat == 1:
            flow = create_anat_datasource()
            flow.inputs.inputnode.subject = subject_id
            flow.inputs.inputnode.anat = sub_dict['anat']

            anat_flow = flow.clone('anat_gather_%d' % num_strat)

            strat_initial.set_leaf_properties(anat_flow, 'outputspec.anat')

        num_strat += 1

        strat_list.append(strat_initial)



    '''
    Inserting Anatomical Preprocessing workflow
    '''
    new_strat_list = []
    num_strat = 0

    if 1 in c.runAnatomicalPreprocessing:

        workflow_bit_id['anat_preproc'] = workflow_counter

        for strat in strat_list:
            # create a new node, Remember to change its name!
            anat_preproc = create_anat_preproc(already_skullstripped).clone('anat_preproc_%d' % num_strat)

            try:
                # connect the new node to the previous leaf
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file, anat_preproc, 'inputspec.anat')

            except:
                logConnectionError('Anatomical Preprocessing No valid Previous for strat', num_strat, strat.get_resource_pool(), '0001')
                num_strat += 1
                continue

            if 0 in c.runAnatomicalPreprocessing:
                # we are forking so create a new node
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(anat_preproc.name)


            strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
            # add stuff to resource pool if we need it

            strat.update_resource_pool({'anatomical_brain':(anat_preproc, 'outputspec.brain')})
            strat.update_resource_pool({'anatomical_reorient':(anat_preproc, 'outputspec.reorient')})
            
            #write to log
            create_log_node(anat_preproc, 'outputspec.brain', num_strat)

            num_strat += 1

    strat_list += new_strat_list



    '''
    T1 -> Template, Non-linear registration (FNIRT or ANTS)
    '''
    
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    
    # either run FSL anatomical-to-MNI registration, or...
    
    if 1 in c.runRegistrationPreprocessing:

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

                fnirt_reg_anat_mni = create_nonlinear_register('anat_mni_fnirt_register_%d' % num_strat)

                try:
                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_mni, 'inputspec.input_brain')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_mni, 'inputspec.input_skull')

                    # pass the reference files                
                    fnirt_reg_anat_mni.inputs.inputspec.reference_brain = c.template_brain_only_for_anat
                    fnirt_reg_anat_mni.inputs.inputspec.reference_skull = c.template_skull_for_anat
                    fnirt_reg_anat_mni.inputs.inputspec.ref_mask = c.ref_mask

                    # assign the FSL FNIRT config file specified in pipeline config.yml
                    fnirt_reg_anat_mni.inputs.inputspec.fnirt_config = c.fnirtConfig
                

                except:
                    logConnectionError('Anatomical Registration (FSL)', num_strat, strat.get_resource_pool(), '0002')
                    raise

                if (0 in c.runRegistrationPreprocessing) or ('ANTS' in c.regOption):
                    tmp = strategy()
                    tmp.resource_pool = dict(strat.resource_pool)
                    tmp.leaf_node = (strat.leaf_node)
                    tmp.leaf_out_file = str(strat.leaf_out_file)
                    tmp.name = list(strat.name)
                    strat = tmp
                    new_strat_list.append(strat)

                strat.append_name(fnirt_reg_anat_mni.name)
                strat.set_leaf_properties(fnirt_reg_anat_mni, 'outputspec.output_brain')

                strat.update_resource_pool({'anatomical_to_mni_linear_xfm':(fnirt_reg_anat_mni, 'outputspec.linear_xfm'),
                                            'anatomical_to_mni_nonlinear_xfm':(fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm'),
                                            'mni_to_anatomical_linear_xfm':(fnirt_reg_anat_mni, 'outputspec.invlinear_xfm'),
                                            'mni_normalized_anatomical':(fnirt_reg_anat_mni, 'outputspec.output_brain')})


                create_log_node(fnirt_reg_anat_mni, 'outputspec.output_brain', num_strat)
            
            
                num_strat += 1
                
        strat_list += new_strat_list


        
        new_strat_list = []
            
        for strat in strat_list:
            
            nodes = getNodeList(strat)
            
            # or run ANTS anatomical-to-MNI registration instead
            if ('ANTS' in c.regOption) and \
                    ('anat_mni_fnirt_register' not in nodes):

                ants_reg_anat_mni = create_wf_calculate_ants_warp('anat_mni' \
                        '_ants_register_%d' % num_strat, c.regWithSkull[0])

                try:

                    # calculating the transform with the skullstripped is
                    # reported to be better, but it requires very high
                    # quality skullstripping. If skullstripping is imprecise
                    # registration with skull is preferred
                    if (1 in c.regWithSkull):

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
                        node, out_file = strat.get_node_from_resource_pool('anatomical_brain')

                        # pass the anatomical to the workflow
                        workflow.connect(node, out_file,
                            ants_reg_anat_mni,'inputspec.anatomical_brain')

                        # pass the reference file
                        ants_reg_anat_mni.inputs.inputspec.reference_brain = \
                            c.template_brain_only_for_anat

                        # get the reorient skull-on anatomical from resource pool
                        node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')

                        # pass the anatomical to the workflow
                        workflow.connect(node, out_file,
                            ants_reg_anat_mni,'inputspec.anatomical_skull')

                        # pass the reference file
                        ants_reg_anat_mni.inputs.inputspec.reference_skull = \
                            c.template_skull_for_anat


                    else:

                        node, out_file = strat.get_node_from_resource_pool('anatomical_brain')

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
                        metric = ['MI','MI','CC']
                    ants_reg_anat_mni.inputs.inputspec.metric_weight = [1,1,1]
                    ants_reg_anat_mni.inputs.inputspec. \
                        radius_or_number_of_bins = [32,32,4]
                    ants_reg_anat_mni.inputs.inputspec. \
                        sampling_strategy = ['Regular','Regular',None]
                    ants_reg_anat_mni.inputs.inputspec. \
                        sampling_percentage = [0.25,0.25,None]
                    ants_reg_anat_mni.inputs.inputspec. \
                        number_of_iterations = [[1000,500,250,100], \
                        [1000,500,250,100], [100,100,70,20]]
                    ants_reg_anat_mni.inputs.inputspec. \
                        convergence_threshold = [1e-8,1e-8,1e-9]
                    ants_reg_anat_mni.inputs.inputspec. \
                        convergence_window_size = [10,10,15]
                    ants_reg_anat_mni.inputs.inputspec. \
                        transforms = ['Rigid','Affine','SyN']
                    ants_reg_anat_mni.inputs.inputspec. \
                        transform_parameters = [[0.1],[0.1],[0.1,3,0]]
                    ants_reg_anat_mni.inputs.inputspec. \
                        shrink_factors = [[8,4,2,1],[8,4,2,1],[6,4,2,1]]
                    ants_reg_anat_mni.inputs.inputspec. \
                        smoothing_sigmas = [[3,2,1,0],[3,2,1,0],[3,2,1,0]]
                

                except:
                    logConnectionError('Anatomical Registration (ANTS)', num_strat, strat.get_resource_pool(), '0003')
                    raise

                if 0 in c.runRegistrationPreprocessing:
                    tmp = strategy()
                    tmp.resource_pool = dict(strat.resource_pool)
                    tmp.leaf_node = (strat.leaf_node)
                    tmp.leaf_out_file = str(strat.leaf_out_file)
                    tmp.name = list(strat.name)
                    strat = tmp
                    new_strat_list.append(strat)

                strat.append_name(ants_reg_anat_mni.name)
                strat.set_leaf_properties(ants_reg_anat_mni, 'outputspec.normalized_output_brain')

                strat.update_resource_pool({'ants_initial_xfm':(ants_reg_anat_mni, 'outputspec.ants_initial_xfm'),
                                            'ants_rigid_xfm':(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm'),
                                            'ants_affine_xfm':(ants_reg_anat_mni, 'outputspec.ants_affine_xfm'),
                                            'anatomical_to_mni_nonlinear_xfm':(ants_reg_anat_mni, 'outputspec.warp_field'),
                                            'mni_to_anatomical_nonlinear_xfm':(ants_reg_anat_mni, 'outputspec.inverse_warp_field'),
                                            'anat_to_mni_ants_composite_xfm':(ants_reg_anat_mni, 'outputspec.composite_transform'),
                                            'mni_normalized_anatomical':(ants_reg_anat_mni, 'outputspec.normalized_output_brain')})

                create_log_node(ants_reg_anat_mni, 'outputspec.normalized_output_brain', num_strat)
          
                num_strat += 1
            
    strat_list += new_strat_list




    '''
    [SYMMETRIC] T1 -> Symmetric Template, Non-linear registration (FNIRT/ANTS)
    '''
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    
    # either run FSL anatomical-to-MNI registration, or...
    
    if (1 in c.runRegistrationPreprocessing) and (1 in c.runVMHC):

        if not os.path.exists(c.template_symmetric_brain_only):
            logger.info("\n\n" + ("ERROR: Missing file - %s" % c.template_symmetric_brain_only) + "\n\n" + \
                        "Error name: cpac_pipeline_0017" + "\n\n")
            raise Exception
        
        if not os.path.exists(c.template_symmetric_skull):
            logger.info("\n\n" + ("ERROR: Missing file - %s" % c.template_symmetric_skull) + "\n\n" + \
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

                fnirt_reg_anat_symm_mni = create_nonlinear_register('anat_symmetric_mni_fnirt_register_%d' % num_strat)

                try:
                
                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni, 'inputspec.input_brain')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni, 'inputspec.input_skull')

                    # pass the reference files                
                    fnirt_reg_anat_symm_mni.inputs.inputspec.reference_brain = c.template_symmetric_brain_only
                    fnirt_reg_anat_symm_mni.inputs.inputspec.reference_skull = c.template_symmetric_skull
                    fnirt_reg_anat_symm_mni.inputs.inputspec.ref_mask = c.dilated_symmetric_brain_mask

                    # assign the FSL FNIRT config file specified in pipeline config.yml
                    fnirt_reg_anat_symm_mni.inputs.inputspec.fnirt_config = c.configFileTwomm

                    #node, out_file = strat.get_node_from_resource_pool('mni_normalized_anatomical')
                    #workflow.connect(node, out_file,
                    #                 fnirt_reg_anat_symm_mni, 'inputspec.wait')

                except:
                    logConnectionError('Symmetric Anatomical Registration (FSL)', num_strat, strat.get_resource_pool(), '0002')
                    raise

                if (0 in c.runRegistrationPreprocessing):
                    tmp = strategy()
                    tmp.resource_pool = dict(strat.resource_pool)
                    tmp.leaf_node = (strat.leaf_node)
                    tmp.leaf_out_file = str(strat.leaf_out_file)
                    tmp.name = list(strat.name)
                    strat = tmp
                    new_strat_list.append(strat)

                strat.append_name(fnirt_reg_anat_symm_mni.name)
                strat.set_leaf_properties(fnirt_reg_anat_symm_mni, 'outputspec.output_brain')

                strat.update_resource_pool({'anatomical_to_symmetric_mni_linear_xfm':(fnirt_reg_anat_symm_mni, 'outputspec.linear_xfm'),
                                            'anatomical_to_symmetric_mni_nonlinear_xfm':(fnirt_reg_anat_symm_mni, 'outputspec.nonlinear_xfm'),
                                            'symmetric_mni_to_anatomical_linear_xfm':(fnirt_reg_anat_symm_mni, 'outputspec.invlinear_xfm'),
                                            'symmetric_mni_normalized_anatomical':(fnirt_reg_anat_symm_mni, 'outputspec.output_brain')})#,
                                            #'mni_normalized_anatomical':(ants_reg_anat_symm_mni, 'outputspec.wait')})


                create_log_node(fnirt_reg_anat_symm_mni, 'outputspec.output_brain', num_strat)
            
            
                num_strat += 1
                
        strat_list += new_strat_list


        
        new_strat_list = []
            
        for strat in strat_list:
            
            nodes = getNodeList(strat)
            
            # or run ANTS anatomical-to-MNI registration instead
            if ('ANTS' in c.regOption) and \
                    ('anat_mni_fnirt_register' not in nodes) and \
                    ('anat_symmetric_mni_fnirt_register' not in nodes):

                ants_reg_anat_symm_mni = create_wf_calculate_ants_warp('anat' \
                        '_symmetric_mni_ants_register_%d' % num_strat, \
                        c.regWithSkull[0])

                try:

                    # calculating the transform with the skullstripped is
                    # reported to be better, but it requires very high
                    # quality skullstripping. If skullstripping is imprecise
                    # registration with skull is preferred
                    if (1 in c.regWithSkull):

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
                        node, out_file = strat.get_node_from_resource_pool('anatomical_brain')

                        # pass the anatomical to the workflow
                        workflow.connect(node, out_file,
                            ants_reg_anat_symm_mni,'inputspec.anatomical_brain')

                        # pass the reference file
                        ants_reg_anat_symm_mni.inputs.inputspec.reference_brain = \
                            c.template_symmetric_brain_only

                        # get the reorient skull-on anatomical from resource pool
                        node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')

                        # pass the anatomical to the workflow
                        workflow.connect(node, out_file,
                            ants_reg_anat_symm_mni,'inputspec.anatomical_skull')

                        # pass the reference file
                        ants_reg_anat_symm_mni.inputs.inputspec.reference_skull = \
                            c.template_symmetric_skull


                    else:

                        # get the skullstripped anatomical from resource pool
                        node, out_file = strat.get_node_from_resource_pool('anatomical_brain')

                        workflow.connect(node, out_file, ants_reg_anat_symm_mni,
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
                        metric = ['MI','MI','CC']
                    ants_reg_anat_symm_mni.inputs.inputspec.metric_weight = [1,1,1]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        radius_or_number_of_bins = [32,32,4]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        sampling_strategy = ['Regular','Regular',None]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        sampling_percentage = [0.25,0.25,None]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        number_of_iterations = [[1000,500,250,100], \
                        [1000,500,250,100], [100,100,70,20]]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        convergence_threshold = [1e-8,1e-8,1e-9]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        convergence_window_size = [10,10,15]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        transforms = ['Rigid','Affine','SyN']
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        transform_parameters = [[0.1],[0.1],[0.1,3,0]]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        shrink_factors = [[8,4,2,1],[8,4,2,1],[6,4,2,1]]
                    ants_reg_anat_symm_mni.inputs.inputspec. \
                        smoothing_sigmas = [[3,2,1,0],[3,2,1,0],[3,2,1,0]]

                    #node, out_file = strat.get_node_from_resource_pool('mni_normalized_anatomical')
                    #workflow.connect(node, out_file,
                    #                 ants_reg_anat_symm_mni, 'inputspec.wait')
                

                except:
                    logConnectionError('Symmetric Anatomical Registration (ANTS)', num_strat, strat.get_resource_pool(), '0003')
                    raise

                if 0 in c.runRegistrationPreprocessing:
                    tmp = strategy()
                    tmp.resource_pool = dict(strat.resource_pool)
                    tmp.leaf_node = (strat.leaf_node)
                    tmp.leaf_out_file = str(strat.leaf_out_file)
                    tmp.name = list(strat.name)
                    strat = tmp
                    new_strat_list.append(strat)

                strat.append_name(ants_reg_anat_symm_mni.name)
                strat.set_leaf_properties(ants_reg_anat_symm_mni, 'outputspec.normalized_output_brain')

                strat.update_resource_pool({'ants_symmetric_initial_xfm':(ants_reg_anat_symm_mni, 'outputspec.ants_initial_xfm'),
                                            'ants_symmetric_rigid_xfm':(ants_reg_anat_symm_mni, 'outputspec.ants_rigid_xfm'),
                                            'ants_symmetric_affine_xfm':(ants_reg_anat_symm_mni, 'outputspec.ants_affine_xfm'),
                                            'anatomical_to_symmetric_mni_nonlinear_xfm':(ants_reg_anat_symm_mni, 'outputspec.warp_field'),
                                            'symmetric_mni_to_anatomical_nonlinear_xfm':(ants_reg_anat_symm_mni, 'outputspec.inverse_warp_field'),
                                            'anat_to_symmetric_mni_ants_composite_xfm':(ants_reg_anat_symm_mni, 'outputspec.composite_transform'),
                                            'symmetric_mni_normalized_anatomical':(ants_reg_anat_symm_mni, 'outputspec.normalized_output_brain')})#,
                                            #'mni_normalized_anatomical':(ants_reg_anat_symm_mni, 'outputspec.wait')})

                create_log_node(ants_reg_anat_symm_mni, 'outputspec.normalized_output_brain', num_strat)
          
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
                seg_preproc = create_seg_preproc(False, 'seg_preproc_%d' % num_strat)
            elif 'anat_mni_ants_register' in nodes:
                seg_preproc = create_seg_preproc(True, 'seg_preproc_%d' % num_strat)

            try:
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 seg_preproc, 'inputspec.brain')

                if 'anat_mni_fnirt_register' in nodes:
                    node, out_file = strat.get_node_from_resource_pool('mni_to_anatomical_linear_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc, 'inputspec.standard2highres_mat')
                elif 'anat_mni_ants_register' in nodes:
                    node, out_file = strat.get_node_from_resource_pool('ants_initial_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc, 'inputspec.standard2highres_init')
                    node, out_file = strat.get_node_from_resource_pool('ants_rigid_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc, 'inputspec.standard2highres_rig')

                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc, 'inputspec.standard2highres_mat')


                seg_preproc.inputs.inputspec.PRIOR_CSF = c.PRIORS_CSF
                seg_preproc.inputs.inputspec.PRIOR_GRAY = c.PRIORS_GRAY
                seg_preproc.inputs.inputspec.PRIOR_WHITE = c.PRIORS_WHITE

                seg_preproc.inputs.csf_threshold.csf_threshold = \
                                        c.cerebralSpinalFluidThreshold
                seg_preproc.inputs.wm_threshold.wm_threshold = \
                                        c.whiteMatterThreshold
                seg_preproc.inputs.gm_threshold.gm_threshold = \
                                        c.grayMatterThreshold
                seg_preproc.get_node('csf_threshold').iterables = ('csf_threshold',
                                        c.cerebralSpinalFluidThreshold)
                seg_preproc.get_node('wm_threshold').iterables = ('wm_threshold',
                                        c.whiteMatterThreshold)
                seg_preproc.get_node('gm_threshold').iterables = ('gm_threshold',
                                        c.grayMatterThreshold)


            except:
                logConnectionError('Segmentation Preprocessing', num_strat, strat.get_resource_pool(), '0004')
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
            strat.update_resource_pool({'anatomical_gm_mask' : (seg_preproc, 'outputspec.gm_mask'),
                                        'anatomical_csf_mask': (seg_preproc, 'outputspec.csf_mask'),
                                        'anatomical_wm_mask' : (seg_preproc, 'outputspec.wm_mask'),
                                        'seg_probability_maps': (seg_preproc, 'outputspec.probability_maps'),
                                        'seg_mixeltype': (seg_preproc, 'outputspec.mixeltype'),
                                        'seg_partial_volume_map': (seg_preproc, 'outputspec.partial_volume_map'),
                                        'seg_partial_volume_files': (seg_preproc, 'outputspec.partial_volume_files')})

            create_log_node(seg_preproc, 'outputspec.partial_volume_map', num_strat)
            num_strat += 1

    strat_list += new_strat_list




    '''
    Inserting Functional Data workflow
    '''

    num_strat = 0

    if 1 in c.runFunctionalDataGathering:
        for strat in strat_list:
            # create a new node, Remember to change its name!
            # Flow = create_func_datasource(sub_dict['rest'])
            # Flow.inputs.inputnode.subject = subject_id
            try: 
                funcFlow = create_func_datasource(sub_dict['rest'], 'func_gather_%d' % num_strat)
                funcFlow.inputs.inputnode.subject = subject_id
            except Exception as xxx:
                logger.info( "Error create_func_datasource failed."+\
                      " (%s:%d)" % dbg_file_lineno() )
                raise

            """
            Add in nodes to get parameters from configuration file
            """
            try:
                # a node which checks if scan _parameters are present for each scan
                scan_params = pe.Node(util.Function(input_names=['subject',
                                                                 'scan',
                                                                 'subject_map',
                                                                 'start_indx',
                                                                 'stop_indx',
                                                                 'tr',
                                                                 'tpattern'],
                                                   output_names=['tr',
                                                                 'tpattern',
                                                                 'ref_slice',
                                                                 'start_indx',
                                                                 'stop_indx'],
                                                   function=get_scan_params),
                                       name='scan_params_%d' % num_strat)
            except Exception as xxx:
                logger.info( "Error creating scan_params node."+\
                      " (%s:%d)" % dbg_file_lineno() )
                raise


            if "Selected Functional Volume" in c.func_reg_input:
               
                try:
                    get_func_volume = pe.Node(interface=preprocess.Calc(),
                        name = 'get_func_volume_%d' % num_strat)
         
                    get_func_volume.inputs.expr = 'a'
                    get_func_volume.inputs.single_idx = c.func_reg_input_volume
                    get_func_volume.inputs.outputtype = 'NIFTI_GZ' 
 
                except Exception as xxx:
                    logger.info( "Error creating get_func_volume node."+\
                          " (%s:%d)" % dbg_file_lineno() )
                    raise

                try:
                    workflow.connect(funcFlow, 'outputspec.rest',
                        get_func_volume, 'in_file_a')

                except Exception as xxx:
                    logger.info( "Error connecting get_func_volume node."+\
                          " (%s:%d)" % dbg_file_lineno() )
                    raise


            # wire in the scan parameter workflow
            try:
                workflow.connect(funcFlow, 'outputspec.subject',
                                 scan_params, 'subject')
            except Exception as xxx:
                logger.info( "Error connecting scan_params 'subject' input."+\
                      " (%s:%d)" % dbg_file_lineno() )
                raise


            try:
                workflow.connect(funcFlow, 'outputspec.scan',
                                 scan_params, 'scan')
            except Exception as xxx:
                logger.info( "Error connecting scan_params 'scan' input."+\
                      " (%s:%d)" % dbg_file_lineno() )
                raise

            # connect in constants
            scan_params.inputs.subject_map = sub_dict
            scan_params.inputs.start_indx = c.startIdx
            scan_params.inputs.stop_indx = c.stopIdx
            scan_params.inputs.tr = c.TR
            scan_params.inputs.tpattern = c.slice_timing_pattern[0]
    
            # node to convert TR between seconds and milliseconds
            try:
                convert_tr = pe.Node(util.Function(input_names=['tr'],
                                                   output_names=['tr'],
                                                   function=get_tr),
                                    name='convert_tr_%d' % num_strat)
            except Exception as xxx:
                logger.info( "Error creating convert_tr node."+\
                      " (%s:%d)" % dbg_file_lineno() )
                raise
    
            try:    
                workflow.connect(scan_params, 'tr',
                                  convert_tr, 'tr')
            except Exception as xxx:
                logger.info( "Error connecting convert_tr 'tr' input."+\
                      " (%s:%d)" % dbg_file_lineno() )
                raise
    

            strat.set_leaf_properties(funcFlow, 'outputspec.rest')
            strat.update_resource_pool({'raw_functional' : (funcFlow, 'outputspec.rest')})

            if "Selected Functional Volume" in c.func_reg_input:
                strat.update_resource_pool({'selected_func_volume': (get_func_volume, 'out_file')})

            num_strat += 1



        """
        Truncate scan length based on configuration information
        """

        num_strat = 0

        for strat in strat_list:
            try:
                trunc_wf=create_wf_edit_func( wf_name = "edit_func_%d"%(num_strat))
            except Exception as xxx:
                logger.info( "Error create_wf_edit_func failed."+\
                      " (%s:%d)" %(dbg_file_lineno()))
                raise

            # find the output data on the leaf node
            try:
                node, out_file = strat.get_leaf_properties()
            except Exception as xxx:
                logger.info( "Error  get_leaf_properties failed."+\
                      " (%s:%d)" % dbg_file_lineno() )
                raise

            # connect the functional data from the leaf node into the wf
            try:
                workflow.connect(node, out_file, trunc_wf, 'inputspec.func')
            except Exception as xxx:
                logger.info( "Error connecting input 'func' to trunc_wf."+\
                      " (%s:%d)" % dbg_file_lineno() )
                print xxx
                raise

            # connect the other input parameters
            try: 
                workflow.connect(scan_params, 'start_indx',
                                 trunc_wf, 'inputspec.start_idx')
            except Exception as xxx:
                logger.info( "Error connecting input 'start_indx' to trunc_wf."+\
                      " (%s:%d)" % dbg_file_lineno() )
                raise

            try:
                workflow.connect(scan_params, 'stop_indx',
                                 trunc_wf, 'inputspec.stop_idx')
            except Exception as xxx:
                logger.info( "Error connecting input 'stop_idx' to trunc_wf."+\
                      " (%s:%d)" % dbg_file_lineno() )
                raise

   
            # replace the leaf node with the output from the recently added workflow 
            strat.set_leaf_properties(trunc_wf, 'outputspec.edited_func')
            num_strat = num_strat+1



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
                    func_slice_timing_correction = pe.Node(interface=preprocess.TShift(),
                        name = 'func_slice_timing_correction_%d'%(num_strat))
                    func_slice_timing_correction.inputs.outputtype = 'NIFTI_GZ'
                except Exception as xxx:
                    logger.info( "Error connecting input 'stop_idx' to trunc_wf."+\
                          " (%s:%d)" % dbg_file_lineno() )
                    raise

                # find the output data on the leaf node
                try:
                    node, out_file = strat.get_leaf_properties()
                except Exception as xxx:
                    logger.info( "Error  get_leaf_properties failed."+\
                          " (%s:%d)" % dbg_file_lineno() )
                    raise

   
                # connect the output of the leaf node as the in_file
                try:
                    workflow.connect(node, out_file,
                        func_slice_timing_correction,'in_file')
                except Exception as xxx:
                    logger.info( "Error connecting input 'infile' to func_slice_timing_correction afni node."+\
                          " (%s:%d)" % dbg_file_lineno() )
                    raise

                logger.info("connected input to slc")
                # we might prefer to use the TR stored in the NIFTI header
                # if not, use the value in the scan_params node
                logger.info( "TR %s" %c.TR)
                if c.TR:
                    try:
                        workflow.connect(scan_params, 'tr',
                            func_slice_timing_correction, 'tr')
                    except Exception as xxx:
                        logger.info( "Error connecting input 'tr' to func_slice_timing_correction afni node."+\
                             " (%s:%d)" % dbg_file_lineno() )
                        print xxx
                        raise
                    logger.info("connected TR")

                # we might prefer to use the slice timing information stored in the NIFTI header
                # if not, use the value in the scan_params node
                logger.info( "slice timing pattern %s"%c.slice_timing_pattern[0])
                try:
                    if not "Use NIFTI Header" in c.slice_timing_pattern[0]:
                        try:
                            logger.info( "connecting slice timing pattern %s"%c.slice_timing_pattern[0])
                            workflow.connect(scan_params, 'tpattern',
                                func_slice_timing_correction, 'tpattern')
                            logger.info( "connected slice timing pattern %s"%c.slice_timing_pattern[0])
                        except Exception as xxx:
                            logger.info( "Error connecting input 'acquisition' to func_slice_timing_correction afni node."+\
                                 " (%s:%d)" % dbg_file_lineno() )
                            print xxx
                            raise
                        logger.info( "connected slice timing pattern %s"%c.slice_timing_pattern[0])
                except Exception as xxx:
                    logger.info( "Error connecting input 'acquisition' to func_slice_timing_correction afni node."+\
                                 " (%s:%d)" % dbg_file_lineno() )
                    print xxx
                    raise

                if (0 in c.runFunctionalPreprocessing):
                    # we are forking so create a new node
                    tmp = strategy()
                    tmp.resource_pool = dict(strat.resource_pool)
                    tmp.leaf_node = (strat.leaf_node)
                    tmp.leaf_out_file = str(strat.leaf_out_file)
                    tmp.name = list(strat.name)
                    strat = tmp
                    new_strat_list.append(strat)
   
                # add the name of the node to the strat name
                strat.append_name(func_slice_timing_correction.name)
   
                # set the leaf node 
                strat.set_leaf_properties(func_slice_timing_correction, 'out_file')

                # add the outputs to the resource pool
                strat.update_resource_pool({'slice_time_corrected':(func_slice_timing_correction, 'out_file')})
                num_strat += 1
   
            # add new strats (if forked) 
            strat_list += new_strat_list
    
            logger.info( " finished connecting slice timing pattern")



        """
        Inserting Functional Image Preprocessing
        Workflow
        """
        new_strat_list = []
        num_strat = 0
    
        workflow_counter += 1
    
        if 1 in c.runFunctionalPreprocessing:
    
            workflow_bit_id['func_preproc'] = workflow_counter
    
            for strat in strat_list:
    
                if '3dAutoMask' in c.functionalMasking:
    
                    try:
                        func_preproc = create_func_preproc(use_bet=False, wf_name='func_preproc_automask_%d' % num_strat)
                    except Exception as xxx:
                        logger.info( "Error allocating func_preproc."+\
                             " (%s:%d)" % dbg_file_lineno() )
                        raise

                    node = None
                    out_file = None
                    try:
                        node, out_file = strat.get_leaf_properties()
                        logger.info("%s::%s==>%s"%(node, out_file,func_preproc)) 
                        try:
                            workflow.connect(node, out_file, func_preproc, 'inputspec.func')
                        except Exception as xxx:
                            logger.info( "Error connecting leafnode to func, func_preproc."+\
                                 " (%s:%d)" % (dbg_file_lineno()) )
                            print xxx
                            raise
                        logger.info("infile rest connected") 
                    except Exception as xxx:
                        logConnectionError('Functional Preprocessing', num_strat, strat.get_resource_pool(), '0005_automask')
                        num_strat += 1
                        raise
    
                    if (0 in c.runFunctionalPreprocessing) or ('BET' in c.functionalMasking):
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
                    strat.update_resource_pool({'mean_functional':(func_preproc, 'outputspec.example_func')})
                    strat.update_resource_pool({'functional_preprocessed_mask':(func_preproc, 'outputspec.preprocessed_mask')})
                    strat.update_resource_pool({'movement_parameters':(func_preproc, 'outputspec.movement_parameters')})
                    strat.update_resource_pool({'max_displacement':(func_preproc, 'outputspec.max_displacement')})
                    strat.update_resource_pool({'preprocessed':(func_preproc, 'outputspec.preprocessed')})
                    strat.update_resource_pool({'functional_brain_mask':(func_preproc, 'outputspec.mask')})
                    strat.update_resource_pool({'motion_correct':(func_preproc, 'outputspec.motion_correct')})
                    strat.update_resource_pool({'coordinate_transformation':(func_preproc, 'outputspec.oned_matrix_save')})
    
    
                    create_log_node(func_preproc, 'outputspec.preprocessed', num_strat)
                    num_strat += 1
    
            strat_list += new_strat_list
    
            new_strat_list = []
                
    
            for strat in strat_list:
                
                nodes = getNodeList(strat)
                
                if ('BET' in c.functionalMasking) and ('func_preproc_automask' not in nodes):
    
                    func_preproc = create_func_preproc( use_bet=True, \
                                                        wf_name='func_preproc_bet_%d' % num_strat)
                    node = None
                    out_file = None
                    try:
                        node, out_file = strat.get_leaf_properties()
                        workflow.connect(node, out_file, func_preproc, 'inputspec.func')
    
                    except Exception as xxx:
                        logConnectionError('Functional Preprocessing', num_strat, strat.get_resource_pool(), '0005_bet')
                        num_strat += 1
                        raise
    
                    if 0 in c.runFunctionalPreprocessing:
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
                    strat.update_resource_pool({'mean_functional':(func_preproc, 'outputspec.example_func')})
                    strat.update_resource_pool({'functional_preprocessed_mask':(func_preproc, 'outputspec.preprocessed_mask')})
                    strat.update_resource_pool({'movement_parameters':(func_preproc, 'outputspec.movement_parameters')})
                    strat.update_resource_pool({'max_displacement':(func_preproc, 'outputspec.max_displacement')})
                    strat.update_resource_pool({'preprocessed':(func_preproc, 'outputspec.preprocessed')})
                    strat.update_resource_pool({'functional_brain_mask':(func_preproc, 'outputspec.mask')})
                    strat.update_resource_pool({'motion_correct':(func_preproc, 'outputspec.motion_correct')})
                    strat.update_resource_pool({'coordinate_transformation':(func_preproc, 'outputspec.oned_matrix_save')})
    
    
                    create_log_node(func_preproc, 'outputspec.preprocessed', num_strat)
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
    
                fristons_model = fristons_twenty_four(wf_name='fristons_parameter_model_%d' % num_strat)
    
                try:
    
                    node, out_file = strat.get_node_from_resource_pool('movement_parameters')
                    workflow.connect(node, out_file,
                                     fristons_model, 'inputspec.movement_file')
    
                except Exception as xxx:
                    logConnectionError('Friston\'s Parameter Model', num_strat, strat.get_resource_pool(), '0006')
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
    
                strat.update_resource_pool({'movement_parameters':(fristons_model, 'outputspec.movement_file')})
    
        
                create_log_node(fristons_model, 'outputspec.movement_file', num_strat)
                
                num_strat += 1
                
        strat_list += new_strat_list




    '''
    Func -> T1 Registration (Initial Linear reg)
    '''

    # Depending on configuration, either passes output matrix to Func -> Template ApplyWarp,
    # or feeds into linear reg of BBReg operation (if BBReg is enabled)

    new_strat_list = []
    num_strat = 0
    workflow_counter += 1
    
    if 1 in c.runRegisterFuncToAnat:

        workflow_bit_id['func_to_anat'] = workflow_counter

        for strat in strat_list:
            func_to_anat = create_register_func_to_anat('func_to_anat_FLIRT_%d' % num_strat)
       
            # Input registration parameters
            func_to_anat.inputs.inputspec.interp = 'trilinear'

            try:
                def pick_wm(seg_prob_list):
                    seg_prob_list.sort()
                    return seg_prob_list[-1]

                if 'Mean Functional' in c.func_reg_input:
                    # Input functional image (mean functional)
                    node, out_file = strat.get_node_from_resource_pool('mean_functional')
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.func')

                elif 'Selected Functional Volume' in c.func_reg_input:
                    # Input functional image (specific volume)
                    node, out_file = strat.get_node_from_resource_pool('selected_func_volume')
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.func')


                # Input skull-stripped anatomical (anat.nii.gz)
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.anat')

   

            except:               
                logConnectionError('Register Functional to Anatomical (pre BBReg)', num_strat, strat.get_resource_pool(), '0007')
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

            strat.update_resource_pool({'mean_functional_in_anat':(func_to_anat, 'outputspec.anat_func_nobbreg'),
                                        'functional_to_anat_linear_xfm':(func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg')})

            # Outputs:
            # functional_to_anat_linear_xfm = func-t1.mat, linear, sent to 'premat' of post-FNIRT applywarp,
            #                                 or to the input of the post-ANTS c3d_affine_tool
            
            #create_log_node(func_to_anat, 'outputspec.mni_func', num_strat)
            num_strat += 1

    strat_list += new_strat_list



    '''
    Func -> T1 Registration (BBREG)
    '''

    # Outputs 'functional_to_anat_linear_xfm', a matrix file of the functional-to-anatomical
    # registration warp to be applied LATER in func_mni_warp, which accepts it as input 'premat'

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

                func_to_anat_bbreg = create_bbregister_func_to_anat('func_to_anat_bbreg_%d' % num_strat)
       
                # Input registration parameters
                func_to_anat_bbreg.inputs.inputspec.bbr_schedule = c.boundaryBasedRegistrationSchedule

                try:
                    def pick_wm(seg_prob_list):
                        seg_prob_list.sort()
                        return seg_prob_list[-1]


                    if 'Mean Functional' in c.func_reg_input:
                        # Input functional image (mean functional)
                        node, out_file = strat.get_node_from_resource_pool('mean_functional')
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg, 'inputspec.func')

                    elif 'Selected Functional Volume' in c.func_reg_input:
                        # Input functional image (specific volume)
                        node, out_file = strat.get_node_from_resource_pool('selected_func_volume')
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg, 'inputspec.func')

                    # Input segmentation probability maps for white matter segmentation
                    node, out_file = strat.get_node_from_resource_pool('seg_probability_maps')
                    workflow.connect(node, (out_file, pick_wm),
                                     func_to_anat_bbreg, 'inputspec.anat_wm_segmentation')

                    # Input anatomical whole-head image (reoriented)
                    node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                    workflow.connect(node, out_file,
                                     func_to_anat_bbreg, 'inputspec.anat_skull')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     func_to_anat_bbreg, 'inputspec.linear_reg_matrix')
   

                except:
                    logConnectionError('Register Functional to Anatomical (BBReg)', num_strat, strat.get_resource_pool(), '0008')
                    raise


                if 0 in c.runBBReg:
                    tmp = strategy()
                    tmp.resource_pool = dict(strat.resource_pool)
                    tmp.leaf_node = (strat.leaf_node)
                    tmp.out_file = str(strat.leaf_out_file)
                
                    # This line is needed here for some reason, otherwise the connection
                    # between func_preproc and nuisance will break - even though this
                    # workflow has nothing to do with it - but excluding this line
                    # below removes the leaf node from the new forked strat
                    tmp.leaf_out_file = str(strat.leaf_out_file)
                
                    tmp.name = list(strat.name)
                    strat = tmp
                    new_strat_list.append(strat)
               

                strat.append_name(func_to_anat_bbreg.name)
                # strat.set_leaf_properties(func_mni_warp, 'out_file')

                strat.update_resource_pool({'mean_functional_in_anat':(func_to_anat_bbreg, 'outputspec.anat_func'),
                                            'functional_to_anat_linear_xfm':(func_to_anat_bbreg, 'outputspec.func_to_anat_linear_xfm')})

                # Outputs:
                # functional_to_anat_linear_xfm = func-t1.mat, linear, sent to 'premat' of post-FNIRT applywarp,
                #                                 or to the input of the post-ANTS c3d_affine_tool
            
                #create_log_node(func_to_anat, 'outputspec.mni_func', num_strat)
                num_strat += 1

    strat_list += new_strat_list
     


    '''
    Inserting Generate Motion Statistics Workflow
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runGenerateMotionStatistics:
        workflow_bit_id['gen_motion_stats'] = workflow_counter
        for strat in strat_list:

            gen_motion_stats = motion_power_statistics('gen_motion_stats_%d' % num_strat)
            gen_motion_stats.inputs.scrubbing_input.threshold = c.scrubbingThreshold
            gen_motion_stats.inputs.scrubbing_input.remove_frames_before = c.numRemovePrecedingFrames
            gen_motion_stats.inputs.scrubbing_input.remove_frames_after = c.numRemoveSubsequentFrames
            gen_motion_stats.get_node('scrubbing_input').iterables = ('threshold', c.scrubbingThreshold)

            try:
                # #**special case where the workflow is not getting outputs from resource pool
                # but is connected to functional datasource
                workflow.connect(funcFlow, 'outputspec.subject',
                             gen_motion_stats, 'inputspec.subject_id')

                workflow.connect(funcFlow, 'outputspec.scan',
                             gen_motion_stats, 'inputspec.scan_id')

                node, out_file = strat.get_node_from_resource_pool('motion_correct')
                workflow.connect(node, out_file,
                                 gen_motion_stats, 'inputspec.motion_correct')

                node, out_file = strat.get_node_from_resource_pool('movement_parameters')
                workflow.connect(node, out_file,
                                 gen_motion_stats, 'inputspec.movement_parameters')

                node, out_file = strat.get_node_from_resource_pool('max_displacement')
                workflow.connect(node, out_file,
                                 gen_motion_stats, 'inputspec.max_displacement')

                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 gen_motion_stats, 'inputspec.mask')

                node, out_file = strat.get_node_from_resource_pool('coordinate_transformation')
                workflow.connect(node, out_file,
                                 gen_motion_stats, 'inputspec.oned_matrix_save')

            except:
                logConnectionError('Generate Motion Statistics', num_strat, strat.get_resource_pool(), '0009')
                raise

            if 0 in c.runGenerateMotionStatistics:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(gen_motion_stats.name)

            strat.update_resource_pool({'frame_wise_displacement':(gen_motion_stats, 'outputspec.FD_1D'),
                                        'scrubbing_frames_excluded':(gen_motion_stats, 'outputspec.frames_ex_1D'),
                                        'scrubbing_frames_included':(gen_motion_stats, 'outputspec.frames_in_1D'),
                                        'power_params':(gen_motion_stats, 'outputspec.power_params'),
                                        'motion_params':(gen_motion_stats, 'outputspec.motion_params')})
            
            create_log_node(gen_motion_stats, 'outputspec.motion_params', num_strat)
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
            
                if 'anat_mni_fnirt_register' in nodes:
                    nuisance = create_nuisance(False, 'nuisance_%d' % num_strat)
                else:
                    nuisance = create_nuisance(True, 'nuisance_%d' % num_strat)

                nuisance.get_node('residuals').iterables = ([('selector', c.Corrections),
                                                             ('compcor_ncomponents', c.nComponents)])

                nuisance.inputs.inputspec.lat_ventricles_mask = c.lateral_ventricles_mask

                try:
                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.subject')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_gm_mask')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.gm_mask')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_wm_mask')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.wm_mask')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_csf_mask')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.csf_mask')

                    node, out_file = strat.get_node_from_resource_pool('movement_parameters')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.motion_components')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.func_to_anat_linear_xfm')
                
                    if 'anat_mni_fnirt_register' in nodes:
                        node, out_file = strat.get_node_from_resource_pool('mni_to_anatomical_linear_xfm')
                        workflow.connect(node, out_file,
                                         nuisance, 'inputspec.mni_to_anat_linear_xfm')
                    else:
                        # pass the ants_affine_xfm to the input for the
                        # INVERSE transform, but ants_affine_xfm gets inverted
                        # within the workflow

                        node, out_file = strat.get_node_from_resource_pool('ants_initial_xfm')
                        workflow.connect(node, out_file,
                                         nuisance, 'inputspec.anat_to_mni_initial_xfm')

                        node, out_file = strat.get_node_from_resource_pool('ants_rigid_xfm')
                        workflow.connect(node, out_file,
                                         nuisance, 'inputspec.anat_to_mni_rigid_xfm')

                        node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                        workflow.connect(node, out_file,
                                         nuisance, 'inputspec.anat_to_mni_affine_xfm')


                except:
                    logConnectionError('Nuisance', num_strat, strat.get_resource_pool(), '0010')
                    raise


                if 0 in c.runNuisance:
                    tmp = strategy()
                    tmp.resource_pool = dict(strat.resource_pool)
                    tmp.leaf_node = (strat.leaf_node)
                    tmp.leaf_out_file = str(strat.leaf_out_file)
                    tmp.name = list(strat.name)
                    strat = tmp
                    new_strat_list.append(strat)

                strat.append_name(nuisance.name)

                strat.set_leaf_properties(nuisance, 'outputspec.subject')

                strat.update_resource_pool({'functional_nuisance_residuals':(nuisance, 'outputspec.subject')})

                create_log_node(nuisance, 'outputspec.subject', num_strat)
            
                num_strat += 1

    strat_list += new_strat_list



    '''
    Inserting Median Angle Correction Workflow
    '''
    
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runMedianAngleCorrection:
        workflow_bit_id['median_angle_corr'] = workflow_counter
        for strat in strat_list:
            median_angle_corr = create_median_angle_correction('median_angle_corr_%d' % num_strat)

            median_angle_corr.get_node('median_angle_correct').iterables = ('target_angle_deg', c.targetAngleDeg)
            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 median_angle_corr, 'inputspec.subject')
            except:
                logConnectionError('Median Angle Correction', num_strat, strat.get_resource_pool(), '0011')
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

            strat.update_resource_pool({'functional_median_angle_corrected':(median_angle_corr, 'outputspec.subject')})
            
            create_log_node(median_angle_corr, 'outputspec.subject', num_strat)

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
                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 alff, 'inputspec.rest_mask')

            except:
                logConnectionError('ALFF', num_strat, strat.get_resource_pool(), '0012')
                raise
            strat.append_name(alff.name)
            strat.update_resource_pool({'alff_img':(alff, 'outputspec.alff_img')})
            strat.update_resource_pool({'falff_img':(alff, 'outputspec.falff_img')})
            
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
        for strat in strat_list:
            frequency_filter = pe.Node(util.Function(input_names=['realigned_file',
                                                                  'bandpass_freqs',
                                                                  'sample_period'],
                                                     output_names=['bandpassed_file'],
                                                     function=bandpass_voxels),
                                       name='frequency_filter_%d' % num_strat)

            frequency_filter.iterables = ('bandpass_freqs', c.nuisanceBandpassFreq)
            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 frequency_filter, 'realigned_file')

            except:
                logConnectionError('Frequency Filtering', num_strat, strat.get_resource_pool(), '0013')
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

            strat.update_resource_pool({'functional_freq_filtered':(frequency_filter, 'bandpassed_file')})

            create_log_node(frequency_filter, 'bandpassed_file', num_strat)
            
            num_strat += 1

    strat_list += new_strat_list



    '''
    Inserting Scrubbing Workflow
    '''
    
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1


    if 1 in c.runScrubbing:

        workflow_bit_id['scrubbing'] = workflow_counter

        for strat in strat_list:

            nodes = getNodeList(strat)

            if 'gen_motion_stats' in nodes:

                scrubbing = create_scrubbing_preproc('scrubbing_%d' % num_strat)

                try:

                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file,
                                     scrubbing, 'inputspec.preprocessed')

                    node, out_file = strat.get_node_from_resource_pool('scrubbing_frames_included')
                    workflow.connect(node, out_file,
                                     scrubbing, 'inputspec.frames_in_1D')

                    node, out_file = strat.get_node_from_resource_pool('movement_parameters')
                    workflow.connect(node, out_file,
                                     scrubbing, 'inputspec.movement_parameters')

                except:
                    logConnectionError('Scrubbing Workflow', num_strat, strat.get_resource_pool(), '0014')
                    raise

                if 0 in c.runScrubbing:
                    tmp = strategy()
                    tmp.resource_pool = dict(strat.resource_pool)
                    tmp.leaf_node = (strat.leaf_node)
                    tmp.leaf_out_file = str(strat.leaf_out_file)
                    tmp.name = list(strat.name)
                    strat = tmp
                    new_strat_list.append(strat)

                strat.append_name(scrubbing.name)

                strat.set_leaf_properties(scrubbing, 'outputspec.preprocessed')

                strat.update_resource_pool({'scrubbing_movement_parameters' : (scrubbing, 'outputspec.scrubbed_movement_parameters'),
                                            'scrubbed_preprocessed': (scrubbing, 'outputspec.preprocessed')})
            
                create_log_node(scrubbing, 'outputspec.preprocessed', num_strat)

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
    
                functional_brain_mask_to_standard = pe.Node(interface=fsl.ApplyWarp(),
                                                            name='func_mni_fsl_warp_mask_%d' % num_strat)
                functional_brain_mask_to_standard.inputs.interp = 'nn'
                functional_brain_mask_to_standard.inputs.ref_file = c.template_skull_for_func

                mean_functional_warp = pe.Node(interface=fsl.ApplyWarp(), name='mean_func_fsl_warp_%d' % num_strat)
                mean_functional_warp.inputs.ref_file = c.template_brain_only_for_func
                
                
                motion_correct_warp = pe.Node(interface=fsl.ApplyWarp(), name="motion_correct_fsl_warp_%d" % num_strat)
                motion_correct_warp.inputs.ref_file = c.template_brain_only_for_func
                
    
                try:

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     func_mni_warp, 'field_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     func_mni_warp, 'premat')
  
                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file,
                                     func_mni_warp, 'in_file')
    

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     functional_brain_mask_to_standard, 'field_file')
                    workflow.connect(node, out_file,
                                     mean_functional_warp, 'field_file')
                    workflow.connect(node, out_file,
                                     motion_correct_warp, 'field_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     functional_brain_mask_to_standard, 'premat') 
                    workflow.connect(node, out_file,
                                     mean_functional_warp, 'premat')
                    workflow.connect(node, out_file,
                                     motion_correct_warp, 'premat') 

                    node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                    workflow.connect(node, out_file,
                                     functional_brain_mask_to_standard, 'in_file')


                    node, out_file = strat.get_node_from_resource_pool('mean_functional')
                    workflow.connect(node, out_file, mean_functional_warp, 'in_file')
                    
                    
                    node, out_file = strat.get_node_from_resource_pool('motion_correct')
                    workflow.connect(node, out_file, motion_correct_warp, 'in_file')

                    

                except:
                    logConnectionError('Functional Timeseries Registration to MNI space (FSL)', num_strat, strat.get_resource_pool(), '0015')
                    raise
    
                strat.update_resource_pool({'functional_mni':(func_mni_warp, 'out_file'),
                                            'functional_brain_mask_to_standard':(functional_brain_mask_to_standard, 'out_file'),
                                            'mean_functional_in_mni':(mean_functional_warp, 'out_file'),
                                            'motion_correct_to_standard':(motion_correct_warp, 'out_file')})
                strat.append_name(func_mni_warp.name)
                create_log_node(func_mni_warp, 'out_file', num_strat)
            
                num_strat += 1
                
                
        strat_list += new_strat_list  
            
        for strat in strat_list:
            
            nodes = getNodeList(strat)
             
            if ('ANTS' in c.regOption) and ('anat_mni_fnirt_register' not in nodes):

                # ANTS warp application

                def fsl_to_itk_conversion(source_file, reference, func_name):

                    # converts FSL-format .mat affine xfm into ANTS-format
                    # .txt; .mat affine comes from Func->Anat registration
                    fsl_to_itk_func_mni = create_wf_c3d_fsl_to_itk(0, name=\
                            'fsl_to_itk_%s_%d' % (func_name, \
                            num_strat))

                    try:

                        # convert the .mat from linear Func->Anat to
                        # ANTS format
                        node, out_file = strat.get_node_from_resource_pool(\
                                'functional_to_anat_linear_xfm')
                        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                                'inputspec.affine_file')

                        node, out_file = strat.get_node_from_resource_pool(\
                                reference)
                        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                                'inputspec.reference_file')

                        node, out_file = strat.get_node_from_resource_pool(\
                                source_file)
                        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                                'inputspec.source_file')

                    except:
                        logConnectionError('Functional Timeseries ' \
                            'Registration to MNI space (ANTS)', num_strat, \
                            strat.get_resource_pool(), '0016')
                        raise

                    strat.update_resource_pool({'itk_func_anat_affine_%s' % \
                            (func_name): (fsl_to_itk_func_mni, \
                            'outputspec.itk_transform')})

                    strat.append_name(fsl_to_itk_func_mni.name)
                    create_log_node(fsl_to_itk_func_mni, \
                            'outputspec.itk_transform', num_strat)



                def collect_transforms_func_mni(func_name):

                    # collects series of warps to be applied
                    collect_transforms_func_mni = \
                            create_wf_collect_transforms(0, name=\
                            'collect_transforms_%s_%d' % \
                            (func_name, num_strat))

                    try:

                        # Field file from anatomical nonlinear registration
                        node, out_file = strat.get_node_from_resource_pool(\
                                'anatomical_to_mni_nonlinear_xfm')
                        workflow.connect(node, out_file,
                                collect_transforms_func_mni,
                                'inputspec.warp_file')

                        # initial transformation from anatomical registration
                        node, out_file = strat.get_node_from_resource_pool(\
                                'ants_initial_xfm')
                        workflow.connect(node, out_file,
                                collect_transforms_func_mni,
                                'inputspec.linear_initial')

                        # affine transformation from anatomical registration
                        node, out_file = strat.get_node_from_resource_pool(\
                                'ants_affine_xfm')
                        workflow.connect(node, out_file,
                                collect_transforms_func_mni,
                                'inputspec.linear_affine')

                        # rigid transformation from anatomical registration
                        node, out_file = strat.get_node_from_resource_pool(\
                                'ants_rigid_xfm')
                        workflow.connect(node, out_file,
                                collect_transforms_func_mni,
                                'inputspec.linear_rigid')

                        # Premat from Func->Anat linear reg and bbreg
                        # (if bbreg is enabled)
                        node, out_file = strat.get_node_from_resource_pool(\
                                'itk_func_anat_affine_%s' % func_name)
                        workflow.connect(node, out_file, 
                                collect_transforms_func_mni,
                                'inputspec.fsl_to_itk_affine')

                    except:
                        logConnectionError('Functional Timeseries ' \
                            'Registration to MNI space (ANTS)', num_strat, \
                            strat.get_resource_pool(), '0016')
                        raise

                    strat.update_resource_pool({'itk_collected_warps_%s' % \
                            (func_name): (collect_transforms_func_mni, \
                            'outputspec.transformation_series')})

                    strat.append_name(collect_transforms_func_mni.name)
                    create_log_node(collect_transforms_func_mni, \
                            'outputspec.transformation_series', num_strat)




                def ants_apply_warps_func_mni(input_node, input_outfile, \
                        reference, interp, input_image_type, func_name):

                    # apply ants warps
                    apply_ants_warp_func_mni = create_wf_apply_ants_warp(0, \
                            name='apply_ants_warp_%s_%d' % \
                            (func_name, num_strat))

                    apply_ants_warp_func_mni.inputs.inputspec. \
                            reference_image = reference

                    apply_ants_warp_func_mni.inputs.inputspec.dimension = 3

                    apply_ants_warp_func_mni.inputs.inputspec. \
                            interpolation = interp

                    # input_image_type:
                    # (0 or 1 or 2 or 3)
                    # Option specifying the input image type of scalar
                    # (default), vector, tensor, or time series.
                    apply_ants_warp_func_mni.inputs.inputspec. \
                            input_image_type = input_image_type

                    try:

                        # this <node, out_file> pulls in directly because
                        # it pulls in the leaf in some instances
                        workflow.connect(input_node, input_outfile,
                            apply_ants_warp_func_mni, 'inputspec.input_image')

                        node, out_file = strat.get_node_from_resource_pool(\
                                'itk_collected_warps_%s' % func_name)
                        workflow.connect(node, out_file, 
                                apply_ants_warp_func_mni,
                                'inputspec.transforms')


                    except:
                        logConnectionError('Functional Timeseries ' \
                            'Registration to MNI space (ANTS)', num_strat, \
                            strat.get_resource_pool(), '0016')
                        raise

                    strat.update_resource_pool({func_name: \
                            (apply_ants_warp_func_mni, \
                            'outputspec.output_image')})

                    strat.append_name(apply_ants_warp_func_mni.name)
                    create_log_node(apply_ants_warp_func_mni, \
                            'outputspec.output_image', num_strat)




                # 4D FUNCTIONAL apply warp
                fsl_to_itk_conversion('mean_functional', 'anatomical_brain', 'functional_mni')
                collect_transforms_func_mni('functional_mni')

                node, out_file = strat.get_leaf_properties()
                ants_apply_warps_func_mni(node, out_file, c.template_brain_only_for_func, 'Linear', 3, 'functional_mni')



                # FUNCTIONAL MASK apply warp
                fsl_to_itk_conversion('functional_brain_mask', 'anatomical_brain', 'functional_brain_mask_to_standard')
                collect_transforms_func_mni('functional_brain_mask_to_standard')

                node, out_file = strat.get_node_from_resource_pool('func' \
                        'tional_brain_mask')
                ants_apply_warps_func_mni(node, out_file, c.template_brain_only_for_func, 'NearestNeighbor', 0, 'functional_brain_mask_to_standard')



                # FUNCTIONAL MEAN apply warp
                fsl_to_itk_conversion('mean_functional', 'anatomical_brain', 'mean_functional_in_mni')
                collect_transforms_func_mni('mean_functional_in_mni')

                node, out_file = strat.get_node_from_resource_pool('mean' \
                        '_functional')
                ants_apply_warps_func_mni(node, out_file, c.template_brain_only_for_func, 'Linear', 0, 'mean_functional_in_mni')


                # 4D FUNCTIONAL MOTION-CORRECTED apply warp
                fsl_to_itk_conversion('mean_functional', 'anatomical_brain', 'motion_correct_to_standard')
                collect_transforms_func_mni('motion_correct_to_standard')

                node, out_file = strat.get_node_from_resource_pool('motion_correct')
                ants_apply_warps_func_mni(node, out_file, c.template_brain_only_for_func, 'Linear', 3, 'motion_correct_to_standard')

            
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
                vmhc = create_vmhc(True, 'vmhc_%d' % num_strat)

            vmhc.inputs.inputspec.standard_for_func = c.template_skull_for_func
            vmhc.inputs.fwhm_input.fwhm = c.fwhm
            vmhc.get_node('fwhm_input').iterables = ('fwhm', c.fwhm)


            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.rest_res')
                
                node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.example_func2highres_mat')

                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.rest_mask')

                node, out_file = strat.get_node_from_resource_pool('mean_functional')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.mean_functional')

                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.brain')



                if ('ANTS' in c.regOption) and \
                    ('anat_mni_fnirt_register' not in nodes) and \
                    ('anat_symmetric_mni_fnirt_register' not in nodes):

                    node, out_file = strat.get_node_from_resource_pool('ants_symmetric_initial_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.ants_symm_initial_xfm')

                    node, out_file = strat.get_node_from_resource_pool('ants_symmetric_rigid_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.ants_symm_rigid_xfm')

                    node, out_file = strat.get_node_from_resource_pool('ants_symmetric_affine_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.ants_symm_affine_xfm')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_symmetric_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.ants_symm_warp_field')
                                     
                else:

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_symmetric_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     vmhc, 'inputspec.fnirt_nonlinear_warp')


            except:
                logConnectionError('VMHC', num_strat, strat.get_resource_pool(), '0019')
                raise


            strat.update_resource_pool({'vmhc_raw_score':(vmhc, 'outputspec.VMHC_FWHM_img')})
            strat.update_resource_pool({'vmhc_fisher_zstd':(vmhc, 'outputspec.VMHC_Z_FWHM_img')})
            strat.update_resource_pool({'vmhc_fisher_zstd_zstat_map':(vmhc, 'outputspec.VMHC_Z_stat_FWHM_img')})
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
                err_msg = 'Cluster size specified: %d, is not supported. '\
                          'Change to 7, 19, or 27 and try again' % cluster_size
                raise Exception(err_msg)
            else:
                preproc.inputs.inputspec.cluster_size = cluster_size
                reho = preproc.clone('reho_%d' % num_strat)

            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 reho, 'inputspec.rest_res_filt')

                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 reho, 'inputspec.rest_mask')
            except:
                logConnectionError('ReHo', num_strat, strat.get_resource_pool(), '0020')
                raise


            strat.update_resource_pool({'raw_reho_map':(reho, 'outputspec.raw_reho_map')})
            strat.append_name(reho.name)
            
            create_log_node(reho, 'outputspec.raw_reho_map', num_strat)
            
            num_strat += 1
    strat_list += new_strat_list



    '''
    Spatial Regression Based Time Series
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runSpatialRegression:

        for strat in strat_list:

            resample_spatial_map_to_native_space = pe.Node(interface=fsl.FLIRT(),
                                                         name='resample_spatial_map_to_native_space_%d' % num_strat)
            resample_spatial_map_to_native_space.inputs.interp = 'nearestneighbour'
            resample_spatial_map_to_native_space.inputs.apply_xfm = True
            resample_spatial_map_to_native_space.inputs.in_matrix_file = c.identityMatrix

            spatial_map_dataflow = create_spatial_map_dataflow(c.spatialPatternMaps, 'spatial_map_dataflow_%d' % num_strat)

            spatial_map_timeseries = get_spatial_map_timeseries('spatial_map_timeseries_%d' % num_strat)
            spatial_map_timeseries.inputs.inputspec.demean = c.spatialDemean

            try:

                node, out_file = strat.get_node_from_resource_pool('functional_mni')
                node2, out_file2 = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')

                # resample the input functional file and functional mask to spatial map
                workflow.connect(node, out_file,
                                 resample_spatial_map_to_native_space, 'reference')
                workflow.connect(spatial_map_dataflow, 'select_spatial_map.out_file',
                                 resample_spatial_map_to_native_space, 'in_file')
                               
                # connect it to the spatial_map_timeseries
                workflow.connect(resample_spatial_map_to_native_space, 'out_file',
                                 spatial_map_timeseries, 'inputspec.spatial_map')
                workflow.connect(node2, out_file2,
                                 spatial_map_timeseries, 'inputspec.subject_mask')
                workflow.connect(node, out_file,
                                 spatial_map_timeseries, 'inputspec.subject_rest')
                               
                

            except:
                logConnectionError('Spatial map timeseries extraction', num_strat, strat.get_resource_pool(), '0029')
                raise

            if 0 in c.runSpatialRegression:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(spatial_map_timeseries.name)

            strat.update_resource_pool({'spatial_map_timeseries' : (spatial_map_timeseries, 'outputspec.subject_timeseries')})
            
            create_log_node(spatial_map_timeseries, 'outputspec.subject_timeseries', num_strat)

            num_strat += 1

    strat_list += new_strat_list



    '''
    ROI Based Time Series
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runROITimeseries:

        for strat in strat_list:

            if c.roiSpecificationFile != None:

                resample_functional_to_roi = pe.Node(interface=fsl.FLIRT(),
                                                      name='resample_functional_to_roi_%d' % num_strat)
                resample_functional_to_roi.inputs.interp = 'trilinear'
                resample_functional_to_roi.inputs.apply_xfm = True
                resample_functional_to_roi.inputs.in_matrix_file = c.identityMatrix
    
                roi_dataflow = create_roi_mask_dataflow(c.roiSpecificationFile, 'roi', 'roi_dataflow_%d' % num_strat)
    
                roi_timeseries = get_roi_timeseries('roi_timeseries_%d' % num_strat)
                roi_timeseries.inputs.inputspec.output_type = c.roiTSOutputs
            
            
            if c.roiSpecificationFileForSCA != None:
            
                # same workflow, except to run TSE and send it to the resource pool
                # so that it will not get sent to SCA
                resample_functional_to_roi_for_sca = pe.Node(interface=fsl.FLIRT(),
                                                      name='resample_functional_to_roi_for_sca_%d' % num_strat)
                resample_functional_to_roi_for_sca.inputs.interp = 'trilinear'
                resample_functional_to_roi_for_sca.inputs.apply_xfm = True
                resample_functional_to_roi_for_sca.inputs.in_matrix_file = c.identityMatrix
                
                roi_dataflow_for_sca = create_roi_mask_dataflow(c.roiSpecificationFileForSCA, 'roi', 'roi_dataflow_for_sca_%d' % num_strat)
    
                roi_timeseries_for_sca = get_roi_timeseries('roi_timeseries_for_sca_%d' % num_strat)
                roi_timeseries_for_sca.inputs.inputspec.output_type = c.roiTSOutputs

            try:

                if c.roiSpecificationFile != None:

                    node, out_file = strat.get_node_from_resource_pool('functional_mni')
    
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
                
                
                if c.roiSpecificationFileForSCA != None:
                
                    node, out_file = strat.get_node_from_resource_pool('functional_mni')
                
                    # TSE only, not meant for SCA
                    # resample the input functional file to roi
                    workflow.connect(node, out_file,
                                     resample_functional_to_roi_for_sca, 'in_file')
                    workflow.connect(roi_dataflow_for_sca, 'outputspec.out_file',
                                     resample_functional_to_roi_for_sca, 'reference')
    
                    # connect it to the roi_timeseries
                    workflow.connect(roi_dataflow_for_sca, 'outputspec.out_file',
                                     roi_timeseries_for_sca, 'input_roi.roi')
                    workflow.connect(resample_functional_to_roi_for_sca, 'out_file',
                                     roi_timeseries_for_sca, 'inputspec.rest')

            except:
                logConnectionError('ROI Timeseries analysis', num_strat, strat.get_resource_pool(), '0031')
                raise

            if 0 in c.runROITimeseries:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            if c.roiSpecificationFile != None:
                strat.append_name(roi_timeseries.name)
                strat.update_resource_pool({'roi_timeseries' : (roi_timeseries, 'outputspec.roi_outputs')})
                create_log_node(roi_timeseries, 'outputspec.roi_outputs', num_strat)

            if c.roiSpecificationFileForSCA != None:
                strat.append_name(roi_timeseries_for_sca.name)
                strat.update_resource_pool({'roi_timeseries_for_SCA' : (roi_timeseries_for_sca, 'outputspec.roi_outputs')})
                create_log_node(roi_timeseries_for_sca, 'outputspec.roi_outputs', num_strat)

            if (c.roiSpecificationFile != None) or (c.roiSpecificationFileForSCA != None):
                num_strat += 1


    strat_list += new_strat_list



    '''
    Voxel Based Time Series 
    '''

    new_strat_list = []
    num_strat = 0
    if 1 in c.runVoxelTimeseries:


        for strat in strat_list:

            if c.maskSpecificationFile != None:

                resample_functional_to_mask = pe.Node(interface=fsl.FLIRT(),
                                                      name='resample_functional_to_mask_%d' % num_strat)
                resample_functional_to_mask.inputs.interp = 'trilinear'
                resample_functional_to_mask.inputs.apply_xfm = True
                resample_functional_to_mask.inputs.in_matrix_file = c.identityMatrix
    
                mask_dataflow = create_roi_mask_dataflow(c.maskSpecificationFile, 'voxel', 'mask_dataflow_%d' % num_strat)
    
                voxel_timeseries = get_voxel_timeseries('voxel_timeseries_%d' % num_strat)
                voxel_timeseries.inputs.inputspec.output_type = c.voxelTSOutputs
            
            if c.maskSpecificationFileForSCA != None:
            
                resample_functional_to_mask_for_sca = pe.Node(interface=fsl.FLIRT(),
                                                      name='resample_functional_to_mask_for_sca_%d' % num_strat)
                resample_functional_to_mask_for_sca.inputs.interp = 'trilinear'
                resample_functional_to_mask_for_sca.inputs.apply_xfm = True
                resample_functional_to_mask_for_sca.inputs.in_matrix_file = c.identityMatrix
    
                mask_dataflow_for_sca = create_roi_mask_dataflow(c.maskSpecificationFileForSCA, 'voxel', 'mask_dataflow_for_sca_%d' % num_strat)
    
                voxel_timeseries_for_sca = get_voxel_timeseries('voxel_timeseries_for_sca_%d' % num_strat)
                voxel_timeseries_for_sca.inputs.inputspec.output_type = c.voxelTSOutputs
            

            try:

                if c.maskSpecificationFile != None:

                    node, out_file = strat.get_node_from_resource_pool('functional_mni')
    
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
                
                if c.maskSpecificationFileForSCA != None:
                    
                    node, out_file = strat.get_node_from_resource_pool('functional_mni')
                    
                    # resample the input functional file to mask
                    workflow.connect(node, out_file,
                                     resample_functional_to_mask_for_sca, 'in_file')
                    workflow.connect(mask_dataflow_for_sca, 'outputspec.out_file',
                                     resample_functional_to_mask_for_sca, 'reference')
    
                    # connect it to the voxel_timeseries
                    workflow.connect(mask_dataflow_for_sca, 'outputspec.out_file',
                                     voxel_timeseries_for_sca, 'input_mask.mask')
                    workflow.connect(resample_functional_to_mask_for_sca, 'out_file',
                                     voxel_timeseries_for_sca, 'inputspec.rest')
                

            except:
                logConnectionError('Voxel timeseries analysis', num_strat, strat.get_resource_pool(), '0030')
                raise

            if 0 in c.runVoxelTimeseries:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            if c.maskSpecificationFile != None:
                strat.append_name(voxel_timeseries.name)
                strat.update_resource_pool({'voxel_timeseries': (voxel_timeseries, 'outputspec.mask_outputs')})
                create_log_node(voxel_timeseries, 'outputspec.mask_outputs', num_strat)

            if c.maskSpecificationFileForSCA != None:
                strat.append_name(voxel_timeseries_for_sca.name)
                strat.update_resource_pool({'voxel_timeseries_for_SCA': (voxel_timeseries_for_sca, 'outputspec.mask_outputs')})
                create_log_node(voxel_timeseries_for_sca, 'outputspec.mask_outputs', num_strat)

            if (c.maskSpecificationFile != None) or (c.maskSpecificationFileForSCA != None):
                num_strat += 1


    strat_list += new_strat_list



    '''
    Inserting SCA
    Workflow for ROI INPUT
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runSCA and (1 in c.runROITimeseries):
        for strat in strat_list:

            sca_roi = create_sca('sca_roi_%d' % num_strat)


            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 sca_roi, 'inputspec.functional_file')

                node, out_file = strat.get_node_from_resource_pool('roi_timeseries_for_SCA')
                workflow.connect(node, (out_file, extract_one_d),
                                 sca_roi, 'inputspec.timeseries_one_d')
            except:
                logConnectionError('SCA ROI', num_strat, strat.get_resource_pool(), '0032')
                raise


            strat.update_resource_pool({'sca_roi_correlation_stack':(sca_roi, 'outputspec.correlation_stack'),
                                        'sca_roi_correlation_files':(sca_roi, 'outputspec.correlation_files')})
            
            create_log_node(sca_roi, 'outputspec.correlation_stack', num_strat)
            
            strat.append_name(sca_roi.name)
            num_strat += 1
            
    strat_list += new_strat_list



    '''
    Inserting SCA
    Workflow for Voxel INPUT
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runSCA and (1 in c.runVoxelTimeseries):
        for strat in strat_list:

            sca_seed = create_sca('sca_seed_%d' % num_strat)

            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 sca_seed, 'inputspec.functional_file')

                node, out_file = strat.get_node_from_resource_pool('voxel_timeseries_for_SCA')
                workflow.connect(node, (out_file, extract_one_d),
                                 sca_seed, 'inputspec.timeseries_one_d')
            except:
                logConnectionError('SCA', num_strat, strat.get_resource_pool(), '0036')
                raise


            strat.update_resource_pool({'sca_seed_correlation_files':(sca_seed, 'outputspec.correlation_files')})

            strat.append_name(sca_seed.name)
            num_strat += 1

    strat_list += new_strat_list



    '''
    Temporal Regression for Dual Regression
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runDualReg and (1 in c.runSpatialRegression):

        for strat in strat_list:

            dr_temp_reg = create_temporal_reg('temporal_dual_regression_%d' % num_strat)
            dr_temp_reg.inputs.inputspec.normalize = c.mrsNorm
            dr_temp_reg.inputs.inputspec.demean = c.mrsDemean

            try:
                node, out_file = strat.get_node_from_resource_pool('spatial_map_timeseries')
                
                node2, out_file2 = strat.get_leaf_properties()
                node3, out_file3 = strat.get_node_from_resource_pool('functional_brain_mask')

                workflow.connect(node2, out_file2,
                                 dr_temp_reg, 'inputspec.subject_rest')

                workflow.connect(node, out_file,
                                 dr_temp_reg, 'inputspec.subject_timeseries')

                workflow.connect(node3, out_file3,
                                 dr_temp_reg, 'inputspec.subject_mask')

            except:
                logConnectionError('Temporal multiple regression for dual regression', num_strat, strat.get_resource_pool(), '0033')
                raise

            strat.update_resource_pool({'dr_tempreg_maps_stack':(dr_temp_reg, 'outputspec.temp_reg_map'),
                                        'dr_tempreg_maps_files':(dr_temp_reg, 'outputspec.temp_reg_map_files')})
            strat.update_resource_pool({'dr_tempreg_maps_zstat_stack':(dr_temp_reg, 'outputspec.temp_reg_map_z'),
                                        'dr_tempreg_maps_zstat_files':(dr_temp_reg, 'outputspec.temp_reg_map_z_files')})
            
            strat.append_name(dr_temp_reg.name)
            
            create_log_node(dr_temp_reg, 'outputspec.temp_reg_map', num_strat)
            
            num_strat += 1
            
    elif 1 in c.runDualReg and (0 in c.runSpatialRegression):
        logger.info("\n\n" + "WARNING: Dual Regression - Spatial regression was turned off for at least one of the strategies.")
        logger.info("Spatial regression is required for dual regression." + "\n\n")
            
    strat_list += new_strat_list



    '''
    Temporal Regression for SCA
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runMultRegSCA and (1 in c.runROITimeseries):
        for strat in strat_list:

            sc_temp_reg = create_temporal_reg('temporal_regression_sca_%d' % num_strat, which='RT')
            sc_temp_reg.inputs.inputspec.normalize = c.mrsNorm
            sc_temp_reg.inputs.inputspec.demean = c.mrsDemean

            try:
                node, out_file = strat.get_node_from_resource_pool('functional_mni')
                node2, out_file2 = strat.get_node_from_resource_pool('roi_timeseries_for_SCA')
                node3, out_file3 = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')

                workflow.connect(node, out_file,
                                 sc_temp_reg, 'inputspec.subject_rest')

                workflow.connect(node2, (out_file2, extract_txt),
                                 sc_temp_reg, 'inputspec.subject_timeseries')

                workflow.connect(node3, out_file3,
                                 sc_temp_reg, 'inputspec.subject_mask')

            except:
                logConnectionError('Temporal multiple regression for seed based connectivity', num_strat, strat.get_resource_pool(), '0037')
                raise


            strat.update_resource_pool({'sca_tempreg_maps_stack':(sc_temp_reg, 'outputspec.temp_reg_map'),
                                        'sca_tempreg_maps_files':(sc_temp_reg, 'outputspec.temp_reg_map_files')})
            strat.update_resource_pool({'sca_tempreg_maps_zstat_stack':(sc_temp_reg, 'outputspec.temp_reg_map_z'),
                                        'sca_tempreg_maps_zstat_files':(sc_temp_reg, 'outputspec.temp_reg_map_z_files')})
            
            create_log_node(sc_temp_reg, 'outputspec.temp_reg_map', num_strat)
            
            strat.append_name(sc_temp_reg.name)
            num_strat += 1
    strat_list += new_strat_list



    '''
    Inserting Surface Registration
    '''

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runSurfaceRegistraion:
        workflow_bit_id['surface_registration'] = workflow_counter
        for strat in strat_list:

            surface_reg = create_surface_registration('surface_reg_%d' % num_strat)
            surface_reg.inputs.inputspec.recon_subjects = c.reconSubjectsDirectory
            surface_reg.inputs.inputspec.subject_id = subject_id


            try:

                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 surface_reg, 'inputspec.rest')

                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 surface_reg, 'inputspec.brain')

            except:
                logConnectionError('Surface Registration Workflow', num_strat, strat.get_resource_pool(), '0048')
                raise

            if 0 in c.runSurfaceRegistraion:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(surface_reg.name)

            strat.update_resource_pool({'bbregister_registration' : (surface_reg, 'outputspec.out_reg_file'),
                                        'left_hemisphere_surface' :  (surface_reg, 'outputspec.lh_surface_file'),
                                        'right_hemisphere_surface' : (surface_reg, 'outputspec.rh_surface_file')})

            num_strat += 1

    strat_list += new_strat_list



    '''
    Inserting vertices based timeseries
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runVerticesTimeSeries:
        for strat in strat_list:

            vertices_timeseries = get_vertices_timeseries('vertices_timeseries_%d' % num_strat)

            try:

                node, out_file = strat.get_node_from_resource_pool('left_hemisphere_surface')
                workflow.connect(node, out_file,
                                 vertices_timeseries, 'inputspec.lh_surface_file')

                node, out_file = strat.get_node_from_resource_pool('right_hemisphere_surface')
                workflow.connect(node, out_file,
                                 vertices_timeseries, 'inputspec.rh_surface_file')

            except:
                logConnectionError('Vertices Timeseries Extraction', num_strat, strat.get_resource_pool(), '0049')
                raise

            if 0 in c.runVerticesTimeSeries:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(vertices_timeseries.name)

            strat.update_resource_pool({'vertices_timeseries' : (vertices_timeseries, 'outputspec.surface_outputs')})

            num_strat += 1

    strat_list += new_strat_list



    inputnode_fwhm = None
    if c.fwhm != None:

        inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')
        inputnode_fwhm.iterables = ("fwhm", c.fwhm)


    '''
    Inserting Network centrality
    '''

    new_strat_list = []
    num_strat = 0


    if 1 in c.runNetworkCentrality:
        # For each desired strategy
        for strat in strat_list:
            
            # Resample the functional mni to the centrality mask resolution
            resample_functional_to_template = pe.Node(interface=fsl.FLIRT(),
                                                  name='resample_functional_to_template_%d' % num_strat)
            resample_functional_to_template.inputs.interp = 'trilinear'
            resample_functional_to_template.inputs.apply_xfm = True
            resample_functional_to_template.inputs.in_matrix_file = c.identityMatrix

            template_dataflow = create_roi_mask_dataflow(c.templateSpecificationFile,
                                                         'centrality',
                                                         'template_dataflow_%d' % num_strat)

            # Connect in each workflow for the centrality method of interest
            def connectCentralityWorkflow(methodOption,
                                          thresholdOption,
                                          threshold,
                                          weightOptions,
                                          mList):
                # Create centrality workflow
                network_centrality = create_resting_state_graphs(\
                                     c.memoryAllocatedForDegreeCentrality,
                                     'network_centrality_%d-%d' \
                                     %(num_strat,methodOption))
                # Connect registered function input image to inputspec
                workflow.connect(resample_functional_to_template, 'out_file',
                                 network_centrality, 'inputspec.subject')
                # Subject mask/parcellation image
                workflow.connect(template_dataflow, 'outputspec.out_file',
                                 network_centrality, 'inputspec.template')
                # Give which method we're doing (0 - deg, 1 - eig, 2 - lfcd)
                network_centrality.inputs.inputspec.method_option = \
                methodOption
                # Type of threshold (0 - p-value, 1 - sparsity, 2 - corr)
                network_centrality.inputs.inputspec.threshold_option = \
                thresholdOption
                # Connect threshold value (float)
                network_centrality.inputs.inputspec.threshold = threshold
                # List of two booleans, first for binary, second for weighted
                network_centrality.inputs.inputspec.weight_options = \
                weightOptions
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
                
            # Init merge node for appending method output lists to one another
            merge_node = pe.Node(util.Function(input_names=['deg_list',
                                                            'eig_list',
                                                            'lfcd_list'],
                                          output_names = ['merged_list'],
                                          function = merge_lists),
                            name = 'merge_node_%d' % num_strat)
            
            # If we're calculating degree centrality
            if c.degWeightOptions.count(True) > 0:
                connectCentralityWorkflow(0,
                                          c.degCorrelationThresholdOption,
                                          c.degCorrelationThreshold,
                                          c.degWeightOptions,
                                          'deg_list')

            # If we're calculating eigenvector centrality
            if c.eigWeightOptions.count(True) > 0:
                connectCentralityWorkflow(1,
                                          c.eigCorrelationThresholdOption,
                                          c.eigCorrelationThreshold,
                                          c.eigWeightOptions,
                                          'eig_list')
            
            # If we're calculating lFCD
            if c.lfcdWeightOptions.count(True) > 0:
                connectCentralityWorkflow(2,
                                          2,
                                          c.lfcdCorrelationThreshold,
                                          c.lfcdWeightOptions,
                                          'lfcd_list')

            try:

                node, out_file = strat.get_node_from_resource_pool('functional_mni')

                # resample the input functional file to template(roi/mask)
                workflow.connect(node, out_file,
                                 resample_functional_to_template, 'in_file')
                workflow.connect(template_dataflow, 'outputspec.out_file',
                                 resample_functional_to_template, 'reference')
                strat.update_resource_pool({'centrality_outputs' : (merge_node, 'merged_list')})

                # if smoothing is required
                if c.fwhm != None :

                    z_score = get_cent_zscore('centrality_zscore_%d' % num_strat)

                    smoothing = pe.MapNode(interface=fsl.MultiImageMaths(),
                                       name='network_centrality_smooth_%d' % num_strat,
                                       iterfield=['in_file'])

                    zstd_smoothing = pe.MapNode(interface=fsl.MultiImageMaths(),
                                       name='network_centrality_zstd_smooth_%d' % num_strat,
                                       iterfield=['in_file'])


                    # calculate zscores
                    workflow.connect(template_dataflow, 'outputspec.out_file',
                                     z_score, 'inputspec.mask_file')
# workflow.connect(network_centrality, 'outputspec.centrality_outputs',
# z_score, 'inputspec.input_file')
                    workflow.connect(merge_node, 'merged_list',
                                     z_score, 'inputspec.input_file')


                    # connecting raw centrality outputs to smoothing
                    workflow.connect(template_dataflow, 'outputspec.out_file',
                                     smoothing, 'operand_files')
                    workflow.connect(merge_node, 'merged_list',
                                    smoothing, 'in_file')
                    workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                     smoothing, 'op_string')


                    # connecting zscores to smoothing
                    workflow.connect(template_dataflow, 'outputspec.out_file',
                                     zstd_smoothing, 'operand_files')
                    workflow.connect(z_score, 'outputspec.z_score_img',
                                    zstd_smoothing, 'in_file')
                    workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                     zstd_smoothing, 'op_string')

                    strat.append_name(smoothing.name)
                    strat.update_resource_pool({'centrality_outputs_zstd': (z_score, 'outputspec.z_score_img'),
                                                'centrality_outputs_smoothed': (smoothing, 'out_file'),
                                                'centrality_outputs_smoothed_zstd': (zstd_smoothing, 'out_file')})
                    
                    strat.append_name(smoothing.name)
                    create_log_node(smoothing, 'out_file', num_strat)

            except:
                logConnectionError('Network Centrality', num_strat, strat.get_resource_pool(), '0050')
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
     WARP OUTPUTS TO TEMPLATE
    """""""""""""""""""""""""""""""""""""""""""""""""""

    '''
    OUTPUT TO STANDARD
    '''

    def output_to_standard(output_name, output_resource, strat, num_strat,
                           map_node=0, input_image_type=0):

        nodes = getNodeList(strat)
           
        if 'apply_ants_warp_functional_mni' in nodes:

            # ANTS WARP APPLICATION

            fsl_to_itk_convert = create_wf_c3d_fsl_to_itk(map_node, input_image_type, name= \
                    '%s_fsl_to_itk_%d' % (output_name, num_strat))

            collect_transforms = create_wf_collect_transforms(map_node, \
                    name='%s_collect_transforms_%d' \
                    % (output_name, num_strat))

            apply_ants_warp = create_wf_apply_ants_warp(map_node, name= \
                    '%s_to_standard_%d' % (output_name, num_strat))

            apply_ants_warp.inputs.inputspec.dimension = 3
            apply_ants_warp.inputs.inputspec.interpolation = 'Linear'
            apply_ants_warp.inputs.inputspec. \
                    reference_image = c.template_brain_only_for_func

            apply_ants_warp.inputs.inputspec. \
                    input_image_type = input_image_type
                    
            

            try:

                # affine from FLIRT func->anat linear registration
                node, out_file = strat.get_node_from_resource_pool('func' \
                        'tional_to_anat_linear_xfm')
                workflow.connect(node, out_file, fsl_to_itk_convert,
                        'inputspec.affine_file')

                # reference used in FLIRT func->anat linear registration
                node, out_file = strat.get_node_from_resource_pool('anat' \
                        'omical_brain')
                workflow.connect(node, out_file, fsl_to_itk_convert,
                        'inputspec.reference_file')

                # output file to be converted
                node, out_file = strat. \
                        get_node_from_resource_pool(output_resource)
                workflow.connect(node, out_file, fsl_to_itk_convert,
                        'inputspec.source_file')


                # nonlinear warp from anatomical->template ANTS registration
                node, out_file = strat.get_node_from_resource_pool('anat' \
                        'omical_to_mni_nonlinear_xfm')
                workflow.connect(node, out_file, collect_transforms,
                        'inputspec.warp_file')

                # linear initial from anatomical->template ANTS registration
                node, out_file = strat.get_node_from_resource_pool('ants' \
                        '_initial_xfm')
                workflow.connect(node, out_file, collect_transforms,
                        'inputspec.linear_initial')

                # linear affine from anatomical->template ANTS registration
                node, out_file = strat.get_node_from_resource_pool('ants' \
                        '_affine_xfm')
                workflow.connect(node, out_file, collect_transforms,
                        'inputspec.linear_affine')

                # rigid affine from anatomical->template ANTS registration
                node, out_file = strat.get_node_from_resource_pool('ants' \
                        '_rigid_xfm')
                workflow.connect(node, out_file, collect_transforms,
                        'inputspec.linear_rigid')

                # converted FLIRT func->anat affine, now in ITK (ANTS) format
                workflow.connect(fsl_to_itk_convert,
                        'outputspec.itk_transform', collect_transforms,
                        'inputspec.fsl_to_itk_affine')


                # output file to be converted
                node, out_file = strat. \
                        get_node_from_resource_pool(output_resource)
                workflow.connect(node, out_file, apply_ants_warp,
                        'inputspec.input_image')

                # collection of warps to be applied to the output file
                workflow.connect(collect_transforms,
                        'outputspec.transformation_series', apply_ants_warp,
                        'inputspec.transforms')
                        


            except:
                logConnectionError('%s to MNI (ANTS)' % (output_name),
                        num_strat, strat.get_resource_pool(), '0022')
                raise

            strat.update_resource_pool({'%s_to_standard' % (output_name): \
                        (apply_ants_warp, 'outputspec.output_image')})
            strat.append_name(apply_ants_warp.name)
            
            num_strat += 1



        else:

            # FSL WARP APPLICATION

            if map_node == 0:
            
                apply_fsl_warp = pe.Node(interface=fsl.ApplyWarp(),
                        name='%s_to_standard_%d' % (output_name, num_strat))
                          

            elif map_node == 1:
            
                apply_fsl_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                        name='%s_to_standard_%d' % (output_name, num_strat), \
                        iterfield=['in_file'])
                        

            apply_fsl_warp.inputs.ref_file = c.template_skull_for_func


            try:

                # output file to be warped
                node, out_file = strat. \
                        get_node_from_resource_pool(output_resource)
                workflow.connect(node, out_file, apply_fsl_warp, 'in_file')

                # linear affine from func->anat linear FLIRT registration
                node, out_file = strat.get_node_from_resource_pool('func' \
                        'tional_to_anat_linear_xfm')
                workflow.connect(node, out_file, apply_fsl_warp, 'premat')

                # nonlinear warp from anatomical->template FNIRT registration
                node, out_file = strat.get_node_from_resource_pool('anat' \
                        'omical_to_mni_nonlinear_xfm')
                workflow.connect(node, out_file, apply_fsl_warp, 'field_file')
                


            except:
                logConnectionError('%s to MNI (FSL)' % (output_name), \
                        num_strat, strat.get_resource_pool(), '0021')
                raise Exception
                

            strat.update_resource_pool({'%s_to_standard' % (output_name): \
                    (apply_fsl_warp, 'out_file')})
            strat.append_name(apply_fsl_warp.name)
            
            num_strat += 1




    '''
    OUTPUT TO SMOOTH
    '''

    def output_smooth(output_name, output_resource, strat, num_strat, map_node=0):
 
        output_to_standard_smooth = None

        if map_node == 0:
            output_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                    name='%s_smooth_%d' % (output_name, num_strat))


        elif map_node == 1:
            output_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                    name='%s_smooth_%d' % (output_name, num_strat), \
                    iterfield=['in_file'])


        try:

            node, out_file = strat. \
                    get_node_from_resource_pool(output_resource)

            workflow.connect(node, out_file, output_smooth, 'in_file')

            workflow.connect(inputnode_fwhm, ('fwhm', set_gauss), 
                    output_smooth, 'op_string')

            node, out_file = strat. \
                    get_node_from_resource_pool('functional_brain_mask')
            workflow.connect(node, out_file, output_smooth, 'operand_files')


        except:
            logConnectionError('%s smooth' % output_name, num_strat, \
                    strat.get_resource_pool(), '0027')
            raise

        strat.append_name(output_smooth.name)
        strat.update_resource_pool({'%s_smooth' % (output_name): \
                (output_smooth, 'out_file')})


        if 1 in c.runRegisterFuncToMNI:

            if map_node == 0:
                output_to_standard_smooth = pe.Node(interface= \
                        fsl.MultiImageMaths(), name='%s_to_standard_' \
                        'smooth_%d' % (output_name, num_strat))


            elif map_node == 1:
                output_to_standard_smooth = pe.MapNode(interface= \
                        fsl.MultiImageMaths(), name='%s_to_standard_' \
                        'smooth_%d' % (output_name, num_strat), \
                        iterfield=['in_file'])



            try:

                node, out_file = strat.get_node_from_resource_pool('%s_to_' \
                        'standard' % output_name)
                
                workflow.connect(node, out_file, output_to_standard_smooth,
                        'in_file')
                
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                        output_to_standard_smooth, 'op_string')

                node, out_file = strat.get_node_from_resource_pool('func' \
                        'tional_brain_mask_to_standard')
                workflow.connect(node, out_file, output_to_standard_smooth,
                        'operand_files')


            except:
                logConnectionError('%s smooth in MNI' % output_name, \
                        num_strat, strat.get_resource_pool(), '0028')
                raise Exception


            strat.append_name(output_to_standard_smooth.name)
            strat.update_resource_pool({'%s_to_standard_smooth' % \
                    (output_name):(output_to_standard_smooth, 'out_file')})
            create_log_node(output_to_standard_smooth, 'out_file', num_strat)

            
        num_strat += 1




    '''
    z-score standardization functions
    '''

    def z_score_standardize(output_name, output_resource, strat, num_strat, map_node=0):

        z_score_std = get_zscore(output_resource, 'z_score_std_%s_%d' % (output_name, num_strat))

        try:

            node, out_file = strat.get_node_from_resource_pool(output_resource)

            workflow.connect(node, out_file, z_score_std, 'inputspec.input_file')

            # needs the template-space functional mask because we are z-score
            # standardizing outputs that have already been warped to template
            node, out_file = strat. \
                    get_node_from_resource_pool('functional_brain_mask_to_standard')
            workflow.connect(node, out_file, z_score_std, 'inputspec.mask_file')


        except:

            logConnectionError('%s z-score standardize' % output_name, num_strat, \
                    strat.get_resource_pool(), '0127')
            raise

        strat.append_name(z_score_std.name)
        strat.update_resource_pool({'%s_zstd' % (output_resource): \
                (z_score_std, 'outputspec.z_score_img')})



    def fisher_z_score_standardize(output_name, output_resource, timeseries_oned_file, strat, num_strat, map_node=0):

        fisher_z_score_std = get_fisher_zscore(output_resource, map_node, 'fisher_z_score_std_%s_%d' % (output_name, num_strat))

        try:

            node, out_file = strat. \
                    get_node_from_resource_pool(output_resource)

            workflow.connect(node, out_file, fisher_z_score_std, 'inputspec.correlation_file')


            node, out_file = strat. \
                    get_node_from_resource_pool(timeseries_oned_file)
            workflow.connect(node, out_file, fisher_z_score_std, 'inputspec.timeseries_one_d')


        except:

            logConnectionError('%s fisher z-score standardize' % output_name, num_strat, \
                    strat.get_resource_pool(), '0128')
            raise

        strat.append_name(fisher_z_score_std.name)
        strat.update_resource_pool({'%s_fisher_zstd' % (output_resource): \
                (fisher_z_score_std, 'outputspec.fisher_z_score_img')})
                
                
                
                
    '''
    calculate output averages via individual-level mask
    '''
    
    def calc_avg(output_resource, strat, num_strat, map_node=0):
    
        if map_node == 0:
                    
            calc_average = pe.Node(interface=preprocess.Maskave(),
                name='%s_mean_%d' % (output_resource, num_strat))

            mean_to_csv = pe.Node(util.Function(input_names=\
                    ['in_file', 'output_name'],
                    output_names=['output_mean'],
                    function=extract_output_mean),
                    name='%s_mean_to_txt_%d' % (output_resource, \
                    num_strat))
                        
        elif map_node == 1:
            
            calc_average = pe.MapNode(interface=preprocess.Maskave(),
                name='%s_mean_%d' % (output_resource, num_strat), \
                iterfield=['in_file'])

            mean_to_csv = pe.MapNode(util.Function(input_names=\
                    ['in_file', 'output_name'],
                    output_names=['output_mean'],
                    function=extract_output_mean),
                    name='%s_mean_to_txt_%d' % (output_resource, \
                    num_strat), iterfield=['in_file'])
            
        
        mean_to_csv.inputs.output_name = output_resource
        
        
        try:
        
            node, out_file = strat. \
                    get_node_from_resource_pool(output_resource)

            workflow.connect(node, out_file, calc_average, 'in_file')
            
            workflow.connect(calc_average, 'out_file', \
                mean_to_csv, 'in_file')      
        
        
        except:
        
            logConnectionError('%s calc average' % \
                output_name, num_strat, strat.get_resource_pool(), '0128')
            raise
        
        
        strat.append_name(calc_average.name)
        strat.update_resource_pool({'output_means.@%s_average' % (output_resource): (mean_to_csv, 'output_mean')})




    '''
    Transforming Dual Regression outputs to MNI
    '''

    new_strat_list = []
    num_strat = 0

    if (1 in c.runRegisterFuncToMNI) and (1 in c.runDualReg) and (1 in c.runSpatialRegression):
        for strat in strat_list:

            output_to_standard('dr_tempreg_maps_stack', 'dr_tempreg_maps_stack', strat, num_strat, input_image_type=3)

            output_to_standard('dr_tempreg_maps_zstat_stack', 'dr_tempreg_maps_zstat_stack', strat, num_strat, input_image_type=3)

            # dual reg 'files', too
            output_to_standard('dr_tempreg_maps_files', 'dr_tempreg_maps_files', strat, num_strat, 1)
            output_to_standard('dr_tempreg_maps_zstat_files', 'dr_tempreg_maps_zstat_files', strat, num_strat, 1)
             
            num_strat += 1
    
    strat_list += new_strat_list



    '''
    Transforming alff/falff outputs to MNI
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runRegisterFuncToMNI and (1 in c.runALFF):

        for strat in strat_list:

            output_to_standard('alff', 'alff_img', strat, num_strat)
            output_to_standard('falff', 'falff_img', strat, num_strat)
               
            num_strat += 1
    
    strat_list += new_strat_list




    '''
    Transforming ReHo outputs to MNI
    '''
    
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runReHo):
        for strat in strat_list:

            output_to_standard('reho', 'raw_reho_map', strat, num_strat)

            num_strat += 1

    strat_list += new_strat_list



    '''
    Transforming SCA ROI outputs to MNI
    '''
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runSCA) and (1 in c.runROITimeseries):
        for strat in strat_list:

            output_to_standard('sca_roi_stack', 'sca_roi_correlation_stack', strat, num_strat, input_image_type=3)
            output_to_standard('sca_roi_files', 'sca_roi_correlation_files', strat, num_strat, 1)
            
            num_strat += 1

    strat_list += new_strat_list



    '''
    Transforming SCA Voxel outputs to MNI
    '''
    new_strat_list = []
    num_strat = 0

    if 1 in c.runRegisterFuncToMNI and (1 in c.runSCA) and (1 in c.runVoxelTimeseries):
        for strat in strat_list:

            output_to_standard('sca_seed', 'sca_seed_correlation_files', strat, num_strat)
            
            num_strat += 1
    
    strat_list += new_strat_list





    """""""""""""""""""""""""""""""""""""""""""""""""""
     SMOOTHING NORMALIZED OUTPUTS
    """""""""""""""""""""""""""""""""""""""""""""""""""

    '''
    Smoothing Temporal Regression for SCA scores
    '''
    new_strat_list = []
    num_strat = 0

    if (1 in c.runMultRegSCA) and (1 in c.runROITimeseries) and c.fwhm != None:
        for strat in strat_list:

            sc_temp_reg_maps_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='sca_tempreg_maps_stack_smooth_%d' % num_strat, iterfield=['in_file'])
            sc_temp_reg_maps_files_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='sca_tempreg_maps_files_smooth_%d' % num_strat, iterfield=['in_file'])
            sc_temp_reg_maps_Z_stack_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='sca_tempreg_maps_zstat_stack_smooth_%d' % num_strat, iterfield=['in_file'])
            sc_temp_reg_maps_Z_files_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='sca_tempreg_maps_zstat_files_smooth_%d' % num_strat, iterfield=['in_file'])
                             
            '''                 
            sc_temp_output_smooth_average = pe.MapNode(interface=preprocess.Maskave(),
                    name='sca_tempreg_maps_smooth_mean_%d' % num_strat, \
                    iterfield=['in_file'])

            mean_to_csv = pe.MapNode(util.Function(input_names=['in_file', 'output_name'],
                    output_names=['output_mean'],
                    function=extract_output_mean),
                    name='%s_smooth_mean_to_txt_%d' % (output_name, \
                    num_strat), iterfield=['in_file'])
            '''           

            try:
                node, out_file = strat.get_node_from_resource_pool('sca_tempreg_maps_stack')
                node5, out_file5 = strat.get_node_from_resource_pool('sca_tempreg_maps_files')
                node2, out_file2 = strat.get_node_from_resource_pool('sca_tempreg_maps_zstat_stack')
                node3, out_file3 = strat.get_node_from_resource_pool('sca_tempreg_maps_zstat_files')
                node4, out_file4 = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')

                # non-normalized stack
                workflow.connect(node, out_file,
                                 sc_temp_reg_maps_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 sc_temp_reg_maps_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 sc_temp_reg_maps_smooth, 'operand_files')

                # non-normalized files
                workflow.connect(node5, out_file5,
                                 sc_temp_reg_maps_files_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 sc_temp_reg_maps_files_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 sc_temp_reg_maps_files_smooth, 'operand_files')

                # normalized stack
                workflow.connect(node2, out_file2,
                                 sc_temp_reg_maps_Z_stack_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 sc_temp_reg_maps_Z_stack_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 sc_temp_reg_maps_Z_stack_smooth, 'operand_files')

                # normalized files
                workflow.connect(node3, out_file3,
                                 sc_temp_reg_maps_Z_files_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 sc_temp_reg_maps_Z_files_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 sc_temp_reg_maps_Z_files_smooth, 'operand_files')

            except:
                logConnectionError('SCA Temporal regression smooth', num_strat, strat.get_resource_pool(), '0038')
                raise
            strat.append_name(sc_temp_reg_maps_smooth.name)
            strat.update_resource_pool({'sca_tempreg_maps_stack_smooth':(sc_temp_reg_maps_smooth, 'out_file'),
                                        'sca_tempreg_maps_files_smooth':(sc_temp_reg_maps_files_smooth, 'out_file'),
                                        'sca_tempreg_maps_zstat_stack_smooth':(sc_temp_reg_maps_Z_stack_smooth, 'out_file'),
                                        'sca_tempreg_maps_zstat_files_smooth':(sc_temp_reg_maps_Z_files_smooth, 'out_file')})

            create_log_node(sc_temp_reg_maps_smooth, 'out_file', num_strat)
            num_strat += 1
    strat_list += new_strat_list
    
    
    
    '''
    calc averages of sca_tempreg outputs
    '''
    
    new_strat_list = []
    num_strat = 0

    if (1 in c.runMultRegSCA) and (1 in c.runROITimeseries):
    
        for strat in strat_list:
                          
            if 1 in c.runRegisterFuncToMNI:
            
                calc_avg("sca_tempreg_maps_files", strat, num_strat, 1)
                
                if c.fwhm != None:
                
                    calc_avg("sca_tempreg_maps_files_smooth", strat, num_strat, 1)              
            
            num_strat += 1
            
    strat_list += new_strat_list



    '''
    Smoothing Temporal Regression for Dual Regression
    '''
    new_strat_list = []
    num_strat = 0

    if (1 in c.runDualReg) and (1 in c.runSpatialRegression) and c.fwhm != None:
        for strat in strat_list:

            dr_temp_reg_maps_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                                              name='dr_tempreg_maps_stack_smooth_%d' % num_strat)
            dr_temp_reg_maps_Z_stack_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                                              name='dr_tempreg_maps_zstat_stack_smooth_%d' % num_strat)
            dr_temp_reg_maps_files_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='dr_tempreg_maps_files_smooth_%d' % num_strat, iterfield=['in_file'])
            dr_temp_reg_maps_Z_files_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='dr_tempreg_maps_zstat_files_smooth_%d' % num_strat, iterfield=['in_file'])

            try:
                node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_stack_to_standard')
                node2, out_file2 = strat.get_node_from_resource_pool('dr_tempreg_maps_zstat_stack_to_standard')
                node5, out_file5 = strat.get_node_from_resource_pool('dr_tempreg_maps_files_to_standard')
                node3, out_file3 = strat.get_node_from_resource_pool('dr_tempreg_maps_zstat_files_to_standard')
                node4, out_file4 = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')

                # non-normalized stack
                workflow.connect(node, out_file,
                                 dr_temp_reg_maps_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 dr_temp_reg_maps_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 dr_temp_reg_maps_smooth, 'operand_files')

                # normalized stack
                workflow.connect(node2, out_file2,
                                 dr_temp_reg_maps_Z_stack_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 dr_temp_reg_maps_Z_stack_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 dr_temp_reg_maps_Z_stack_smooth, 'operand_files')

                # normalized files
                workflow.connect(node5, out_file5,
                                 dr_temp_reg_maps_files_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 dr_temp_reg_maps_files_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 dr_temp_reg_maps_files_smooth, 'operand_files')

                # normalized z-stat files
                workflow.connect(node3, out_file3,
                                 dr_temp_reg_maps_Z_files_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 dr_temp_reg_maps_Z_files_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 dr_temp_reg_maps_Z_files_smooth, 'operand_files')

            except:
                logConnectionError('Dual regression temp reg smooth', num_strat, strat.get_resource_pool(), '0039')
                raise
            strat.append_name(dr_temp_reg_maps_smooth.name)
            strat.update_resource_pool({'dr_tempreg_maps_stack_to_standard_smooth':(dr_temp_reg_maps_smooth, 'out_file'),
                                        'dr_tempreg_maps_zstat_stack_to_standard_smooth':(dr_temp_reg_maps_Z_stack_smooth, 'out_file'),
                                        'dr_tempreg_maps_files_to_standard_smooth':(dr_temp_reg_maps_files_smooth, 'out_file'),
                                        'dr_tempreg_maps_zstat_files_to_standard_smooth':(dr_temp_reg_maps_Z_files_smooth, 'out_file')})
            create_log_node(dr_temp_reg_maps_smooth, 'out_file', num_strat)
            num_strat += 1
    strat_list += new_strat_list



    '''
    calc averages of dr_tempreg outputs
    '''
    
    new_strat_list = []
    num_strat = 0

    if (1 in c.runDualReg) and (1 in c.runSpatialRegression):
    
        for strat in strat_list:
    
            calc_avg("dr_tempreg_maps_files", strat, num_strat, 1)
            
            #if c.fwhm != None:
                       
            if 1 in c.runRegisterFuncToMNI:
            
                calc_avg("dr_tempreg_maps_files_to_standard", strat, num_strat, 1)
                
                if c.fwhm != None:
                
                    calc_avg("dr_tempreg_maps_files_to_standard_smooth", strat, num_strat, 1)              
            
            num_strat += 1
            
    strat_list += new_strat_list





    '''
    Smoothing motion-corrected functional to MNI output
    '''

    new_strat_list = []
    num_strat = 0
    if (1 in c.runRegisterFuncToMNI) and (1 in c.runFunctionalPreprocessing) and (c.fwhm != None):
        for strat in strat_list:

            output_smooth('motion_correct', 'motion_correct', strat, num_strat)

            num_strat += 1

    strat_list += new_strat_list



    '''    
    Smoothing ALFF fALFF Z scores and or possibly Z scores in MNI 
    '''
    
    new_strat_list = []
    num_strat = 0
    if (1 in c.runALFF) and c.fwhm != None:
        for strat in strat_list:

            output_smooth('alff', 'alff_img', strat, num_strat)
            output_smooth('falff', 'falff_img', strat, num_strat)

            num_strat += 1

    strat_list += new_strat_list



    '''
    calc averages of alff/falff outputs
    '''
    
    new_strat_list = []
    num_strat = 0

    if 1 in c.runALFF:
    
        for strat in strat_list:
    
            calc_avg("alff_img", strat, num_strat)
            calc_avg("falff_img", strat, num_strat)
            
            if c.fwhm != None:
            
                calc_avg("alff_smooth", strat, num_strat)
                calc_avg("falff_smooth", strat, num_strat)
            
            if 1 in c.runRegisterFuncToMNI:
            
                calc_avg("alff_to_standard", strat, num_strat)
                calc_avg("falff_to_standard", strat, num_strat)
                
                if c.fwhm != None:
                
                    calc_avg("alff_to_standard_smooth", strat, num_strat)
                    calc_avg("falff_to_standard_smooth", strat, num_strat)                
            
            num_strat += 1
            
    strat_list += new_strat_list



    '''
    z-standardize alff/falff MNI-standardized outputs
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runZScoring and (1 in c.runALFF):

        for strat in strat_list:

            if c.fwhm != None:
                z_score_standardize('alff_to_standard_smooth', 'alff_to_standard_smooth', strat, num_strat)
                z_score_standardize('falff_to_standard_smooth', 'falff_to_standard_smooth', strat, num_strat)
            
            z_score_standardize('alff_to_standard', 'alff_to_standard', strat, num_strat)
            z_score_standardize('falff_to_standard', 'falff_to_standard', strat, num_strat)

            num_strat += 1

    strat_list += new_strat_list




    '''
    Smoothing ReHo outputs and or possibly ReHo outputs in MNI 
    '''
    
    new_strat_list = []
    num_strat = 0

    if (1 in c.runReHo) and c.fwhm != None:
        for strat in strat_list:

            output_smooth('reho', 'raw_reho_map', strat, num_strat)

            num_strat += 1

    strat_list += new_strat_list



    '''
    calc averages of reho outputs
    '''
    
    new_strat_list = []
    num_strat = 0

    if 1 in c.runReHo:
    
        for strat in strat_list:
    
            calc_avg("raw_reho_map", strat, num_strat)
            
            if c.fwhm != None:
            
                calc_avg("reho_smooth", strat, num_strat)
            
            if 1 in c.runRegisterFuncToMNI:
            
                calc_avg("reho_to_standard", strat, num_strat)
                
                if c.fwhm != None:
                
                    calc_avg("reho_to_standard_smooth", strat, num_strat)              
            
            num_strat += 1
            
    strat_list += new_strat_list




    '''
    z-standardize ReHo MNI-standardized outputs
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runZScoring and (1 in c.runReHo):

        for strat in strat_list:

            if c.fwhm != None:
                z_score_standardize('reho_zstd', 'reho_to_standard_smooth', strat, num_strat)
            
            z_score_standardize('reho', 'reho_to_standard', strat, num_strat)

            num_strat += 1

    strat_list += new_strat_list


    

    '''
    Smoothing SCA roi based Z scores and or possibly Z scores in MNI 
    '''
    if (1 in c.runSCA) and (1 in c.runROITimeseries) and c.fwhm != None:
        for strat in strat_list:

            output_smooth('sca_roi_stack', 'sca_roi_correlation_stack', strat, num_strat)
            output_smooth('sca_roi_files', 'sca_roi_correlation_files', strat, num_strat, 1)
            
            num_strat += 1

    strat_list += new_strat_list
    
    
    
    '''
    calc averages of SCA files outputs
    '''
    
    new_strat_list = []
    num_strat = 0

    if (1 in c.runSCA) and (1 in c.runROITimeseries):
    
        for strat in strat_list:
    
            calc_avg("sca_roi_correlation_files", strat, num_strat, 1)
            
            if c.fwhm != None:
            
                calc_avg("sca_roi_files_smooth", strat, num_strat, 1)
            
            if 1 in c.runRegisterFuncToMNI:
            
                calc_avg("sca_roi_files_to_standard", strat, num_strat, 1)
                
                if c.fwhm != None:
                
                    calc_avg("sca_roi_files_to_standard_smooth", strat, num_strat, 1)              
            
            num_strat += 1
            
    strat_list += new_strat_list
    



    '''
    fisher-z-standardize SCA ROI MNI-standardized outputs
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runZScoring and (1 in c.runSCA) and (1 in c.runROITimeseries):

        for strat in strat_list:

            if c.fwhm != None:
                #fisher_z_score_standardize('sca_roi_stack', 'sca_roi_stack_to_standard_smooth', 'roi_timeseries_for_SCA', strat, num_strat)
                fisher_z_score_standardize('sca_roi_files_fisher_zstd', 'sca_roi_files_to_standard_smooth', 'roi_timeseries_for_SCA', strat, num_strat, 1)
            
            #fisher_z_score_standardize('sca_roi_stack', 'sca_roi_stack_to_standard', 'roi_timeseries_for_SCA', strat, num_strat)
            fisher_z_score_standardize('sca_roi_files', 'sca_roi_files_to_standard', 'roi_timeseries_for_SCA', strat, num_strat, 1)

            num_strat += 1

    strat_list += new_strat_list



    '''
    Smoothing SCA seed based Z scores and or possibly Z scores in MNI 
    '''
    new_strat_list = []
    num_strat = 0

    if (1 in c.runSCA) and (1 in c.runVoxelTimeseries) and c.fwhm != None:
        for strat in strat_list:

            output_smooth('sca_seed', 'sca_seed_correlation_files', strat, num_strat)

            num_strat += 1

    strat_list += new_strat_list
    
    
    
    '''
    calc averages of SCA seed outputs
    '''
    
    new_strat_list = []
    num_strat = 0

    if (1 in c.runSCA) and (1 in c.runVoxelTimeseries):
    
        for strat in strat_list:
    
            calc_avg("sca_seed_correlation_files", strat, num_strat)
            
            if c.fwhm != None:
            
                calc_avg("sca_seed_smooth", strat, num_strat)
            
            if 1 in c.runRegisterFuncToMNI:
            
                calc_avg("sca_seed_to_standard", strat, num_strat)
                
                if c.fwhm != None:
                
                    calc_avg("sca_seed_to_standard_smooth", strat, num_strat)              
            
            num_strat += 1
            
    strat_list += new_strat_list



    '''
    fisher-z-standardize SCA seed MNI-standardized outputs
    '''

    new_strat_list = []
    num_strat = 0

    if 1 in c.runZScoring and (1 in c.runSCA) and (1 in c.runVoxelTimeseries):

        for strat in strat_list:

            if c.fwhm != None:
                fisher_z_score_standardize('sca_seed_fisher_zstd', 'sca_seed_to_standard_smooth', 'voxel_timeseries_for_SCA', strat, num_strat)
            
            fisher_z_score_standardize('sca_seed', 'sca_seed_to_standard', 'voxel_timeseries_for_SCA', strat, num_strat)

            num_strat += 1

    strat_list += new_strat_list





    """""""""""""""""""""""""""""""""""""""""""""""""""
     QUALITY CONTROL
    """""""""""""""""""""""""""""""""""""""""""""""""""


    if 1 in c.generateQualityControlImages:

        #register color palettes
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
    
        hist = pe.Node(util.Function(input_names=['measure_file',
                                                   'measure'],
                                     output_names=['hist_path'],
                                     function=gen_histogram),
                        name='histogram')

        for strat in strat_list:

            nodes = getNodeList(strat)

            #make SNR plot

            if 1 in c.runFunctionalPreprocessing:

                try:

                    hist_ = hist.clone('hist_snr_%d' % num_strat)
                    hist_.inputs.measure = 'snr'

                    drop_percent = pe.Node(util.Function(input_names=['measure_file',
                                                     'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_snr_%d' % num_strat)
                    drop_percent.inputs.percent_ = 99

                    preproc, out_file = strat.get_node_from_resource_pool('preprocessed')
                    brain_mask, mask_file = strat.get_node_from_resource_pool('functional_brain_mask')
                    func_to_anat_xfm, xfm_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    anat_ref, ref_file = strat.get_node_from_resource_pool('anatomical_brain')
                    mfa, mfa_file = strat.get_node_from_resource_pool('mean_functional_in_anat')

                    std_dev = pe.Node(util.Function(input_names=['mask_', 'func_'],
                                                    output_names=['new_fname'],
                                                      function=gen_std_dev),
                                        name='std_dev_%d' % num_strat)

                    std_dev_anat = pe.Node(util.Function(input_names=['func_',
                                                                      'ref_',
                                                                      'xfm_',
                                                                      'interp_'],
                                                         output_names=['new_fname'],
                                                         function=gen_func_anat_xfm),
                                           name='std_dev_anat_%d' % num_strat)

                    snr = pe.Node(util.Function(input_names=['std_dev', 'mean_func_anat'],
                                                output_names=['new_fname'],
                                                function=gen_snr),
                                  name='snr_%d' % num_strat)

                    ###
                    snr_val = pe.Node(util.Function(input_names=['measure_file'],
                                                output_names=['snr_storefl'],
                                                function=cal_snr_val),
                                  name='snr_val%d' % num_strat)


                    std_dev_anat.inputs.interp_ = 'trilinear'

                    montage_snr = create_montage('montage_snr_%d' % num_strat,
                                    'red_to_blue', 'snr')


                    workflow.connect(preproc, out_file,
                                     std_dev, 'func_')

                    workflow.connect(brain_mask, mask_file,
                                     std_dev, 'mask_')

                    workflow.connect(std_dev, 'new_fname',
                                     std_dev_anat, 'func_')

                    workflow.connect(func_to_anat_xfm, xfm_file,
                                     std_dev_anat, 'xfm_')

                    workflow.connect(anat_ref, ref_file,
                                     std_dev_anat, 'ref_')

                    workflow.connect(std_dev_anat, 'new_fname',
                                     snr, 'std_dev')

                    workflow.connect(mfa, mfa_file,
                                     snr, 'mean_func_anat')

                    workflow.connect(snr, 'new_fname',
                                     hist_, 'measure_file')

                    workflow.connect(snr, 'new_fname',
                                     drop_percent, 'measure_file')

                    workflow.connect(snr, 'new_fname',
                                     snr_val, 'measure_file')   ###


                    workflow.connect(drop_percent, 'modified_measure_file',
                                     montage_snr, 'inputspec.overlay')

                    workflow.connect(anat_ref, ref_file,
                                     montage_snr, 'inputspec.underlay')


                    strat.update_resource_pool({'qc___snr_a': (montage_snr, 'outputspec.axial_png'),
                                                'qc___snr_s': (montage_snr, 'outputspec.sagittal_png'),
                                                'qc___snr_hist': (hist_, 'hist_path'),
                                                'qc___snr_val': (snr_val, 'snr_storefl')})   ###
                    if not 3 in qc_montage_id_a:
                        qc_montage_id_a[3] = 'snr_a'
                        qc_montage_id_s[3] = 'snr_s'
                        qc_hist_id[3] = 'snr_hist'

                except:
                    logStandardError('QC', 'unable to get resources for SNR plot', '0051')
                    raise


            #make motion parameters plot

            if 1 in c.runFunctionalPreprocessing:

                try:

                    mov_param, out_file = strat.get_node_from_resource_pool('movement_parameters')
                    mov_plot = pe.Node(util.Function(input_names=['motion_parameters'],
                                                     output_names=['translation_plot',
                                                                   'rotation_plot'],
                                                     function=gen_motion_plt),
                                       name='motion_plt_%d' % num_strat)

                    workflow.connect(mov_param, out_file,
                                     mov_plot, 'motion_parameters')
                    strat.update_resource_pool({'qc___movement_trans_plot': (mov_plot, 'translation_plot'),
                                                'qc___movement_rot_plot': (mov_plot, 'rotation_plot')})

                    if not 6 in qc_plot_id:
                        qc_plot_id[6] = 'movement_trans_plot'

                    if not 7 in qc_plot_id:
                        qc_plot_id[7] = 'movement_rot_plot'


                except:
                    logStandardError('QC', 'unable to get resources for Motion Parameters plot', '0052')
                    raise


            # make FD plot and volumes removed
            if (1 in c.runGenerateMotionStatistics) and ('gen_motion_stats' in nodes):

                try:

                    fd, out_file = strat.get_node_from_resource_pool('frame_wise_displacement')
                    excluded, out_file_ex = strat.get_node_from_resource_pool('scrubbing_frames_excluded')

                    fd_plot = pe.Node(util.Function(input_names=['arr',
                                                                 'ex_vol',
                                                                 'measure'],
                                                    output_names=['hist_path'],
                                                    function=gen_plot_png),
                                      name='fd_plot_%d' % num_strat)
                    fd_plot.inputs.measure = 'FD'
                    workflow.connect(fd, out_file,
                                     fd_plot, 'arr')
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


                montage_skull = create_montage('montage_skull_%d' % num_strat,
                                    'red', 'skull_vis')   ###

                skull_edge = pe.Node(util.Function(input_names=['file_'],
                                                   output_names=['new_fname'],
                                                   function=make_edge),
                                     name='skull_edge_%d' % num_strat)


                workflow.connect(skull, out_file_s,
                                 skull_edge, 'file_')

                workflow.connect(anat_underlay, out_file,
                                 montage_skull, 'inputspec.underlay')

                workflow.connect(skull_edge, 'new_fname',
                                 montage_skull, 'inputspec.overlay')

                strat.update_resource_pool({'qc___skullstrip_vis_a': (montage_skull, 'outputspec.axial_png'),
                                            'qc___skullstrip_vis_s': (montage_skull, 'outputspec.sagittal_png')})

                if not 1 in qc_montage_id_a:
                        qc_montage_id_a[1] = 'skullstrip_vis_a'
                        qc_montage_id_s[1] = 'skullstrip_vis_s'

            except:
                logStandardError('QC', 'Cannot generate QC montages for Skull Stripping: Resources Not Found', '0054')
                raise


            ### make QC montages for mni normalized anatomical image

            try:
                mni_anat_underlay, out_file = strat.get_node_from_resource_pool('mni_normalized_anatomical')

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

                montage_anat = create_montage('montage_anat_%d' % num_strat,
                                    'red', 't1_edge_on_mean_func_in_t1')   ###

                anat_edge = pe.Node(util.Function(input_names=['file_'],
                                                   output_names=['new_fname'],
                                                   function=make_edge),
                                     name='anat_edge_%d' % num_strat)

                workflow.connect(anat, out_file,
                                 anat_edge, 'file_')


                workflow.connect(m_f_a, out_file_mfa,
                                 montage_anat, 'inputspec.underlay')

                workflow.connect(anat_edge, 'new_fname',
                                 montage_anat, 'inputspec.overlay')

                strat.update_resource_pool({'qc___mean_func_with_t1_edge_a': (montage_anat, 'outputspec.axial_png'),
                                            'qc___mean_func_with_t1_edge_s': (montage_anat, 'outputspec.sagittal_png')})

                if not 4 in qc_montage_id_a:
                        qc_montage_id_a[4] = 'mean_func_with_t1_edge_a'
                        qc_montage_id_s[4] = 'mean_func_with_t1_edge_s'


            except:
                logStandardError('QC', 'Cannot generate QC montages for Mean Functional in T1 with T1 edge: Resources Not Found', '0056')
                raise

            # make QC montage for Mean Functional in MNI with MNI edge

            try:
                m_f_i, out_file = strat.get_node_from_resource_pool('mean_functional_in_mni')

                montage_mfi = create_montage('montage_mfi_%d' % num_strat,
                                    'red', 'MNI_edge_on_mean_func_mni')   ###

#                  MNI_edge = pe.Node(util.Function(input_names=['file_'],
#                                                     output_names=['new_fname'],
#                                                     function=make_edge),
#                                       name='MNI_edge_%d' % num_strat)
#                  #MNI_edge.inputs.file_ = c.template_brain_only_for_func
#                 workflow.connect(MNI_edge, 'new_fname',
#                                  montage_mfi, 'inputspec.overlay')

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

                    histogram = hist.clone('hist_%s_%d' % (measure, num_strat))
                    histogram.inputs.measure = measure

                    drop_percent = pe.MapNode(util.Function(input_names=['measure_file',
                                                         'percent_'],
                                           output_names=['modified_measure_file'],
                                           function=drop_percent_),
                                           name='dp_%s_%d' % (measure, num_strat), iterfield=['measure_file'])
                    drop_percent.inputs.percent_ = 99.999

                    overlay, out_file = strat.get_node_from_resource_pool(measure)

                    montage = create_montage('montage_%s_%d' % (measure, num_strat),
                                    'cyan_to_yellow', measure)
                    montage.inputs.inputspec.underlay = c.template_brain_only_for_func

                    workflow.connect(overlay, out_file,
                                     drop_percent, 'measure_file')

                    workflow.connect(drop_percent, 'modified_measure_file',
                                     montage, 'inputspec.overlay')

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
                    print "[!] Creation of QA montages for %s has failed.\n" % measure
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
                            QA_montages('alff_to_standard_smooth_zstd', 11)
                            QA_montages('falff_to_standard_smooth_zstd', 12)

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
                            QA_montages('reho_to_standard_smooth_fisher_zstd', 17)

                        else:
                            QA_montages('reho_to_standard_fisher_zstd', 18)


            '''
            # SCA ROI QA montages
            if (1 in c.runSCA) and (1 in c.runROITimeseries):

                if 1 in c.runRegisterFuncToMNI:
                    QA_montages('sca_roi_to_standard', 19)

                    if c.fwhm != None:
                        QA_montages('sca_roi_to_standard_smooth', 20)

                    if 1 in c.runZScoring:

                        if c.fwhm != None:
                            QA_montages('sca_roi_to_standard_smooth_fisher_zstd', 22)

                        else:
                            QA_montages('sca_roi_to_standard_fisher_zstd', 21)
            '''


            # SCA Seed QA montages
            if (1 in c.runSCA) and (1 in c.runVoxelTimeseries):

                if 1 in c.runRegisterFuncToMNI:
                    QA_montages('sca_seed_to_standard', 23)

                    if c.fwhm != None:
                        QA_montages('sca_seed_to_standard_smooth', 24)

                    if 1 in c.runZScoring:

                        if c.fwhm != None:
                            QA_montages('sca_seed_to_standard_smooth_fisher_zstd', 26)

                        else:
                            QA_montages('sca_seed_to_standard_fisher_zstd', 25)


            # SCA Multiple Regression
            if (1 in c.runMultRegSCA) and (1 in c.runROITimeseries):

                if 1 in c.runRegisterFuncToMNI:
                    QA_montages('sca_tempreg_maps_files', 27)
                    QA_montages('sca_tempreg_maps_zstat_files', 28)

                    if c.fwhm != None:
                        QA_montages('sca_tempreg_maps_files_smooth', 29)
                        QA_montages('sca_tempreg_maps_zstat_files_smooth', 30)


            '''
            # Dual Regression QA montages
            if (1 in c.runDualReg) and (1 in c.runSpatialRegression):

                QA_montages('dr_tempreg_maps_files', 31)
                QA_montages('dr_tempreg_maps_zstat_files', 32)

                if 1 in c.runRegisterFuncToMNI:
                    QA_montages('dr_tempreg_maps_files_to_standard', 33)
                    QA_montages('dr_tempreg_maps_zstat_files_to_standard', 34)

                    if c.fwhm != None:
                        QA_montages('dr_tempreg_maps_files_to_standard_smooth', 35)
                        QA_montages('dr_tempreg_maps_zstat_files_to_standard_smooth', 36)
            '''


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





            '''
            # make QC montages for SCA ROI Smoothed Derivative
            if (1 in c.runSCA) and (1 in c.runROITimeseries):

                hist_sca_roi = hist.clone('hist_sca_roi_%d' % num_strat)
                hist_sca_roi.inputs.measure = 'sca_roi'

                drop_percent_sca_roi = pe.MapNode(util.Function(input_names=['measure_file',
                                                     'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_sca_roi_%d' % num_strat, iterfield=['measure_file'])
                drop_percent_sca_roi.inputs.percent_ = 99.999

                if 1 in c.runZScoring:

                    hist_sca_roi_zstd = hist.clone('hist_sca_roi_zstd_%d' % num_strat)
                    hist_sca_roi_zstd.inputs.measure = 'sca_roi'

                    drop_percent_sca_roi_zstd = pe.MapNode(util.Function(input_names=['measure_file',
                                                         'percent_'],
                                           output_names=['modified_measure_file'],
                                           function=drop_percent_),
                                           name='dp_sca_roi_zstd_%d' % num_strat, iterfield=['measure_file'])
                    drop_percent_sca_roi_zstd.inputs.percent_ = 99.999

                    if c.fwhm != None:

                        sca_roi_smooth_zstd_overlay, out_file = strat.get_node_from_resource_pool('sca_roi_to_standard_smooth_fisher_zstd')
                        montage_sca_roi_smooth_zstd = create_montage('montage_sca_roi_standard_smooth_zstd_%d' % num_strat,
                                        'cyan_to_yellow', 'sca_roi_smooth')

                        montage_sca_roi_smooth_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func

                        workflow.connect(sca_roi_smooth_zstd_overlay, out_file,
                                         drop_percent_sca_roi_zstd, 'measure_file')

                        workflow.connect(drop_percent_sca_roi_zstd, 'modified_measure_file',
                                         montage_sca_roi_smooth_zstd, 'inputspec.overlay')

                        workflow.connect(sca_roi_smooth_zstd_overlay, out_file,
                                         hist_sca_roi_zstd, 'measure_file')
                        strat.update_resource_pool({'qc___sca_roi_smooth_a': (montage_sca_roi_smooth_zstd, 'outputspec.axial_png'),
                                                'qc___sca_roi_smooth_s': (montage_sca_roi_smooth_zstd, 'outputspec.sagittal_png'),
                                                'qc___sca_roi_smooth_hist': (hist_sca_roi_zstd, 'hist_path')})

                        if not 9 in qc_montage_id_a:
                            qc_montage_id_a[9] = 'sca_roi_smooth_a'
                            qc_montage_id_s[9] = 'sca_roi_smooth_s'
                            qc_hist_id[9] = 'sca_roi_smooth_hist'


                    else:

                        sca_roi_zstd_overlay, out_file = strat.get_node_from_resource_pool('sca_roi_to_standard_fisher_zstd')
                        montage_sca_roi_zstd = create_montage('montage_sca_roi_zstd_standard_%d' % num_strat,
                                        'cyan_to_yellow', 'sca_roi')

                        montage_sca_roi_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func
                        workflow.connect(sca_roi_zstd_overlay, out_file,
                                         drop_percent_sca_roi_zstd, 'measure_file')

                        workflow.connect(drop_percent_sca_roi_zstd, 'modified_measure_file',
                                         montage_sca_roi_zstd, 'inputspec.overlay')

                        workflow.connect(sca_roi_zstd_overlay, out_file,
                                         hist_sca_roi_zstd, 'measure_file')

                        strat.update_resource_pool({'qc___sca_roi_a': (montage_sca_roi_zstd, 'outputspec.axial_png'),
                                                'qc___sca_roi_s': (montage_sca_roi_zstd, 'outputspec.sagittal_png'),
                                                'qc___sca_roi_hist': (hist_sca_roi_zstd, 'hist_path')})

                        if not 9 in qc_montage_id_a:
                            qc_montage_id_a[9] = 'sca_roi_a'
                            qc_montage_id_s[9] = 'sca_roi_s'
                            qc_hist_id[9] = 'sca_roi_hist'



                if c.fwhm != None:

                    sca_roi_smooth_overlay, out_file = strat.get_node_from_resource_pool('sca_roi_to_standard_smooth')
                    montage_sca_roi_smooth = create_montage('montage_sca_roi_standard_smooth_%d' % num_strat,
                                    'cyan_to_yellow', 'sca_roi_smooth')

                    montage_sca_roi_smooth.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(sca_roi_smooth_overlay, out_file,
                                     drop_percent_sca_roi, 'measure_file')

                    workflow.connect(drop_percent_sca_roi, 'modified_measure_file',
                                     montage_sca_roi_smooth, 'inputspec.overlay')

                    workflow.connect(sca_roi_smooth_overlay, out_file,
                                     hist_sca_roi, 'measure_file')
                    strat.update_resource_pool({'qc___sca_roi_smooth_a': (montage_sca_roi_smooth, 'outputspec.axial_png'),
                                            'qc___sca_roi_smooth_s': (montage_sca_roi_smooth, 'outputspec.sagittal_png'),
                                            'qc___sca_roi_smooth_hist': (hist_sca_roi, 'hist_path')})

                    if not 9 in qc_montage_id_a:
                        qc_montage_id_a[9] = 'sca_roi_smooth_a'
                        qc_montage_id_s[9] = 'sca_roi_smooth_s'
                        qc_hist_id[9] = 'sca_roi_smooth_hist'


                else:

                    sca_roi_overlay, out_file = strat.get_node_from_resource_pool('sca_roi_to_standard')
                    montage_sca_roi = create_montage('montage_sca_roi_standard_%d' % num_strat,
                                    'cyan_to_yellow', 'sca_roi')

                    montage_sca_roi.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(sca_roi_overlay, out_file,
                                     drop_percent_sca_roi, 'measure_file')

                    workflow.connect(drop_percent_sca_roi, 'modified_measure_file',
                                     montage_sca_roi, 'inputspec.overlay')

                    workflow.connect(sca_roi_overlay, out_file,
                                     hist_sca_roi, 'measure_file')

                    strat.update_resource_pool({'qc___sca_roi_a': (montage_sca_roi, 'outputspec.axial_png'),
                                            'qc___sca_roi_s': (montage_sca_roi, 'outputspec.sagittal_png'),
                                            'qc___sca_roi_hist': (hist_sca_roi, 'hist_path')})

                    if not 9 in qc_montage_id_a:
                        qc_montage_id_a[9] = 'sca_roi_a'
                        qc_montage_id_s[9] = 'sca_roi_s'
                        qc_hist_id[9] = 'sca_roi_hist'




            
            # make QC montages for SCA Smoothed Derivative
            if (1 in c.runSCA) and (1 in c.runVoxelTimeseries):

                hist_sca_seed = hist.clone('hist_sca_seeds_%d' % num_strat)
                hist_sca_seed.inputs.measure = 'sca_seeds'

                drop_percent_sca_seed = pe.MapNode(util.Function(input_names=['measure_file',
                                                     'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_sca_seed_%d' % num_strat, iterfield=['measure_file'])
                drop_percent_sca_seed.inputs.percent_ = 99.999


                if 1 in c.runZScoring:

                    hist_sca_seed_zstd = hist.clone('hist_sca_seeds_zstd_%d' % num_strat)
                    hist_sca_seed_zstd.inputs.measure = 'sca_seeds'

                    drop_percent_sca_seed_zstd = pe.MapNode(util.Function(input_names=['measure_file',
                                                         'percent_'],
                                           output_names=['modified_measure_file'],
                                           function=drop_percent_),
                                           name='dp_sca_seed_zstd_%d' % num_strat, iterfield=['measure_file'])
                    drop_percent_sca_seed_zstd.inputs.percent_ = 99.999

                    if c.fwhm != None:

                        sca_seed_smooth_zstd_overlay, out_file = strat.get_node_from_resource_pool('sca_seed_to_standard_smooth_fisher_zstd')
                        montage_sca_seeds_smooth_zstd = create_montage('montage_seed_standard_smooth_zstd_%d' % num_strat,
                                        'cyan_to_yellow', 'sca_seed_smooth')

                        montage_sca_seeds_smooth_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func
                        workflow.connect(sca_seed_smooth_zstd_overlay, out_file,
                                         drop_percent_sca_seed_zstd, 'measure_file')

                        workflow.connect(drop_percent_sca_seed_zstd, 'modified_measure_file',
                                         montage_sca_seeds_smooth_zstd, 'inputspec.overlay')

                        workflow.connect(sca_seed_smooth_zstd_overlay, out_file,
                                         hist_sca_seed_zstd, 'measure_file')

                        strat.update_resource_pool({'qc___sca_seeds_smooth_a': (montage_sca_seeds_smooth_zstd, 'outputspec.axial_png'),
                                                'qc___sca_seeds_smooth_s': (montage_sca_seeds_smooth_zstd, 'outputspec.sagittal_png'),
                                                'qc___sca_seeds_smooth_hist': (hist_sca_seed_zstd, 'hist_path')})

                        if not 10 in qc_montage_id_a:
                            qc_montage_id_a[10] = 'sca_seeds_smooth_a'
                            qc_montage_id_s[10] = 'sca_seeds_smooth_s'
                            qc_hist_id[10] = 'sca_seeds_smooth_hist'

                    else:
                
                        sca_seed_zstd_overlay, out_file = strat.get_node_from_resource_pool('sca_seed_to_standard_fisher_zstd')
                        montage_sca_seeds_zstd = create_montage('montage_sca_seed_standard_zstd_%d' % num_strat,
                                        'cyan_to_yellow', 'sca_seed')

                        montage_sca_seeds_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func
                        workflow.connect(sca_seed_zstd_overlay, out_file,
                                         drop_percent_sca_seed_zstd, 'measure_file')

                        workflow.connect(drop_percent_sca_seed_zstd, 'modified_measure_file',
                                         montage_sca_seeds_zstd, 'inputspec.overlay')

                        workflow.connect(sca_seed_zstd_overlay, out_file,
                                         hist_sca_seed_zstd, 'measure_file')
                        strat.update_resource_pool({'qc___sca_seeds_a': (montage_sca_seeds_zstd, 'outputspec.axial_png'),
                                                'qc___sca_seeds_s': (montage_sca_seeds_zstd, 'outputspec.sagittal_png'),
                                                'qc___sca_seeds_hist': (hist_sca_seed_zstd, 'hist_path')})

                        if not 10 in qc_montage_id_a:
                            qc_montage_id_a[10] = 'sca_seeds_a'
                            qc_montage_id_s[10] = 'sca_seeds_s'
                            qc_hist_id[10] = 'sca_seeds_hist'



                if c.fwhm != None:

                    sca_seed_smooth_overlay, out_file = strat.get_node_from_resource_pool('sca_seed_to_standard_smooth')
                    montage_sca_seeds_smooth = create_montage('montage_seed_standard_smooth_%d' % num_strat,
                                    'cyan_to_yellow', 'sca_seed_smooth')

                    montage_sca_seeds_smooth.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(sca_seed_smooth_overlay, out_file,
                                     drop_percent_sca_seed, 'measure_file')

                    workflow.connect(drop_percent_sca_seed, 'modified_measure_file',
                                     montage_sca_seeds_smooth, 'inputspec.overlay')

                    workflow.connect(sca_seed_smooth_overlay, out_file,
                                     hist_sca_seed, 'measure_file')

                    strat.update_resource_pool({'qc___sca_seeds_smooth_a': (montage_sca_seeds_smooth, 'outputspec.axial_png'),
                                            'qc___sca_seeds_smooth_s': (montage_sca_seeds_smooth, 'outputspec.sagittal_png'),
                                            'qc___sca_seeds_smooth_hist': (hist_sca_seed, 'hist_path')})

                    if not 10 in qc_montage_id_a:
                        qc_montage_id_a[10] = 'sca_seeds_smooth_a'
                        qc_montage_id_s[10] = 'sca_seeds_smooth_s'
                        qc_hist_id[10] = 'sca_seeds_smooth_hist'


                else:
                
                    sca_seed_overlay, out_file = strat.get_node_from_resource_pool('sca_seed_to_standard')
                    montage_sca_seeds = create_montage('montage_sca_seed_standard_%d' % num_strat,
                                    'cyan_to_yellow', 'sca_seed')

                    montage_sca_seeds.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(sca_seed_overlay, out_file,
                                     drop_percent_sca_seed, 'measure_file')

                    workflow.connect(drop_percent_sca_seed, 'modified_measure_file',
                                     montage_sca_seeds, 'inputspec.overlay')

                    workflow.connect(sca_seed_overlay, out_file,
                                     hist_sca_seed, 'measure_file')
                    strat.update_resource_pool({'qc___sca_seeds_a': (montage_sca_seeds, 'outputspec.axial_png'),
                                            'qc___sca_seeds_s': (montage_sca_seeds, 'outputspec.sagittal_png'),
                                            'qc___sca_seeds_hist': (hist_sca_seed, 'hist_path')})

                    if not 10 in qc_montage_id_a:
                        qc_montage_id_a[10] = 'sca_seeds_a'
                        qc_montage_id_s[10] = 'sca_seeds_s'
                        qc_hist_id[10] = 'sca_seeds_hist'

            


            # make QC montages for Network Centrality
            if 1 in c.runNetworkCentrality:

                hist_ = hist.clone('hist_centrality_%d' % num_strat)
                hist_.inputs.measure = 'centrality'

                drop_percent = pe.MapNode(util.Function(input_names=['measure_file',
                                                     'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_centrality_%d' % num_strat, iterfield=['measure_file'])
                drop_percent.inputs.percent_ = 99.999
                if c.fwhm != None:

                    centrality_overlay, out_file = strat.get_node_from_resource_pool('centrality_outputs_smoothed')
                    montage_centrality = create_montage('montage_centrality_%d' % num_strat,
                                    'cyan_to_yellow', 'centrality')

                    montage_centrality.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(centrality_overlay, out_file,
                                     drop_percent, 'measure_file')

                    workflow.connect(drop_percent, 'modified_measure_file',
                                     montage_centrality, 'inputspec.overlay')

                    workflow.connect(centrality_overlay, out_file,
                                     hist_, 'measure_file')
                    strat.update_resource_pool({'qc___centrality_smooth_a': (montage_centrality, 'outputspec.axial_png'),
                                            'qc___centrality_smooth_s': (montage_centrality, 'outputspec.sagittal_png'),
                                            'qc___centrality_smooth_hist': (hist_, 'hist_path')})
                    if not 11 in qc_montage_id_a:
                        qc_montage_id_a[11] = 'centrality_smooth_a'
                        qc_montage_id_s[11] = 'centrality_smooth_s'
                        qc_hist_id[11] = 'centrality_smooth_hist'



                else:

                    centrality_overlay, out_file = strat.get_node_from_resource_pool('centrality_outputs')
                    montage_centrality = create_montage('montage_centrality_standard_%d' % num_strat,
                                    'cyan_to_yellow', 'centrality')

                    montage_centrality.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(centrality_overlay, out_file,
                                     drop_percent, 'measure_file')

                    workflow.connect(drop_percent, 'modified_measure_file',
                                     montage_centrality, 'inputspec.overlay')

                    workflow.connect(centrality_overlay, out_file,
                                     hist_, 'measure_file')
                    strat.update_resource_pool({'qc___centrality_a': (montage_centrality, 'outputspec.axial_png'),
                                            'qc___centrality_s': (montage_centrality, 'outputspec.sagittal_png'),
                                            'qc___centrality_hist': (hist_, 'hist_path')})
                    if not 11 in qc_montage_id_a:
                        qc_montage_id_a[11] = 'centrality_a'
                        qc_montage_id_s[11] = 'centrality_s'
                        qc_hist_id[11] = 'centrality_hist'


            #QC Montages for MultiReg SCA
            if (1 in c.runMultRegSCA) and (1 in c.runROITimeseries):


                hist_ = hist.clone('hist_dr_sca_%d' % num_strat)
                hist_.inputs.measure = 'temporal_regression_sca'

                drop_percent = pe.MapNode(util.Function(input_names=['measure_file',
                                                      'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_temporal_regression_sca_%d' % num_strat, iterfield=['measure_file'])
                drop_percent.inputs.percent_ = 99.98

                if c.fwhm != None:

                    temporal_regression_sca_overlay, out_file = strat.get_node_from_resource_pool('sca_tempreg_maps_zstat_files_smooth')
                    montage_temporal_regression_sca = create_montage('montage_temporal_regression_sca_%d' % num_strat,
                                      'cyan_to_yellow', 'temporal_regression_sca_smooth')

                    montage_temporal_regression_sca.inputs.inputspec.underlay = c.template_brain_only_for_func
                    strat.update_resource_pool({'qc___temporal_regression_sca_smooth_a': (montage_temporal_regression_sca, 'outputspec.axial_png'),
                                            'qc___temporal_regression_sca_smooth_s': (montage_temporal_regression_sca, 'outputspec.sagittal_png'),
                                            'qc___temporal_regression_sca_smooth_hist': (hist_, 'hist_path')})

                    if not 12 in qc_montage_id_a:
                        qc_montage_id_a[12] = 'temporal_regression_sca_smooth_a'
                        qc_montage_id_s[12] = 'temporal_regression_sca_smooth_s'
                        qc_hist_id[12] = 'temporal_regression_sca_smooth_hist'

                else:
                    temporal_regression_sca_overlay, out_file = strat.get_node_from_resource_pool('sca_tempreg_maps_zstat_files')
                    montage_temporal_regression_sca = create_montage('montage_temporal_regression_sca_%d' % num_strat,
                                      'cyan_to_yellow', 'temporal_regression_sca')

                    montage_temporal_regression_sca.inputs.inputspec.underlay = c.template_brain_only_for_func
                    strat.update_resource_pool({'qc___temporal_regression_sca_a': (montage_temporal_regression_sca, 'outputspec.axial_png'),
                                            'qc___temporal_regression_sca_s': (montage_temporal_regression_sca, 'outputspec.sagittal_png'),
                                            'qc___temporal_regression_sca_hist': (hist_, 'hist_path')})

                    if not 12 in qc_montage_id_a:
                        qc_montage_id_a[12] = 'temporal_regression_sca_a'
                        qc_montage_id_s[12] = 'temporal_regression_sca_s'
                        qc_hist_id[12] = 'temporal_regression_sca_hist'




                workflow.connect(temporal_regression_sca_overlay, out_file,
                                 drop_percent, 'measure_file')

                workflow.connect(drop_percent, 'modified_measure_file',
                                 montage_temporal_regression_sca, 'inputspec.overlay')
                workflow.connect(temporal_regression_sca_overlay, out_file,
                                     hist_, 'measure_file')

            #QC Montages for MultiReg DR
            if (1 in c.runDualReg) and (1 in c.runSpatialRegression):


                hist_ = hist.clone('hist_temp_dr_%d' % num_strat)
                hist_.inputs.measure = 'temporal_dual_regression'

                drop_percent = pe.MapNode(util.Function(input_names=['measure_file',
                                                      'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_temporal_dual_regression_%d' % num_strat, iterfield=['measure_file'])
                drop_percent.inputs.percent_ = 99.98

                if c.fwhm != None:

                    temporal_dual_regression_overlay, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_zstat_files_to_standard_smooth')
                    montage_temporal_dual_regression = create_montage('montage_temporal_dual_regression_%d' % num_strat,
                                      'cyan_to_yellow', 'temporal_dual_regression_smooth')

                    montage_temporal_dual_regression.inputs.inputspec.underlay = c.template_brain_only_for_func
                    strat.update_resource_pool({'qc___temporal_dual_regression_smooth_a': (montage_temporal_dual_regression, 'outputspec.axial_png'),
                                            'qc___temporal_dual_regression_smooth_s': (montage_temporal_dual_regression, 'outputspec.sagittal_png'),
                                            'qc___temporal_dual_regression_smooth_hist': (hist_, 'hist_path')})
                    if not 13 in qc_montage_id_a:
                        qc_montage_id_a[13] = 'temporal_dual_regression_smooth_a'
                        qc_montage_id_s[13] = 'temporal_dual_regression_smooth_s'
                        qc_hist_id[13] = 'temporal_dual_regression_smooth_hist'


                else:
                    temporal_dual_regression_overlay, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_zstat_files_to_standard')
                    montage_temporal_dual_regression = create_montage('montage_temporal_dual_regression_%d' % num_strat,
                                      'cyan_to_yellow', 'temporal_dual_regression')

                    montage_temporal_dual_regression.inputs.inputspec.underlay = c.template_brain_only_for_func
                    strat.update_resource_pool({'qc___temporal_dual_regression_a': (montage_temporal_dual_regression, 'outputspec.axial_png'),
                                            'qc___temporal_dual_regression_s': (montage_temporal_dual_regression, 'outputspec.sagittal_png'),
                                            'qc___temporal_dual_regression_hist': (hist_, 'hist_path')})
                    if not 13 in qc_montage_id_a:
                        qc_montage_id_a[13] = 'temporal_dual_regression_a'
                        qc_montage_id_s[13] = 'temporal_dual_regression_s'
                        qc_hist_id[13] = 'temporal_dual_regression_hist'






                workflow.connect(temporal_dual_regression_overlay, out_file,
                                 drop_percent, 'measure_file')

                workflow.connect(drop_percent, 'modified_measure_file',
                                 montage_temporal_dual_regression, 'inputspec.overlay')
                workflow.connect(temporal_dual_regression_overlay, out_file,
                                     hist_, 'measure_file')

            


            if 1 in c.runVMHC:
                hist_ = hist.clone('hist_vmhc_%d' % num_strat)
                hist_.inputs.measure = 'vmhc'

                drop_percent = pe.Node(util.Function(input_names=['measure_file',
                                                     'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_vmhc%d' % num_strat)
                drop_percent.inputs.percent_ = 99.98

                vmhc_overlay, out_file = strat.get_node_from_resource_pool('vmhc_fisher_zstd_zstat_map')
                montage_vmhc = create_montage('montage_vmhc_%d' % num_strat,
                                  'cyan_to_yellow', 'vmhc_smooth')

                montage_vmhc.inputs.inputspec.underlay = c.template_brain_only_for_func
                workflow.connect(vmhc_overlay, out_file,
                                 drop_percent, 'measure_file')

                workflow.connect(drop_percent, 'modified_measure_file',
                                 montage_vmhc, 'inputspec.overlay')
                workflow.connect(vmhc_overlay, out_file,
                                     hist_, 'measure_file')
                strat.update_resource_pool({'qc___vmhc_smooth_a': (montage_vmhc, 'outputspec.axial_png'),
                                            'qc___vmhc_smooth_s': (montage_vmhc, 'outputspec.sagittal_png'),
                                            'qc___vmhc_smooth_hist': (hist_, 'hist_path')})

                if not 14 in qc_montage_id_a:
                    qc_montage_id_a[14] = 'vmhc_smooth_a'
                    qc_montage_id_s[14] = 'vmhc_smooth_s'
                    qc_hist_id[14] = 'vmhc_smooth_hist'




            if 1 in c.runReHo:
                hist_ = hist.clone('hist_reho_%d' % num_strat)
                hist_.inputs.measure = 'reho'

                drop_percent = pe.Node(util.Function(input_names=['measure_file',
                                                     'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_reho_%d' % num_strat)
                drop_percent.inputs.percent_ = 99.999

                if 1 in c.runZScoring:

                    hist_reho_zstd = hist.clone('hist_reho_zstd_%d' % num_strat)
                    hist_reho_zstd.inputs.measure = 'reho_zstd'

                    drop_percent_zstd = pe.Node(util.Function(input_names=['measure_file',
                                                         'percent_'],
                                           output_names=['modified_measure_file'],
                                           function=drop_percent_),
                                           name='dp_reho_zstd_%d' % num_strat)
                    drop_percent_zstd.inputs.percent_ = 99.999

                    if c.fwhm != None:
                        reho_zstd_overlay, out_file = strat.get_node_from_resource_pool('reho_to_standard_smooth_zstd')
                        montage_reho_zstd = create_montage('montage_reho_zstd_%d' % num_strat,
                                      'cyan_to_yellow', 'reho_standard_smooth_zstd')
                        montage_reho_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func
                        workflow.connect(reho_zstd_overlay, out_file,
                                         hist_reho_zstd, 'measure_file')
                        strat.update_resource_pool({'qc___reho_zstd_smooth_a': (montage_reho_zstd, 'outputspec.axial_png'),
                                                'qc___reho_zstd_smooth_s': (montage_reho_zstd, 'outputspec.sagittal_png'),
                                                'qc___reho_zstd_smooth_hist': (hist_reho_zstd, 'hist_path')})

                        if not 15 in qc_montage_id_a:
                            qc_montage_id_a[15] = 'reho_zstd_smooth_a'
                            qc_montage_id_s[15] = 'reho_zstd_smooth_s'
                            qc_hist_id[15] = 'reho_zstd_smooth_hist'


                    else:
                        reho_zstd_overlay, out_file = strat.get_node_from_resource_pool('reho_to_standard_zstd')
                        montage_reho_zstd = create_montage('montage_reho_zstd_%d' % num_strat,
                                      'cyan_to_yellow', 'reho_standard_zstd')
                        montage_reho_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func
                        workflow.connect(reho_zstd_overlay, out_file,
                                         hist_reho_zstd, 'measure_file')
                        strat.update_resource_pool({'qc___reho_zstd_a': (montage_reho_zstd, 'outputspec.axial_png'),
                                                'qc___reho_zstd_s': (montage_reho_zstd, 'outputspec.sagittal_png'),
                                                'qc___reho_zstd_hist': (hist_reho_zstd, 'hist_path')})

                        if not 15 in qc_montage_id_a:
                            qc_montage_id_a[15] = 'reho_zstd_a'
                            qc_montage_id_s[15] = 'reho_zstd_s'
                            qc_hist_id[15] = 'reho_zstd_hist'


                    workflow.connect(reho_zstd_overlay, out_file,
                                     drop_percent_zstd, 'measure_file')

                    workflow.connect(drop_percent_zstd, 'modified_measure_file',
                                     montage_reho_zstd, 'inputspec.overlay')


                

                if c.fwhm != None:
                    reho_overlay, out_file = strat.get_node_from_resource_pool('reho_to_standard_smooth')
                    montage_reho = create_montage('montage_reho_%d' % num_strat,
                                  'cyan_to_yellow', 'reho_standard_smooth')
                    montage_reho.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(reho_overlay, out_file,
                                     hist_, 'measure_file')
                    strat.update_resource_pool({'qc___reho_smooth_a': (montage_reho, 'outputspec.axial_png'),
                                                'qc___reho_smooth_s': (montage_reho, 'outputspec.sagittal_png'),
                                                'qc___reho_smooth_hist': (hist_, 'hist_path')})

                    if not 15 in qc_montage_id_a:
                        qc_montage_id_a[15] = 'reho_smooth_a'
                        qc_montage_id_s[15] = 'reho_smooth_s'
                        qc_hist_id[15] = 'reho_smooth_hist'


                else:
                    reho_overlay, out_file = strat.get_node_from_resource_pool('reho_to_standard')
                    montage_reho = create_montage('montage_reho_%d' % num_strat,
                                  'cyan_to_yellow', 'reho_standard')
                    montage_reho.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(reho_overlay, out_file,
                                     hist_, 'measure_file')
                    strat.update_resource_pool({'qc___reho_a': (montage_reho, 'outputspec.axial_png'),
                                                'qc___reho_s': (montage_reho, 'outputspec.sagittal_png'),
                                                'qc___reho_hist': (hist_, 'hist_path')})

                    if not 15 in qc_montage_id_a:
                        qc_montage_id_a[15] = 'reho_a'
                        qc_montage_id_s[15] = 'reho_s'
                        qc_hist_id[15] = 'reho_hist'


                workflow.connect(reho_overlay, out_file,
                                 drop_percent, 'measure_file')

                workflow.connect(drop_percent, 'modified_measure_file',
                                 montage_reho, 'inputspec.overlay')



            if 1 in c.runALFF:
                hist_alff = hist.clone('hist_alff_%d' % num_strat)
                hist_alff.inputs.measure = 'alff'

                hist_falff = hist.clone('hist_falff_%d' % num_strat)
                hist_falff.inputs.measure = 'falff'


                drop_percent = pe.Node(util.Function(input_names=['measure_file',
                                                     'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_alff_%d' % num_strat)
                drop_percent.inputs.percent_ = 99.7

                drop_percent_falff = drop_percent.clone('dp_falff_%d' % num_strat)
                drop_percent_falff.inputs.percent_ = 99.999

                if 1 in c.runZScoring:

                    hist_alff_zstd = hist.clone('hist_alff_zstd_%d' % num_strat)
                    hist_alff_zstd.inputs.measure = 'alff_zstd'

                    hist_falff_zstd = hist.clone('hist_falff_zstd_%d' % num_strat)
                    hist_falff_zstd.inputs.measure = 'falff_zstd'


                    drop_percent_zstd = pe.Node(util.Function(input_names=['measure_file',
                                                         'percent_'],
                                           output_names=['modified_measure_file'],
                                           function=drop_percent_),
                                           name='dp_alff_zstd_%d' % num_strat)
                    drop_percent_zstd.inputs.percent_ = 99.7

                    drop_percent_falff_zstd = drop_percent.clone('dp_falff_zstd_%d' % num_strat)
                    drop_percent_falff_zstd.inputs.percent_ = 99.999


                    if c.fwhm != None:
                        alff_zstd_overlay, out_file = strat.get_node_from_resource_pool('alff_to_standard_smooth_zstd')
                        falff_zstd_overlay, out_file_f = strat.get_node_from_resource_pool('falff_to_standard_smooth_zstd')
                        montage_alff_zstd = create_montage('montage_alff_zstd_%d' % num_strat,
                                      'cyan_to_yellow', 'alff_standard_smooth_zstd')
                        montage_alff_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func
                        montage_falff_zstd = create_montage('montage_falff_zstd_%d' % num_strat,
                                      'cyan_to_yellow', 'falff_standard_smooth_zstd')
                        montage_falff_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func
                        workflow.connect(alff_zstd_overlay, out_file,
                                         hist_alff_zstd, 'measure_file')

                        workflow.connect(falff_zstd_overlay, out_file_f,
                                         hist_falff_zstd, 'measure_file')
                        strat.update_resource_pool({'qc___alff_smooth_a': (montage_alff_zstd, 'outputspec.axial_png'),
                                                'qc___alff_smooth_s': (montage_alff_zstd, 'outputspec.sagittal_png'),
                                                'qc___falff_smooth_a': (montage_falff_zstd, 'outputspec.axial_png'),
                                                'qc___falff_smooth_s': (montage_falff_zstd, 'outputspec.sagittal_png'),
                                                'qc___alff_smooth_hist': (hist_alff_zstd, 'hist_path'),
                                                'qc___falff_smooth_hist': (hist_falff_zstd, 'hist_path')})

                        if not 16 in qc_montage_id_a:
                            qc_montage_id_a[16] = 'alff_smooth_a'
                            qc_montage_id_s[16] = 'alff_smooth_s'
                            qc_hist_id[16] = 'alff_smooth_hist'

                        if not 17 in qc_montage_id_a:
                            qc_montage_id_a[17] = 'falff_smooth_a'
                            qc_montage_id_s[17] = 'falff_smooth_s'
                            qc_hist_id[17] = 'falff_smooth_hist'



                    else:
                        alff_zstd_overlay, out_file = strat.get_node_from_resource_pool('alff_to_standard_zstd')
                        falff_zstd_overlay, out_file = strat.get_node_from_resource_pool('falff_to_standard_zstd')
                        montage_alff_zstd = create_montage('montage_alff_zstd_%d' % num_strat,
                                      'cyan_to_yellow', 'alff_standard_zstd')
                        montage_alff_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func
                        montage_falff_zstd = create_montage('montage_falff_zstd_%d' % num_strat,
                                      'cyan_to_yellow', 'falff_standard_zstd')
                        montage_falff_zstd.inputs.inputspec.underlay = c.template_brain_only_for_func
                        workflow.connect(alff_zstd_overlay, out_file,
                                         hist_alff_zstd, 'measure_file')

                        workflow.connect(falff_zstd_overlay, out_file_f,
                                         hist_falff_zstd, 'measure_file')
                        strat.update_resource_pool({'qc___alff_zstd_a': (montage_alff_zstd, 'outputspec.axial_png'),
                                                'qc___alff_zstd_s': (montage_alff_zstd, 'outputspec.sagittal_png'),
                                                'qc___falff_zstd_a': (montage_falff_zstd, 'outputspec.axial_png'),
                                                'qc___falff_zstd_s': (montage_falff_zstd, 'outputspec.sagittal_png'),
                                                'qc___alff_zstd_hist': (hist_alff_zstd, 'hist_path'),
                                                'qc___falff_zstd_hist': (hist_falff_zstd, 'hist_path')})

                        if not 16 in qc_montage_id_a:
                            qc_montage_id_a[16] = 'alff_a'
                            qc_montage_id_s[16] = 'alff_smooth_s'
                            qc_hist_id[16] = 'alff_smooth_hist'

                        if not 16 in qc_montage_id_a:
                            qc_montage_id_a[17] = 'falff_a'
                            qc_montage_id_s[17] = 'falff_s'
                            qc_hist_id[17] = 'falff_hist'


                    workflow.connect(alff_zstd_overlay, out_file,
                                     drop_percent_zstd, 'measure_file')

                    workflow.connect(drop_percent_zstd, 'modified_measure_file',
                                     montage_alff_zstd, 'inputspec.overlay')

                    workflow.connect(falff_zstd_overlay, out_file,
                                     drop_percent_falff_zstd, 'measure_file')

                    workflow.connect(drop_percent_falff_zstd, 'modified_measure_file',
                                     montage_falff_zstd, 'inputspec.overlay')




                if c.fwhm != None:
                    alff_overlay, out_file = strat.get_node_from_resource_pool('alff_to_standard_smooth')
                    falff_overlay, out_file_f = strat.get_node_from_resource_pool('falff_to_standard_smooth')
                    montage_alff = create_montage('montage_alff_%d' % num_strat,
                                  'cyan_to_yellow', 'alff_standard_smooth')
                    montage_alff.inputs.inputspec.underlay = c.template_brain_only_for_func
                    montage_falff = create_montage('montage_falff_%d' % num_strat,
                                  'cyan_to_yellow', 'falff_standard_smooth')
                    montage_falff.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(alff_overlay, out_file,
                                     hist_alff, 'measure_file')

                    workflow.connect(falff_overlay, out_file_f,
                                     hist_falff, 'measure_file')
                    strat.update_resource_pool({'qc___alff_smooth_a': (montage_alff, 'outputspec.axial_png'),
                                            'qc___alff_smooth_s': (montage_alff, 'outputspec.sagittal_png'),
                                            'qc___falff_smooth_a': (montage_falff, 'outputspec.axial_png'),
                                            'qc___falff_smooth_s': (montage_falff, 'outputspec.sagittal_png'),
                                            'qc___alff_smooth_hist': (hist_alff, 'hist_path'),
                                            'qc___falff_smooth_hist': (hist_falff, 'hist_path')})

                    if not 16 in qc_montage_id_a:
                        qc_montage_id_a[16] = 'alff_smooth_a'
                        qc_montage_id_s[16] = 'alff_smooth_s'
                        qc_hist_id[16] = 'alff_smooth_hist'

                    if not 17 in qc_montage_id_a:
                        qc_montage_id_a[17] = 'falff_smooth_a'
                        qc_montage_id_s[17] = 'falff_smooth_s'
                        qc_hist_id[17] = 'falff_smooth_hist'


                else:
                    alff_overlay, out_file = strat.get_node_from_resource_pool('alff_to_standard')
                    falff_overlay, out_file = strat.get_node_from_resource_pool('falff_to_standard')
                    montage_alff = create_montage('montage_alff_%d' % num_strat,
                                  'cyan_to_yellow', 'alff_standard')
                    montage_alff.inputs.inputspec.underlay = c.template_brain_only_for_func
                    montage_falff = create_montage('montage_falff_%d' % num_strat,
                                  'cyan_to_yellow', 'falff_standard')
                    montage_falff.inputs.inputspec.underlay = c.template_brain_only_for_func
                    workflow.connect(alff_overlay, out_file,
                                     hist_alff, 'measure_file')

                    workflow.connect(falff_overlay, out_file_f,
                                     hist_falff, 'measure_file')
                    strat.update_resource_pool({'qc___alff_a': (montage_alff, 'outputspec.axial_png'),
                                            'qc___alff_s': (montage_alff, 'outputspec.sagittal_png'),
                                            'qc___falff_a': (montage_falff, 'outputspec.axial_png'),
                                            'qc___falff_s': (montage_falff, 'outputspec.sagittal_png'),
                                            'qc___alff_hist': (hist_alff, 'hist_path'),
                                            'qc___falff_hist': (hist_falff, 'hist_path')})

                    if not 16 in qc_montage_id_a:
                        qc_montage_id_a[16] = 'alff_a'
                        qc_montage_id_s[16] = 'alff_smooth_s'
                        qc_hist_id[16] = 'alff_smooth_hist'

                    if not 16 in qc_montage_id_a:
                        qc_montage_id_a[17] = 'falff_a'
                        qc_montage_id_s[17] = 'falff_s'
                        qc_hist_id[17] = 'falff_hist'



                workflow.connect(alff_overlay, out_file,
                                 drop_percent, 'measure_file')

                workflow.connect(drop_percent, 'modified_measure_file',
                                 montage_alff, 'inputspec.overlay')

                workflow.connect(falff_overlay, out_file,
                                 drop_percent_falff, 'measure_file')

                workflow.connect(drop_percent_falff, 'modified_measure_file',
                                 montage_falff, 'inputspec.overlay')

            '''


            num_strat += 1
            
                
    logger.info('\n\n' + 'Pipeline building completed.' + '\n\n')



    ###################### end of workflow ###########

    # Run the pipeline only if the user signifies.
    # otherwise, only construct the pipeline (above)
    if run == 1:

        try:
            workflow.write_graph(graph2use='orig')
        except:
            pass
   
   
   
        ## this section creates names for the different branched strategies.
        ## it identifies where the pipeline has forked and then appends the
        ## name of the forked nodes to the branch name in the output directory
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
                lastNodeChar = node[len(node)-1]

                while lastNodeChar == '_' or lastNodeChar == '-' or is_number(lastNodeChar):
                    # make 'renamedNode' the node name with the last character
                    # stripped off, continue this until the _# at the end
                    # of it is gone - does it this way instead of just cutting
                    # off the last two characters in case of a large amount of
                    # strats which can reach double digits
                    renamedNode = renamedNode[:-1]
                    lastNodeChar = renamedNode[len(renamedNode)-1]

                   
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
           
            for fork in forkPoint:

                if 'ants' in fork:
                    forklabel = 'ANTS'
                if 'fnirt' in fork:
                    forklabel = 'FNIRT'
                if 'automask' in fork:
                    forklabel = '3dAutoMask(func)'
                if 'bet' in fork:
                    forklabel = 'BET(func)'
                if 'bbreg' in fork:
                    forklabel = 'bbreg'
                if 'frequency' in fork:
                    forklabel = 'freq-filter'
                if 'nuisance' in fork:
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
        for scanID in sub_dict['rest']:
            scan_ids.append('scan_'+ str(scanID))
        
        pipes = []
        origStrat = 0
        
        for strat in strat_list:
            rp = strat.get_resource_pool()
    
            # build helper dictionary to assist with a clean strategy label for symlinks
    
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
    
    
            if any('median_angle_corr'in name for name in strat.get_name()):
                strategy_tag_helper_symlinks['_target_angle_deg'] = 1
            else:
                strategy_tag_helper_symlinks['_target_angle_deg'] = 0
    
    
            if any('nuisance'in name for name in strat.get_name()):
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
                
                if workflow_bit_id.get(name) != None:
                        strat_tag += name + '_'
                        
                        print name, ' --- ', 2 ** workflow_bit_id[name]
                        hash_val += 2 ** workflow_bit_id[name]

    
            if p_name == None or p_name == 'None':
                
                if forkPointsDict[strat]:
                    pipeline_id = c.pipelineName + forkPointsDict[strat]
                else:
                    pipeline_id = c.pipelineName
                    #if running multiple pipelines with gui, need to change this in future
                    p_name = None

            else:

                if forkPointsDict[strat]:
                    pipeline_id = c.pipelineName + forkPointsDict[strat]
                else:
                    pipeline_id = p_name
                    #if running multiple pipelines with gui, need to change this in future
                    p_name = None
    
            logger.info('strat_tag,  ---- , hash_val,  ---- , pipeline_id: %s, ---- %s, ---- %s' % (strat_tag, hash_val, pipeline_id))
            pip_ids.append(pipeline_id)
            wf_names.append(strat.get_name())
    
            for key in sorted(rp.keys()):
    
                ds = pe.Node(nio.DataSink(), name='sinker_%d' % sink_idx)
                ds.inputs.base_directory = c.outputDirectory
                ds.inputs.container = os.path.join('pipeline_%s' % pipeline_id, subject_id)
                ds.inputs.regexp_substitutions = [(r"/_sca_roi(.)*[/]", '/'),
                                                  (r"/_smooth_centrality_(\d)+[/]", '/'),
                                                  (r"/_z_score(\d)+[/]", "/"),
                                                  (r"/_dr_tempreg_maps_zstat_files_smooth_(\d)+[/]", "/"),
                                                  (r"/_sca_tempreg_maps_zstat_files_smooth_(\d)+[/]", "/"),
                                                  (r"/qc___", '/qc/')]
                node, out_file = rp[key]
                workflow.connect(node, out_file,
                                 ds, key)
                logger.info('node, out_file, key: %s, %s, %s' % (node, out_file, key))
    
    
                link_node = pe.Node(interface=util.Function(input_names=['in_file', 'strategies',
                                        'subject_id', 'pipeline_id', 'helper', 'create_sym_links'],
                                        output_names=[],
                                        function=process_outputs),
                                        name='process_outputs_%d' % sink_idx)

                link_node.inputs.strategies = strategies
                link_node.inputs.subject_id = subject_id
                link_node.inputs.pipeline_id = 'pipeline_%s' % (pipeline_id)
                link_node.inputs.helper = dict(strategy_tag_helper_symlinks)


                if 1 in c.runSymbolicLinks:             
                    link_node.inputs.create_sym_links = True
                else:
                    link_node.inputs.create_sym_links = False

    
                workflow.connect(ds, 'out_file', link_node, 'in_file')

                sink_idx += 1
                logger.info('sink index: %s' % sink_idx)
    

            d_name = os.path.join(c.outputDirectory, ds.inputs.container)
            if not os.path.exists(d_name):
                os.makedirs(d_name)
            
    
            try:
                G = nx.DiGraph()
                strat_name = strat.get_name()
                G.add_edges_from([(strat_name[s], strat_name[s + 1]) for s in range(len(strat_name) - 1)])
                dotfilename = os.path.join(d_name, 'strategy.dot')
                nx.write_dot(G, dotfilename)
                format_dot(dotfilename, 'png')
            except:
                logStandardWarning('Datasink', 'Cannot Create the strategy and pipeline graph, dot or/and pygraphviz is not installed')
                pass
    
    
            logger.info('%s*' % d_name)
            num_strat += 1
            
            pipes.append(pipeline_id)
    
    
        # creates the HTML files used to represent the logging-based status
        create_log_template(pip_ids, wf_names, scan_ids, subject_id, log_dir)
    
    


    
        logger.info('\n\n' + ('Strategy forks: %s' % pipes) + '\n\n')


        pipeline_start_date = strftime("%Y-%m-%d")
        pipeline_start_datetime = strftime("%Y-%m-%d %H:%M:%S")
        pipeline_starttime_string = pipeline_start_datetime.replace(' ','_')
        pipeline_starttime_string = pipeline_starttime_string.replace(':','-')
        
        strat_no = 0
       
        subject_info['resource_pool'] = []

        for strat in strat_list:

            strat_label = 'strat_%d' % strat_no

            subject_info[strat_label] = strat.get_name()

            subject_info['resource_pool'].append(strat.get_resource_pool())

            strat_no += 1


        subject_info['status'] = 'Running'

        subject_info_pickle = open(os.getcwd() + '/subject_info.p', 'wb')

        pickle.dump(subject_info, subject_info_pickle)

        subject_info_pickle.close()
        

        # Actually run the pipeline now, for the current subject
        workflow.run(plugin='MultiProc', plugin_args={'n_procs': c.numCoresPerSubject})
        

        subject_info['status'] = 'Completed'

        subject_info_pickle = open(os.getcwd() + '/subject_info_%s.p' % subject_id , 'wb')

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
    
        subject_dir = os.path.join(c.outputDirectory, 'pipeline_' + pipeline_id, subject_id)

        create_output_mean_csv(subject_dir)


        for count, scanID in enumerate(pip_ids):
            for scan in scan_ids:
                create_log_node(None, None, count, scan).run()
            
            
    
        if 1 in c.generateQualityControlImages:
    
            for pip_id in pip_ids:
    
                f_path = os.path.join(os.path.join(c.outputDirectory, 'pipeline_' + pip_id), subject_id)
    
                f_path = os.path.join(f_path, 'qc_files_here')
    
                generateQCPages(f_path, qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id)
    
    
            ### Automatically generate QC index page
            create_all_qc.run(os.path.join(c.outputDirectory, 'pipeline_' + pip_id))       
        


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

            elapsed_time_data.append(int(((time.time() - pipeline_start_time)/60)))


            # elapsedTimeBin list:
            #  [0] - cumulative elapsed time (minutes) across all subjects
            #  [1] - number of times the elapsed time has been appended
            #        (effectively a measure of how many subjects have run)



            # needs to happen:
                 # write more doc for all this
                 # warning in .csv that some runs may be partial
                 # code to delete .tmp file


            timing_temp_file_path = os.path.join(c.outputDirectory, '%s_pipeline_timing.tmp' % unique_pipeline_id)

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
                pipelineTimeDict['Cores_Per_Subject'] = c.numCoresPerSubject
                pipelineTimeDict['Simultaneous_Subjects'] = c.numSubjectsAtOnce
                pipelineTimeDict['Number_of_Subjects'] = num_subjects
                pipelineTimeDict['Start_Time'] = pipeline_start_stamp
                pipelineTimeDict['End_Time'] = strftime("%Y-%m-%d_%H:%M:%S")
                pipelineTimeDict['Elapsed_Time_(minutes)'] = elapsedTimeBin[0]
                pipelineTimeDict['Status'] = 'Complete'
                
                gpaTimeFields= ['Pipeline', 'Cores_Per_Subject', 'Simultaneous_Subjects', 'Number_of_Subjects', 'Start_Time', 'End_Time', 'Elapsed_Time_(minutes)', 'Status']
                timeHeader = dict((n, n) for n in gpaTimeFields)
                
                timeCSV = open(os.path.join(c.outputDirectory, 'cpac_individual_timing_%s.csv' % c.pipelineName), 'a')
                readTimeCSV = open(os.path.join(c.outputDirectory, 'cpac_individual_timing_%s.csv' % c.pipelineName), 'rb')
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
        
        # Remove working directory when done
        sub_w_path = os.path.join(c.workingDirectory, wfname)
        if c.removeWorkingDir:
            try:
                if os.path.exists(sub_w_path):
                    import shutil
                    logger.info("removing dir -> %s" % sub_w_path)
                    shutil.rmtree(sub_w_path)
            except:
                logStandardWarning('Datasink', ('Couldn\'t remove subjects %s working directory' % wfname))
                pass
        
        endString = ("End of subject workflow %s \n\n" % wfname) + "CPAC run complete:\n" + ("pipeline configuration- %s \n" % c.pipelineName) + \
        ("subject workflow- %s \n\n" % wfname) + ("Elapsed run time (minutes): %s \n\n" % ((time.time() - pipeline_start_time)/60)) + \
        ("Timing information saved in %s/cpac_individual_timing_%s.csv \n" % (c.outputDirectory, c.pipelineName)) + \
        ("System time of start:      %s \n" % pipeline_start_datetime) + ("System time of completion: %s" % strftime("%Y-%m-%d %H:%M:%S"))
    
        logger.info(endString)
    

    return workflow


# Run the prep_workflow function with specific arguments
def run(config, subject_list_file, indx, strategies,
        maskSpecificationFile, roiSpecificationFile, templateSpecificationFile,
        p_name=None, **kwargs):
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
    creds_path : string (optional)
        filepath to the AWS keys credentials file
    bucket_name : string (optional)
        name of the S3 bucket to pull data from
    bucket_prefix : string (optional)
        base directory where S3 inputs are stored and downloaded from
    bucket_upload_prefix : string (optional)
        base directory where local outputs are sent to in S3
    local_prefix : string (optional)
        base directory where the local subject list files were built
    '''

    # Import packages
    import commands
    from CPAC.AWS import fetch_creds
    from CPAC.AWS import aws_utils
    commands.getoutput('source ~/.bashrc')
    import pickle
    import yaml

    # Init variables
    creds_path = kwargs.get('creds_path')
    bucket_name = kwargs.get('bucket_name')
    bucket_prefix = kwargs.get('bucket_prefix')
    bucket_upload_prefix = kwargs.get('bucket_upload_prefix')
    local_prefix = kwargs.get('local_prefix')

    # Import configuration file
    c = Configuration(yaml.load(open(os.path.realpath(config), 'r')))

    # Try and load in the subject list
    try:
        sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        raise Exception ("Subject list is not in proper YAML format. Please check your file")

    # Grab the subject of interest
    sub_dict = sublist[int(indx)-1]
    sub_id = sub_dict['subject_id']

    # Build and download subject's list
    # If we're using AWS
    if creds_path:
        bucket = fetch_creds.return_bucket(creds_path, bucket_name)
        print 'Using data from S3 bucket: %s' % bucket_name
        # Check to see if outputs are already uploaded
        upl_files = [str(k.name) for k in bucket.list(prefix=bucket_upload_prefix)]

        aws_utils.build_download_sublist(bucket,
                                         bucket_prefix,
                                         local_prefix, [sub_dict])
    # Otherwise, state use of local disk and move on
    else:
        print 'Using local disk for input/output'

    # Load in the different spec files to Configuration object
    c.maskSpecificationFile = maskSpecificationFile
    c.roiSpecificationFile = roiSpecificationFile
    c.templateSpecificationFile = templateSpecificationFile

    try:
        # Build and run the pipeline
        prep_workflow(sub_dict, c, pickle.load(open(strategies, 'r')), 1, p_name)
    except Exception as e:
        print 'Could not complete cpac run for subject: %s!' % sub_id
        print 'Error: %s' % e

    # Now upload results to S3
    if creds_path:
        sub_output_dir = os.path.join(c.outputDirectory, 'pipeline_*')
        sub_work_dir = os.path.join(c.workingDirectory, '*_' + sub_id + '_*')
        output_list = aws_utils.collect_subject_files(sub_output_dir,
                                                      sub_id)
        working_list = aws_utils.collect_subject_files(sub_work_dir,
                                                       sub_id)
        dst_list = [o.replace(c.outputDirectory, bucket_upload_prefix)
                    for o in output_list]
        aws_utils.s3_upload(bucket, output_list, dst_list, make_public=True)

        # Delete subject working/output directories
        for wfile in working_list:
            os.remove(wfile)
        for ofile in output_list:
            os.remove(ofile)
