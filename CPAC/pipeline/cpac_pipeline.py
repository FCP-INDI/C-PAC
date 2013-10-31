import os
import sys
import copy
import argparse
import time
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
from   nipype.pipeline.utils import format_dot
import nipype.interfaces.ants as ants
from nipype.interfaces.ants import ApplyTransforms
import nipype.interfaces.c3 as c3
from nipype import config
from nipype import logging
from multiprocessing import Process
import pkg_resources as p
import CPAC
import shutil
from CPAC.anat_preproc.anat_preproc import create_anat_preproc
from CPAC.func_preproc.func_preproc import create_func_preproc
from CPAC.seg_preproc.seg_preproc import create_seg_preproc

from CPAC.registration import create_nonlinear_register, create_register_func_to_anat, create_bbregister_func_to_anat, \
                              create_ants_nonlinear_xfm, create_apply_ants_xfm
from CPAC.nuisance import create_nuisance, bandpass_voxels

from CPAC.median_angle import create_median_angle_correction
from CPAC.generate_motion_statistics import motion_power_statistics
from CPAC.generate_motion_statistics import fristons_twenty_four
from CPAC.scrubbing import create_scrubbing_preproc
from CPAC.timeseries import create_surface_registration, get_voxel_timeseries, \
                            get_roi_timeseries, get_vertices_timeseries, \
                            get_spatial_map_timeseries
from CPAC.network_centrality import create_resting_state_graphs, get_zscore
from CPAC.utils.datasource import *
from CPAC.utils import Configuration, create_all_qc   ### no create_log_template here, move in CPAC/utils/utils.py
from CPAC.qc.qc import create_montage, create_montage_gm_wm_csf
from CPAC.qc.utils import register_pallete, make_edge, drop_percent_, \
                          gen_histogram, gen_plot_png, gen_motion_plt, \
                          gen_std_dev, gen_func_anat_xfm, gen_snr, generateQCPages, cal_snr_val   ###
from CPAC.utils.utils import extract_one_d, set_gauss, \
                             prepare_symbolic_links, \
                             get_scan_params, get_tr, extract_txt, create_log, create_log_template   ###
from CPAC.vmhc.vmhc import create_vmhc
from CPAC.reho.reho import create_reho
from CPAC.alff.alff import create_alff
from CPAC.sca.sca import create_sca, create_temporal_reg
from CPAC.interfaces.afni import preprocess 
import zlib
import linecache
from string import Template

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
            if resource_key in self.resource_pool:
                return self.resource_pool[resource_key]
        except:
            print 'no node for output: ', resource_key
            raise

    def update_resource_pool(self, resources):
        for key, value in resources.items():
            if key in self.resource_pool:
                print 'Warning key %s already exists in resource pool, replacing with %s ' % (key, value)

            self.resource_pool[key] = value

    
def prep_workflow(sub_dict, c, strategies, p_name=None):


    timing = open(os.path.join(c.outputDirectory, 'cpac_pipeline_timing.txt'), 'wt')

    # Start timing here
    pipeline_start_time = time.time()



    print '********************', c.standardResolutionBrain

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


    wfname = 'resting_preproc_' + str(subject_id)
    workflow = pe.Workflow(name=wfname)
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {'hash_method': 'timestamp', 'crashdump_dir': os.path.abspath(c.crashLogDirectory)}
    config.update_config({'logging': {'log_directory': log_dir, 'log_to_file': True}})
    logging.update_logging(config)

    if c.reGenerateOutputs is True:

        import commands
        cmd = "find %s -name \'*sink*\' -exec rm -rf {} \\;" % os.path.join(c.workingDirectory, wfname)
        print cmd
        commands.getoutput(cmd)
        cmd = "find %s -name \'*link*\' -exec rm -rf {} \\;" % os.path.join(c.workingDirectory, wfname)
        print cmd
        commands.getoutput(cmd)
        cmd = "find %s -name \'*log*\' -exec rm -rf {} \\;" % os.path.join(c.workingDirectory, wfname)
        print cmd
        commands.getoutput(cmd)

    mflow = None
    pflow = None

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

    strat_list = []

    workflow_bit_id = {}
    workflow_counter = 0

    """
    Initialize Anatomical Input Data Flow
    """
    new_strat_list = []
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

    print strat_list

    """
    Inserting Anatomical Preprocessing workflow
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runAnatomicalPreprocessing:

        workflow_bit_id['anat_preproc'] = workflow_counter

        for strat in strat_list:
            # create a new node, Remember to change its name!
            anat_preproc = create_anat_preproc().clone('anat_preproc_%d' % num_strat)


            try:
                # connect the new node to the previous leaf
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file, anat_preproc, 'inputspec.anat')

            except:
                print 'Invalid Connection: Anat Preprocessing No valid Previous for strat : ', num_strat, ' resource_pool: ', strat.get_resource_pool()
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



    """
    T1 -> Template, Non-linear registration (FNIRT or ANTS)
    """
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    
    # either run FSL anatomical-to-MNI registration, or...
    
    if 1 in c.runRegistrationPreprocessing:

        workflow_bit_id['anat_mni_register'] = workflow_counter
        for strat in strat_list:

            if 'FSL' in c.regOption:

                fnirt_reg_anat_mni = create_nonlinear_register('anat_mni_fnirt_register_%d' % num_strat)

                try:
                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_mni, 'inputspec.input_brain')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_mni, 'inputspec.input_skull')

                    # pass the reference files                
                    fnirt_reg_anat_mni.inputs.inputspec.reference_brain = c.standardResolutionBrainAnat
                    fnirt_reg_anat_mni.inputs.inputspec.reference_skull = c.standardAnat

                    # assign the FSL FNIRT config file specified in pipeline config.yml
                    fnirt_reg_anat_mni.inputs.inputspec.fnirt_config = c.fnirtConfig
                

                except:
                    print 'Invalid Connection: Anatomical Registration:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                if 0 in c.runRegistrationPreprocessing:
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


                """
                ### add mni normalized anatomical to qc output html page
                strat.update_resource_pool({'qc___mni_normalized_anatomical': (reg_anat_mni, 'outputspec.output_brain')})

                if not 6 in qc_montage_id_a:
                    qc_montage_id_a[6] = 'mni_normalized_anatomical'
                    ###qc_montage_id_s[6] = 'skullstrip_vis_s'
                """


                create_log_node(fnirt_reg_anat_mni, 'outputspec.output_brain', num_strat)
            
            
                num_strat += 1
            
            
            # or run ANTS anatomical-to-MNI registration instead        
            elif 'ANTS' in c.regOption:

                ants_reg_anat_mni = create_ants_nonlinear_xfm('anat_mni_ants_register_%d' % num_strat)

                try:
                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file,
                                     ants_reg_anat_mni, 'inputspec.anatomical_brain')

                    # pass the reference file           
                    ants_reg_anat_mni.inputs.inputspec.reference_brain = c.standardResolutionBrainAnat
                

                except:
                    print 'Invalid Connection: Anatomical Registration:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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
                strat.set_leaf_properties(ants_reg_anat_mni, 'outputspec.output_brain')

                strat.update_resource_pool({#'anatomical_to_mni_linear_xfm':(reg_anat_mni, 'outputspec.linear_xfm'),
                                            #'anatomical_to_mni_nonlinear_xfm':(reg_anat_mni, 'outputspec.nonlinear_xfm'),
                                            #'mni_to_anatomical_linear_xfm':(reg_anat_mni, 'outputspec.invlinear_xfm'),
                                            'ants_affine_xfm':(ants_reg_anat_mni, 'outputspec.affine_transformation'),
                                            #'ants_warp_field':(ants_reg_anat_mni, 'outputspec.warp_field'),
                                            # changed 'ants_warp_field' to 'anatomical_to_mni_nonlinear_xfm' because
                                            # it is the warp field output of ANTS, and FLIRT/FNIRT also outputs a warp
                                            # field named this. need to verify that these work interchangeably
                                            'anatomical_to_mni_nonlinear_xfm':(ants_reg_anat_mni, 'outputspec.warp_field'),
                                            'mni_to_anatomical_linear_xfm':(ants_reg_anat_mni, 'outputspec.inverse_warp'), #<---- this is the mni to anatomical NONLINEAR xfm
                                            'mni_normalized_anatomical':(ants_reg_anat_mni, 'outputspec.output_brain')})

                create_log_node(ants_reg_anat_mni, 'outputspec.output_brain', num_strat)  
          
                num_strat += 1
            
    strat_list += new_strat_list
    
    
 
    """
    Inserting Segmentation Preprocessing
    Workflow
    """

    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runSegmentationPreprocessing:
        workflow_bit_id['seg_preproc'] = workflow_counter
        for strat in strat_list:

            if 'FSL' in c.regOption:
                seg_preproc = create_seg_preproc(False, 'seg_preproc_%d' % num_strat)
            elif 'ANTS' in c.regOption:
                seg_preproc = create_seg_preproc(True, 'seg_preproc_%d' % num_strat)

            try:
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 seg_preproc, 'inputspec.brain')

                if 'FSL' in c.regOption:
                    node, out_file = strat.get_node_from_resource_pool('mni_to_anatomical_linear_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc, 'inputspec.standard2highres_mat')
                elif 'ANTS' in c.regOption:
                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     seg_preproc, 'inputspec.standard2highres_mat')                


                seg_preproc.inputs.inputspec.PRIOR_CSF = c.PRIOR_CSF
                seg_preproc.inputs.inputspec.PRIOR_GRAY = c.PRIOR_GRAY
                seg_preproc.inputs.inputspec.PRIOR_WHITE = c.PRIOR_WHITE

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
                print 'Invalid Connection: Segmentation Preprocessing:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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

    """
    Inserting Functional Input Data workflow
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runFunctionalDataGathering:
        for strat in strat_list:
            # create a new node, Remember to change its name!
            # Flow = create_func_datasource(sub_dict['rest'])
            # Flow.inputs.inputnode.subject = subject_id
            funcFlow = create_func_datasource(sub_dict['rest'], 'func_gather_%d' % num_strat)
            funcFlow.inputs.inputnode.subject = subject_id

            if 0 in c.runFunctionalDataGathering:
                # we are forking so create a new node
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.set_leaf_properties(funcFlow, 'outputspec.rest')

            num_strat += 1

    strat_list += new_strat_list


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

            slice_timing = sub_dict.get('scan_parameters')
            # a node which checks if scan _parameters are present for each scan
            scan_params = pe.Node(util.Function(input_names=['subject',
                                                             'scan',
                                                             'subject_map',
                                                             'start_indx',
                                                             'stop_indx'],
                                               output_names=['tr',
                                                             'tpattern',
                                                             'ref_slice',
                                                             'start_indx',
                                                             'stop_indx'],
                                               function=get_scan_params),
                                 name='scan_params_%d' % num_strat)

            convert_tr = pe.Node(util.Function(input_names=['tr'],
                                               output_names=['tr'],
                                               function=get_tr),
                                 name='convert_tr_%d' % num_strat)

            # if scan parameters are available slice timing correction is
            # turned on
            if slice_timing:

                func_preproc = create_func_preproc(slice_timing_correction=True, wf_name='func_preproc_%d' % num_strat)

                # getting the scan parameters
                workflow.connect(funcFlow, 'outputspec.subject',
                                 scan_params, 'subject')
                workflow.connect(funcFlow, 'outputspec.scan',
                                 scan_params, 'scan')
                scan_params.inputs.subject_map = sub_dict
                scan_params.inputs.start_indx = c.startIdx
                scan_params.inputs.stop_indx = c.stopIdx

                # passing the slice timing correction parameters
                workflow.connect(scan_params, 'tr',
                                 func_preproc, 'scan_params.tr')
                workflow.connect(scan_params, 'ref_slice',
                                 func_preproc, 'scan_params.ref_slice')
                workflow.connect(scan_params, 'tpattern',
                                 func_preproc, 'scan_params.acquisition')

                workflow.connect(scan_params, 'start_indx',
                                 func_preproc, 'inputspec.start_idx')
                workflow.connect(scan_params, 'stop_indx',
                                 func_preproc, 'inputspec.stop_idx')

                workflow.connect(scan_params, 'tr',
                                 convert_tr, 'tr')
            else:
                func_preproc = create_func_preproc(slice_timing_correction=False, wf_name='func_preproc_%d' % num_strat)
                func_preproc.inputs.inputspec.start_idx = c.startIdx
                func_preproc.inputs.inputspec.stop_idx = c.stopIdx

                convert_tr.inputs.tr = c.TR

            node = None
            out_file = None
            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file, func_preproc, 'inputspec.rest')

            except:
                print 'Invalid Connection: Functional Preprocessing No valid Previous for strat : ', num_strat, ' resource_pool: ', strat.get_resource_pool()
                num_strat += 1
                continue

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
            if slice_timing:
                strat.update_resource_pool({'slice_timing_corrected': (func_preproc, 'outputspec.slice_time_corrected')})
            strat.update_resource_pool({'mean_functional':(func_preproc, 'outputspec.example_func')})
            strat.update_resource_pool({'functional_preprocessed_mask':(func_preproc, 'outputspec.preprocessed_mask')})
            strat.update_resource_pool({'movement_parameters':(func_preproc, 'outputspec.movement_parameters')})
            strat.update_resource_pool({'max_displacement':(func_preproc, 'outputspec.max_displacement')})
            #strat.update_resource_pool({'xform_matrix':(func_preproc, 'outputspec.xform_matrix')})
            strat.update_resource_pool({'preprocessed':(func_preproc, 'outputspec.preprocessed')})
            strat.update_resource_pool({'functional_brain_mask':(func_preproc, 'outputspec.mask')})
            strat.update_resource_pool({'motion_correct':(func_preproc, 'outputspec.motion_correct')})


            create_log_node(func_preproc, 'outputspec.preprocessed', num_strat)
            num_strat += 1

            
    strat_list += new_strat_list



    """
    Inserting Friston's 24 parameter  Workflow
    Incase this workflow runs , it overwrites the movement_parameters file
    So the file contains 24 parameters for motion and that gets wired to all the workflows
    that depend on. The effect should be seen when regressing out nuisance signals and motion
    is used as one of the regressors
    """
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

            except:
                print 'Invalid Connection: fristons_parameter_model ', num_strat, ' resource_pool: ', strat.get_resource_pool()
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



    """
    Inserting Anatomical to Functional Registration
    """
    """
    new_strat_list = []
    num_strat = 0

    ### REDUNDANT REGISTRATION!!!  ---  DOUBLE CHECK

    workflow_counter += 1
    if 1 in c.runAnatomicalToFunctionalRegistration:
        workflow_bit_id['anat_to_func_register'] = workflow_counter
        for strat in strat_list:
            anat_to_func_reg = pe.Node(interface=fsl.FLIRT(),
                               name='anat_to_func_register_%d' % num_strat)
            anat_to_func_reg.inputs.cost = 'corratio'
            anat_to_func_reg.inputs.dof = 6

            inv_anat_to_func = pe.Node(interface=fsl.ConvertXFM(),
                                       name='inv_anat_to_func_register_%d' % num_strat)
            inv_anat_to_func.inputs.invert_xfm = True

            func_gm = pe.Node(interface=fsl.ApplyXfm(),
                               name='func_gm_%d' % num_strat)
            func_gm.inputs.apply_xfm = True
            func_gm.inputs.interp = 'nearestneighbour'

            func_csf = pe.Node(interface=fsl.ApplyXfm(),
                               name='func_csf_%d' % num_strat)
            func_csf.inputs.apply_xfm = True
            func_csf.inputs.interp = 'nearestneighbour'

            func_wm = pe.Node(interface=fsl.ApplyXfm(),
                               name='func_wm_%d' % num_strat)
            func_wm.inputs.apply_xfm = True
            func_wm.inputs.interp = 'nearestneighbour'

            try:
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 anat_to_func_reg, 'in_file')
                workflow.connect(anat_to_func_reg, 'out_matrix_file',
                                 inv_anat_to_func, 'in_file')

                node, out_file = strat.get_node_from_resource_pool('mean_functional')
                workflow.connect(node, out_file,
                                 anat_to_func_reg, 'reference')
                workflow.connect(node, out_file,
                                 func_gm, 'reference')
                workflow.connect(node, out_file,
                                 func_csf, 'reference')
                workflow.connect(node, out_file,
                                 func_wm, 'reference')

                node, out_file = strat.get_node_from_resource_pool('anatomical_gm_mask')
                workflow.connect(node, out_file,
                                 func_gm, 'in_file')
                workflow.connect(anat_to_func_reg, 'out_matrix_file',
                                 func_gm, 'in_matrix_file')

                node, out_file = strat.get_node_from_resource_pool('anatomical_csf_mask')
                workflow.connect(node, out_file,
                                 func_csf, 'in_file')
                workflow.connect(anat_to_func_reg, 'out_matrix_file',
                                 func_csf, 'in_matrix_file')

                node, out_file = strat.get_node_from_resource_pool('anatomical_wm_mask')
                workflow.connect(node, out_file,
                                 func_wm, 'in_file')
                workflow.connect(anat_to_func_reg, 'out_matrix_file',
                                 func_wm, 'in_matrix_file')

            except:
                print 'Invalid Connection: Anatomical to Functional Registration:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise

            if 0 in c.runAnatomicalToFunctionalRegistration:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(anat_to_func_reg.name)

            strat.update_resource_pool({'anatomical_to_functional_xfm':(anat_to_func_reg, 'out_matrix_file'),
                                        'inverse_anatomical_to_functional_xfm':(inv_anat_to_func, 'out_file'),
                                        'functional_gm_mask':(func_gm, 'out_file'),
                                        'functional_wm_mask':(func_wm, 'out_file'),
                                        'functional_csf_mask':(func_csf, 'out_file')})
            
            create_log_node(anat_to_func_reg, 'out_matrix_file', num_strat)

            num_strat += 1

    strat_list += new_strat_list
    """

    """
    Func -> T1 Registration (Initial Linear reg)
    """

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

                # Input functional image (func.nii.gz)
                node, out_file = strat.get_node_from_resource_pool('mean_functional')
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.func')

                # Input skull-stripped anatomical (anat.nii.gz)
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.anat')

   

            except:
                print 'Invalid Connection: Register Functional to Anatomical (pre BBReg):', num_strat, ' resource_pool: ', strat.get_resource_pool()
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
                                        #'anatomical_wm_edge':(func_to_anat, 'outputspec.anat_wm_edge'),
                                        'functional_to_anat_linear_xfm':(func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg')})
                                        #'functional_to_mni_linear_xfm':(func_to_anat, 'outputspec.func_to_anat_linear_xfm'),
                                        #'mni_to_functional_linear_xfm':(func_to_anat, 'outputspec.mni_to_func_linear_xfm')})

            # Outputs:
            # functional_to_anat_linear_xfm = func-t1.mat, linear, sent to 'premat' of post-FNIRT applywarp,
            #                                 or to the input of the post-ANTS c3d_affine_tool
            
            #create_log_node(func_to_anat, 'outputspec.mni_func', num_strat)
            num_strat += 1

    strat_list += new_strat_list




    """
    Func -> T1 Registration (BBREG)
    """

    # Outputs 'functional_to_anat_linear_xfm', a matrix file of the functional-to-anatomical
    # registration warp to be applied LATER in func_mni_warp, which accepts it as input 'premat'

    new_strat_list = []
    num_strat = 0
    workflow_counter += 1
    
    if 1 in c.runBBReg:
        workflow_bit_id['func_to_anat_bbreg'] = workflow_counter
        for strat in strat_list:
            func_to_anat_bbreg = create_bbregister_func_to_anat('func_to_anat_bbreg_%d' % num_strat)
       
            # Input registration parameters
            func_to_anat_bbreg.inputs.inputspec.bbr_schedule = c.boundaryBasedRegistrationSchedule

            try:
                def pick_wm(seg_prob_list):
                    seg_prob_list.sort()
                    return seg_prob_list[-1]

                # Input functional image (func.nii.gz)
                node, out_file = strat.get_node_from_resource_pool('mean_functional')
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
                print 'Invalid Connection: Register Functional to Anatomical (BBReg):', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise

            if 0 in c.runBBReg:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(func_to_anat.name)
            # strat.set_leaf_properties(func_mni_warp, 'out_file')

            strat.update_resource_pool({#'mean_functional_in_mni':(copy_tr, 'out_file'),
                                        'mean_functional_in_anat':(func_to_anat_bbreg, 'outputspec.anat_func'),
                                        #'anatomical_wm_edge':(func_to_anat, 'outputspec.anat_wm_edge'),
                                        'functional_to_anat_linear_xfm':(func_to_anat_bbreg, 'outputspec.func_to_anat_linear_xfm')})
                                        #'functional_to_mni_linear_xfm':(func_to_anat, 'outputspec.func_to_anat_linear_xfm'),
                                        #'mni_to_functional_linear_xfm':(func_to_anat, 'outputspec.mni_to_func_linear_xfm')})

            # Outputs:
            # functional_to_anat_linear_xfm = func-t1.mat, linear, sent to 'premat' of post-FNIRT applywarp,
            #                                 or to the input of the post-ANTS c3d_affine_tool
            
            #create_log_node(func_to_anat, 'outputspec.mni_func', num_strat)
            num_strat += 1

    strat_list += new_strat_list


    """
    Func -> T1 Registration (BBREG), this has an ApplyWarp (mni warp) which
    takes in the premat, for some reason - may not actually go anywhere
    """
    """
    # Outputs 'functional_to_anat_linear_xfm', a matrix file of the functional-to-anatomical
    # registration warp to be applied LATER in func_mni_warp, which accepts it as input 'premat'

    new_strat_list = []
    num_strat = 0
    workflow_counter += 1
    
    if 1 in c.runRegisterFuncToAnat:
        workflow_bit_id['func_to_anat'] = workflow_counter
        for strat in strat_list:
            func_to_anat = create_bbregister_func_to_anat('func_to_anat_%d' % num_strat)
       
            # Input registration parameters
            func_to_anat.inputs.inputspec.mni = c.standardResolutionBrain
            func_to_anat.inputs.inputspec.interp = 'trilinear'
            func_to_anat.inputs.inputspec.bbr_schedule = c.boundaryBasedRegistrationSchedule

            try:
                def pick_wm(seg_prob_list):
                    seg_prob_list.sort()
                    return seg_prob_list[-1]

                # Input functional image (func.nii.gz)
                node, out_file = strat.get_node_from_resource_pool('mean_functional')
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.func')

                # Input segmentation probability maps for white matter segmentation
                node, out_file = strat.get_node_from_resource_pool('seg_probability_maps')
                workflow.connect(node, (out_file, pick_wm),
                                 func_to_anat, 'inputspec.anat_wm_segmentation')

                # Input anatomical whole-head image
                node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.anat_skull')

                # Input skull-stripped anatomical (anat.nii.gz)
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.anat')


                # Necessary???
                #node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_linear_xfm')
                #workflow.connect(node, out_file,
                # func_to_anat, 'inputspec.anat_to_mni_linear_xfm')

                # Nonlinear transform - WHY IS THIS HERE?
                node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.anat_to_mni_nonlinear_xfm')
   

            except:
                print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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

            strat.update_resource_pool({#'mean_functional_in_mni':(copy_tr, 'out_file'),
                                        'mean_functional_in_mni':(func_to_anat, 'outputspec.mni_func'),
                                        'mean_functional_in_anat':(func_to_anat, 'outputspec.anat_func'),
                                        #'anatomical_wm_edge':(func_to_anat, 'outputspec.anat_wm_edge'),
                                        'functional_to_anat_linear_xfm':(func_to_anat, 'outputspec.func_to_anat_linear_xfm')})
                                        #'functional_to_mni_linear_xfm':(func_to_anat, 'outputspec.func_to_anat_linear_xfm'),
                                        #'mni_to_functional_linear_xfm':(func_to_anat, 'outputspec.mni_to_func_linear_xfm')})

            # Outputs:
            # functional_to_anat_linear_xfm = func-t1.mat, linear, sent to 'premat' of post-FNIRT applywarp,
            # or to the input of the post-ANTS c3d_affine_tool
            
            #create_log_node(func_to_anat, 'outputspec.mni_func', num_strat)
            num_strat += 1

    strat_list += new_strat_list
    """


    """
    Inserting Generate Motion Statistics Workflow
    """
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

            except:
                print 'Invalid Connection: Generate Motion Statistics:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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


    """
    Inserting Nuisance Workflow
    """
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runNuisance:
        workflow_bit_id['nuisance'] = workflow_counter
        for strat in strat_list:

            if 'FSL' in c.regOption:
                nuisance = create_nuisance(False, 'nuisance_%d' % num_strat)
            elif 'ANTS' in c.regOption:
                nuisance = create_nuisance(True, 'nuisance_%d' % num_strat)


            nuisance.get_node('residuals').iterables = ([('selector', c.Corrections),
                                                         ('compcor_ncomponents', c.nComponents)])

            nuisance.inputs.inputspec.harvard_oxford_mask = c.harvardOxfordMask

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

                if 'FSL' in c.regOption:
                    node, out_file = strat.get_node_from_resource_pool('mni_to_anatomical_linear_xfm')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.mni_to_anat_linear_xfm')
                elif 'ANTS' in c.regOption:
                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     nuisance, 'inputspec.mni_to_anat_linear_xfm')


            except:
                print 'Invalid Connection: Nuisance:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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

    """
    Inserting Median Angle Correction Workflow
    """
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
                print 'Invalid Connection: Median Angle Correction:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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

    """
    Inserting ALFF/fALFF
    Workflow
    """
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
                print 'Invalid Connection: ALFF:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise
            strat.append_name(alff.name)
            strat.update_resource_pool({'alff_img':(alff, 'outputspec.alff_img')})
            strat.update_resource_pool({'falff_img':(alff, 'outputspec.falff_img')})
            strat.update_resource_pool({'alff_Z_img':(alff, 'outputspec.alff_Z_img')})
            strat.update_resource_pool({'falff_Z_img':(alff, 'outputspec.falff_Z_img')})
            
            create_log_node(alff, 'outputspec.falff_img', num_strat)

            num_strat += 1
    strat_list += new_strat_list




    """
    Inserting Frequency Filtering Node
    """
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
                print 'Invalid Connection: Frequency Filtering:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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




    """
    Inserting Scrubbing Workflow
    """
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1

    ###
    if c.runScrubbing == 10:
        c.runScrubbing = [1, 0]

    if 1 in c.runScrubbing:
        workflow_bit_id['scrubbing'] = workflow_counter
        for strat in strat_list:

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
                print 'Invalid Connection: Scrubbing Workflow:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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





    """
    Func -> Template, uses ApplyWarp (FSL) or Merge + WarpImageMultiTransform (ANTS) to apply
    also includes mean functional warp
    """
    new_strat_list = []
    num_strat = 0
    if 1 in c.runRegisterFuncToMNI:

        for strat in strat_list:

            # Run FSL ApplyWarp
            if 'FSL' in c.regOption:

                func_mni_warp = pe.Node(interface=fsl.ApplyWarp(),
                                        name='func_mni_fsl_warp_%d' % num_strat)
                func_mni_warp.inputs.ref_file = c.standardResolutionBrain
    
                functional_brain_mask_to_standard = pe.Node(interface=fsl.ApplyWarp(),
                                                            name='func_mni_fsl_warp_mask_%d' % num_strat)
                functional_brain_mask_to_standard.inputs.interp = 'nn'
                functional_brain_mask_to_standard.inputs.ref_file = c.standard

                mean_functional_warp = pe.Node(interface=fsl.ApplyWarp(), name='mean_func_fsl_warp_%d' % num_strat)
                mean_functional_warp.inputs.ref_file = c.standardResolutionBrain
    
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

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     functional_brain_mask_to_standard, 'premat') 
                    workflow.connect(node, out_file,
                                     mean_functional_warp, 'premat') 

                    node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                    workflow.connect(node, out_file,
                                     functional_brain_mask_to_standard, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool('mean_functional')
                    workflow.connect(node, out_file, mean_functional_warp, 'in_file')

                    

                except:
                    print 'Invalid Connection: Register Functional timeseries to MNI space:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise
    
                strat.update_resource_pool({'functional_mni':(func_mni_warp, 'out_file'),
                                            'functional_brain_mask_to_standard':(functional_brain_mask_to_standard, 'out_file'),
                                            'mean_functional_in_mni':(mean_functional_warp, 'out_file')})
                strat.append_name(func_mni_warp.name)
                create_log_node(func_mni_warp, 'out_file', num_strat)
            
                num_strat += 1


            # Run ANTS apply (WarpImageMultiTransform) instead
            elif 'ANTS' in c.regOption:

                # THIS SHOULD ALL BE COMBINED INTO ITS OWN WORKFLOW IN REGISTRATION.PY

                # converts FSL-format .mat affine xfm into ANTS-format .txt
                # .mat affine comes from Func->Anat registration
                fsl_reg_2_itk = pe.Node(c3.C3dAffineTool(), name='fsl_reg_2_itk_%d' % num_strat)
                fsl_reg_2_itk.inputs.itk_transform = True
                fsl_reg_2_itk.inputs.fsl2ras = True

                # converts FSL-format .mat affine xfm into ANTS-format .txt
                # .mat affine comes from Func->Anat registration
                fsl_reg_2_itk_mask = pe.Node(c3.C3dAffineTool(), name='fsl_reg_2_itk_mask_%d' % num_strat)
                fsl_reg_2_itk_mask.inputs.itk_transform = True
                fsl_reg_2_itk_mask.inputs.fsl2ras = True

                # converts FSL-format .mat affine xfm into ANTS-format .txt
                # .mat affine comes from Func->Anat registration
                fsl_reg_2_itk_mean = pe.Node(c3.C3dAffineTool(), name='fsl_reg_2_itk_mean_%d' % num_strat)
                fsl_reg_2_itk_mean.inputs.itk_transform = True
                fsl_reg_2_itk_mean.inputs.fsl2ras = True



                #collects series of transformations to be applied to the moving images
                collect_transforms = pe.Node(util.Merge(3), name='collect_transforms_%d' % num_strat)

                #performs series of transformations on moving images
                warp_images = pe.Node(ants.WarpTimeSeriesImageMultiTransform(), name='func_mni_ants_warp_images_%d' % num_strat)
                warp_images.inputs.reference_image = c.standardResolutionBrain
                warp_images.inputs.dimension = 4
 

                #collects series of transformations to be applied to the moving images
                collect_transforms_mask = pe.Node(util.Merge(3), name='collect_transforms_mask_%d' % num_strat)

                #performs series of transformations on moving images
                warp_images_mask = pe.Node(ants.WarpImageMultiTransform(), name='func_mni_ants_warp_images_mask_%d' % num_strat)
                warp_images_mask.inputs.reference_image = c.standard
                warp_images_mask.inputs.use_nearest = True
                warp_images_mask.inputs.dimension = 3
                
                
                #collects series of transformations to be applied to the moving images
                collect_transforms_mean = pe.Node(util.Merge(3), name='collect_transforms_mean_%d' % num_strat)

                # mean functional warp
                mean_functional_ants_warp = pe.Node(ants.WarpImageMultiTransform(), name='mean_func_ants_warp_%d' % num_strat)
                mean_functional_ants_warp.inputs.reference_image = c.standardResolutionBrain
                mean_functional_ants_warp.inputs.dimension = 3

    
                try:

                    # FUNCTIONAL

                    # convert the .mat from linear Func->Anat to ANTS format
                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file, fsl_reg_2_itk, 'transform_file')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file, fsl_reg_2_itk, 'reference_file')

                    node, out_file = strat.get_node_from_resource_pool('mean_functional')
                    workflow.connect(node, out_file, fsl_reg_2_itk, 'source_file')


                    # Premat from Func->Anat linear reg and bbreg (if bbreg is enabled)
                    workflow.connect(fsl_reg_2_itk, 'itk_transform', collect_transforms, 'in3')
    
                    # Field file from anatomical nonlinear registration
                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file, collect_transforms, 'in1')

                    # affine transformation from anatomical registration
                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file, collect_transforms, 'in2')

                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file, warp_images, 'input_image')

                    workflow.connect(collect_transforms, 'out', warp_images, 'transformation_series')



                    # FUNCTIONAL MASK

                    # convert the .mat from linear Func->Anat to ANTS format
                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file, fsl_reg_2_itk_mask, 'transform_file')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file, fsl_reg_2_itk_mask, 'reference_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                    workflow.connect(node, out_file, fsl_reg_2_itk_mask, 'source_file')


                    # Premat from Func->Anat linear reg and bbreg (if bbreg is enabled)
                    workflow.connect(fsl_reg_2_itk_mask, 'itk_transform', collect_transforms_mask, 'in3')
    
                    # Field file from anatomical nonlinear registration
                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file, collect_transforms_mask, 'in1')

                    # affine transformation from anatomical registration
                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file, collect_transforms_mask, 'in2')

                    node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                    workflow.connect(node, out_file, warp_images_mask, 'input_image')

                    workflow.connect(collect_transforms_mask, 'out', warp_images_mask, 'transformation_series')
                    
                    
                    
                    # MEAN FUNCTIONAL

                    # convert the .mat from linear Func->Anat to ANTS format
                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file, fsl_reg_2_itk_mean, 'transform_file')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file, fsl_reg_2_itk_mean, 'reference_file')

                    node, out_file = strat.get_node_from_resource_pool('mean_functional')
                    workflow.connect(node, out_file, fsl_reg_2_itk_mean, 'source_file')


                    # Premat from Func->Anat linear reg and bbreg (if bbreg is enabled)
                    workflow.connect(fsl_reg_2_itk_mean, 'itk_transform', collect_transforms_mean, 'in3')
    
                    # Field file from anatomical nonlinear registration
                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file, collect_transforms_mean, 'in1')

                    # affine transformation from anatomical registration
                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file, collect_transforms_mean, 'in2')

                    node, out_file = strat.get_node_from_resource_pool('mean_functional')
                    workflow.connect(node, out_file, mean_functional_ants_warp, 'input_image')

                    workflow.connect(collect_transforms_mean, 'out', mean_functional_ants_warp, 'transformation_series')



    

                except:
                    print 'Invalid Connection: Register Functional timeseries to MNI space:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise
    
                strat.update_resource_pool({'functional_mni':(warp_images, 'output_image'),
                                            'functional_brain_mask_to_standard':(warp_images_mask, 'output_image'),
                                            'mean_functional_in_mni':(mean_functional_ants_warp, 'output_image')})

                strat.append_name(warp_images.name)
                create_log_node(warp_images, 'output_image', num_strat)
            
                num_strat += 1



    strat_list += new_strat_list
    



    """
    Inserting VMHC
    Workflow
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runVMHC:
        for strat in strat_list:

            if 'FSL' in c.regOption:
                preproc = create_vmhc(False)
            elif 'ANTS' in c.regOption:
                preproc = create_vmhc(True)

            preproc.inputs.inputspec.brain_symmetric = \
                                            c.brainSymmetric
            preproc.inputs.inputspec.symm_standard = \
                                            c.symmStandard
            preproc.inputs.inputspec.twomm_brain_mask_dil = \
                                            c.twommBrainMaskDiluted
            preproc.inputs.inputspec.config_file_twomm = \
                                            c.configFileTwomm
            preproc.inputs.inputspec.standard = \
                                            c.standard
            preproc.inputs.fwhm_input.fwhm = c.fwhm
            preproc.get_node('fwhm_input').iterables = ('fwhm',
                                                        c.fwhm)


            vmhc = preproc.clone('vmhc_%d' % num_strat)



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

                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.brain')

                node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.reorient')

                node, out_file = strat.get_node_from_resource_pool('mean_functional')
                workflow.connect(node, out_file,
                                 vmhc, 'inputspec.mean_functional')

            except:
                print 'Invalid Connection: VMHC:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise


            strat.update_resource_pool({'vmhc_raw_score':(vmhc, 'outputspec.VMHC_FWHM_img')})
            strat.update_resource_pool({'vmhc_z_score':(vmhc, 'outputspec.VMHC_Z_FWHM_img')})
            strat.update_resource_pool({'vmhc_z_score_stat_map':(vmhc, 'outputspec.VMHC_Z_stat_FWHM_img')})
            strat.append_name(vmhc.name)
            
            create_log_node(vmhc, 'outputspec.VMHC_FWHM_img', num_strat)
            
            num_strat += 1

    strat_list += new_strat_list




    """
    Inserting REHO
    Workflow
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runReHo:
        for strat in strat_list:

            preproc = create_reho()
            preproc.inputs.inputspec.cluster_size = c.clusterSize
            reho = preproc.clone('reho_%d' % num_strat)


            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 reho, 'inputspec.rest_res_filt')

                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 reho, 'inputspec.rest_mask')
            except:
                print 'Invalid Connection: REHO:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise


            strat.update_resource_pool({'raw_reho_map':(reho, 'outputspec.raw_reho_map')})
            strat.update_resource_pool({'reho_Z_img':(reho, 'outputspec.z_score')})
            strat.append_name(reho.name)
            
            create_log_node(reho, 'outputspec.raw_reho_map', num_strat)
            
            num_strat += 1
    strat_list += new_strat_list



    """
    Transforming ALFF Z scores and fAlff Z scores to MNI
    """
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runALFF):
        for strat in strat_list:

            if 'FSL' in c.regOption:

                alff_Z_to_standard = pe.Node(interface=fsl.ApplyWarp(),
                               name='alff_Z_to_standard_%d' % num_strat)

                alff_Z_to_standard.inputs.ref_file = c.standard

                falff_Z_to_standard = alff_Z_to_standard.clone('falff_Z_to_standard_%d' % num_strat)
                falff_Z_to_standard.inputs.ref_file = c.standard

                try:

                    node, out_file = strat.get_node_from_resource_pool('alff_Z_img')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard, 'premat')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard, 'field_file')


                    node, out_file = strat.get_node_from_resource_pool('falff_Z_img')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard, 'premat')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard, 'field_file')



                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'alff_Z_to_standard':(alff_Z_to_standard, 'out_file')})
                strat.update_resource_pool({'falff_Z_to_standard':(falff_Z_to_standard, 'out_file')})
                strat.append_name(alff_Z_to_standard.name)
            
                num_strat += 1


            elif 'ANTS' in c.regOption:

                alff_Z_to_standard = create_apply_ants_xfm(3, 0, name='alff_Z_to_standard_%d' % num_strat)

                alff_Z_to_standard.inputs.inputspec.warp_reference = c.standard


                falff_Z_to_standard = alff_Z_to_standard.clone('falff_Z_to_standard_%d' % num_strat)
               
                falff_Z_to_standard.inputs.inputspec.warp_reference = c.standard


                try:

                    node, out_file = strat.get_node_from_resource_pool('alff_Z_img')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard, 'inputspec.in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard, 'inputspec.func_anat_affine')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard, 'inputspec.conversion_reference')

                    node, out_file = strat.get_node_from_resource_pool('alff_Z_img')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard, 'inputspec.conversion_source')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard, 'inputspec.nonlinear_field')

                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard, 'inputspec.ants_affine')


                    node, out_file = strat.get_node_from_resource_pool('falff_Z_img')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard, 'inputspec.in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard, 'inputspec.func_anat_affine')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard, 'inputspec.conversion_reference')

                    node, out_file = strat.get_node_from_resource_pool('falff_Z_img')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard, 'inputspec.conversion_source')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard, 'inputspec.nonlinear_field')

                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard, 'inputspec.ants_affine')



                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'alff_Z_to_standard':(alff_Z_to_standard, 'outputspec.out_file')})
                strat.update_resource_pool({'falff_Z_to_standard':(falff_Z_to_standard, 'outputspec.out_file')})
                strat.append_name(alff_Z_to_standard.name)
            
                num_strat += 1
    
    strat_list += new_strat_list


    inputnode_fwhm = None
    if c.fwhm != None:

        inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')
        inputnode_fwhm.iterables = ("fwhm", c.fwhm)



    """    
    Smoothing ALFF fALFF Z scores and or possibly Z scores in MNI 
    """
    new_strat_list = []
    num_strat = 0
    if (1 in c.runALFF) and c.fwhm != None:
        for strat in strat_list:


            alff_Z_to_standard_smooth = None
            falff_Z_to_standard_smooth = None

            alff_Z_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                        name='alff_Z_smooth_%d' % num_strat)

            falff_Z_smooth = alff_Z_smooth.clone('falff_Z_smooth_%d' % num_strat)


            try:
                node, out_file = strat.get_node_from_resource_pool('alff_Z_img')
                workflow.connect(node, out_file,
                                 alff_Z_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                alff_Z_smooth, 'op_string')
                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 alff_Z_smooth, 'operand_files')

                node, out_file = strat.get_node_from_resource_pool('falff_Z_img')
                workflow.connect(node, out_file,
                                 falff_Z_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                falff_Z_smooth, 'op_string')
                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 falff_Z_smooth, 'operand_files')



            except:
                print 'Invalid Connection: ALFF smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise
            strat.append_name(alff_Z_smooth.name)
            strat.update_resource_pool({'alff_Z_smooth':(alff_Z_smooth, 'out_file')})
            strat.update_resource_pool({'falff_Z_smooth':(falff_Z_smooth, 'out_file')})
            
  
            if c.runRegisterFuncToMNI:

                alff_Z_to_standard_smooth = alff_Z_smooth.clone('alff_Z_to_standard_smooth_%d' % num_strat)
                falff_Z_to_standard_smooth = alff_Z_smooth.clone('falff_Z_to_standard_smooth_%d' % num_strat)

                try:


                    node, out_file = strat.get_node_from_resource_pool('alff_Z_to_standard')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard_smooth, 'in_file')
                    workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                    alff_Z_to_standard_smooth, 'op_string')

                    node, out_file = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')
                    workflow.connect(node, out_file,
                                     alff_Z_to_standard_smooth, 'operand_files')

                    node, out_file = strat.get_node_from_resource_pool('falff_Z_to_standard')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard_smooth, 'in_file')
                    workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                    falff_Z_to_standard_smooth, 'op_string')
                    node, out_file = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')
                    workflow.connect(node, out_file,
                                     falff_Z_to_standard_smooth, 'operand_files')




                except:

                    print 'Invalid Connection: ALFF smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.append_name(alff_Z_to_standard_smooth.name)
                strat.append_name(falff_Z_to_standard_smooth.name)
                strat.update_resource_pool({'alff_Z_to_standard_smooth':(alff_Z_to_standard_smooth, 'out_file')})
                strat.update_resource_pool({'falff_Z_to_standard_smooth':(falff_Z_to_standard_smooth, 'out_file')})

                create_log_node(alff_Z_to_standard_smooth, 'out_file', num_strat)
                create_log_node(falff_Z_to_standard_smooth, 'out_file', num_strat)

            num_strat += 1
    strat_list += new_strat_list



    """
    Transforming ReHo Z scores to MNI
    """
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runReHo):
        for strat in strat_list:

            if 'FSL' in c.regOption:

                reho_Z_to_standard = pe.Node(interface=fsl.ApplyWarp(),
                               name='reho_Z_to_standard_%d' % num_strat)

                reho_Z_to_standard.inputs.ref_file = c.standard


                try:

                    node, out_file = strat.get_node_from_resource_pool('reho_Z_img')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard, 'premat')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard, 'field_file')

                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'reho_Z_to_standard':(reho_Z_to_standard, 'out_file')})
                strat.append_name(reho_Z_to_standard.name)
                num_strat += 1


            elif 'ANTS' in c.regOption:
                
                reho_Z_to_standard = create_apply_ants_xfm(3, 0, name='reho_Z_to_standard_%d' % num_strat)

                reho_Z_to_standard.inputs.inputspec.warp_reference = c.standard
                '''

                # converts FSL-format .mat affine xfm into ANTS-format .txt
                # .mat affine comes from Func->Anat registration
                fsl_reg_2_itk_reho = pe.Node(c3.C3dAffineTool(), name='fsl_reg_2_itk_reho_%d' % num_strat)
                fsl_reg_2_itk_reho.inputs.itk_transform = True
                fsl_reg_2_itk_reho.inputs.fsl2ras = True



                #collects series of transformations to be applied to the moving images
                collect_transforms_reho = pe.Node(util.Merge(3), name='collect_transforms_reho_%d' % num_strat)

                #performs series of transformations on moving images
                warp_images_reho = pe.Node(ants.WarpImageMultiTransform(), name='reho_ants_warp_images_%d' % num_strat)
                warp_images_reho.inputs.reference_image = c.standard
                warp_images_reho.inputs.dimension = 3
                '''
                try:
                    '''
                    # MEAN FUNCTIONAL

                    # convert the .mat from linear Func->Anat to ANTS format
                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file, fsl_reg_2_itk_reho, 'transform_file')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file, fsl_reg_2_itk_reho, 'reference_file')

                    node, out_file = strat.get_node_from_resource_pool('reho_Z_img')
                    workflow.connect(node, out_file, fsl_reg_2_itk_reho, 'source_file')


                    # Premat from Func->Anat linear reg and bbreg (if bbreg is enabled)
                    workflow.connect(fsl_reg_2_itk_reho, 'itk_transform', collect_transforms_reho, 'in3')
    
                    # Field file from anatomical nonlinear registration
                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file, collect_transforms_reho, 'in1')

                    # affine transformation from anatomical registration
                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file, collect_transforms_reho, 'in2')

                    node, out_file = strat.get_node_from_resource_pool('reho_Z_img')
                    workflow.connect(node, out_file, warp_images_reho, 'input_image')

                    workflow.connect(collect_transforms_reho, 'out', warp_images_reho, 'transformation_series')

                    '''
                    node, out_file = strat.get_node_from_resource_pool('reho_Z_img')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard, 'inputspec.in_file')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard, 'inputspec.nonlinear_field')

                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard, 'inputspec.ants_affine')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard, 'inputspec.func_anat_affine')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard, 'inputspec.conversion_reference')

                    node, out_file = strat.get_node_from_resource_pool('reho_Z_img')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard, 'inputspec.conversion_source')
                    

                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'reho_Z_to_standard':(reho_Z_to_standard, 'outputspec.out_file')})#warp_images_reho, 'output_image')})#reho_Z_to_standard, 'outputspec.out_file')})
                strat.append_name(reho_Z_to_standard.name)#warp_images_reho.name)#reho_Z_to_standard.name)
                num_strat += 1


    strat_list += new_strat_list


    """
    Smoothing ReHo Z scores and or possibly Z scores in MNI 
    """
    new_strat_list = []
    num_strat = 0

    if (1 in c.runReHo) and c.fwhm != None:
        for strat in strat_list:


            reho_Z_to_standard_smooth = None

            reho_Z_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                        name='reho_Z_smooth_%d' % num_strat)

            try:
                node, out_file = strat.get_node_from_resource_pool('reho_Z_img')
                workflow.connect(node, out_file,
                                 reho_Z_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                reho_Z_smooth, 'op_string')
                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 reho_Z_smooth, 'operand_files')

            except:
                print 'Invalid Connection: reho_Z smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise
            strat.append_name(reho_Z_smooth.name)
            strat.update_resource_pool({'reho_Z_smooth':(reho_Z_smooth, 'out_file')})

            if 1 in c.runRegisterFuncToMNI:

                reho_Z_to_standard_smooth = reho_Z_smooth.clone('reho_Z_to_standard_smooth_%d' % num_strat)

                try:

                    node, out_file = strat.get_node_from_resource_pool('reho_Z_to_standard')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard_smooth, 'in_file')
                    workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                    reho_Z_to_standard_smooth, 'op_string')
                    node, out_file = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')
                    workflow.connect(node, out_file,
                                     reho_Z_to_standard_smooth, 'operand_files')



                except:

                    print 'Invalid Connection: reho_Z_to_standard smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.append_name(reho_Z_to_standard_smooth.name)
                strat.update_resource_pool({'reho_Z_to_standard_smooth':(reho_Z_to_standard_smooth, 'out_file')})
                create_log_node(reho_Z_to_standard_smooth, 'out_file', num_strat)
            
            num_strat += 1
    strat_list += new_strat_list


    """
    Spatial Regression Based Time Series
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runSpatialRegression:

        for strat in strat_list:

            resample_functional_to_spatial_map = pe.Node(interface=fsl.FLIRT(),
                                                         name='resample_functional_to_spatial_map_%d' % num_strat)
            resample_functional_to_spatial_map.inputs.interp = 'trilinear'
            resample_functional_to_spatial_map.inputs.apply_xfm = True
            resample_functional_to_spatial_map.inputs.in_matrix_file = c.identityMatrix
            
            resample_functional_mask_to_spatial_map = pe.Node(interface=fsl.FLIRT(),
                                                         name='resample_functional_mask_to_spatial_map_%d' % num_strat)
            resample_functional_mask_to_spatial_map.inputs.interp = 'nearestneighbour'
            resample_functional_mask_to_spatial_map.inputs.apply_xfm = True
            resample_functional_mask_to_spatial_map.inputs.in_matrix_file = c.identityMatrix

            spatial_map_dataflow = create_spatial_map_dataflow(c.spatialPatternMaps, 'spatial_map_dataflow_%d' % num_strat)

            spatial_map_timeseries = get_spatial_map_timeseries('spatial_map_timeseries_%d' % num_strat)
            spatial_map_timeseries.inputs.inputspec.demean = c.spatialDemean

            try:

                node, out_file = strat.get_node_from_resource_pool('functional_mni')
                node2, out_file2 = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')

                # resample the input functional file and functional mask to spatial map
                workflow.connect(node, out_file,
                                 resample_functional_to_spatial_map, 'in_file')
                workflow.connect(spatial_map_dataflow, 'select_spatial_map.out_file',
                                 resample_functional_to_spatial_map, 'reference')
                workflow.connect(node2, out_file2,
                                 resample_functional_mask_to_spatial_map, 'in_file')
                workflow.connect(spatial_map_dataflow, 'select_spatial_map.out_file',
                                 resample_functional_mask_to_spatial_map, 'reference')

                # connect it to the spatial_map_timeseries
                workflow.connect(spatial_map_dataflow, 'select_spatial_map.out_file',
                                 spatial_map_timeseries, 'inputspec.spatial_map')
                workflow.connect(resample_functional_mask_to_spatial_map, 'out_file',
                                 spatial_map_timeseries, 'inputspec.subject_mask')
                workflow.connect(resample_functional_to_spatial_map, 'out_file',
                                 spatial_map_timeseries, 'inputspec.subject_rest')
                               
                

            except:
                print 'Invalid Connection: Spatial map timeSeries extraction workflow:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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

            strat.update_resource_pool({'spatial_map_timeseries' : (spatial_map_timeseries, 'outputspec.subject_timeseries'),
                                       'functional_to_spatial_map' : (resample_functional_to_spatial_map, 'out_file'),
                                       'functional_mask_to_spatial_map' : (resample_functional_mask_to_spatial_map, 'out_file')})
            
            create_log_node(spatial_map_timeseries, 'outputspec.subject_timeseries', num_strat)

            num_strat += 1

    strat_list += new_strat_list


    """
    Voxel Based Time Series 
    """
    new_strat_list = []
    num_strat = 0
    if 1 in c.runVoxelTimeseries:


        for strat in strat_list:

            resample_functional_to_mask = pe.Node(interface=fsl.FLIRT(),
                                                  name='resample_functional_to_mask_%d' % num_strat)
            resample_functional_to_mask.inputs.interp = 'trilinear'
            resample_functional_to_mask.inputs.apply_xfm = True
            resample_functional_to_mask.inputs.in_matrix_file = c.identityMatrix

            mask_dataflow = create_mask_dataflow(c.maskSpecificationFile, 'mask_dataflow_%d' % num_strat)

            voxel_timeseries = get_voxel_timeseries('voxel_timeseries_%d' % num_strat)
            voxel_timeseries.inputs.inputspec.output_type = c.voxelTSOutputs

            try:

                node, out_file = strat.get_node_from_resource_pool('functional_mni')

                # resample the input functional file to mask
                workflow.connect(node, out_file,
                                 resample_functional_to_mask, 'in_file')
                workflow.connect(mask_dataflow, 'select_mask.out_file',
                                 resample_functional_to_mask, 'reference')

                # connect it to the voxel_timeseries
                workflow.connect(mask_dataflow, 'select_mask.out_file',
                                 voxel_timeseries, 'input_mask.mask')
                workflow.connect(resample_functional_to_mask, 'out_file',
                                 voxel_timeseries, 'inputspec.rest')

            except:
                print 'Invalid Connection: Voxel TimeSeries Analysis Workflow:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise

            if 0 in c.runVoxelTimeseries:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(voxel_timeseries.name)

            strat.update_resource_pool({'voxel_timeseries': (voxel_timeseries, 'outputspec.mask_outputs')})
            
            create_log_node(voxel_timeseries, 'outputspec.mask_outputs', num_strat)

            num_strat += 1

    strat_list += new_strat_list




    """
    ROI Based Time Series
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runROITimeseries:

        for strat in strat_list:

            resample_functional_to_roi = pe.Node(interface=fsl.FLIRT(),
                                                  name='resample_functional_to_roi_%d' % num_strat)
            resample_functional_to_roi.inputs.interp = 'trilinear'
            resample_functional_to_roi.inputs.apply_xfm = True
            resample_functional_to_roi.inputs.in_matrix_file = c.identityMatrix

            roi_dataflow = create_roi_dataflow(c.roiSpecificationFile, 'roi_dataflow_%d' % num_strat)

            roi_timeseries = get_roi_timeseries('roi_timeseries_%d' % num_strat)
            roi_timeseries.inputs.inputspec.output_type = c.roiTSOutputs

            try:

                node, out_file = strat.get_node_from_resource_pool('functional_mni')

                # resample the input functional file to roi
                workflow.connect(node, out_file,
                                 resample_functional_to_roi, 'in_file')
                workflow.connect(roi_dataflow, 'select_roi.out_file',
                                 resample_functional_to_roi, 'reference')

                # connect it to the roi_timeseries
                workflow.connect(roi_dataflow, 'select_roi.out_file',
                                 roi_timeseries, 'input_roi.roi')
                workflow.connect(resample_functional_to_roi, 'out_file',
                                 roi_timeseries, 'inputspec.rest')

            except:
                print 'Invalid Connection: ROI TimeSeries Analysis Workflow:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise

            if 0 in c.runROITimeseries:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.leaf_out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name(roi_timeseries.name)

            strat.update_resource_pool({'roi_timeseries' : (roi_timeseries, 'outputspec.roi_outputs')})

            create_log_node(roi_timeseries, 'outputspec.roi_outputs', num_strat)
            num_strat += 1

    strat_list += new_strat_list


    """
    Inserting SCA
    Workflow for ROI INPUT
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runSCA and (1 in c.runROITimeseries):
        for strat in strat_list:

            sca_roi = create_sca('sca_roi_%d' % num_strat)


            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 sca_roi, 'inputspec.functional_file')

                node, out_file = strat.get_node_from_resource_pool('roi_timeseries')
                workflow.connect(node, (out_file, extract_one_d),
                                 sca_roi, 'inputspec.timeseries_one_d')
            except:
                print 'Invalid Connection: SCA ROI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise


            strat.update_resource_pool({'sca_roi_correlations':(sca_roi, 'outputspec.correlation_file')})
            strat.update_resource_pool({'sca_roi_Z':(sca_roi, 'outputspec.Z_score')})
            
            create_log_node(sca_roi, 'outputspec.correlation_file', num_strat)
            
            strat.append_name(sca_roi.name)
            num_strat += 1
    strat_list += new_strat_list


    """
    Temporal Regression for Dual Regression
    """
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
                print 'Invalid Connection: Temporal multiple regression for dual regression:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise


            strat.update_resource_pool({'dr_tempreg_maps_stack':(dr_temp_reg, 'outputspec.temp_reg_map')})
            strat.update_resource_pool({'dr_tempreg_maps_z_stack':(dr_temp_reg, 'outputspec.temp_reg_map_z'),
                                        'dr_tempreg_maps_z_files':(dr_temp_reg, 'outputspec.temp_reg_map_z_stack')})
            strat.append_name(dr_temp_reg.name)
            
            create_log_node(dr_temp_reg, 'outputspec.temp_reg_map', num_strat)
            
            num_strat += 1
    strat_list += new_strat_list

    """
    Transforming Dual Regression Z stats to MNI
    """
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runDualReg) and (1 in c.runSpatialRegression):
        for strat in strat_list:

            if 'FSL' in c.regOption:

                dr_tempreg_stack_Z_to_standard = pe.Node(interface=fsl.ApplyWarp(),
                               name='dr_tempreg_stack_Z_to_standard_%d' % num_strat)
                dr_tempreg_stack_Z_to_standard.inputs.ref_file = c.standard
            
                dr_tempreg_files_Z_to_standard = pe.MapNode(interface=fsl.ApplyWarp(),
                                                            name = 'dr_tempreg_files_Z_to_standard_%d' % num_strat,
                                                            iterfield=['in_file'])
                dr_tempreg_files_Z_to_standard.inputs.ref_file = c.standard
            
                dr_tempreg_stack_to_standard = dr_tempreg_stack_Z_to_standard.clone(name= 'dr_tempreg_stack_to_standard_%d' % num_strat)


                try:
                
                    ##dual tempreg z stack
                    node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_z_stack')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_Z_to_standard, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_Z_to_standard, 'premat')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_Z_to_standard, 'field_file')
                
                    ####dual tempreg z files
                    node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_z_files')
                    workflow.connect(node, out_file,
                                     dr_tempreg_files_Z_to_standard, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_files_Z_to_standard, 'premat')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_files_Z_to_standard, 'field_file')
                
                    ##dual tempreg stack
                    node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_stack')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_to_standard, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_to_standard, 'premat')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_to_standard, 'field_file')                



                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'dr_tempreg_maps_z_stack_to_standard':(dr_tempreg_stack_Z_to_standard, 'out_file')})
                strat.update_resource_pool({'dr_tempreg_maps_z_files_to_standard':(dr_tempreg_files_Z_to_standard, 'out_file')})
                strat.update_resource_pool({'dr_tempreg_maps_stack_to_standard':(dr_tempreg_stack_to_standard, 'out_file')})
                strat.append_name(dr_tempreg_stack_to_standard.name)
            
                num_strat += 1


            elif 'ANTS' in c.regOption:

                dr_tempreg_stack_Z_to_standard = create_apply_ants_xfm(3, 0, 
                               name='dr_tempreg_stack_Z_to_standard_%d' % num_strat)
                dr_tempreg_stack_Z_to_standard.inputs.inputspec.warp_reference = c.standard
      
                dr_tempreg_files_Z_to_standard = create_apply_ants_xfm(3, 1,
                                                            name = 'dr_tempreg_files_Z_to_standard_%d' % num_strat)
                dr_tempreg_files_Z_to_standard.inputs.inputspec.warp_reference = c.standard
            
                dr_tempreg_stack_to_standard = dr_tempreg_stack_Z_to_standard.clone(name= 'dr_tempreg_stack_to_standard_%d' % num_strat)


                try:

                    ##dual tempreg z stack
                    node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_z_stack')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_Z_to_standard, 'inputspec.in_file')

                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_Z_to_standard, 'inputspec.ants_affine')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_Z_to_standard, 'inputspec.func_anat_affine')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_Z_to_standard, 'inputspec.nonlinear_field')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_Z_to_standard, 'inputspec.conversion_reference')

                    node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_z_stack')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_Z_to_standard, 'inputspec.conversion_source')
                
                    ####dual tempreg z files
                    node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_z_files')
                    workflow.connect(node, out_file,
                                     dr_tempreg_files_Z_to_standard, 'inputspec.in_file')

                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_files_Z_to_standard, 'inputspec.ants_affine')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_files_Z_to_standard, 'inputspec.func_anat_affine')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_files_Z_to_standard, 'inputspec.nonlinear_field')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     dr_tempreg_files_Z_to_standard, 'inputspec.conversion_reference')

                    node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_z_files')
                    workflow.connect(node, out_file,
                                     dr_tempreg_files_Z_to_standard, 'inputspec.conversion_source')
                
                    ##dual tempreg stack
                    node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_stack')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_to_standard, 'inputspec.in_file')

                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_to_standard, 'inputspec.ants_affine')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_to_standard, 'inputspec.func_anat_affine')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_to_standard, 'inputspec.nonlinear_field') 

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_to_standard, 'inputspec.conversion_reference')

                    node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_stack')
                    workflow.connect(node, out_file,
                                     dr_tempreg_stack_to_standard, 'inputspec.conversion_source')


                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'dr_tempreg_maps_z_stack_to_standard':(dr_tempreg_stack_Z_to_standard, 'outputspec.out_file')})
                strat.update_resource_pool({'dr_tempreg_maps_z_files_to_standard':(dr_tempreg_files_Z_to_standard, 'outputspec.out_file')})
                strat.update_resource_pool({'dr_tempreg_maps_stack_to_standard':(dr_tempreg_stack_to_standard, 'outputspec.out_file')})
                strat.append_name(dr_tempreg_stack_to_standard.name)
            
                num_strat += 1



    strat_list += new_strat_list
    


    """
    Inserting SCA
    Workflow for Voxel INPUT
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runSCA and (1 in c.runVoxelTimeseries):
        for strat in strat_list:

            sca_seed = create_sca('sca_seed_%d' % num_strat)


            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 sca_seed, 'inputspec.functional_file')

                node, out_file = strat.get_node_from_resource_pool('voxel_timeseries')
                workflow.connect(node, (out_file, extract_one_d),
                                 sca_seed, 'inputspec.timeseries_one_d')
            except:
                print 'Invalid Connection: SCA:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise


            strat.update_resource_pool({'sca_seed_correlations':(sca_seed, 'outputspec.correlation_file')})
            strat.update_resource_pool({'sca_seed_Z':(sca_seed, 'outputspec.Z_score')})
            strat.append_name(sca_seed.name)
            num_strat += 1
    strat_list += new_strat_list


    """
    Temporal Regression for SCA
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runMultRegSCA and (1 in c.runROITimeseries):
        for strat in strat_list:

            sc_temp_reg = create_temporal_reg('temporal_regression_sca_%d' % num_strat, which='RT')
            sc_temp_reg.inputs.inputspec.normalize = c.mrsNorm
            sc_temp_reg.inputs.inputspec.demean = c.mrsDemean

            try:
                node, out_file = strat.get_node_from_resource_pool('functional_mni')
                node2, out_file2 = strat.get_node_from_resource_pool('roi_timeseries')
                node3, out_file3 = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')

                workflow.connect(node, out_file,
                                 sc_temp_reg, 'inputspec.subject_rest')

                workflow.connect(node2, (out_file2, extract_txt),
                                 sc_temp_reg, 'inputspec.subject_timeseries')

                workflow.connect(node3, out_file3,
                                 sc_temp_reg, 'inputspec.subject_mask')

            except:
                print 'Invalid Connection: Temporal multiple regression for seed based connectivity:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise


            strat.update_resource_pool({'sca_tempreg_maps_stack':(sc_temp_reg, 'outputspec.temp_reg_map')})
            strat.update_resource_pool({'sca_tempreg_maps_z_stack':(sc_temp_reg, 'outputspec.temp_reg_map_z'),
                                        'sca_tempreg_maps_z_files':(sc_temp_reg, 'outputspec.temp_reg_map_z_stack')})
            
            create_log_node(sc_temp_reg, 'outputspec.temp_reg_map', num_strat)
            
            strat.append_name(sc_temp_reg.name)
            num_strat += 1
    strat_list += new_strat_list


    """
    Smoothing Temporal Regression for SCA scores
    """
    new_strat_list = []
    num_strat = 0

    if (1 in c.runMultRegSCA) and (1 in c.runROITimeseries) and c.fwhm != None:
        for strat in strat_list:

            sc_temp_reg_maps_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='sca_tempreg_maps_stack_smooth_%d' % num_strat, iterfield=['in_file'])
            sc_temp_reg_maps_Z_stack_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='sca_tempreg_maps_Z_stack_smooth_%d' % num_strat, iterfield=['in_file'])
            sc_temp_reg_maps_Z_files_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='sca_tempreg_maps_Z_files_smooth_%d' % num_strat, iterfield=['in_file'])

            try:
                node, out_file = strat.get_node_from_resource_pool('sca_tempreg_maps_stack')
                node2, out_file2 = strat.get_node_from_resource_pool('sca_tempreg_maps_z_stack')
                node3, out_file3 = strat.get_node_from_resource_pool('sca_tempreg_maps_z_files')
                node4, out_file4 = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')

                # non-normalized stack
                workflow.connect(node, out_file,
                                 sc_temp_reg_maps_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 sc_temp_reg_maps_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 sc_temp_reg_maps_smooth, 'operand_files')

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
                print 'Invalid Connection: sca_tempreg smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise
            strat.append_name(sc_temp_reg_maps_smooth.name)
            strat.update_resource_pool({'sca_tempreg_maps_stack_smooth':(sc_temp_reg_maps_smooth, 'out_file'),
                                       'sca_tempreg_maps_z_stack_smooth':(sc_temp_reg_maps_Z_stack_smooth, 'out_file'),
                                       'sca_tempreg_maps_z_files_smooth':(sc_temp_reg_maps_Z_files_smooth, 'out_file')})

            create_log_node(sc_temp_reg_maps_smooth, 'out_file', num_strat)
            num_strat += 1
    strat_list += new_strat_list


    """
    Smoothing Temporal Regression for Dual Regression
    """
    new_strat_list = []
    num_strat = 0

    if (1 in c.runDualReg) and (1 in c.runSpatialRegression) and c.fwhm != None:
        for strat in strat_list:

            dr_temp_reg_maps_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                                              name='dr_tempreg_maps_stack_smooth_%d' % num_strat)
            dr_temp_reg_maps_Z_stack_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                                              name='dr_tempreg_maps_Z_stack_smooth_%d' % num_strat)
            dr_temp_reg_maps_Z_files_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                              name='dr_tempreg_maps_Z_files_smooth_%d' % num_strat, iterfield=['in_file'])

            try:
                node, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_stack_to_standard')
                node2, out_file2 = strat.get_node_from_resource_pool('dr_tempreg_maps_z_stack_to_standard')
                node3, out_file3 = strat.get_node_from_resource_pool('dr_tempreg_maps_z_files_to_standard')
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
                workflow.connect(node3, out_file3,
                                 dr_temp_reg_maps_Z_files_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                 dr_temp_reg_maps_Z_files_smooth, 'op_string')

                workflow.connect(node4, out_file4,
                                 dr_temp_reg_maps_Z_files_smooth, 'operand_files')

            except:
                print 'Invalid Connection: dr_tempreg smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise
            strat.append_name(dr_temp_reg_maps_smooth.name)
            strat.update_resource_pool({'dr_tempreg_maps_stack_smooth':(dr_temp_reg_maps_smooth, 'out_file'),
                                       'dr_tempreg_maps_z_stack_smooth':(dr_temp_reg_maps_Z_stack_smooth, 'out_file'),
                                       'dr_tempreg_maps_z_files_smooth':(dr_temp_reg_maps_Z_files_smooth, 'out_file')})
            create_log_node(dr_temp_reg_maps_smooth, 'out_file', num_strat)
            num_strat += 1
    strat_list += new_strat_list


    """
    Transforming SCA Voxel Z scores to MNI
    """
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runSCA) and (1 in c.runVoxelTimeseries):
        for strat in strat_list:

            if 'FSL' in c.regOption:

                sca_seed_Z_to_standard = pe.MapNode(interface=fsl.ApplyWarp(),
                               name='sca_seed_Z_to_standard_%d' % num_strat, iterfield=['in_file'])

                sca_seed_Z_to_standard.inputs.ref_file = c.standard


                try:

                    node, out_file = strat.get_node_from_resource_pool('sca_seed_Z')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard, 'premat')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard, 'field_file')

                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'sca_seed_Z_to_standard':(sca_seed_Z_to_standard, 'out_file')})
                strat.append_name(sca_seed_Z_to_standard.name)
                num_strat += 1


            elif 'ANTS' in c.regOption:

                sca_seed_Z_to_standard = create_apply_ants_xfm(3, 1,
                               name='sca_seed_Z_to_standard_%d' % num_strat)

                sca_seed_Z_to_standard.inputs.inputspec.warp_reference = c.standard


                try:

                    node, out_file = strat.get_node_from_resource_pool('sca_seed_Z')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard, 'inputspec.in_file')

                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard, 'inputspec.ants_affine')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard, 'inputspec.func_anat_affine')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard, 'inputspec.conversion_reference')

                    node, out_file = strat.get_node_from_resource_pool('sca_seed_Z')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard, 'inputspec.conversion_source')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard, 'inputspec.nonlinear_field')


                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'sca_seed_Z_to_standard':(sca_seed_Z_to_standard, 'outputspec.out_file')})
                strat.append_name(sca_seed_Z_to_standard.name)
                num_strat += 1



    strat_list += new_strat_list



    """
    Transforming SCA ROI Z scores to MNI
    """
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runSCA) and (1 in c.runROITimeseries):
        for strat in strat_list:

            if 'FSL' in c.regOption:

                sca_roi_Z_to_standard = pe.MapNode(interface=fsl.ApplyWarp(),
                               name='sca_roi_Z_to_standard_%d' % num_strat, iterfield=['in_file'])

                sca_roi_Z_to_standard.inputs.ref_file = c.standard


                try:

                    node, out_file = strat.get_node_from_resource_pool('sca_roi_Z')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard, 'in_file')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard, 'premat')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard, 'field_file')

                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'sca_roi_Z_to_standard':(sca_roi_Z_to_standard, 'out_file')})
                strat.append_name(sca_roi_Z_to_standard.name)
                num_strat += 1


            elif 'ANTS' in c.regOption:

                sca_roi_Z_to_standard = create_apply_ants_xfm(3, 1,
                               name='sca_roi_Z_to_standard_%d' % num_strat)

                sca_roi_Z_to_standard.inputs.inputspec.warp_reference = c.standard


                try:

                    node, out_file = strat.get_node_from_resource_pool('sca_roi_Z')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard, 'inputspec.in_file')

                    node, out_file = strat.get_node_from_resource_pool('ants_affine_xfm')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard, 'inputspec.ants_affine')

                    node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard, 'inputspec.func_anat_affine')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard, 'inputspec.nonlinear_field')

                    node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard, 'inputspec.conversion_reference')

                    node, out_file = strat.get_node_from_resource_pool('sca_roi_Z')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard, 'inputspec.conversion_source')

                except:
                    print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.update_resource_pool({'sca_roi_Z_to_standard':(sca_roi_Z_to_standard, 'outputspec.out_file')})
                strat.append_name(sca_roi_Z_to_standard.name)
                num_strat += 1


    strat_list += new_strat_list



    """
    
    Smoothing SCA seed based Z scores and or possibly Z scores in MNI 
    """
    new_strat_list = []
    num_strat = 0

    if (1 in c.runSCA) and (1 in c.runVoxelTimeseries) and c.fwhm != None:
        for strat in strat_list:


            sca_seed_Z_to_standard_smooth = None

            sca_seed_Z_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='sca_seed_Z_smooth_%d' % num_strat, iterfield=['in_file'])

            try:
                node, out_file = strat.get_node_from_resource_pool('sca_seed_Z')
                workflow.connect(node, out_file,
                                 sca_seed_Z_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                sca_seed_Z_smooth, 'op_string')
                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 sca_seed_Z_smooth, 'operand_files')

            except:
                print 'Invalid Connection: sca_Z smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise
            strat.append_name(sca_seed_Z_smooth.name)
            strat.update_resource_pool({'sca_seed_Z_smooth':(sca_seed_Z_smooth, 'out_file')})

            if 1 in c.runRegisterFuncToMNI:

                sca_seed_Z_to_standard_smooth = sca_seed_Z_smooth.clone('sca_seed_Z_to_standard_smooth_%d' % num_strat)

                try:

                    node, out_file = strat.get_node_from_resource_pool('sca_seed_Z_to_standard')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard_smooth, 'in_file')
                    workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                    sca_seed_Z_to_standard_smooth, 'op_string')
                    node, out_file = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')
                    workflow.connect(node, out_file,
                                     sca_seed_Z_to_standard_smooth, 'operand_files')



                except:

                    print 'Invalid Connection: sca_seed_Z_to_standard smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.append_name(sca_seed_Z_to_standard_smooth.name)
                strat.update_resource_pool({'sca_seed_Z_to_standard_smooth':(sca_seed_Z_to_standard_smooth, 'out_file')})
                create_log_node(sca_seed_Z_to_standard_smooth, 'out_file', num_strat)
            num_strat += 1
    strat_list += new_strat_list


    """
    
    Smoothing SCA roi based Z scores and or possibly Z scores in MNI 
    """
    if (1 in c.runSCA) and (1 in c.runROITimeseries) and c.fwhm != None:
        for strat in strat_list:


            sca_roi_Z_to_standard_smooth = None

            sca_roi_Z_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='sca_roi_Z_smooth_%d' % num_strat, iterfield=['in_file'])

            try:
                node, out_file = strat.get_node_from_resource_pool('sca_roi_Z')
                workflow.connect(node, out_file,
                                 sca_roi_Z_smooth, 'in_file')
                workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                sca_roi_Z_smooth, 'op_string')
                node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
                workflow.connect(node, out_file,
                                 sca_roi_Z_smooth, 'operand_files')

            except:
                print 'Invalid Connection: sca_Z smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise
            strat.append_name(sca_roi_Z_smooth.name)
            strat.update_resource_pool({'sca_roi_Z_smooth':(sca_roi_Z_smooth, 'out_file')})

            if 1 in c.runRegisterFuncToMNI:

                sca_roi_Z_to_standard_smooth = sca_roi_Z_smooth.clone('sca_roi_Z_to_standard_smooth_%d' % num_strat)

                try:


                    node, out_file = strat.get_node_from_resource_pool('sca_roi_Z_to_standard')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard_smooth, 'in_file')
                    workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                    sca_roi_Z_to_standard_smooth, 'op_string')
                    node, out_file = strat.get_node_from_resource_pool('functional_brain_mask_to_standard')
                    workflow.connect(node, out_file,
                                     sca_roi_Z_to_standard_smooth, 'operand_files')



                except:

                    print 'Invalid Connection: sca_roi_Z_to_standard smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.append_name(sca_roi_Z_smooth.name)
                strat.update_resource_pool({'sca_roi_Z_to_standard_smooth':(sca_roi_Z_to_standard_smooth, 'out_file')})
            num_strat += 1
    strat_list += new_strat_list



    """
    Inserting Surface Registration
    """
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
                print 'Invalid Connection: Surface Registration Workflow:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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

    """
    Inserting vertices based timeseries
    """
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
                print 'Invalid Connection: Vertices Time Series Extraction Workflow:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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

    """
    Inserting Network centrality
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runNetworkCentrality:

        for strat in strat_list:



            resample_functional_to_template = pe.Node(interface=fsl.FLIRT(),
                                                  name='resample_functional_to_template_%d' % num_strat)
            resample_functional_to_template.inputs.interp = 'trilinear'
            resample_functional_to_template.inputs.apply_xfm = True
            resample_functional_to_template.inputs.in_matrix_file = c.identityMatrix

            template_dataflow = create_mask_dataflow(c.templateSpecificationFile, 'template_dataflow_%d' % num_strat)

            network_centrality = create_resting_state_graphs(c.memoryAllocatedForDegreeCentrality, 'network_centrality_%d' % num_strat)
            network_centrality.inputs.inputspec.threshold_option = c.correlationThresholdOption
            network_centrality.inputs.inputspec.threshold = c.correlationThreshold
            network_centrality.inputs.centrality_options.weight_options = c.centralityWeightOptions
            network_centrality.inputs.centrality_options.method_options = c.centralityMethodOptions



            try:

                node, out_file = strat.get_node_from_resource_pool('functional_mni')

                # resample the input functional file to template(roi/mask)
                workflow.connect(node, out_file,
                                 resample_functional_to_template, 'in_file')
                workflow.connect(template_dataflow, 'outputspec.out_file',
                                 resample_functional_to_template, 'reference')

                workflow.connect(resample_functional_to_template, 'out_file',
                                 network_centrality, 'inputspec.subject')
                workflow.connect(template_dataflow, 'outputspec.out_file',
                                 network_centrality, 'inputspec.template')


                strat.append_name(network_centrality.name)

                strat.update_resource_pool({'centrality_outputs' : (network_centrality, 'outputspec.centrality_outputs')})
                
                create_log_node(network_centrality, 'outputspec.centrality_outputs', num_strat)

                # if smoothing is required
                if c.fwhm != None :

                    z_score = get_zscore('centrality_zscore_%d' % num_strat)

                    smoothing = pe.MapNode(interface=fsl.MultiImageMaths(),
                                       name='network_centrality_smooth_%d' % num_strat,
                                       iterfield=['in_file'])


                    # calculate zscores
                    workflow.connect(template_dataflow, 'outputspec.out_file',
                                     z_score, 'inputspec.mask_file')
                    workflow.connect(network_centrality, 'outputspec.centrality_outputs',
                                     z_score, 'inputspec.input_file')

                    # connecting zscores to smoothing
                    workflow.connect(template_dataflow, 'outputspec.out_file',
                                     smoothing, 'operand_files')
                    workflow.connect(z_score, 'outputspec.z_score_img',
                                    smoothing, 'in_file')
                    workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                                     smoothing, 'op_string')

                    strat.append_name(smoothing.name)
                    strat.update_resource_pool({'centrality_outputs_smoothed' : (smoothing, 'out_file'),
                                                'centrality_outputs_zscore' :   (z_score, 'outputspec.z_score_img')})
                    
                    strat.append_name(smoothing.name)
                    create_log_node(smoothing, 'out_file', num_strat)

            except:
                print 'Invalid Connection: Network Centrality Workflow:', num_strat, ' resource_pool: ', strat.get_resource_pool()
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



    """
    Quality Control
    """


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
                    print 'unable to get resources for SNR plot'
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
                    print 'unable to get resources for Motion Parameters Plot'
                    raise


            # make FD plot and volumes removed
            if 1 in c.runGenerateMotionStatistics:

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
                    print 'unable to get resources for FD Plot'
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

                print 'Cannot generate QC montages for Skull Stripping: Resources Not Found'
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

                print 'Cannot generate QC montages for mni normalized anatomical: Resources Not Found'
                raise



            # make QC montages for CSF WM GM

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
                print 'Cannot generate QC montages for WM GM CSF masks: Resources Not Found'
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
                print 'Cannot generate QC montages for Mean Functional in T1 with T1 edge: Resources Not Found'
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
#                  #MNI_edge.inputs.file_ = c.standardResolutionBrain
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
                print 'Cannot generate QC montages for Mean Functional in MNI with MNI edge: Resources Not Found'
                raise


            # make QC montages for SCA ROI Smoothed Derivative
            if (1 in c.runSCA) and (1 in c.runROITimeseries):

                hist_ = hist.clone('hist_sca_roi_%d' % num_strat)
                hist_.inputs.measure = 'sca_roi'

                drop_percent = pe.MapNode(util.Function(input_names=['measure_file',
                                                     'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_sca_roi_%d' % num_strat, iterfield=['measure_file'])
                drop_percent.inputs.percent_ = 99.999
                if c.fwhm != None:

                    sca_overlay, out_file = strat.get_node_from_resource_pool('sca_roi_Z_to_standard_smooth')
                    montage_sca_roi = create_montage('montage_sca_roi_standard_smooth_%d' % num_strat,
                                    'cyan_to_yellow', 'sca_roi_smooth')

                    montage_sca_roi.inputs.inputspec.underlay = c.standardResolutionBrain

                    workflow.connect(sca_overlay, out_file,
                                     drop_percent, 'measure_file')

                    workflow.connect(drop_percent, 'modified_measure_file',
                                     montage_sca_roi, 'inputspec.overlay')

                    workflow.connect(sca_overlay, out_file,
                                     hist_, 'measure_file')
                    strat.update_resource_pool({'qc___sca_roi_smooth_a': (montage_sca_roi, 'outputspec.axial_png'),
                                            'qc___sca_roi_smooth_s': (montage_sca_roi, 'outputspec.sagittal_png'),
                                            'qc___sca_roi_smooth_hist': (hist_, 'hist_path')})

                    if not 9 in qc_montage_id_a:
                        qc_montage_id_a[9] = 'sca_roi_smooth_a'
                        qc_montage_id_s[9] = 'sca_roi_smooth_s'
                        qc_hist_id[9] = 'sca_roi_smooth_hist'


                else:

                    sca_overlay, out_file = strat.get_node_from_resource_pool('sca_roi_Z_to_standard')
                    montage_sca_roi = create_montage('montage_sca_roi_standard_%d' % num_strat,
                                    'cyan_to_yellow', 'sca_roi')

                    montage_sca_roi.inputs.inputspec.underlay = c.standardResolutionBrain
                    workflow.connect(sca_overlay, out_file,
                                     drop_percent, 'measure_file')

                    workflow.connect(drop_percent, 'modified_measure_file',
                                     montage_sca_roi, 'inputspec.overlay')

                    workflow.connect(sca_overlay, out_file,
                                     hist_, 'measure_file')

                    strat.update_resource_pool({'qc___sca_roi_a': (montage_sca_roi, 'outputspec.axial_png'),
                                            'qc___sca_roi_s': (montage_sca_roi, 'outputspec.sagittal_png'),
                                            'qc___sca_roi_hist': (hist_, 'hist_path')})

                    if not 9 in qc_montage_id_a:
                        qc_montage_id_a[9] = 'sca_roi_a'
                        qc_montage_id_s[9] = 'sca_roi_s'
                        qc_hist_id[9] = 'sca_roi_hist'



            # make QC montages for SCA Smoothed Derivative
            if (1 in c.runSCA) and (1 in c.runVoxelTimeseries):
                hist_ = hist.clone('hist_sca_seeds_%d' % num_strat)
                hist_.inputs.measure = 'sca_seeds'

                drop_percent = pe.MapNode(util.Function(input_names=['measure_file',
                                                     'percent_'],
                                       output_names=['modified_measure_file'],
                                       function=drop_percent_),
                                       name='dp_sca_seed_%d' % num_strat, iterfield=['measure_file'])
                drop_percent.inputs.percent_ = 99.999
                if c.fwhm != None:

                    sca_overlay, out_file = strat.get_node_from_resource_pool('sca_seed_Z_to_standard_smooth')
                    montage_sca_seeds = create_montage('montage_seed_standard_smooth_%d' % num_strat,
                                    'cyan_to_yellow', 'sca_seed_smooth')

                    montage_sca_seeds.inputs.inputspec.underlay = c.standardResolutionBrain
                    workflow.connect(sca_overlay, out_file,
                                     drop_percent, 'measure_file')

                    workflow.connect(drop_percent, 'modified_measure_file',
                                     montage_sca_seeds, 'inputspec.overlay')

                    workflow.connect(sca_overlay, out_file,
                                     hist_, 'measure_file')

                    strat.update_resource_pool({'qc___sca_seeds_smooth_a': (montage_sca_seeds, 'outputspec.axial_png'),
                                            'qc___sca_seeds_smooth_s': (montage_sca_seeds, 'outputspec.sagittal_png'),
                                            'qc___sca_seeds_smooth_hist': (hist_, 'hist_path')})

                    if not 10 in qc_montage_id_a:
                        qc_montage_id_a[10] = 'sca_seeds_smooth_a'
                        qc_montage_id_s[10] = 'sca_seeds_smooth_s'
                        qc_hist_id[10] = 'sca_seeds_smooth_hist'

                else:
                
                    sca_overlay, out_file = strat.get_node_from_resource_pool('sca_seed_Z_to_standard')
                    montage_sca_seeds = create_montage('montage_sca_seed_standard_%d' % num_strat,
                                    'cyan_to_yellow', 'sca_seed')

                    montage_sca_seeds.inputs.inputspec.underlay = c.standardResolutionBrain
                    workflow.connect(sca_overlay, out_file,
                                     drop_percent, 'measure_file')

                    workflow.connect(drop_percent, 'modified_measure_file',
                                     montage_sca_seeds, 'inputspec.overlay')

                    workflow.connect(sca_overlay, out_file,
                                     hist_, 'measure_file')
                    strat.update_resource_pool({'qc___sca_seeds_a': (montage_sca_seeds, 'outputspec.axial_png'),
                                            'qc___sca_seeds_s': (montage_sca_seeds, 'outputspec.sagittal_png'),
                                            'qc___sca_seeds_hist': (hist_, 'hist_path')})

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

                    montage_centrality.inputs.inputspec.underlay = c.standardResolutionBrain
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

                    montage_centrality.inputs.inputspec.underlay = c.standardResolutionBrain
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

                    temporal_regression_sca_overlay, out_file = strat.get_node_from_resource_pool('sca_tempreg_maps_z_files_smooth')
                    montage_temporal_regression_sca = create_montage('montage_temporal_regression_sca_%d' % num_strat,
                                      'cyan_to_yellow', 'temporal_regression_sca_smooth')

                    montage_temporal_regression_sca.inputs.inputspec.underlay = c.standardResolutionBrain
                    strat.update_resource_pool({'qc___temporal_regression_sca_smooth_a': (montage_temporal_regression_sca, 'outputspec.axial_png'),
                                            'qc___temporal_regression_sca_smooth_s': (montage_temporal_regression_sca, 'outputspec.sagittal_png'),
                                            'qc___temporal_regression_sca_smooth_hist': (hist_, 'hist_path')})

                    if not 12 in qc_montage_id_a:
                        qc_montage_id_a[12] = 'temporal_regression_sca_smooth_a'
                        qc_montage_id_s[12] = 'temporal_regression_sca_smooth_s'
                        qc_hist_id[12] = 'temporal_regression_sca_smooth_hist'

                else:
                    temporal_regression_sca_overlay, out_file = strat.get_node_from_resource_pool('sca_tempreg_maps_z_files')
                    montage_temporal_regression_sca = create_montage('montage_temporal_regression_sca_%d' % num_strat,
                                      'cyan_to_yellow', 'temporal_regression_sca')

                    montage_temporal_regression_sca.inputs.inputspec.underlay = c.standardResolutionBrain
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

                    temporal_dual_regression_overlay, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_z_files_smooth')
                    montage_temporal_dual_regression = create_montage('montage_temporal_dual_regression_%d' % num_strat,
                                      'cyan_to_yellow', 'temporal_dual_regression_smooth')

                    montage_temporal_dual_regression.inputs.inputspec.underlay = c.standardResolutionBrain
                    strat.update_resource_pool({'qc___temporal_dual_regression_smooth_a': (montage_temporal_dual_regression, 'outputspec.axial_png'),
                                            'qc___temporal_dual_regression_smooth_s': (montage_temporal_dual_regression, 'outputspec.sagittal_png'),
                                            'qc___temporal_dual_regression_smooth_hist': (hist_, 'hist_path')})
                    if not 13 in qc_montage_id_a:
                        qc_montage_id_a[13] = 'temporal_dual_regression_smooth_a'
                        qc_montage_id_s[13] = 'temporal_dual_regression_smooth_s'
                        qc_hist_id[13] = 'temporal_dual_regression_smooth_hist'


                else:
                    temporal_dual_regression_overlay, out_file = strat.get_node_from_resource_pool('dr_tempreg_maps_z_files')
                    montage_temporal_dual_regression = create_montage('montage_temporal_dual_regression_%d' % num_strat,
                                      'cyan_to_yellow', 'temporal_dual_regression')

                    montage_temporal_dual_regression.inputs.inputspec.underlay = c.standardResolutionBrain
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

                vmhc_overlay, out_file = strat.get_node_from_resource_pool('vmhc_z_score_stat_map')
                montage_vmhc = create_montage('montage_vmhc_%d' % num_strat,
                                  'cyan_to_yellow', 'vmhc_smooth')

                montage_vmhc.inputs.inputspec.underlay = c.standardResolutionBrain
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
                                       name='dp_reho%d' % num_strat)
                drop_percent.inputs.percent_ = 99.999

                if c.fwhm != None:
                    reho_overlay, out_file = strat.get_node_from_resource_pool('reho_Z_to_standard_smooth')
                    montage_reho = create_montage('montage_reho_%d' % num_strat,
                                  'cyan_to_yellow', 'reho_standard_smooth')
                    montage_reho.inputs.inputspec.underlay = c.standardResolutionBrain
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
                    reho_overlay, out_file = strat.get_node_from_resource_pool('reho_Z_to_standard')
                    montage_reho = create_montage('montage_reho_%d' % num_strat,
                                  'cyan_to_yellow', 'reho_standard')
                    montage_reho.inputs.inputspec.underlay = c.standardResolutionBrain
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
                                       name='dp_alff%d' % num_strat)
                drop_percent.inputs.percent_ = 99.7

                drop_percent_falff = drop_percent.clone('dp_falff%d' % num_strat)
                drop_percent_falff.inputs.percent_ = 99.999

                if c.fwhm != None:
                    alff_overlay, out_file = strat.get_node_from_resource_pool('alff_Z_to_standard_smooth')
                    falff_overlay, out_file_f = strat.get_node_from_resource_pool('falff_Z_to_standard_smooth')
                    montage_alff = create_montage('montage_alff_%d' % num_strat,
                                  'cyan_to_yellow', 'alff_standard_smooth')
                    montage_alff.inputs.inputspec.underlay = c.standardResolutionBrain
                    montage_falff = create_montage('montage_falff_%d' % num_strat,
                                  'cyan_to_yellow', 'falff_standard_smooth')
                    montage_falff.inputs.inputspec.underlay = c.standardResolutionBrain
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
                    alff_overlay, out_file = strat.get_node_from_resource_pool('alff_Z_to_standard')
                    falff_overlay, out_file = strat.get_node_from_resource_pool('falff_Z_to_standard')
                    montage_alff = create_montage('montage_alff_%d' % num_strat,
                                  'cyan_to_yellow', 'alff_standard')
                    montage_alff.inputs.inputspec.underlay = c.standardResolutionBrain
                    montage_falff = create_montage('montage_falff_%d' % num_strat,
                                  'cyan_to_yellow', 'falff_standard')
                    montage_falff.inputs.inputspec.underlay = c.standardResolutionBrain
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




            num_strat += 1



    ###################### end of workflow ###########


    try:
        workflow.write_graph(graph2use='orig')
    except:
        pass


    """
    Datasink
    """
    import networkx as nx
    num_strat = 0
    sink_idx = 0
    pip_ids = []
    
    wf_names = []
    scan_ids = ['scan_anat']
    for id in sub_dict['rest']:
        scan_ids.append('scan_'+ str(id))
    
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
                    
                    print name, ' ~~~ ', 2 ** workflow_bit_id[name]
                    hash_val += 2 ** workflow_bit_id[name]

        if p_name == None or p_name == 'None':
            pipeline_id = ''
            pipeline_id = linecache.getline(os.path.realpath(os.path.join(CPAC.__path__[0], 'utils', 'pipeline_names.py')), hash_val)
            pipeline_id = pipeline_id.rstrip('\r\n')
            if pipeline_id == '':
                print 'hash value ', hash_val, ' is greater than the number of words'
                print 'resorting to crc32 value as pipeline_id'
                pipeline_id = zlib.crc32(strat_tag)
        else:
            pipeline_id = p_name
            #if running multiple pipelines with gui, need to change this in future
            p_name = None

        print 'strat_tag,  ~~~~~ , hash_val,  ~~~~~~ , pipeline_id: ', strat_tag, ' ~~~~~ ', hash_val, ' ~~~~~~ ', pipeline_id
        pip_ids.append(pipeline_id)
        wf_names.append(strat.get_name())

        for key in sorted(rp.keys()):

            ds = pe.Node(nio.DataSink(), name='sinker_%d' % sink_idx)
            ds.inputs.base_directory = c.outputDirectory
            ds.inputs.container = os.path.join('pipeline_%s' % pipeline_id, subject_id)
            ds.inputs.regexp_substitutions = [(r"/_sca_roi(.)*[/]", '/'),
                                              (r"/_smooth_centrality_(\d)+[/]", '/'),
                                              (r"/_z_score(\d)+[/]", "/"),
                                              (r"/_dr_tempreg_maps_Z_files_smooth_(\d)+[/]", "/"),
                                              (r"/_sca_tempreg_maps_Z_files_smooth_(\d)+[/]", "/"),
                                              (r"/qc___", '/qc/')]
            node, out_file = rp[key]
            workflow.connect(node, out_file,
                             ds, key)
            print 'node, out_file, key: ', node, out_file, key

            if 1 in c.runSymbolicLinks:

                link_node = pe.Node(interface=util.Function(input_names=['in_file', 'strategies',
                                        'subject_id', 'pipeline_id', 'helper'],
                                        output_names=[],
                                        function=prepare_symbolic_links),
                                        name='link_%d' % sink_idx)
               
                link_node.inputs.strategies = strategies
                link_node.inputs.subject_id = subject_id
                link_node.inputs.pipeline_id = 'pipeline_%s' % (pipeline_id)
                link_node.inputs.helper = dict(strategy_tag_helper_symlinks)

                workflow.connect(ds, 'out_file', link_node, 'in_file')
            sink_idx += 1
            print 'sink index: ', sink_idx

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
            print "Cannot Create the strategy and pipeline graph, dot or/and pygraphviz is not installed"
            pass


        print d_name, '*'
        num_strat += 1

    create_log_template(pip_ids, wf_names, scan_ids, subject_id, log_dir)



    workflow.run(plugin='MultiProc', plugin_args={'n_procs': c.numCoresPerSubject})

    ###workflow.run(plugin='Linear')

    """
    try:

        workflow.run(plugin='MultiProc', plugin_args={'n_procs': c.numCoresPerSubject})

    except Exception as e:

        print "Error: CPAC Pipeline has failed."
        print ""
        print e
        print type(e)
        ###raise Exception
    """

    for count, id in enumerate(pip_ids):
        for scan in scan_ids:
            create_log_node(None, None, count, scan).run()
        
    sub_w_path = os.path.join(c.workingDirectory, wfname)
        

    if 1 in c.generateQualityControlImages:

        for pip_id in pip_ids:

            f_path = os.path.join(os.path.join(c.outputDirectory, 'pipeline_' + pip_id), subject_id)

            f_path = os.path.join(f_path, 'qc_files_here')

            generateQCPages(f_path, qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id)


        ### Automatically generate QC index page
        create_all_qc.run(c.outputDirectory)



    if c.removeWorkingDir:
        try:
            if os.path.exists(sub_w_path):
                import shutil
                print "removing dir -> ", sub_w_path
                shutil.rmtree(sub_w_path)
        except:
            print "Couldn't remove subjects %s working directory" % (wfname)
            pass

    print "End of subject workflow ", wfname

    print >>timing, "CPAC run complete."
    print >>timing, "Elapsed run time (minutes): ", ((time.time() - pipeline_start_time)/60)
    print >>timing, ""

    timing.close()

    return workflow




def run(config, subject_list_file, indx, strategies, \
     maskSpecificationFile, roiSpecificationFile, templateSpecificationFile, p_name = None):
    import commands
    commands.getoutput('source ~/.bashrc')
    import os
    import sys
    import argparse
    import pickle
    import yaml


    c = Configuration(yaml.load(open(os.path.realpath(config), 'r')))

    try:
        sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    except:
        raise Exception ("Subject list is not in proper YAML format. Please check your file")

    sub_dict = sublist[int(indx) - 1]


    c.maskSpecificationFile = maskSpecificationFile
    c.roiSpecificationFile = roiSpecificationFile
    c.templateSpecificationFile = templateSpecificationFile

    prep_workflow(sub_dict, c, pickle.load(open(strategies, 'r')), p_name)
