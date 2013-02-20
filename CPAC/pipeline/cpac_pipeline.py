import os
import sys
import copy
import argparse
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
from   nipype.pipeline.utils import format_dot

from multiprocessing import Process
import CPAC
from CPAC.anat_preproc.anat_preproc import create_anat_preproc
from CPAC.func_preproc.func_preproc import create_func_preproc
from CPAC.seg_preproc.seg_preproc import create_seg_preproc

from CPAC.registration import create_nonlinear_register, create_bbregister_func_to_mni
from CPAC.nuisance import create_nuisance, bandpass_voxels

from CPAC.median_angle import create_median_angle_correction
from CPAC.generate_motion_statistics import motion_power_statistics
from CPAC.generate_motion_statistics import fristons_twenty_four
from CPAC.scrubbing import create_scrubbing_preproc
from CPAC.timeseries import create_surface_registration, get_voxel_timeseries,\
                            get_roi_timeseries, get_vertices_timeseries
from CPAC.network_centrality import create_resting_state_graphs, get_zscore
from CPAC.utils.datasource import *
from CPAC.utils import Configuration
from CPAC.utils.utils import extract_one_d, set_gauss, \
                             prepare_symbolic_links, \
                             get_scan_params, get_tr
from CPAC.vmhc.vmhc import create_vmhc
from CPAC.reho.reho import create_reho
from CPAC.alff.alff import create_alff
from CPAC.sca.sca import create_sca
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
            


def prep_workflow(sub_dict, c, strategies):
    
    print '********************',c.standardResolutionBrain
    
    if sub_dict['unique_id']:
        subject_id = sub_dict['subject_id'] +"_"+ sub_dict['unique_id']
    else:
        subject_id = sub_dict['subject_id']
        
        
    wfname = 'resting_preproc_' + str(subject_id)
    workflow = pe.Workflow(name=wfname)
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {'hash_method': 'timestamp', 'crashdump_dir': os.path.abspath(c.crashLogDirectory)}

    mflow = None
    pflow = None

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

            strat.append_name('anat_preproc')



            strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
            # add stuff to resource pool if we need it 

            strat.update_resource_pool({'anatomical_brain':(anat_preproc, 'outputspec.brain')})
            strat.update_resource_pool({'anatomical_reorient':(anat_preproc, 'outputspec.reorient')})

            num_strat += 1

    strat_list += new_strat_list


    """
    Inserting Registration Preprocessing
    Workflow
    """
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runRegistrationPreprocessing:
        workflow_bit_id['anat_mni_register'] = workflow_counter
        for strat in strat_list:
            reg_anat_mni = create_nonlinear_register('anat_mni_register_%d' % num_strat)

            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 reg_anat_mni, 'inputspec.input_brain')
                node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                workflow.connect(node, out_file,
                                 reg_anat_mni, 'inputspec.input_skull')

                reg_anat_mni.inputs.inputspec.reference_brain = c.standardResolutionBrainAnat
                reg_anat_mni.inputs.inputspec.reference_skull = c.standardAnat
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

            strat.append_name('anat_mni_register')
            strat.set_leaf_properties(reg_anat_mni, 'outputspec.output_brain')

            strat.update_resource_pool({'anatomical_to_mni_linear_xfm':(reg_anat_mni, 'outputspec.linear_xfm'),
                                        'anatomical_to_mni_nonlinear_xfm':(reg_anat_mni, 'outputspec.nonlinear_xfm'),
                                        'mni_to_anatomical_linear_xfm':(reg_anat_mni, 'outputspec.invlinear_xfm'),
                                        'mni_normalized_anatomical':(reg_anat_mni, 'outputspec.output_brain')})

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
            seg_preproc = create_seg_preproc('seg_preproc_%d'%num_strat)
            
            try:
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 seg_preproc, 'inputspec.brain')

                node, out_file = strat.get_node_from_resource_pool('mni_to_anatomical_linear_xfm')
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
                
            strat.append_name('seg_preproc')
            strat.update_resource_pool({'anatomical_gm_mask' : (seg_preproc,'outputspec.gm_mask'),
                                        'anatomical_csf_mask': (seg_preproc, 'outputspec.csf_mask'),
                                        'anatomical_wm_mask' : (seg_preproc, 'outputspec.wm_mask'),
                                        'seg_probability_maps': (seg_preproc,'outputspec.probability_maps'),
                                        'seg_mixeltype': (seg_preproc, 'outputspec.mixeltype'),
                                        'seg_partial_volume_map': (seg_preproc, 'outputspec.partial_volume_map'),
                                        'seg_partial_volume_files': (seg_preproc, 'outputspec.partial_volume_files')})

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
            #Flow = create_func_datasource(sub_dict['rest'])
            #Flow.inputs.inputnode.subject = subject_id
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
            #a node which checks if scan _parameters are present for each scan
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
                                               function = get_scan_params ),
                                 name  = 'scan_params_%d' % num_strat)
            
            convert_tr = pe.Node(util.Function(input_names = ['tr'],
                                               output_names= ['tr'],
                                               function = get_tr),
                                 name = 'convert_tr_%d' % num_strat)
            
            #if scan parameters are available slice timing correction is 
            #turned on 
            if slice_timing:
                
                func_preproc = create_func_preproc(slice_timing_correction=True, wf_name = 'func_preproc_%d' % num_strat)
                
                #getting the scan parameters
                workflow.connect(funcFlow, 'outputspec.subject',
                                 scan_params, 'subject')
                workflow.connect(funcFlow, 'outputspec.scan',
                                 scan_params, 'scan')
                scan_params.inputs.subject_map = sub_dict
                scan_params.inputs.start_indx = c.startIdx
                scan_params.inputs.stop_indx = c.stopIdx
                
                #passing the slice timing correction parameters
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
                func_preproc = create_func_preproc(slice_timing_correction=False, wf_name ='func_preproc_%d' % num_strat)
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
                
            strat.append_name('func_preproc')

            strat.set_leaf_properties(func_preproc, 'outputspec.preprocessed')

            # add stuff to resource pool if we need it 
            if slice_timing:
                strat.update_resource_pool({'slice_timing_corrected': (func_preproc, 'outputspec.slice_time_corrected')})
            strat.update_resource_pool({'mean_functional':(func_preproc, 'outputspec.example_func')})
            strat.update_resource_pool({'functional_preprocessed_mask':(func_preproc, 'outputspec.preprocessed_mask')})
            strat.update_resource_pool({'movement_parameters':(func_preproc, 'outputspec.movement_parameters')})
            strat.update_resource_pool({'max_displacement':(func_preproc, 'outputspec.max_displacement')})
            strat.update_resource_pool({'preprocessed':(func_preproc, 'outputspec.preprocessed')})
            strat.update_resource_pool({'functional_brain_mask':(func_preproc, 'outputspec.mask')})
            strat.update_resource_pool({'motion_correct':( func_preproc, 'outputspec.motion_correct')})
            
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
            strat.append_name('fristons_parameter_model')

            strat.update_resource_pool({'movement_parameters':(fristons_model, 'outputspec.movement_file')})

            num_strat += 1
    strat_list += new_strat_list



    """
    Inserting Anatomical to Functional Registration
    """
    new_strat_list = []
    num_strat = 0

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
                
            strat.append_name('anat_to_func_register')

            strat.update_resource_pool({'anatomical_to_functional_xfm':(anat_to_func_reg, 'out_matrix_file'),
                                        'inverse_anatomical_to_functional_xfm':(inv_anat_to_func, 'out_file'),
                                        'functional_gm_mask':(func_gm, 'out_file'),
                                        'functional_wm_mask':(func_wm, 'out_file'),
                                        'functional_csf_mask':(func_csf, 'out_file')})

            num_strat += 1

    strat_list += new_strat_list
    
    """
    Inserting Register Functional to MNI Workflow
    """
    new_strat_list = []
    num_strat = 0
    workflow_counter += 1
    if 1 in c.runRegisterFuncToMNI:
        workflow_bit_id['func_to_mni'] = workflow_counter
        for strat in strat_list:
            func_to_mni = create_bbregister_func_to_mni('func_to_mni_%d' % num_strat)
            func_to_mni.inputs.inputspec.mni = c.standardResolutionBrain
            func_to_mni.inputs.inputspec.interp = 'trilinear'
            func_to_mni.inputs.inputspec.bbr_schedule = c.boundaryBasedRegistrationSchedule

            try:
                def pick_wm(seg_prob_list):
                    seg_prob_list.sort()
                    return seg_prob_list[-1]
                
                node, out_file = strat.get_node_from_resource_pool('mean_functional')
                workflow.connect(node, out_file,
                                 func_to_mni, 'inputspec.func')

                node, out_file = strat.get_node_from_resource_pool('seg_probability_maps')
                workflow.connect(node, (out_file, pick_wm),
                                 func_to_mni, 'inputspec.anat_wm_segmentation')
                
                node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                workflow.connect(node, out_file,
                                 func_to_mni, 'inputspec.anat_skull')
                
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 func_to_mni, 'inputspec.anat')
                
                node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_linear_xfm')
                workflow.connect(node, out_file,
                                 func_to_mni, 'inputspec.anat_to_mni_linear_xfm')
                
                node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
                workflow.connect(node, out_file,
                                 func_to_mni, 'inputspec.anat_to_mni_nonlinear_xfm')
            except:
                print 'Invalid Connection: Register Functional to MNI:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise
            
            if 0 in c.runRegisterFuncToMNI:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                tmp.name = list(strat.name)
                strat = tmp
                new_strat_list.append(strat)

            strat.append_name('func_to_mni')
            #strat.set_leaf_properties(func_mni_warp, 'out_file')

            strat.update_resource_pool({'mean_functional_in_mni':(func_to_mni, 'outputspec.mni_func'),
                                        'mean_functional_in_anat':(func_to_mni, 'outputspec.anat_func'),
                                        'anatomical_wm_edge':(func_to_mni, 'outputspec.anat_wm_edge'),
                                        'functional_to_anat_linear_xfm':(func_to_mni, 'outputspec.func_to_anat_linear_xfm'),
                                        'functional_to_mni_linear_xfm':(func_to_mni, 'outputspec.func_to_mni_linear_xfm'),
                                        'mni_to_functional_linear_xfm':(func_to_mni, 'outputspec.mni_to_func_linear_xfm')})

            num_strat += 1

    strat_list += new_strat_list

    """
    Inserting Generate Motion Statistics Workflow
    """
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runGenerateMotionStatistics:
        workflow_bit_id['gen_motion_stats'] = workflow_counter
        for strat in strat_list:
            
            gen_motion_stats = motion_power_statistics('gen_motion_stats_%d'% num_strat)
            gen_motion_stats.inputs.scrubbing_input.threshold = c.scrubbingThreshold
            gen_motion_stats.inputs.scrubbing_input.remove_frames_before = c.numRemovePrecedingFrames
            gen_motion_stats.inputs.scrubbing_input.remove_frames_after = c.numRemoveSubsequentFrames
            gen_motion_stats.get_node('scrubbing_input').iterables = ('threshold', c.scrubbingThreshold)
            
            try:
                ##**special case where the workflow is not getting outputs from resource pool
                #but is connected to functional datasource
                workflow.connect(funcFlow, 'outputspec.subject',
                             gen_motion_stats, 'inputspec.subject_id')
            
                workflow.connect(funcFlow, 'outputspec.scan',
                             gen_motion_stats, 'inputspec.scan_id')
                
                node,out_file = strat.get_node_from_resource_pool('motion_correct')
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
                
            strat.append_name('gen_motion_stats')

            strat.update_resource_pool({'frame_wise_displacement':(gen_motion_stats, 'outputspec.FD_1D'),
                                        'scrubbing_frames_excluded':(gen_motion_stats, 'outputspec.frames_ex_1D'),
                                        'scrubbing_frames_included':(gen_motion_stats, 'outputspec.frames_in_1D'),
                                        'power_params':(gen_motion_stats,'outputspec.power_params'),
                                        'motion_params':(gen_motion_stats, 'outputspec.motion_params')})

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
            nuisance = create_nuisance('nuisance_%d' % num_strat)
            
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
                
                node, out_file = strat.get_node_from_resource_pool('mni_to_anatomical_linear_xfm')
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

            strat.append_name('nuisance')

            strat.set_leaf_properties(nuisance, 'outputspec.subject')

            strat.update_resource_pool({'functional_nuisance_residuals':(nuisance, 'outputspec.subject')})

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

            strat.append_name('median_angle_corr')

            strat.set_leaf_properties(median_angle_corr, 'outputspec.subject')

            strat.update_resource_pool({'functional_median_angle_corrected':(median_angle_corr, 'outputspec.subject')})

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
            alff = create_alff('alff_%d' % num_strat)
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
            strat.append_name('alff_falff')
            strat.update_resource_pool({'alff_img':(alff, 'outputspec.alff_img')})
            strat.update_resource_pool({'falff_img':(alff, 'outputspec.falff_img')})
            strat.update_resource_pool({'alff_Z_img':(alff, 'outputspec.alff_Z_img')})
            strat.update_resource_pool({'falff_Z_img':(alff, 'outputspec.falff_Z_img')})

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

            strat.append_name('frequency_filter')

            strat.set_leaf_properties(frequency_filter, 'bandpassed_file')

            strat.update_resource_pool({'functional_freq_filtered':(frequency_filter, 'bandpassed_file')})

            num_strat += 1

    strat_list += new_strat_list




    """
    Inserting Scrubbing Workflow
    """
    new_strat_list = []
    num_strat = 0

    workflow_counter += 1
    if 1 in c.runScrubbing:
        workflow_bit_id['scrubbing'] = workflow_counter
        for strat in strat_list:
            
            scrubbing = create_scrubbing_preproc('scrubbing_%d'%num_strat)
            
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

            strat.append_name('scrubbing')

            strat.set_leaf_properties(scrubbing, 'outputspec.preprocessed')

            strat.update_resource_pool({'scrubbing_movement_parameters' : (scrubbing, 'outputspec.scrubbed_movement_parameters'),
                                        'scrubbed_preprocessed': (scrubbing, 'outputspec.preprocessed')})

            num_strat += 1

    strat_list += new_strat_list

    """
    Register Functional timeseries to MNI space
    """
    new_strat_list = []
    num_strat = 0

    for strat in strat_list:
        func_mni_warp = pe.Node(interface=fsl.ApplyWarp(),
                                name='func_mni_warp_%d' % num_strat)
        func_mni_warp.inputs.ref_file = c.standardResolutionBrain

        functional_brain_mask_to_standard = pe.Node(interface=fsl.ApplyWarp(),
                                                    name='functional_brain_mask_to_standard1_%d' % num_strat)
        functional_brain_mask_to_standard.inputs.interp = 'nn'
        functional_brain_mask_to_standard.inputs.ref_file = c.standard

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

            node, out_file = strat.get_node_from_resource_pool('functional_brain_mask')
            workflow.connect(node, out_file,
                             functional_brain_mask_to_standard, 'in_file')


            node, out_file = strat.get_node_from_resource_pool('functional_to_anat_linear_xfm')
            workflow.connect(node, out_file,
                             functional_brain_mask_to_standard, 'premat')

            node, out_file = strat.get_node_from_resource_pool('anatomical_to_mni_nonlinear_xfm')
            workflow.connect(node, out_file,
                             functional_brain_mask_to_standard, 'field_file')
        except:
            print 'Invalid Connection: Register Functional timeseries to MNI space:', num_strat, ' resource_pool: ', strat.get_resource_pool()
            raise

        strat.update_resource_pool({'functional_mni':(func_mni_warp, 'out_file'),
                                    'functional_brain_mask_to_standard':(functional_brain_mask_to_standard, 'out_file')})

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

            preproc = create_vmhc()
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

            except:
                print 'Invalid Connection: VMHC:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise


            strat.update_resource_pool({'vmhc_raw_score':(vmhc, 'outputspec.VMHC_FWHM_img')})
            strat.update_resource_pool({'vmhc_z_score':(vmhc, 'outputspec.VMHC_Z_FWHM_img')})
            strat.update_resource_pool({'vmhc_z_score_stat_map':(vmhc, 'outputspec.VMHC_Z_stat_FWHM_img')})
            strat.append_name('vmhc')
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
            strat.append_name('reho')
            num_strat += 1
    strat_list += new_strat_list



    """
    Transforming ALFF Z scores and fAlff Z scores to MNI
    """
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runALFF):
        for strat in strat_list:

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
            strat.append_name('alff_falff_to_standard')
            num_strat += 1
    strat_list += new_strat_list


    inputnode_fwhm = None
    if len(c.fwhm) > 0:

        inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')
        inputnode_fwhm.iterables = ("fwhm", c.fwhm)



    """    
    Smoothing ALFF fALFF Z scores and or possibly Z scores in MNI 
    """
    new_strat_list = []
    num_strat = 0
    if (1 in c.runALFF) and len(c.fwhm) > 0:
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
            strat.append_name('alff_falff_Z_smooth')
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

                strat.append_name('alff_falff_Z_to_standard_smooth')
                strat.update_resource_pool({'alff_Z_to_standard_smooth':(alff_Z_to_standard_smooth, 'out_file')})
                strat.update_resource_pool({'falff_Z_to_standard_smooth':(falff_Z_to_standard_smooth, 'out_file')})



            num_strat += 1
    strat_list += new_strat_list



    """
    Transforming ReHo Z scores to MNI
    """
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runReHo):
        for strat in strat_list:

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
            strat.append_name('reho_Z_transform')
            num_strat += 1
    strat_list += new_strat_list


    """
    
    Smoothing ReHo Z scores and or possibly Z scores in MNI 
    """
    new_strat_list = []
    num_strat = 0

    if (1 in c.runReHo) and len(c.fwhm) > 0:
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
            strat.append_name('reho_Z_smooth')
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
                    workflow.connect(functional_brain_mask_to_standard, 'out_file',
                                     reho_Z_to_standard_smooth, 'operand_files')



                except:

                    print 'Invalid Connection: reho_Z_to_standard smooth:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                    raise

                strat.append_name('reho_Z_to_standard_smooth')
                strat.update_resource_pool({'reho_Z_to_standard_smooth':(reho_Z_to_standard_smooth, 'out_file')})
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
            resample_functional_to_mask.inputs.interp = 'nearestneighbour'
            resample_functional_to_mask.inputs.apply_xfm = True
            resample_functional_to_mask.inputs.in_matrix_file = c.identityMatrix

            mask_dataflow = create_mask_dataflow(c.maskSpecificationFile, 'mask_dataflow_%d' % num_strat)

            voxel_timeseries = get_voxel_timeseries('voxel_timeseries_%d' % num_strat)
            voxel_timeseries.inputs.inputspec.output_type = c.voxelTSOutputs

            try:

                node, out_file = strat.get_node_from_resource_pool('functional_mni')

                #resample the input functional file to mask 
                workflow.connect(node, out_file,
                                 resample_functional_to_mask, 'in_file' )
                workflow.connect(mask_dataflow, 'select_mask.out_file',
                                 resample_functional_to_mask, 'reference')

                #connect it to the voxel_timeseries
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

            strat.append_name('voxel_timeseries')

            strat.update_resource_pool({'voxel_timeseries': (voxel_timeseries, 'outputspec.mask_outputs')})

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
                                                  name='resample_functional_to_roi_%d'%num_strat)
            resample_functional_to_roi.inputs.interp = 'nearestneighbour'
            resample_functional_to_roi.inputs.apply_xfm = True
            resample_functional_to_roi.inputs.in_matrix_file = c.identityMatrix

            roi_dataflow = create_roi_dataflow(c.roiSpecificationFile, 'roi_dataflow_%d'%num_strat)

            roi_timeseries = get_roi_timeseries('roi_timeseries_%d'%num_strat) 
            roi_timeseries.inputs.inputspec.output_type = c.roiTSOutputs

            try:

                node, out_file = strat.get_node_from_resource_pool('functional_mni')

                #resample the input functional file to roi
                workflow.connect(node, out_file,
                                 resample_functional_to_roi,'in_file' )
                workflow.connect(roi_dataflow, 'select_roi.out_file',
                                 resample_functional_to_roi, 'reference')

                #connect it to the roi_timeseries
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

            strat.append_name('roi_timeseries')

            strat.update_resource_pool({'roi_timeseries' : (roi_timeseries, 'outputspec.roi_outputs')})

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
            strat.append_name('sca_rois')
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
            strat.append_name('sca_seeds')
            num_strat += 1
    strat_list += new_strat_list


    """
    Transforming SCA Voxel Z scores to MNI
    """
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runSCA) and (1 in c.runVoxelTimeseries):
        for strat in strat_list:

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
            strat.append_name('sca_seed_based_to_mni')
            num_strat += 1
    strat_list += new_strat_list



    """
    Transforming SCA ROI Z scores to MNI
    """
    new_strat_list = []
    num_strat = 0


    if 1 in c.runRegisterFuncToMNI and (1 in c.runSCA) and (1 in c.runROITimeseries):
        for strat in strat_list:

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
            strat.append_name('sca_roi_based_to_mni')
            num_strat += 1
    strat_list += new_strat_list



    """
    
    Smoothing SCA seed based Z scores and or possibly Z scores in MNI 
    """
    new_strat_list = []
    num_strat = 0

    if (1 in c.runSCA) and (1 in c.runVoxelTimeseries) and len(c.fwhm) > 0:
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
            strat.append_name('sca_seed_Z_smooth')
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

                strat.append_name('sca_seed_Z_to_standard_smooth')
                strat.update_resource_pool({'sca_seed_Z_to_standard_smooth':(sca_seed_Z_to_standard_smooth, 'out_file')})
            num_strat += 1
    strat_list += new_strat_list


    """
    
    Smoothing SCA roi based Z scores and or possibly Z scores in MNI 
    """
    if (1 in c.runSCA) and (1 in c.runROITimeseries) and len(c.fwhm) > 0:
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
            strat.append_name('sca_roi_Z_to_smooth')
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

                strat.append_name('sca_roi_Z_to_standard_smooth')
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
            
            surface_reg = create_surface_registration('surface_reg_%d'%num_strat)
            surface_reg.inputs.inputspec.recon_subjects = c.reconSubjectsDirectory
            surface_reg.inputs.inputspec.subject_id = subject_id
            
            try:
                
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 surface_reg,'inputspec.rest' )
                
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
                
            strat.append_name('surface_registration')
            
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
            
            vertices_timeseries = get_vertices_timeseries('vertices_timeseries_%d'%num_strat)
            
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
                
            strat.append_name('surface_registration')
            
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
                                                  name='resample_functional_to_template_%d'%num_strat)
            resample_functional_to_template.inputs.interp = 'nearestneighbour'
            resample_functional_to_template.inputs.apply_xfm = True
            resample_functional_to_template.inputs.in_matrix_file = c.identityMatrix
            
            template_dataflow = create_mask_dataflow(c.templateSpecificationFile, 'template_dataflow_%d'%num_strat)
            
            network_centrality = create_resting_state_graphs(c.memoryAllocatedForDegreeCentrality, 'network_centrality_%d'%num_strat)
            network_centrality.inputs.inputspec.threshold_option = c.correlationThresholdOption
            network_centrality.inputs.inputspec.threshold = c.correlationThreshold
            network_centrality.inputs.centrality_options.weight_options = c.centralityWeightOptions
            network_centrality.inputs.centrality_options.method_options = c.centralityMethodOptions
            
            
            
            try:
                
                node, out_file = strat.get_node_from_resource_pool('functional_mni')
                
                #resample the input functional file to template(roi/mask) 
                workflow.connect(node, out_file,
                                 resample_functional_to_template, 'in_file' )
                workflow.connect(template_dataflow, 'outputspec.out_file',
                                 resample_functional_to_template, 'reference')
                
                workflow.connect(resample_functional_to_template, 'out_file',
                                 network_centrality, 'inputspec.subject')
                workflow.connect(template_dataflow, 'outputspec.out_file',
                                 network_centrality, 'inputspec.template')


                strat.append_name('network_centrality')
    
                strat.update_resource_pool({'centrality_outputs' : (network_centrality, 'outputspec.centrality_outputs')})
                
        
                #if smoothing is required
                if len(c.fwhm) > 0 :
                    
                    z_score = get_zscore('centrality_zscore_%d'%num_strat)
                    
                    smoothing = pe.MapNode(interface=fsl.MultiImageMaths(),
                                       name='smooth_centrality_%d'% num_strat, 
                                       iterfield=['in_file'])

                    
                    #calculate zscores
                    workflow.connect(template_dataflow, 'outputspec.out_file',
                                     z_score, 'inputspec.mask_file')
                    workflow.connect(network_centrality, 'outputspec.centrality_outputs',
                                     z_score, 'inputspec.input_file')
                    
                    #connecting zscores to smoothing
                    workflow.connect(template_dataflow, 'outputspec.out_file',
                                     smoothing, 'operand_files')
                    workflow.connect(z_score, 'outputspec.z_score_img',
                                    smoothing, 'in_file')
                    workflow.connect(inputnode_fwhm,('fwhm', set_gauss),
                                     smoothing, 'op_string')

                    
                    strat.update_resource_pool({'centrality_outputs_smoothed' : (smoothing, 'out_file'),
                                                'centrality_outputs_zscore' :   (z_score, 'outputspec.z_score_img')})
                    
                    
                
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
    for strat in strat_list:
        rp = strat.get_resource_pool()

        # build helper dictionary to assist with a clean strategy label for symlinks

        strategy_tag_helper_symlinks = {}

        if 'scrubbing' in strat.get_name():

            strategy_tag_helper_symlinks['_threshold'] = 1

        else:

            strategy_tag_helper_symlinks['_threshold'] = 0



        if 'seg_preproc' in strat.get_name():

            strategy_tag_helper_symlinks['_csf_threshold'] = 1
            strategy_tag_helper_symlinks['_wm_threshold'] = 1
            strategy_tag_helper_symlinks['_gm_threshold'] = 1

        else:
            strategy_tag_helper_symlinks['_csf_threshold'] = 0
            strategy_tag_helper_symlinks['_wm_threshold'] = 0
            strategy_tag_helper_symlinks['_gm_threshold'] = 0


        if 'median_angle_corr' in strat.get_name():

            strategy_tag_helper_symlinks['_target_angle_deg'] = 1

        else:
            strategy_tag_helper_symlinks['_target_angle_deg'] = 0


        if 'nuisance' in strat.get_name():

            strategy_tag_helper_symlinks['nuisance'] = 1

        else:
            strategy_tag_helper_symlinks['nuisance'] = 0

        strat_tag = ""

        hash_val = 0

        for name in strat.get_name():
            if not ('alff' in name.lower()) and not ('vmhc' in name.lower()) \
            and not ('reho' in name.lower()) and not ('sca' in name.lower()) \
            and not ('network_centrality' in name.lower()) and not ('timeseries' in name.lower()):

                strat_tag += name + '_'

                print name, ' ~~~ ', 2 ** workflow_bit_id[name]
                hash_val += 2 ** workflow_bit_id[name]


        pipeline_id = ''
        pipeline_id = linecache.getline(os.path.realpath(os.path.join(CPAC.__path__[0], 'utils', 'pipeline_names.py')), hash_val)
        pipeline_id = pipeline_id.rstrip('\r\n')
        if pipeline_id == '':
            print 'hash value ', hash_val, ' is greater than the number of words'
            print 'resorting to crc32 value as pipeline_id'
            pipeline_id = zlib.crc32(strat_tag)

        print strat_tag, ' ~~~~~ ', hash_val, ' ~~~~~~ ', pipeline_id

        for key in rp.keys():
            ds = pe.Node(nio.DataSink(), name='sinker_%d' % sink_idx)
            ds.inputs.base_directory = c.sinkDirectory
            ds.inputs.container = os.path.join('pipeline_%s' % pipeline_id, subject_id)
            ds.inputs.regexp_substitutions = [(r"/_sca_roi(.)*[/]", '/'),
                                              (r"/_smooth_centrality_(\d)+[/]", '/'),
                                              (r"/_z_score(\d)+[/]", "/")]
            node, out_file = rp[key]
            workflow.connect(node, out_file,
                             ds, key)

            if 1 in c.runSymbolicLinks:

                link_node = pe.Node(interface=util.Function(input_names=['in_file', 'strategies',
                                        'subject_id', 'pipeline_id', 'helper'],
                                        output_names=[],
                                        function=prepare_symbolic_links),
                                        name='link_%d' % sink_idx, iterfield=['in_file'])

                link_node.inputs.strategies = strategies
                link_node.inputs.subject_id = subject_id
                link_node.inputs.pipeline_id = 'pipeline_%s' % (pipeline_id)
                link_node.inputs.helper = dict(strategy_tag_helper_symlinks)

                workflow.connect(ds, 'out_file', link_node, 'in_file')
            sink_idx += 1

        d_name = os.path.join(c.sinkDirectory, ds.inputs.container)
        if not os.path.exists(d_name):
            os.makedirs(d_name)



        try:
            G = nx.DiGraph()
            strat_name = strat.get_name()
            G.add_edges_from([(strat_name[s], strat_name[s+1]) for s in range(len(strat_name)-1)])
            dotfilename = os.path.join(d_name, 'strategy.dot')
            nx.write_dot(G, dotfilename)
            format_dot(dotfilename, 'png')
        except:
            print "Cannot Create the strategy and pipeline graph, dot or/and pygraphviz is not installed"
            pass


        print d_name, '*'
        num_strat += 1

    workflow.run(plugin='MultiProc',
                         plugin_args={'n_procs': c.numCoresPerSubject})
#    workflow.run(updatehash=True)
    sub_w_path = os.path.join(c.workingDirectory, wfname)
    
    if c.removeWorkingDir:
        try:
            if os.path.exists(sub_w_path):
                import shutil
                print "removing dir -> ", sub_w_path
                shutil.rmtree(sub_w_path)
        except:
            print "Couldn't remove subjects %s working directory"%(wfname)
            pass
    
    print "End of subject workflow ", wfname
    
    return workflow




def run(config, subject_list_file, indx, strategies):
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

    prep_workflow(sub_dict, c, pickle.load(open(strategies, 'r')))
