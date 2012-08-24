import os
import sys
import copy
import argparse
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util

from multiprocessing import Process


from CPAC.anat_preproc.anat_preproc import create_anat_preproc
from CPAC.func_preproc.func_preproc import create_func_preproc
from CPAC.reg_preproc.reg_preproc import create_reg_preproc
from CPAC.seg_preproc.seg_preproc import create_seg_preproc

from CPAC.registration import create_nonlinear_register
from CPAC.nuisance import create_nuisance

from CPAC.utils.func_datasource import *
from CPAC.utils.anat_datasource import *

class strategy:

    def __init__(self):
        self.resource_pool = {}
        self.leaf_node = None
        self.leaf_out_file = None

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

def prep_workflow(sub_dict, seed_list, c):
    print "running for subject ", sub_dict
    subject_id = sub_dict['Subject_id'] +"_"+ sub_dict['Unique_id']
    wfname = 'resting_preproc_' + str(subject_id)
    workflow = pe.Workflow(name=wfname)
    workflow.base_dir = c.workingDirectory
    workflow.crash_dir = c.crashLogDirectory
    workflow.config['execution'] = {'hash_method': 'timestamp'}

    mflow = None
    pflow = None

    strat_list = []


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
                tmp.out_file = str(strat.leaf_out_file)
                strat = tmp
                new_strat_list.append(strat)

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
    
    if 1 in c.runRegistrationPreprocessing:
        for strat in strat_list:
            reg_anat_mni = create_nonlinear_register('anat_mni_register_%d' % num_strat)
            
            try:
                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                 reg_anat_mni, 'inputspec.input_brain')
                node, out_file = strat.get_node_from_resource_pool('anatomical_reorient')
                workflow.connect(node, out_file,
                                 reg_anat_mni, 'inputspec.input_skull')
                
                reg_anat_mni.inputs.inputspec.reference_brain = c.standardResolutionBrain
                reg_anat_mni.inputs.inputspec.reference_skull = c.standard
            except:
                print 'Invalid Connection: Anatomical Registration:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise
            
            if 0 in c.runRegistrationPreprocessing:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                strat = tmp
                new_strat_list.append(strat)
            
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

    """
    Inserting Functional Input Data workflow
    """
    new_strat_list = []
    num_strat = 0

    if 1 in c.runFunctionalDataGathering:
        for strat in strat_list:
            # create a new node, Remember to change its name! 
            Flow = create_func_datasource(sub_dict['rest'])
            Flow.inputs.inputnode.subject = subject_id

            funcFlow = Flow.clone('func_gather_%d' % num_strat)

            if 0 in c.runFunctionalDataGathering:
                # we are forking so create a new node
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
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

    if 1 in c.runFunctionalPreprocessing:
        for strat in strat_list:
            # create a new node, Remember to change its name! 
            preproc = create_func_preproc()
            preproc.inputs.inputspec.start_idx = c.startIdx
            preproc.inputs.inputspec.stop_idx = c.stopIdx
            func_preproc = preproc.clone('func_preproc_%d' % num_strat)

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
                tmp.out_file = str(strat.leaf_out_file)
                strat = tmp
                new_strat_list.append(strat)

            strat.set_leaf_properties(func_preproc, 'outputspec.example_func')

            # add stuff to resource pool if we need it 

            strat.update_resource_pool({'mean_functional':(func_preproc, 'outputspec.example_func')})
            strat.update_resource_pool({'functional_preprocessed_mask':(func_preproc, 'outputspec.preprocessed_mask')})
            strat.update_resource_pool({'movement_parameters':(func_preproc, 'outputspec.movement_parameters')})
            strat.update_resource_pool({'max_displacement':(func_preproc, 'outputspec.max_displacement')})
            strat.update_resource_pool({'preprocessed':(func_preproc, 'outputspec.preprocessed')})
            strat.update_resource_pool({'functional_dilated_mask':(func_preproc, 'outputspec.mask')})

            num_strat += 1

    strat_list += new_strat_list

    """
    Inserting Functional to Anatomical Registration
    """
    new_strat_list = []
    num_strat = 0
    
    if 1 in c.runAnatomicalToFunctionalRegistration:
        for strat in strat_list:
            anat_to_func_reg = pe.Node(interface=fsl.FLIRT(),
                               name = 'anat_to_func_register_%d' % num_strat)
            anat_to_func_reg.inputs.cost = 'corratio'
            anat_to_func_reg.inputs.dof = 6
            
            func_gm = pe.Node(interface=fsl.ApplyXfm(),
                               name = 'func_gm_%d' % num_strat)
            func_gm.inputs.reference = c.standardResolutionBrain
            func_gm.inputs.apply_xfm = True

            func_csf = pe.Node(interface=fsl.ApplyXfm(),
                               name = 'func_csf_%d' % num_strat)
            func_csf.inputs.reference = c.standardResolutionBrain
            func_csf.inputs.apply_xfm = True
            
            func_wm = pe.Node(interface=fsl.ApplyXfm(),
                               name = 'func_wm_%d' % num_strat)
            func_wm.inputs.reference = c.standardResolutionBrain
            func_wm.inputs.apply_xfm = True
            
            try:
                node, out_file = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, out_file,
                                 anat_to_func_reg, 'in_file')

                node, out_file = strat.get_node_from_resource_pool('mean_functional')
                workflow.connect(node, out_file,
                                 anat_to_func_reg, 'reference')
                
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
                tmp.out_file = str(strat.leaf_out_file)
                strat = tmp
                new_strat_list.append(strat)
            
            strat.update_resource_pool({'anatomical_to_functional_xfm':(anat_to_func_reg, 'out_matrix_file'),
                                        'functional_gm_mask':(func_gm, 'out_file'),
                                        'functional_wm_mask':(func_wm, 'out_file'),
                                        'functional_csf_mask':(func_csf, 'out_file')})

            num_strat += 1
    
    strat_list += new_strat_list

    """
    Inserting Nuisance Workflow
    """
    new_strat_list = []
    num_strat = 0
    
    if 1 in c.runNuisance:
        for strat in strat_list:
            nuisance = create_nuisance('nuisance_%d' % num_strat)
            nuisance.get_node('residuals').iterables = ('selector', c.Corrections,
                                                        'compcor_ncomponents', c.nComponents)
            try:
                node, out_file = strat.get_node_from_resource_pool('functional_gm_mask')
                workflow.connect(node, out_file,
                                 nuisance, 'inputspec.gm_mask')
                
                node, out_file = strat.get_node_from_resource_pool('functional_wm_mask')
                workflow.connect(node, out_file,
                                 nuisance, 'inputspec.wm_mask')
                
                node, out_file = strat.get_node_from_resource_pool('functional_csf_mask')
                workflow.connect(node, out_file,
                                 nuisance, 'inputspec.csf_mask')
                
                node, out_file = strat.get_node_from_resource_pool('movement_parameters')
                workflow.connect(node, out_file,
                                 nuisance, 'inputspec.motion_components')
            except:
                print 'Invalid Connection: Nuisance:', num_strat, ' resource_pool: ', strat.get_resource_pool()
                raise

            if 0 in c.runNuisance:
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                strat = tmp
                new_strat_list.append(strat)
                
            strat.set_leaf_properties(nuisance, 'outputspec.subject')
            
            num_strat += 1
            
    strat_list += new_strat_list




    idx = 0
    for strat in strat_list:

        r_pool = strat.get_resource_pool()

        print '-----------------------------------------'

        for key in r_pool.keys():


            node, out_file = r_pool[key]

            print idx, ' ', key, ' ', node.name, ' ', out_file


        idx += 1


    workflow.write_graph(graph2use="orig")
#    workflow.run(plugin='MultiProc',
#                         plugin_args={'n_procs': c.numCoresPerSubject})





def run(config_file, subject_list_file):

    path, fname = os.path.split(os.path.realpath(config_file))
    sys.path.append(path)
    c = __import__(fname.split('.')[0])

    path, fname = os.path.split(os.path.realpath(subject_list_file))
    sys.path.append(path)
    s = __import__(fname.split('.')[0])

    sublist = s.subject_list

    flines = open(c.seedFile, 'r').readlines()
    seed_list = [fline.rstrip('\r\n') for fline in flines]

    if not c.runOnGrid:

        procss = [Process(target=prep_workflow, args=(sub_dict, seed_list, c)) for sub_dict in sublist]

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
        for sub_dict in sublist:
            prep_workflow(sub_dict, seed_list, c)

#    symlink_creator(c.sinkDirectory)

#    if c.derivatives[6]:
#        start_group_analysis(c)


if __name__ == "__main__":
    import os
    import commands
    cmd = '/bin/bash ~/.bashrc'
    print cmd
    sys.stderr.write(commands.getoutput(cmd))

    sys.exit(main())


