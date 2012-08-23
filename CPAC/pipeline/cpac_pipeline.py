import os
import sys
import copy
import argparse
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from multiprocessing import Process


from CPAC.anat_preproc.anat_preproc import create_anat_preproc
from CPAC.func_preproc.func_preproc import create_func_preproc
from CPAC.reg_preproc.reg_preproc import create_reg_preproc
from CPAC.seg_preproc.seg_preproc import create_seg_preproc
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



def getSubjectsandScansList(c, s):

    def get_list(fname):
        flines = open(fname, 'r').readlines()
        return [fline.rstrip('\r\n') for fline in flines]

    subjects_list = s.subjects_list

    if c.derivatives[1]:
        return subjects_list, get_list(c.seedFile)
    else:
        return subjects_list, []


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

    run_gather_anat = c.runAnatomicalDataGathering

    if not isinstance(run_gather_anat, list):

        run_gather_anat = [run_gather_anat]

    new_strat_list = []
    num_strat = 0


    strat_initial = None

    for gather_anat in run_gather_anat:


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

    run_anat = c.runAnatomicalPreprocessing

    if not isinstance(run_anat, list):

        run_anat = [run_anat]

    new_strat_list = []
    num_strat = 0

    tmp_run_anat = list(run_anat)

    for el in tmp_run_anat:

        if el == 0:
            tmp_run_anat.remove(el)



    if len(tmp_run_anat) > 0:


        for strat in strat_list:
            # create a new node, Remember to change its name! 
            anat_preproc = create_anat_preproc().clone('anat_preproc_%d' % num_strat)

            try:
                # connect the new node to the previous leaf
                node, out_file = strat.get_leaf_properties()
#                if('outputspec.anat' in out_file) and ('anat' in node.name):
                workflow.connect(node, out_file, anat_preproc, 'inputspec.anat')
#                else:
#                    print num_strat, ' cannot connect ', out_file, ' to Anat Preproc: inputspec.anat '

#                    continue


            except:

                print 'Invalid Connection: Anat Preprocessing No valid Previous for strat : ', num_strat, ' resource_pool: ', strat.get_resource_pool()
                num_strat += 1
                continue


            if 0 in run_anat:

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

            num_strat = num_strat+1

    strat_list = list(strat_list) +(list(new_strat_list))


    """
        Inserting Functional Input Data workflow

    """
    new_strat_list = []
    num_strat = 0

    gather_func = c.runFunctionalDataGathering

    if not isinstance(gather_func, list):

        gather_func = [gather_func]


    tmp_gather_func = list(gather_func)

    for el in tmp_gather_func:

        if el == 0:
            tmp_gather_func.remove(el)


    if len(tmp_gather_func) > 0:


        tmp_strat = None
        for strat in strat_list:


            # create a new node, Remember to change its name! 
            Flow = create_func_datasource(sub_dict['rest'])
            Flow.inputs.inputnode.subject = subject_id

            funcFlow = Flow.clone('func_gather_%d' % num_strat)

            if 0 in gather_func:

                # we are forking so create a new node
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                strat = tmp
                new_strat_list.append(strat)

            strat.set_leaf_properties(funcFlow, 'outputspec.rest')

            num_strat = num_strat+1

    strat_list = strat_list+new_strat_list


    """
        Inserting Functional Image Preprocessing
        Workflow
    """



    new_strat_list = []
    num_strat = 0

    run_func = c.runFunctionalPreprocessing

    if not isinstance(run_func, list):

        run_func = [run_func]

    tmp_run_func = list(run_func)

    for el in tmp_run_func:

        if el == 0:
            tmp_run_func.remove(el)



    if len(tmp_run_func) > 0:

        tmp_strat = None
        for strat in strat_list:

            # create a new node, Remember to change its name! 
            preproc = create_func_preproc()
            preproc.inputs.inputspec.start_idx = c.startIdx
            preproc.inputs.inputspec.stop_idx = c.stopIdx
            func_preproc = preproc.clone('func_preproc_%d' % num_strat)

            node = None
            out_file = None
            try:
                # connect the new node to the previous leaf
                node, out_file = strat.get_leaf_properties()

#                if('outputspec.rest' ==  out_file) and ('func_gather' in node.name):
                workflow.connect(node, out_file, func_preproc, 'inputspec.rest')
#                else:

#                    print num_strat, ' cannot connect ', out_file, ' to Func Preproc: inputspec.rest '
#
#                    num_strat += 1
#                    continue

            except:

                print 'Invalid Connection: Functional Preprocessing No valid Previous for strat : ', num_strat, ' resource_pool: ', strat.get_resource_pool()
                num_strat += 1
                continue

            if 0 in run_func:
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

            num_strat = num_strat+1

    strat_list = strat_list+new_strat_list


    """
        Inserting Registration Preprocessing
        Workflow
    """

    new_strat_list = []
    num_strat = 0

    run_reg = c.runRegistrationPreprocessing

    if not isinstance(run_reg, list):

        run_reg = [run_reg]

    tmp_run_reg = list(run_reg)

    for el in tmp_run_reg:

        if el == 0:
            tmp_run_reg.remove(el)


    if len(tmp_run_reg) > 0:

        for strat in strat_list:

            # create a new node, Remember to change its name! 
            preproc = create_reg_preproc()

            preproc.inputs.inputspec.standard_res_brain = \
                                                c.standardResolutionBrain
            preproc.inputs.inputspec.standard = c.standard
            preproc.inputs.inputspec.config_file = c.configFile
            preproc.inputs.inputspec.standard_brain_mask_dil = \
                                                c.standardBrainMaskDiluted

            reg_preproc = preproc.clone('reg_preproc_%d' % num_strat)

            node = None
            out_file = None
            try:

                # connect the new node to the previous leaf
                node, out_file = strat.get_leaf_properties()

                r_pool_keys = strat.get_resource_pool().keys()

#                if ('outputspec.example_func' == out_file) and ('anatomical_brain' in r_pool_keys) and ('anatomical_reorient' in r_pool_keys):
                workflow.connect(node, out_file, reg_preproc, 'inputspec.example_func')

                # connect from resource pool
                node, var = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, var, reg_preproc, 'inputspec.brain')

                # connect from resource pool
                node, var = strat.get_node_from_resource_pool('anatomical_reorient')
                workflow.connect(node, var, reg_preproc, 'inputspec.reorient')
#                else:
#                    num_strat += 1
#                    continue


            except:
                print 'Invalid Connection: Registration Preprocessing No valid Previous for strat : ', num_strat, ' resource_pool: ', strat.get_resource_pool()
                num_strat += 1
                continue


            if 0 in run_reg:
                # we are forking so create a new node
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                strat = tmp
                new_strat_list.append(strat)

            strat.set_leaf_properties(reg_preproc, 'outputspec.stand2highres_warp')

            # add stuff to resource pool if we need it 
            strat.update_resource_pool({'highres2example_func_mat':(reg_preproc, 'outputspec.highres2example_func_mat')})
            strat.update_resource_pool({'highres2standard_warp':(reg_preproc, 'outputspec.highres2standard_warp')})
            strat.update_resource_pool({'example_func2highres_mat':(reg_preproc, 'outputspec.example_func2highres_mat')})

            num_strat = num_strat+1

    strat_list = strat_list+new_strat_list

    """
        Inserting Segmentation Preprocessing
        Workflow
    """

    new_strat_list = []
    num_strat = 0

    run_seg = c.runSegmentationPreprocessing

    if not isinstance(run_seg, list):

        run_seg = [run_seg]

    tmp_run_seg = list(run_seg)

    for el in tmp_run_seg:

        if el == 0:
            tmp_run_seg.remove(el)


    if len(tmp_run_seg) > 0:

        for strat in strat_list:

            # create a new node, Remember to change its name! 
            preproc = create_seg_preproc()

            preproc.inputs.inputspec.standard_res_brain = c.standardResolutionBrain
            preproc.inputs.inputspec.PRIOR_CSF = c.PRIOR_CSF
            preproc.inputs.inputspec.PRIOR_WHITE = c.PRIOR_WHITE
            preproc.inputs.inputspec.PRIOR_GRAY = c.PRIOR_GRAY
            preproc.inputs.inputspec.standard_res_brain = c.standardResolutionBrain
            preproc.inputs.csf_threshold.csf_threshold = c.cerebralSpinalFluidThreshold
            preproc.inputs.wm_threshold.wm_threshold = c.whiteMatterThreshold
            preproc.inputs.gm_threshold.gm_threshold = c.grayMatterThreshold
            preproc.get_node('csf_threshold').iterables = ('csf_threshold', c.cerebralSpinalFluidThreshold)
            preproc.get_node('wm_threshold').iterables = ('wm_threshold', c.whiteMatterThreshold)
            preproc.get_node('gm_threshold').iterables = ('gm_threshold', c.grayMatterThreshold)
            seg_preproc = preproc.clone('seg_preproc_%d' % num_strat)


            node = None
            out_file = None
            try:
                # connect the new node to the previous leaf
                node, out_file = strat.get_leaf_properties()

                r_pool_keys = strat.get_resource_pool().keys()
#                if 'outputspec.stand2highres_warp' == out_file and ('highres2example_func_mat' in r_pool_keys) and ('functional_preprocessed_mask' in r_pool_keys) and ('anatomical_brain' in r_pool_keys) and ('mean_functional' in r_pool_keys):
                workflow.connect(node, out_file, seg_preproc, 'inputspec.stand2highres_warp')

                # connect from resource pool
                node, var = strat.get_node_from_resource_pool('highres2example_func_mat')
                workflow.connect(node, var, seg_preproc, 'inputspec.highres2example_func_mat')

                # connect from resource pool
                node, var = strat.get_node_from_resource_pool('functional_preprocessed_mask')
                workflow.connect(node, var, seg_preproc, 'inputspec.preprocessed_mask')

                # connect from resource pool
                node, var = strat.get_node_from_resource_pool('anatomical_brain')
                workflow.connect(node, var, seg_preproc, 'inputspec.brain')

                # connect from resource pool
                node, var = strat.get_node_from_resource_pool('mean_functional')
                workflow.connect(node, var, seg_preproc, 'inputspec.example_func')

#                else:
#                    num_strat += 1
#                    continue

            except:
                print 'Invalid Connection: Segmentation Preprocessing No valid Previous for strat : ', num_strat, ' resource_pool: ', strat.get_resource_pool()
                num_strat += 1
                continue


            if 0 in run_seg:
                # we are forking so create a new node
                tmp = strategy()
                tmp.resource_pool = dict(strat.resource_pool)
                tmp.leaf_node = (strat.leaf_node)
                tmp.out_file = str(strat.leaf_out_file)
                strat = tmp
                new_strat_list.append(strat)

            strat.set_leaf_properties(seg_preproc, 'outputspec.wm_mask')

            # add stuff to resource pool if we need it 

            strat.update_resource_pool({'wm_mask':(seg_preproc, 'outputspec.wm_mask')})
            strat.update_resource_pool({'csf_mask':(seg_preproc, 'outputspec.csf_mask')})
            strat.update_resource_pool({'gm_mask':(seg_preproc, 'outputspec.gm_mask')})
            strat.update_resource_pool({'probability_maps':(seg_preproc, 'outputspec.probability_maps')})
            strat.update_resource_pool({'mixeltype':(seg_preproc, 'outputspec.mixeltype')})
            strat.update_resource_pool({'partial_volume_map':(seg_preproc, 'outputspec.partial_volume_map')})
            strat.update_resource_pool({'partial_volume_files':(seg_preproc, 'outputspec.partial_volume_files')})

            num_strat = num_strat+1

    strat_list = strat_list+new_strat_list

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


    sublist, seed_list = getSubjectsandScansList(c, s)



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


