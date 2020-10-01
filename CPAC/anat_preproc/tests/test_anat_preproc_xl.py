import os
import time
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from CPAC.utils import Configuration, Strategy
import CPAC.utils.test_init as test_utils
from CPAC.utils.interfaces.function import Function
from CPAC.utils.datasource import resolve_resolution
from CPAC.anat_preproc.anat_preproc import create_anat_preproc

def file_node(path, file_node_num=0):
    '''
    Create an identity node given a file path

    Parameters
    ----------
    path : string
        file path
    file_node_num : int
        number of file node
    
    Returns
    -------
    input_node : identity node
        an identity node of the input file path
    '''

    input_node = pe.Node(
        util.IdentityInterface(fields=['file']), name='file_node_{0}'.format(file_node_num)
    )
    input_node.inputs.file = path

    return input_node, 'file'


def test_anat_preproc_afni(working_path, input_path, benchmark_path=None, test_wf_name='test_anat_preproc_afni'):
    '''
    Test create_anat_preproc() with AFNI

    Parameters
    ----------
    working_path : string
        nipype working directory
    input_path : string
        input file path
    benchmark_path : string
        benchmark file path
    test_wf_name : string
        name of test workflow
    
    Returns
    -------
    None
    '''

    # create a configuration object
    config = Configuration({
        "num_ants_threads": 4,
        "workingDirectory": os.path.join(working_path, "working"),
        "crashLogDirectory": os.path.join(working_path, "crash"),
        "outputDirectory": working_path,
        "non_local_means_filtering": False,
        "n4_bias_field_correction": False,
        "skullstrip_mask_vol": False,
        "skullstrip_shrink_factor": 0.6,
        "skullstrip_var_shrink_fac": True,
        "skullstrip_shrink_factor_bot_lim": 0.4,
        "skullstrip_avoid_vent": True,
        "skullstrip_n_iterations": 250,
        "skullstrip_pushout": True,
        "skullstrip_touchup": True,
        "skullstrip_fill_hole": 10,
        "skullstrip_NN_smooth": 72,
        "skullstrip_smooth_final": 20,
        "skullstrip_avoid_eyes": True,
        "skullstrip_use_edge": True,
        "skullstrip_exp_frac": 0.1,
        "skullstrip_push_to_edge": False,
        "skullstrip_use_skull": False,
        "skullstrip_perc_int": 0,
        "skullstrip_max_inter_iter": 4,
        "skullstrip_fac": 1,
        "skullstrip_blur_fwhm": 0,
        "skullstrip_monkey": False,
    })

    # mock the strategy
    strat = Strategy()

    resource_dict = {
        "anatomical": input_path
    }

    file_node_num = 0
    for resource, filepath in resource_dict.items():
        strat.update_resource_pool({
            resource: file_node(filepath, file_node_num)
        })
        strat.append_name(resource+'_0')
        file_node_num += 1

    # build the workflow
    workflow = pe.Workflow(name=test_wf_name)
    workflow.base_dir = config.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(config.crashLogDirectory)
    }

    # call create_anat_preproc
    anat_preproc = create_anat_preproc(method='afni', 
                        already_skullstripped=False, 
                        config=config,
                        wf_name='anat_preproc', 
                        sub_dir=None)

    # connect AFNI options
    anat_preproc.inputs.AFNI_options.set(
        mask_vol=config.skullstrip_mask_vol,
        shrink_factor=config.skullstrip_shrink_factor,
        var_shrink_fac=config.skullstrip_var_shrink_fac,
        shrink_fac_bot_lim=config.skullstrip_shrink_factor_bot_lim,
        avoid_vent=config.skullstrip_avoid_vent,
        niter=config.skullstrip_n_iterations,
        pushout=config.skullstrip_pushout,
        touchup=config.skullstrip_touchup,
        fill_hole=config.skullstrip_fill_hole,
        avoid_eyes=config.skullstrip_avoid_eyes,
        use_edge=config.skullstrip_use_edge,
        exp_frac=config.skullstrip_exp_frac,
        smooth_final=config.skullstrip_smooth_final,
        push_to_edge=config.skullstrip_push_to_edge,
        use_skull=config.skullstrip_use_skull,
        perc_int=config.skullstrip_perc_int,
        max_inter_iter=config.skullstrip_max_inter_iter,
        blur_fwhm=config.skullstrip_blur_fwhm,
        fac=config.skullstrip_fac,
        monkey=config.skullstrip_monkey
    )

    node, out_file = strat['anatomical']
    workflow.connect(node, out_file,
        anat_preproc, 'inputspec.anat')

    # run workflow
    workflow.run()

    if benchmark_path is not None:
        out_path = os.path.join(config.workingDirectory, test_wf_name, 'anat_preproc', 'anat_skullstrip',
            input_path[input_path.rindex('/')+1:input_path.rindex('.nii.gz')]+'_resample_skullstrip.nii.gz')

        # calculate corretion between function output and benchmark output
        corr = test_utils.pearson_correlation(out_path, benchmark_path)
        print(f'\nCorrelation = {round(corr,3)}\n')

        assert(corr > .99)



def test_anat_preproc_fsl(working_path, input_path, benchmark_path=None, test_wf_name='test_anat_preproc_fsl'):
    '''
    Test create_anat_preproc() with FSL

    Parameters
    ----------
    working_path : string
        nipype working directory
    input_path : string
        input file path
    benchmark_path : string
        benchmark file path
    test_wf_name : string
        name of test workflow
    
    Returns
    -------
    None
    '''

    config = Configuration({
        "num_ants_threads": 4,
        "workingDirectory": os.path.join(working_path, "working"),
        "crashLogDirectory": os.path.join(working_path, "crash"),
        "outputDirectory": working_path,
        "non_local_means_filtering": False,
        "n4_bias_field_correction": False,
        "bet_frac":  0.5,
        "bet_mask_boolean": True, 
        "bet_mesh_boolean" : False, 
        "bet_outline" : False, 
        "bet_padding" : False, 
        "bet_radius" : 0,
        "bet_reduce_bias" : False, 
        "bet_remove_eyes" : False, 
        "bet_robust" : False,
        "bet_skull" : False,
        "bet_surfaces" : False,
        "bet_threshold" : False,
        "bet_vertical_gradient" : 0.0,
    })

    strat = Strategy()

    resource_dict = {
        "anatomical": input_path
    }

    file_node_num = 0
    for resource, filepath in resource_dict.items():
        strat.update_resource_pool({
            resource: file_node(filepath, file_node_num)
        })
        strat.append_name(resource+'_0')
        file_node_num += 1     

    workflow = pe.Workflow(name=test_wf_name)
    workflow.base_dir = config.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(config.crashLogDirectory)
    }

    anat_preproc = create_anat_preproc(method='fsl', 
                        already_skullstripped=False, 
                        config=config,
                        wf_name='anat_preproc', 
                        sub_dir=None)

    anat_preproc.inputs.BET_options.set(
        frac=config.bet_frac,
        mask_boolean=config.bet_mask_boolean,
        mesh_boolean=config.bet_mesh_boolean,
        outline=config.bet_outline,
        padding=config.bet_padding,
        radius=config.bet_radius,
        reduce_bias=config.bet_reduce_bias,
        remove_eyes=config.bet_remove_eyes,
        robust=config.bet_robust,
        skull=config.bet_skull,
        surfaces=config.bet_surfaces,
        threshold=config.bet_threshold,
        vertical_gradient=config.bet_vertical_gradient,
    )

    node, out_file = strat['anatomical']
    workflow.connect(node, out_file,
        anat_preproc, 'inputspec.anat')

    workflow.run()

    out_path = os.path.join(config.workingDirectory, test_wf_name, 'anat_preproc', 'anat_skullstrip',
        input_path[input_path.rindex('/')+1:input_path.rindex('.nii.gz')]+'_resample_brain.nii.gz')

    if benchmark_path is not None:

        corr = test_utils.pearson_correlation(out_path, benchmark_path)
        print(f'\nCorrelation = {round(corr,3)}\n')

        assert(corr > .99)


if __name__ == '__main__':

    working_path = '/Users/xinhui.li/C-PAC/cpac_runs/test_anat_preproc'
    input_path = '/Users/xinhui.li/C-PAC/cpac_runs/test_anat_preproc/sub-0025427/ses-1/anat/sub-0025427_ses-1_run-1_T1w.nii.gz'

    # test AFNI
    benchmark_path = '/Users/xinhui.li/C-PAC/cpac_runs/test_anat_preproc/afni_output/sub-0025427_ses-1_run-1_T1w_resample_afni_skullstrip.nii.gz'
    start_time = time.time()
    test_anat_preproc_afni(working_path, input_path, benchmark_path, test_wf_name='test_anat_preproc_afni')
    end_time = time.time()
    print(f'\nRun time = {round(end_time - start_time, 3)} seconds\n')

    # test AFNI using low resolution data (10mm, 11mm, 12mm, 15mm, 20mm)
    # start_time = time.time()
    # test_anat_preproc_afni(working_path, input_path, benchmark_path=None, test_wf_name='test_anat_preproc_afni_12mm')
    # end_time = time.time()
    # print(f'\nRun time = {round(end_time - start_time, 3)} seconds\n')

    # test FSL
    # benchmark_path = '/Users/xinhui.li/C-PAC/cpac_runs/test_anat_preproc/fsl_output/sub-0025427_ses-1_run-1_T1w_resample_fsl_skullstrip.nii.gz'
    # test_anat_preproc_fsl(working_path, input_path, benchmark_path, test_wf_name='test_anat_preproc_fsl')

    '''
    Run Time Comparison

    3.4mm: 87.406 seconds
    10mm: 56.596 seconds
    11mm: 73.643 seconds
    12mm or lower resolution: ERROR! AFNI requires # of slices >= 16

    build the graph/don't run the workflow: 0.134 seconds
    '''