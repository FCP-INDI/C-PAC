# test/unit/network_centrality/resting_state_centrality_test.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module performs testing on the functions in
CPAC/network_centrality/resting_state_centrality.py
'''

# Get centrality workflow parameters
def get_wflow_params(reg_type):
    '''
    Function to get basic running parameters for the centrality
    workflows

    Parameters
    ----------
    reg_type : string
        registration sub-folder to set as the base directory of the
        centrality workflows

    Returns
    -------
    fwhm : float
        the FWHM for the Gaussian smoothing kernel
    func_mni : string
        filepath to the resampled functional mni input to workflow
    mask_template : string
        filepath to the centrality mask template
    mem_gb : float
        the amount of memory (GB) allocated to computing centrality
    base_dir : string
        filepath to the base directory where centrality wflow working
        directories will sit
    '''

    # Import packages
    import os
    import yaml
    from CPAC.network_centrality import resting_state_centrality
    from CPAC.utils import test_init

    # Init variables
    config_path = test_init.populate_template_config('pipline_config')
    pipeline_config = yaml.load(open(config_path, 'r'))
    test_wflows = {}

    # Get workflow configuration parameters
    fwhm = pipeline_config['fwhm'][0]
    mask_txt_path = pipeline_config['templateSpecificationFile']
    memory_gb = pipeline_config['memoryAllocatedForDegreeCentrality']

    # Get centrality mask path from txt file
    mask_template = []
    with open(mask_txt_path, 'r') as txt_f:
        mask_template.extend([line for line in txt_f])
    mask_template = mask_template[0].rstrip('\n')

    # Get precomputed centrality files directory
    base_dirs = test_init.return_subj_measure_dirs('network_centrality')
    # Only grab registration strategy of interest
    base_dir = [reg_dir for reg_dir in base_dirs \
                if reg_dir.split('/')[-4] == reg_type][0]

    # Grab functional mni as subject input for that strategy
    func_mni_dir = base_dir.replace('network_centrality', 'functional_mni')
    func_mni = os.path.join(func_mni_dir, 'functional_mni_centrality.nii.gz')

    # Return paramters
    return fwhm, func_mni, mask_template, mem_gb, base_dir


# Set up workflow
def init_workflow(func_mni, mask_template, mem_gb):
    '''
    Function which inits a nipype workflow using network_centrality's
    create_resting_state_graphs() function to be used for testing

    Parameters
    ----------
    self : CentralityWorkflowTestCase
        a unittest.TestCase-inherited class

    Returns
    -------
    None
        this function does not return any values, but populates the
        instance attributes for:

        self.mask_template : string
            filepath to the centrality mask template
        self.output_dirs : list
            a list of filepath strings of output base directories
        self.test_wflows : list
            a list of nipype.pipeline.engine.Workflow objects
    '''


    # Init workflows
    wflow = resting_state_centrality.\
            create_resting_state_graphs(allocated_memory=memory_gb)


    # Set up workflow parameters
    wflow.base_dir = out_dir.replace('output', 'tests')
    wflow.inputs.inputspec.subject = func_mni
    wflow.inputs.inputspec.template = mask_template

    # Make the key the strategy being used (last folder)
    wflow_strat = out_dir.split('/')[-1]
    test_wflows[wflow_strat] = wflow

    # Set centrality TestCase instance attributes
    self.mask_template = mask_template
    self.output_dirs = output_dirs
    self.test_wflows = test_wflows


# Test the ants registration strategy
def run_p_value_thresh(self):
    '''
    Function to run the centrality workflows for the ANTS
    registration strategy
    '''

    # Import packages
    import os
    from CPAC.utils import test_init

    # Init variables
    ants_wflow = self.test_wflows['ants']
    smooth_dict = {}

    # Set up workflows and run each
    ants_wflow.base_dir = os.path.join(ants_wflow.base_dir, 'deg')
    ants_wflow.inputs.inputspec.method_option = 0
    ants_wflow.inputs.inputspec.weight_options = [True, True]
    ants_wflow.inputs.inputspec.threshold_option = 0
    ants_wflow.inputs.inputspec.threshold = 0.001
    print 'running degree centrality...'
    ants_wflow.run()

    # Set up workflows and run each
    ants_wflow.base_dir = ants_wflow.base_dir.replace('deg', 'eig')
    ants_wflow.inputs.inputspec.method_option = 1
    ants_wflow.inputs.inputspec.weight_options = [True, True]
    ants_wflow.inputs.inputspec.threshold_option = 0
    ants_wflow.inputs.inputspec.threshold = 0.001
    print 'running eigenvector centrality...'
    ants_wflow.run()

    # Set up workflows and run each
    ants_wflow.base_dir = ants_wflow.base_dir.replace('eig', 'lfcd')
    ants_wflow.inputs.inputspec.method_option = 0
    ants_wflow.inputs.inputspec.weight_options = [True, False]
    ants_wflow.inputs.inputspec.threshold_option = 2
    ants_wflow.inputs.inputspec.threshold = 0.6
    print 'running lfcd...'
    ants_wflow.run()

# Collect test outputs and compare
def test_collect_and_compare(self):
    '''
    Function to collect the precomputed and test outputs and
    compare the images
    '''

    # Import packages
    import glob
    import os
    import nibabel as nb
    import numpy as np

    # Init variables
    outputs_to_test = {}

    # Grab precomputed outputs and corresponding test outputs
    # For each (strategy) precomputed output directory
    for out_dir in self.output_dirs:
        niis = glob.glob(os.path.join(out_dir, '*.nii.gz'))
        test_dir = out_dir.replace('output', 'tests')
        # For each precomputed output nifti
        for nii in niis:
            nii_file = os.path.basename(nii)
            f_list = []
            for root, dirs, files in os.walk(test_dir):
                if files:
                    f_list.extend([os.path.join(root, file) for file in files \
                              if file == nii_file])
            if len(f_list) > 1:
                err_msg = 'More than one file found for %s in %s; '\
                          'please use only one' % (nii_file, str(f_list))
                raise Exception(err_msg)
            elif len(f_list) == 0:
                print 'No test outputs found for %s, skipping comparison' \
                      % nii_file
            else:
                strat = os.path.basename(out_dir)
                if not outputs_to_test.has_key(strat):
                    outputs_to_test[strat] = {nii : f_list[0]}
                else:
                    outputs_to_test[strat][nii] = f_list[0]

    # Iterate through dictionary and assert correlations\
    pass_thr = 0.98
    err_msg = 'Test failed: correlation < %.3f' % pass_thr
    for strat, golden_vs_test in outputs_to_test.items():
        for golden, test in golden_vs_test.items():
            # Load in golden and test images for comparison
            img1 = nb.load(golden).get_data()
            img2 = nb.load(test).get_data()

            # Compute pearson correlation on flattened 4D images
            print 'Comparing %s outputs...' % \
                  (os.path.basename(test.rstrip('.nii.gz')))
            corr = np.corrcoef(img1.flatten(), img2.flatten())[0,1]
            print 'Correlation = %.3f' % corr

            # Assert the correlation is >= pass_threshold
            self.assertGreaterEqual(corr, pass_thr, err_msg)


# Run and test centrality
def run_and_test_centrality(reg_type):
    '''
    Function to init, run, and test the outputs of the network
    centrality workflow

    Parameters
    ----------
    reg_type : string
        the type of registration used for the functional image (in MNI)
    '''

    # Import packages

    # Init variables
    fwhm, func_mni, mask_template, mem_gb, base_dir = \
            get_wflow_params(reg_type)

    # Initialize common workflow
    common_wflow = init_wflow(func_mni, mask_template, mem_gb)

    # 


# Command-line run-able unittest module
if __name__ == '__main__':

    # Init variables
    registration = 'ants'

    # Run and test centrality
    run_and_test_centrality(registration)
