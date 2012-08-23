import numpy as np
from ..cwas import nifti_cwas

def test_nifti_cwas():

def test_cwas():
    from CPAC.cwas import create_cwas
    import numpy as np
    import os, glob
    import time
    
    g_string = '/home/data/Projects/ADHD200/Manuel_sink/*/func/scan_id_rest_1/scan_id_anat_1/csf_threshold_0.4/gm_threshold_0.2/wm_threshold_0.66/nc_5/selector_0.2.3.7/bandpass_freqs_0.009.0.1/rest_residual_filtered_4mm.nii.gz'
    roi_file = '/usr/share/fsl/4.1/data/standard/MNI152_T1_4mm_brain_mask.nii.gz'
    subjects_list = glob.glob(g_string)
    
    c = create_cwas()
    c.inputs.inputspec.roi = roi_file
    c.inputs.inputspec.subjects = subjects_list[:5]
    c.inputs.inputspec.regressor = np.arange(len(subjects_list[:5]))[:,np.newaxis]
    c.inputs.inputspec.f_samples = 5000
    c.inputs.inputspec.parallel_nodes = 10
    c.base_dir = '/home/data/Projects/cwas_tests/results_fs%i_pn%i' % (c.inputs.inputspec.f_samples,
                                                                       c.inputs.inputspec.parallel_nodes)
    start = time.clock()
    c.run(plugin='MultiProc', plugin_args={'n_procs' : 10})
    end = time.clock()
    print "%.2gs" % (end-start)