import numpy as np
from ..cwas import nifti_cwas

def test_nifti_cwas():

def test_cwas():
    import os, glob
    c = create_cwas()
    c.base_dir = os.getcwd()
    c.inputs.inputspec.roi = roi_file
    c.inputs.inputspec.subjects = subjects_list[:4]
    c.inputs.inputspec.regressor = np.arange(len(subjects_list[:4]))[:,np.newaxis]
    c.inputs.inputspec.parallel_nodes = 3
    c.run()