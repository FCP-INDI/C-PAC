import os
import numpy as np
import pandas as pd
import nilearn.datasets
import nibabel as nb

from CPAC.cwas.pipeline import create_cwas

def test_pipeline():

    try:
        os.mkdir('/tmp/cwas')
    except:
        pass

    abide_data = nilearn.datasets.fetch_abide_pcp(n_subjects=10)
    pheno = pd.DataFrame.from_records(abide_data['phenotypic'])
    images = abide_data['func_preproc']

    pheno = pheno[['FILE_ID', 'AGE_AT_SCAN', 'FIQ']]

    pheno.FILE_ID = pheno.FILE_ID.apply(lambda b: b.decode())

    pheno.to_csv('/tmp/cwas/pheno.csv')

    # Sanity ordering check
    assert all(
        FID in images[i]
        for i, FID in enumerate(pheno.FILE_ID)
    )

    cc = nilearn.datasets.fetch_atlas_craddock_2012()
    img = nb.load(cc['scorr_mean'])
    img_data = np.copy(img.get_data()[:,:,:,10])
    img_data[img_data != 2] = 0.0
    img = nb.Nifti1Image(img_data, img.affine)
    nb.save(img, '/tmp/cwas/roi.nii.gz')
    
    workflow = create_cwas('cwas', '/tmp/cwas', '/tmp/cwas')

    subjects = {
        FID: images[i]
        for i, FID in enumerate(pheno.FILE_ID)
    }

    roi = '/tmp/cwas/roi.nii.gz'
    regressor_file = '/tmp/cwas/pheno.csv'
    participant_column = 'FILE_ID'
    columns = 'AGE_AT_SCAN'
    permutations = 50
    parallel_nodes = 3

    workflow.inputs.inputspec.roi = roi
    workflow.inputs.inputspec.subjects = subjects
    workflow.inputs.inputspec.regressor = regressor_file
    workflow.inputs.inputspec.participant_column = participant_column
    workflow.inputs.inputspec.columns = columns
    workflow.inputs.inputspec.permutations = permutations
    workflow.inputs.inputspec.parallel_nodes = parallel_nodes

    workflow.run(plugin='Linear')
