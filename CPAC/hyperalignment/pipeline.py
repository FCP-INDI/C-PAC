import nibabel
import numpy as np

from mvpa2.base.hdf5 import h5save, h5load
from mvpa2.algorithms.searchlight_hyperalignment import SearchlightHyperalignment
from mvpa2.datasets.mri import fmri_dataset
from mvpa2.mappers.zscore import zscore


def prepare_subject_for_hyperalignment(subject_label, bold_fname, mask_fname, out_dir):
    ds = fmri_dataset(samples=bold_fname, mask=mask_fname)
    zscore(ds, chunks_attr=None)
    out_fname = os.path.join(out_dir, 'sub-%s_data.hdf5' % subject_label)
    h5save(out_fname, ds)


def run_hyperalignment(subjects_to_analyze, out_dir):
    ds_all = [
        h5load('%s/sub-%s_data.hdf5' % (out_dir, subject_label))
        for subject_label in subjects_to_analyze
    ]
    slhyper = SearchlightHyperalignment(radius=2,
                                        nblocks=10,
                                        sparse_radius=5,
                                        dtype='float16')
    hmappers = slhyper(ds_all)
    return hmappers


def create_hyperalignment(name='hyperalignment_wf', working_dir=None, crash_dir=None):
    """
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow.
        
    Returns
    -------
    workflow : nipype.pipeline.engine.Workflow
        
    Notes
    -----
    
    Workflow Inputs::
    
        
    Workflow Outputs::

    
    References
    ----------
    """

    if not working_dir:
        working_dir = os.path.join(os.getcwd(), 'Hyperalignment_work_dir')
    if not crash_dir:
        crash_dir = os.path.join(os.getcwd(), 'Hyperalignment_crash_dir')

    wf = pe.Workflow(name=name)
    wf.base_dir = working_dir
    wf.config['execution'] = {'hash_method': 'timestamp',
                              'crashdump_dir': os.path.abspath(crash_dir)}

    inputspec = pe.Node(
        util.IdentityInterface(fields=[
            'subjects',
            'radius',
        ]),
        name='inputspec'
    )

    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'correlations',
            'significance'
        ]),
        name='outputspec'
    )

    hyper = Hyperalignment()
    hypmaps = hyper(ds_train_fs)




        ds = fmri_dataset(samples=bold_fname, mask=mask_fname)
    zscore(ds, chunks_attr=None)
    out_fname = os.path.join(out_dir, 'sub-%s_data.hdf5' % subject_label)
    print('Saving to %s' % out_fname)
    h5save(out_fname, ds)