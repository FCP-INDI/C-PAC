"""Nipype translation of ANTs workflows
------------------------------------

This functionality is adapted from nipreps/niworkflows:
https://www.nipreps.org/niworkflows/api/niworkflows.anat.ants.html
https://github.com/nipreps/niworkflows/blob/master/niworkflows/anat/ants.py
https://www.nipreps.org/niworkflows/
https://fmriprep.org/en/stable/
https://poldracklab.stanford.edu/
"""
import sys

from importlib.util import find_spec, module_from_spec
from CPAC.info import __version__


def get_template(in_template, raise_empty=False, label=None, suffix=None,
                 desc=None, **kwargs):
    '''Override templateflow.api.get to support plugging
    niworkflows-ants into C-PAC.

    Overrides https://www.templateflow.org/python-client/master/_modules/templateflow/api.html#get  # noqa: E501

    Parameters
    ----------
    in_template: dict

    raise_empty: boolean, optional
        Raise exception if no files were matched

    Returns
    -------
    path: str
    '''
    # tpl_mask_path
    if (label == 'brain' and suffix == 'probseg') or (
        desc == 'brain' and suffix == 'mask'
    ):
        return in_template['tpl_mask_path']

    # tpl_regmask_path
    if desc == 'BrainCerebellumExtraction' and suffix == 'mask':
        return in_template['tpl_regmask_path']


def get_template_specs(in_template, template_spec=None, default_resolution=1):
    '''Override niworkflows.utils.misc.get_template_specs to support
    plugging niworkflows-ants into C-PAC.

    Overrides https://github.com/nipreps/niworkflows/blob/df77431a26fa3f3e4b2aff790fa3b4b2c41dd650/niworkflows/utils/misc.py#L17-L77  # noqa: E501

    Parameters
    ----------
    in_template: dict

    template_spec: dict or None

    default_resolution: int, float, or None

    Returns
    -------
    tpl_target_path: str

    common_spec: dict
    '''
    return (in_template['tpl_target_path'], in_template)


def init_brain_extraction_wf(tpl_target_path,
                             tpl_mask_path,
                             tpl_regmask_path,
                             name='brain_extraction_wf',
                             template_spec=None,
                             use_float=True,
                             normalization_quality='precise',
                             omp_nthreads=None,
                             mem_gb=3.0,
                             bids_suffix='T1w',
                             atropos_refine=True,
                             atropos_use_random_seed=True,
                             atropos_model=None,
                             use_laplacian=True,
                             bspline_fitting_distance=200):
    '''We monkeypatch `init_brain_extraction_wf` from niworkflows.anat
    to use `tpl_target_path`, `tpl_mask_path` and `tpl_regmask_path`
    from pipeline config instead of from datalad templates.

    Parameters
    ----------
    tpl_target_path: str
        path to brain extraction template

    tpl_mask_path: str
        path to probabilistic brain mask

    tpl_regmask_path: str
        path to registration mask

    # TODO: import other paramters from niworkflows

    updated

    Returns
    -------
    wf: Workflow
    '''
    # monkeypatch `get_template_specs` and `get_template` to use paths
    # specified in pipeline config rather than datalad templates
    mp_niworkflows_ants = _monkeypatch_in_template('niworkflows.anat.ants')

    in_template = {
        'tpl_target_path': tpl_target_path,
        'tpl_mask_path': tpl_mask_path,
        'tpl_regmask_path': tpl_regmask_path
    }

    return mp_niworkflows_ants.init_brain_extraction_wf(
        in_template=in_template, name=name, template_spec=template_spec,
        use_float=use_float, normalization_quality=normalization_quality,
        omp_nthreads=omp_nthreads, mem_gb=mem_gb, bids_suffix=bids_suffix,
        atropos_refine=atropos_refine,
        atropos_use_random_seed=atropos_use_random_seed,
        atropos_model=atropos_model, use_laplacian=use_laplacian,
        bspline_fitting_distance=bspline_fitting_distance)


def _monkeypatch_in_template(module_name):
    '''Function to monkeypatch `in_template` to work with a dictionary
    of local paths to templates rather than a string name of a
    templateflow template.

    Parameters
    ----------
    module_name: str
        name of module to monkeypatch

    Returns
    -------
    module: module
        monkeypatched module
    '''
    spec = find_spec(module_name)
    patched_source = '\n'.join([
        f'# MONKEYPATCHED for C-PAC {__version__}',
        spec.loader.get_source('niworkflows.anat.ants').replace(
            'from ..utils.misc import get_template_specs',
            'from CPAC.anat_preproc.ants import get_template_specs').replace(
            'from templateflow.api import get as get_template',
            'from CPAC.anat_preproc.ants import get_template'
         ).replace(
            'in_template : str\n        Name of the skull-stripping template '
            '(\'OASIS30ANTs\', \'NKI\', or\n        path).\n        The brain '
            'template from which regions will be projected\n        '
            'Anatomical template created using e.g. LPBA40 data set with\n'
            '        ``buildtemplateparallel.sh`` in ANTs.\n        The '
            'workflow will automatically search for a brain probability\n'
            '        mask created using e.g. LPBA40 data set which have brain '
            'masks\n        defined, and warped to anatomical template '
            'and\n        averaged resulting in a probability image.',
            'in_template : dict\n        Local paths to templates, keyed with '
            '{\n        \'tpl_target_path\', \'tpl_mask_path\', '
            '\'tpl_regmask_path\'\n        }. See CPAC.anat_preproc.ants for '
            'details.'
        )
    ])
    module = module_from_spec(spec)
    patched_code = compile(patched_source, module.__spec__.origin, 'exec')
    exec(patched_code, module.__dict__)
    sys.modules[module_name] = module
    return module
