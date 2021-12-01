"""`Generate eXtensible Connectivity Pipeline-style quality control files <https://fcp-indi.github.io/docs/user/xcpqc>`_

Columns
-------
sub : str
    subject label :cite:`cite-BIDS21`
ses : str
    session label :cite:`cite-BIDS21`
task : str
    task label :cite:`cite-BIDS21`
run : int
    run index :cite:`cite-BIDS21`
desc : str
    description :cite:`cite-BIDS21`
space : str
    space label :cite:`cite-BIDS21`
meanFD : float
    mean Jenkinson framewise displacement :cite:`cite-Jenk02` :func:`CPAC.generate_motion_statistics.calculate_FD_J`
relMeansRMSMotion : float
    "mean value of RMS motion" :cite:`cite-Ciri19`
relMaxRMSMotion : float
    "maximum vaue of RMS motion" :cite:`cite-Ciri19`
meanDVInit : float
    "mean DVARS" :cite:`cite-Ciri19`
meanDVFinal : float
    "mean DVARS" :cite:`cite-Ciri19`
nVolCensored : int
    "total number of volume(s) censored :cite:`cite-Ciri19`
nVolsRemoved : int
    number of volumes in derivative minus number of volumes in original
    functional scan
motionDVCorrInit : float
    "correlation of RMS and DVARS before regresion" :cite:`cite-Ciri19`
motionDVCorrFinal : float
    "correlation of RMS and DVARS after regresion" :cite:`cite-Ciri19`
coregDice : float
    "Coregsitration of Functional and T1w:[…] Dice index" :cite:`cite-Ciri19` :cite:`cite-Penn19`
coregJaccard : float
    "Coregsitration of Functional and T1w:[…] Jaccard index" :cite:`cite-Ciri19` :cite:`cite-Penn19`
coregCrossCorr : float
    "Coregsitration of Functional and T1w:[…] cross correlation" :cite:`cite-Ciri19` :cite:`cite-Penn19`
coregCoverag : float
    "Coregsitration of Functional and T1w:[…] Coverage index" :cite:`cite-Ciri19` :cite:`cite-Penn19`
normDice : float
    "Normalization of T1w/Functional to Template:[…] Dice index" :cite:`cite-Ciri19` :cite:`cite-Penn19`
normJaccard : float
    "Normalization of T1w/Functional to Template:[…] Jaccard index" :cite:`cite-Ciri19` :cite:`cite-Penn19`
normCrossCorr : float
    "Normalization of T1w/Functional to Template:[…] cross correlation" :cite:`cite-Ciri19` :cite:`cite-Penn19`
normCoverage : float
    "Normalization of T1w/Functional to Template:[…] Coverage index" :cite:`cite-Ciri19` :cite:`cite-Penn19`
"""  # noqa E501  # pylint: disable=line-too-long
import json
import os
import re
from io import BufferedReader

import nibabel as nb
import numpy as np
import pandas as pd
from CPAC.generate_motion_statistics.generate_motion_statistics import \
    motion_power_statistics
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.engine import NodeData
from CPAC.utils.interfaces.function import Function
from CPAC.utils.utils import check_prov_for_motion_tool


def dvcorr(dvars, fdj):
    """Function to correlate DVARS and FD-J"""
    dvars = np.loadtxt(dvars)
    fdj = np.loadtxt(fdj)
    if len(dvars) != len(fdj) - 1:
        raise ValueError(
            'len(DVARS) should be 1 less than len(FDJ), but their respective '
            f'lengths are {len(dvars)} and {len(fdj)}.'
        )
    return np.corrcoef(dvars, fdj[1:])[0, 1]


def strings_from_bids(final_func):
    """
    Function to gather BIDS entities into a dictionary

    Parameters
    ----------
    final_func : str

    Returns
    -------
    dict

    Examples
    --------
    >>> fake_path = (
    ...     '/path/to/sub-fakeSubject_ses-fakeSession_task-peer_run-3_'
    ...     'atlas-Schaefer400_space-MNI152NLin6_res-1x1x1_'
    ...     'desc-NilearnPearson_connectome.tsv')
    >>> strings_from_bids(fake_path)['desc']
    'NilearnPearson'
    >>> strings_from_bids(fake_path)['space']
    'MNI152NLin6'
    """
    from_bids = dict(
        tuple(entity.split('-', 1)) if '-' in entity else
        ('suffix', entity) for entity in final_func.split('/')[-1].split('_')
    )
    from_bids = {k: [from_bids[k]] for k in from_bids}
    if 'space' not in from_bids:
        from_bids['space'] = ['native']
    return from_bids


def generate_desc_qc(original_anat, final_anat, original_func, final_func,
                     n_vols_censored, space_T1w_bold, dvars_after=None,
                     fdj_after=None):
    # pylint: disable=too-many-arguments, too-many-locals, invalid-name
    """Function to generate an RBC-style QC CSV

    Parameters
    ----------
    original_anat : str
        path to original 'T1w' image

    final_anat : str
        path to 'desc-preproc_T1w' image

    original_func : str
        path to original 'bold' image

    final_bold : str
        path to 'desc-preproc_bold' image

    n_vols_censored : int

    space_T1w_bold : str
        path to 'space-T1w_desc-mean_bold' image

    dvars_after : str
        path to DVARS on final 'bold' image

    fdj_after : str
        path to framewise displacement (Jenkinson) on final 'bold' image

    Returns
    -------
    str
        path to XCP QC TSV
    """
    columns = (
        'sub,ses,task,run,desc,space,meanFD,relMeansRMSMotion,'
        'relMaxRMSMotion,meanDVInit,meanDVFinal,nVolCensored,nVolsRemoved,'
        'motionDVCorrInit,motionDVCorrFinal,coregDice,coregJaccard,'
        'coregCrossCorr,coregCoverage,normDice,normJaccard,normCrossCorr,'
        'normCoverage'.split(',')
    )

    images = {
        'original_anat': nb.load(original_anat),
        'original_func': nb.load(original_func),
        'final_anat': nb.load(final_anat),
        'final_func': nb.load(final_func),
        'space-T1w_bold': nb.load(space_T1w_bold)
    }

    # `sub` through `space`
    from_bids = strings_from_bids(final_func)

    # `nVolCensored` & `nVolsRemoved`
    shape_params = {'nVolCensored': n_vols_censored,
                    'nVolsRemoved': images['final_func'].shape[3] -
                    images['original_func'].shape[3]}

    # `meanFD (Jenkinson)`
    fdj = get_other_func(final_func, 'framewise-displacement-jenkinson.1D')
    if isinstance(final_func, BufferedReader):
        final_func = final_func.name
    qc_filepath = _generate_filename(final_func)

    desc_span = re.search(r'_desc-.*_', final_func)
    if desc_span:
        desc_span = desc_span.span()
        final_func = '_'.join([
            final_func[:desc_span[0]],
            final_func[desc_span[1]:]
        ])
    del desc_span
    power_params = {'meanFD': np.mean(np.loadtxt(
        get_other_func(final_func,
                       'framewise-displacement-jenkinson.1D')
    ))}

    # `relMeansRMSMotion` & `relMaxRMSMotion`
    mot = np.genfromtxt(get_other_func(final_func,
                                       'movement-parameters.1D')).T
    # Relative RMS of translation
    rms = np.sqrt(mot[3] ** 2 + mot[4] ** 2 + mot[5] ** 2)
    rms_params = {
        'relMeansRMSMotion': [np.mean(rms)],
        'relMaxRMSMotion': [np.max(rms)]
    }

    # `meanDVInit` & `meanDVFinal`
    dvars = get_other_func(final_func, 'dvars.1D')
    meanDV = {'meanDVInit': np.mean(np.loadtxt(dvars))}
    try:
        meanDV['motionDVCorrInit'] = dvcorr(dvars, fdj)
    except ValueError:
        meanDV['motionDVCorrInit'] = 'ValueError'
    if dvars_after:
        if not fdj_after:
            fdj_after = fdj
        dvars_after = np.loadtxt(dvars_after)
        fdj_after = np.loadtxt(fdj_after)
        meanDV['meanDVFinal'] = np.mean(dvars_after)
        try:
            meanDV['motionDVCorrFinal'] = dvcorr(dvars_after, fdj_after)
        except ValueError:
            meanDV['motionDVCorrFinal'] = 'ValueError'

    # Overlap
    overlap_images = {variable: image.get_fdata().ravel() for
                      variable, image in images.items() if
                      variable in ['space-T1w_bold', 'original_anat']}
    intersect = overlap_images['space-T1w_bold'] * overlap_images[
        'original_anat']
    vols = {variable: np.sum(image) for
            variable, image in overlap_images.items()}
    vol_intersect = np.sum(intersect)
    vol_sum = sum(vols.values())
    vol_union = vol_sum - vol_intersect
    overlap_params = {
        'coregDice': 2 * vol_intersect / vol_sum,
        'coregJaccard': vol_intersect / vol_union,
        'coregCrossCorr': np.corrcoef(
            overlap_images['space-T1w_bold'],
            overlap_images['original_anat'])[0, 1],
        'coregCoverage': vol_intersect / min(vols.values()),
        'normDice': 'N/A: native space',
        'normJaccard': 'N/A: native space',
        'normCrossCorr': 'N/A: native space',
        'normCoverage': 'N/A: native space'
    }

    qc_dict = {
        **from_bids,
        **power_params,
        **rms_params,
        **shape_params,
        **overlap_params,
        **meanDV
    }
    df = pd.DataFrame(qc_dict, columns=columns)
    df.to_csv(qc_filepath, sep='\t', index=False)
    return qc_filepath


def _generate_filename(final):
    """Function to generate an XCP QC filename

    Parameters
    ----------
    final : str
        filepath of input

    Returns
    -------
    str
        QC filepath
    """
    if '_desc-' in final:
        delimiter = re.search(r'_desc-.*_', final).group()
        desc_parts = final.split(delimiter)
        desc_parts = desc_parts[1][::-1].split('_')[-1][::-1] if (
                     '_' in desc_parts[1]) else None
        desc_string = '_'.join([part for part in [
            f'desc-{delimiter.split("-", 1)[-1].rstrip("_")}Xcpqc',
            desc_parts,
            'bold.tsv'
        ] if part])
        del delimiter, desc_parts
    else:
        desc_string = 'desc-xcpqc_bold.tsv'
    return os.path.join(os.getcwd(), '_'.join([
        *[part for part in final.split('_')[:-1] if 'desc' not in part],
        desc_string
    ]))


def get_other_func(final_func, suffix_dot_filetype):
    """
    Function to get another filepath given a final functional file
    and a suffix-dot-extension

    Parameters
    ----------
    final_func : str
        full path

    suffix_dot_filetype: str
        e.g., 'timeseries.1D'

    Returns
    -------
    str
        full path
    """
    return '_'.join([
        *final_func.split('_')[:(-2 if 'desc' in final_func else -1)],
        suffix_dot_filetype
    ])


def _get_motion_correct_tool(final_func):
    with open(get_other_func(final_func, 'movement-parameters.json'),
              'r') as mp_json:
        return check_prov_for_motion_tool(
            json.loads(mp_json.read()).get('CpacProvenance'))


def _get_motion_stats(name, motion_correct_tool, strat_pool, final):
    """
    Parameters
    ----------
    name : str

    motion_correct_tool : str

    strat_pool : CPAC.pipeline.engine.ResourePool

    final : CPAC.pipeline.engine.NodeData
        functional
    """
    pipe_num = name.split('_')[-1]
    wf = pe.Workflow(name=f'_{name}_afterNuisance')
    gen_motion_stats = motion_power_statistics(name, motion_correct_tool)

    sub = NodeData(strat_pool, 'subject')
    scan = NodeData(strat_pool, 'scan')

    other_func = Function(input_names=['final_func', 'suffix_dot_filetype'],
                          output_names=['filepath'],
                          function=get_other_func)

    after_funcs = {'mask': pe.Node(other_func,
                                   name=f'func_mask_after_{pipe_num}'),
                   'mp': pe.Node(other_func,
                                 name=f'func_mp_after_{pipe_num}'),
                   'max': pe.Node(other_func,
                                  name=f'func_max_after_{pipe_num}'),
                   'rels': pe.Node(other_func,
                                   name=f'func_rels_after_{pipe_num}'),
                   'tran': pe.Node(other_func,
                                   name=f'func_tran_after_{pipe_num}')}
    after_funcs['mask'].inputs.suffix_dot_filetype = \
        'space-bold_desc-brain_mask.nii.gz'
    after_funcs['mp'].inputs.suffix_dot_filetype = 'movement-parameters.1D'
    after_funcs['max'].inputs.suffix_dot_filetype = 'max-displacement.rms'
    after_funcs['rels'].inputs.suffix_dot_filetype = 'rels-displacement.rms'
    after_funcs['tran'].inputs.suffix_dot_filetype = 'transformations.rms'

    wf.connect([
        (sub.node, gen_motion_stats, [(sub.out, 'inputspec.subject_id')]),
        (scan.node, gen_motion_stats, [(scan.out, 'inputspec.scan_id')]),
        *[(final.node, after_func, [(final.out, 'final_func')])
          for _, after_func in after_funcs.items()],
        (after_funcs['mp'], gen_motion_stats, [
            ('filepath', 'inputspec.movement_parameters')]),
        (after_funcs['max'], gen_motion_stats, [
            ('filepath', 'inputspec.max_displacement')]),
        (after_funcs['rels'], gen_motion_stats, [
            ('filepath', 'inputspec.rels_displacement')]),
        (after_funcs['mask'], gen_motion_stats, [
            ('filepath', 'inputspec.mask')]),
        (after_funcs['tran'], gen_motion_stats, [
            ('filepath', 'inputspec.transformations')]),
    ])

    result = wf.run()

    return result.outputs.DVARS_1D, result.outputs.FDJ_1D


def qc_xcp(wf, cfg, strat_pool, pipe_num, opt=None):
    # pylint: disable=unused-argument, invalid-name
    """
    {'name': 'qc_xcp',
     'config': ['pipeline_setup', 'output_directory', 'quality_control'],
     'switch': ['generate_xcpqc_files'],
     'option_key': 'None',
     'option_val': 'None',
     'inputs': ['bold', 'subject', 'scan', 'desc-preproc_bold', 'desc-preproc_T1w', 'T1w',
                'space-T1w_desc-mean_bold'],
     'outputs': ['xcpqc']}
    """
    original = {}
    final = {}
    original['anat'] = NodeData(strat_pool, 'T1w')
    original['func'] = NodeData(strat_pool, 'bold')
    final['anat'] = NodeData(strat_pool, 'desc-preproc_T1w')
    final['func'] = NodeData(strat_pool, 'desc-preproc_bold')
    t1w_bold = NodeData(strat_pool, 'space-T1w_desc-mean_bold')

    motion_tool = pe.Node(Function(input_names=['final_func'],
                                   output_names=['motion_correct_tool'],
                                   function=_get_motion_correct_tool,
                                   as_module=True),
                          name=f'get_motion_correct_tool_{pipe_num}')

    gen_motion_stats = pe.Node(Function(input_names=['name',
                                                     'motion_correct_tool',
                                                     'strat_pool', 'final'],
                                        output_names=['DVARS_1D', 'FDJ_1D'],
                                        function=_get_motion_stats,
                                        as_module=True),
                               name=f'gen_motion_stats_wf_{pipe_num}')

    gen_motion_stats.inputs.name = f'gen_motion_stats_after_{pipe_num}'
    gen_motion_stats.inputs.strat_pool = strat_pool
    gen_motion_stats.inputs.final = final['func']

    qc_file = pe.Node(Function(input_names=['original_func', 'final_func',
                                            'original_anat', 'final_anat',
                                            'space_T1w_bold',
                                            'n_vols_censored', 'dvars_after',
                                            'fdj_after'],
                               output_names=['qc_file'],
                               function=generate_desc_qc,
                               as_module=True),
                      name=f'xcpqc_{pipe_num}')

    try:
        n_vols_censored = NodeData(strat_pool, 'n_vols_censored')
        wf.connect(n_vols_censored.node, n_vols_censored.out,
                   qc_file, 'n_vols_censored')
    except LookupError:
        qc_file.inputs.n_vols_censored = 'unknown'

    wf.connect([
        (original['anat'].node, qc_file, [
            (original['anat'].out, 'original_anat')]),
        (original['func'].node, qc_file, [
            (original['func'].out, 'original_func')]),
        (final['anat'].node, qc_file, [(final['anat'].out, 'final_anat')]),
        (final['func'].node, qc_file, [(final['func'].out, 'final_func')]),
        (final['func'].node, motion_tool, [
            (final['func'].out, 'final_func')
        ]),
        (t1w_bold.node, qc_file, [(t1w_bold.out, 'space_T1w_bold')]),
        (motion_tool, gen_motion_stats, [
            ('motion_correct_tool', 'motion_correct_tool')]),
        (final['func'].node, gen_motion_stats, [
            (final['func'].out, 'inputspec.motion_correct')]),
        (gen_motion_stats, qc_file, [('DVARS_1D', 'dvars_after'),
                                     ('FDJ_1D', 'fdj_after')])])

    outputs = {
        'xcpqc': (qc_file, 'qc_file'),
    }

    return (wf, outputs)
