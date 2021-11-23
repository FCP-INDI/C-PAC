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


def _from_bids(final_func):
    from_bids = dict(
        tuple(entity.split('-', 1)) if '-' in entity else
        ('suffix', entity) for entity in final_func.split('/')[-1].split('_')
    )
    from_bids = {k: [from_bids[k]] for k in from_bids}
    if 'space' not in from_bids:
        from_bids['space'] = ['native']
    return from_bids


def _from_bids_node(final_func):
    from_bids = _from_bids(final_func)
    return from_bids['sub'], from_bids['run'], from_bids['func']


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
    from_bids = _from_bids(final_func)

    # `nVolCensored` & `nVolsRemoved`
    shape_params = {'nVolCensored': n_vols_censored,
                    'nVolsRemoved': images['final_func'].shape[3] -
                    images['original_func'].shape[3]}

    # `meanFD (Jenkinson)`
    fdj = _get_other_func(final_func, 'framewise-displacement-jenkinson.1D')
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
        _get_other_func(final_func,
                        'framewise-displacement-jenkinson.1D')
    ))}

    # `relMeansRMSMotion` & `relMaxRMSMotion`
    mot = np.genfromtxt(_get_other_func(final_func,
                                        'movement-parameters.1D')).T
    # Relative RMS of translation
    rms = np.sqrt(mot[3] ** 2 + mot[4] ** 2 + mot[5] ** 2)
    rms_params = {
        'relMeansRMSMotion': [np.mean(rms)],
        'relMaxRMSMotion': [np.max(rms)]
    }

    # `meanDVInit` & `meanDVFinal`
    dvars = _get_other_func(final_func, 'dvars.1D')
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
            f'desc-{delimiter.split("-", 1)[-1].rstrip("_")}+xcpqc',
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


def _get_other_func(final_func, suffix_dot_filetype):
    return '_'.join([
        *final_func.split('_')[:(-2 if 'desc' in final_func else -1)],
        suffix_dot_filetype
    ])


def qc_xcp(wf, cfg, strat_pool, pipe_num, opt=None):
    # pylint: disable=unused-argument, invalid-name
    """
    {'name': 'qc_xcp',
     'config': ['pipeline_setup', 'output_directory', 'quality_control'],
     'switch': ['generate_xcpqc_files'],
     'option_key': 'None',
     'option_val': 'None',
     'inputs': ['bold', 'desc-preproc_bold', 'desc-preproc_T1w', 'T1w',
                'space-T1w_desc-mean_bold'],
     'outputs': ['xcpqc']}
    """
    original = {'anat': {}, 'func': {}}
    final = {'anat': {}, 'func': {}}
    t1w_bold = {}
    original['anat']['node'], original['anat']['out'] = strat_pool.get_data(
        'T1w')
    original['func']['node'], original['func']['out'] = strat_pool.get_data(
        'bold')
    final['anat']['node'], final['anat']['out'] = strat_pool.get_data(
        'desc-preproc_T1w')
    final['func']['node'], final['func']['out'] = strat_pool.get_data(
        'desc-preproc_bold')
    t1w_bold['node'], t1w_bold['out'] = strat_pool.get_data(
        'space-T1w_desc-mean_bold')

    with open(_get_other_func(final['func'], 'movement-parameters.json'),
              'r') as mp_json:
        motion_correct_tool = check_prov_for_motion_tool(
            json.loads(mp_json.read()).get('CpacProvenance'))

    from_bids = pe.Node(Function(input_names=['final_func'],
                                 output_names=['sub', 'run', 'func'],
                                 function=_from_bids_node),
                        name=f'motion_stats_from_bids_{pipe_num}')
    gen_motion_stats = motion_power_statistics(
        name=f'gen_motion_stats_after_{pipe_num}',
        motion_correct_tool=motion_correct_tool)

    gen_motion_stats.inputs.inputspec.mask = _get_other_func(
        final['func'], 'space-bold_desc-brain_mask.nii.gz')
    gen_motion_stats.inputs.inputspec.movement_parameters = _get_other_func(
        final['func'], 'movement-parameters.1D')
    gen_motion_stats.inputs.inputspec.max_displacement = _get_other_func(
        final['func'], 'max-displacement.rms')
    gen_motion_stats.inputs.inputspec.rels_displacement = _get_other_func(
        final['func'], 'rels-displacement.rms')
    gen_motion_stats.inputs.inputspec.transformations = _get_other_func(
        final['func'], 'transformations.rms')

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
        n_vols_censored = strat_pool.get_data('n_vols_censored')
        wf.connect(n_vols_censored[0], n_vols_censored[1],
                   qc_file, 'n_vols_censored')
    except LookupError:
        qc_file.inputs.n_vols_censored = 'unknown'

    wf.connect([
        (original['anat']['node'], qc_file, [
            (original['anat']['out'], 'original_anat')]),
        (original['func']['node'], qc_file, [
            (original['func']['out'], 'original_func')]),
        (final['anat']['node'], qc_file, [
            (final['anat']['out'], 'final_anat')]),
        (final['func']['node'], qc_file, [
            (final['func']['out'], 'final_func')]),
        (t1w_bold['node'], qc_file, [(t1w_bold['out'], 'space_T1w_bold')]),
        (from_bids, gen_motion_stats, [('sub', 'inputspec.subject_id'),
                                       ('run', 'inputspec.scan'),
                                       ('func', 'inputspec.motion_correct')]),
        (gen_motion_stats, qc_file, [('DVARS_1D', 'dvars_after'),
                                     ('FDJ_1D', 'fdj_after')])])

    outputs = {
        'xcpqc': (qc_file, 'qc_file'),
    }

    return (wf, outputs)
