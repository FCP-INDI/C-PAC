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
    "Coregsitration of Functional and T1w:[…] Dice index" :cite:`cite-Ciri19`
coregJaccard : float
    "Coregsitration of Functional and T1w:[…] Jaccard index" :cite:`cite-Ciri19`
coregCrossCorr : float
    "Coregsitration of Functional and T1w:[…] cross correlation" :cite:`cite-Ciri19`
coregCoverag : float
    "Coregsitration of Functional and T1w:[…] Coverage index" :cite:`cite-Ciri19`
normDice : float
    "Normalization of T1w/Functional to Template:[…] Dice index" :cite:`cite-Ciri19`
normJaccard : float
    "Normalization of T1w/Functional to Template:[…] Jaccard index" :cite:`cite-Ciri19`
normCrossCorr : float
    "Normalization of T1w/Functional to Template:[…] cross correlation" :cite:`cite-Ciri19`
normCoverage : float
    "Normalization of T1w/Functional to Template:[…] Coverage index" :cite:`cite-Ciri19`
"""  # noqa: E501  # pylint: disable=line-too-long
import os
import re

from io import BufferedReader

import nibabel as nb
import numpy as np
import pandas as pd

from CPAC.utils.interfaces.fixes import FixHeaderApplyTransforms

from CPAC.generate_motion_statistics.generate_motion_statistics import \
    motion_power_statistics
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.qc.qcmetrics import regisQ
from CPAC.utils.interfaces.function import Function
from CPAC.utils.utils import check_prov_for_motion_tool

motion_params = ['movement-parameters', 'dvars',
                 'framewise-displacement-jenkinson']


def _connect_motion(wf, strat_pool, qc_file, brain_mask_key, final, pipe_num):
    # pylint: disable=invalid-name, too-many-arguments
    try:
        nodes = {'censor-indices': strat_pool.node_data('censor-indices')}
        wf.connect(nodes['censor-indices'].node, nodes['censor-indices'].out,
                   qc_file, 'censor_indices')
    except LookupError:
        nodes = {}
        qc_file.inputs.censor_indices = []
    motion_prov = strat_pool.get_cpac_provenance('movement-parameters')
    motion_correct_tool = check_prov_for_motion_tool(motion_prov)
    gen_motion_stats = motion_power_statistics('motion_stats-after_'
                                               f'{pipe_num}',
                                               motion_correct_tool)
    nodes = {
        **nodes,
        **{node_data: strat_pool.node_data(node_data) for node_data in [
            'subject', 'scan', brain_mask_key, 'max-displacement',
            *motion_params]}}
    if motion_correct_tool == '3dvolreg' and strat_pool.check_rpool(
            'coordinate-transformation'):
        nodes['coordinate-transformation'] = strat_pool.node_data(
            'coordinate-transformation')
        wf.connect(nodes['coordinate-transformation'].node,
                   nodes['coordinate-transformation'].out,
                   gen_motion_stats, 'inputspec.transformations')
    elif motion_correct_tool == 'mcflirt' and strat_pool.check_rpool(
            'rels-displacement'):
        nodes['rels-displacement'] = strat_pool.node_data('rels-displacement')
        wf.connect(nodes['rels-displacement'].node,
                   nodes['rels-displacement'].out,
                   gen_motion_stats, 'inputspec.rels_displacement')
    wf.connect([
        (final['func'].node, gen_motion_stats, [
            (final['func'].out, 'inputspec.motion_correct')]),
        (nodes['subject'].node, gen_motion_stats, [
            (nodes['subject'].out, 'inputspec.subject_id')]),
        (nodes['scan'].node, gen_motion_stats, [
            (nodes['scan'].out, 'inputspec.scan_id')]),
        (nodes['movement-parameters'].node, gen_motion_stats, [
            (nodes['movement-parameters'].out,
                'inputspec.movement_parameters')]),
        (nodes['max-displacement'].node, gen_motion_stats, [
            (nodes['max-displacement'].out,
                'inputspec.max_displacement')]),
        (nodes[brain_mask_key].node, gen_motion_stats, [
            (nodes[brain_mask_key].out, 'inputspec.mask')]),
        (gen_motion_stats, qc_file, [
            ('outputspec.DVARS_1D', 'dvars_after'),
            ('outputspec.FDJ_1D', 'fdj_after')]),
        *[(nodes[node].node, qc_file, [
            (nodes[node].out, node.replace('-', '_'))
        ]) for node in motion_params]])
    return wf


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


def generate_xcp_qc(desc, bold2t1w_mask, t1w_mask, bold2template_mask,
                    template_mask, original_func, final_func,
                    movement_parameters, dvars, censor_indices,
                    framewise_displacement_jenkinson, dvars_after, fdj_after,
                    template):
    # pylint: disable=too-many-arguments, too-many-locals, invalid-name
    """Function to generate an RBC-style QC CSV

    Parameters
    ----------
    desc : str
        description string

    original_func : str
        path to original 'bold' image

    final_bold : str
        path to 'desc-preproc_bold' image

    bold2t1w_mask : str
        path to bold-to-T1w transform applied to space-bold_desc-brain_mask
        with space-T1w_desc-brain_mask reference

    t1w_mask : str
        path to space-T1w_desc-brain_mask

    bold2template_mask : str
        path to space-template_desc-bold_mask

    template_mask : str
        path to space-template_desc-T1w_mask

    movement_parameters: str
        path to movement parameters

    dvars : str
        path to DVARS before motion correction

    censor_indices : list
        list of indices of censored volumes

    framewise_displacement_jenkinson : str
        path to framewise displacement (Jenkinson) before motion correction

    dvars_after : str
        path to DVARS on final 'bold' image

    fdj_after : str
        path to framewise displacement (Jenkinson) on final 'bold' image

    template : str
        path to registration template

    Returns
    -------
    str
        path to desc-xcp_quality TSV
    """
    columns = (
        'sub,ses,task,run,desc,space,meanFD,relMeansRMSMotion,'
        'relMaxRMSMotion,meanDVInit,meanDVFinal,nVolCensored,nVolsRemoved,'
        'motionDVCorrInit,motionDVCorrFinal,coregDice,coregJaccard,'
        'coregCrossCorr,coregCoverage,normDice,normJaccard,normCrossCorr,'
        'normCoverage'.split(',')
    )

    images = {
        'original_func': nb.load(original_func),
        'final_func': nb.load(final_func),
    }

    # `sub` through `desc`
    from_bids = {
        **strings_from_bids(original_func),
        'space': os.path.basename(template).split('.', 1)[0].split('_', 1)[0],
        'desc': desc
    }

    # `nVolCensored` & `nVolsRemoved`
    n_vols_censored = len(
        censor_indices) if censor_indices is not None else 'unknown'
    shape_params = {'nVolCensored': n_vols_censored,
                    'nVolsRemoved': images['final_func'].shape[3] -
                    images['original_func'].shape[3]}

    if isinstance(final_func, BufferedReader):
        final_func = final_func.name
    qc_filepath = os.path.join(os.getcwd(), 'xcpqc.tsv')

    desc_span = re.search(r'_desc-.*_', final_func)
    if desc_span:
        desc_span = desc_span.span()
        final_func = '_'.join([
            final_func[:desc_span[0]],
            final_func[desc_span[1]:]
        ])
    del desc_span

    # `meanFD (Jenkinson)`
    power_params = {'meanFD': np.mean(np.loadtxt(
        framewise_displacement_jenkinson))}

    # `relMeansRMSMotion` & `relMaxRMSMotion`
    mot = np.genfromtxt(movement_parameters).T
    # Relative RMS of translation
    rms = np.sqrt(mot[3] ** 2 + mot[4] ** 2 + mot[5] ** 2)
    rms_params = {
        'relMeansRMSMotion': [np.mean(rms)],
        'relMaxRMSMotion': [np.max(rms)]
    }

    # `meanDVInit` & `meanDVFinal`
    meanDV = {'meanDVInit': np.mean(np.loadtxt(dvars))}
    try:
        meanDV['motionDVCorrInit'] = dvcorr(
            dvars, framewise_displacement_jenkinson)
    except ValueError as value_error:
        meanDV['motionDVCorrInit'] = f'ValueError({str(value_error)})'
    if dvars_after:
        if not fdj_after:
            fdj_after = framewise_displacement_jenkinson
        meanDV['meanDVFinal'] = np.mean(np.loadtxt(dvars_after))
        try:
            meanDV['motionDVCorrFinal'] = dvcorr(dvars_after, fdj_after)
        except ValueError as value_error:
            meanDV['motionDVCorrFinal'] = f'ValueError({str(value_error)})'

    # Overlap
    overlap_params = regisQ(bold2t1w_mask=bold2t1w_mask, t1w_mask=t1w_mask,
                            bold2template_mask=bold2template_mask,
                            template_mask=template_mask)

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


def qc_xcp(wf, cfg, strat_pool, pipe_num, opt=None):
    # pylint: disable=invalid-name, unused-argument
    """
    {'name': 'qc_xcp',
     'config': ['pipeline_setup', 'output_directory', 'quality_control'],
     'switch': ['generate_xcpqc_files'],
     'option_key': 'None',
     'option_val': 'None',
     'inputs': ['space-T1w_desc-brain_mask',
                'from-bold_to-T1w_mode-image_desc-linear_xfm',
                'space-bold_desc-brain_mask',
                ['space-template_desc-bold_mask',
                 'space-EPItemplate_desc-bold_mask'],
                'space-template_desc-T1w_mask',
                ['T1w-template-funcreg', 'EPI-template-funcreg']],
     'outputs': ['desc-xcp_quality']}
    """
    qc_file = pe.Node(Function(input_names=['desc', 'bold2t1w_mask',
                                            't1w_mask', 'bold2template_mask',
                                            'template_mask', 'original_func',
                                            'final_func',
                                            'movement_parameters', 'dvars',
                                            'censor_indices',
                                            'framewise_displacement_jenkinson',
                                            'dvars_after', 'fdj_after',
                                            'template'],
                               output_names=['qc_file'],
                               function=generate_xcp_qc,
                               as_module=True),
                      name=f'qcxcp_{pipe_num}')
    qc_file.inputs.desc = 'preproc'
    original = {}
    final = {}
    original['anat'] = strat_pool.node_data('T1w')
    original['func'] = strat_pool.node_data('bold')
    final['anat'] = strat_pool.node_data('desc-preproc_T1w')
    t1w_bold = strat_pool.node_data('space-T1w_desc-mean_bold')

    bold_to_T1w_mask = pe.Node(
        FixHeaderApplyTransforms(dimension=3, interpolation='NearestNeighbor'),
        name='bold_to_T1w_mask')

    nodes = {key: strat_pool.node_data(key) for key in [
        'from-bold_to-T1w_mode-image_desc-linear_xfm',
        'space-bold_desc-brain_mask', 'space-template_desc-T1w_mask']}
    nodes['t1w_mask'] = strat_pool.node_data('space-T1w_desc-brain_mask')
    nodes['bold2template_mask'] = strat_pool.node_data([
        'space-template_desc-bold_mask', 'space-EPItemplate_desc-bold_mask'])

    nodes['template'] = strat_pool.node_data(['T1w-template-funcreg',
                                              'EPI-template-funcreg'])

    wf.connect([
        (nodes['space-bold_desc-brain_mask'].node, bold_to_T1w_mask, [
            (nodes['space-bold_desc-brain_mask'].out, 'input_image')]),
        (nodes['t1w_mask'].node, bold_to_T1w_mask, [
            (nodes['t1w_mask'].out, 'reference_image')]),
        (nodes['t1w_mask'].node, qc_file, [
            (nodes['t1w_mask'].out, 't1w_mask')]),
        (nodes['from-bold_to-T1w_mode-image_desc-linear_xfm'].node,
         bold_to_T1w_mask, [
            (nodes['from-bold_to-T1w_mode-image_desc-linear_xfm'].out,
             'transforms')]),
        (bold_to_T1w_mask, qc_file, [('output_image', 'bold2T1w_mask')]),
        (nodes['bold2template_mask'].node, qc_file, [
            (nodes['bold2template_mask'].out, 'bold2template_mask')]),
        (nodes['space-template_desc-T1w_mask'].node, qc_file, [
            (nodes['space-template_desc-T1w_mask'].out, 'template_mask')]),
        (original['func'].node, qc_file, [
            (original['func'].out, 'original_func')]),
        (final['func'].node, qc_file, [(final['func'].out, 'final_func')]),
        (t1w_bold.node, qc_file, [(t1w_bold.out, 'space_T1w_bold')]),
        (nodes['template'].node, qc_file, [
            (nodes['template'].out, 'template')])])

    return wf, {'desc-xcp_quality': qc_file}


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
    from_bids = {k: from_bids[k] for k in from_bids}
    if 'space' not in from_bids:
        from_bids['space'] = 'native'
    return from_bids
