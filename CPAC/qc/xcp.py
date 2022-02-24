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
"""  # noqa: E501  # pylint: disable=line-too-long
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

motion_params = ['movement-parameters', 'dvars',
                 'framewise-displacement-jenkinson']


def calculate_overlap(image_pair):
    '''
    Function to calculate Dice, Jaccard, CrossCorr and Coverage:cite:`cite-Penn19` from a
    pair of arrays

    Parameters
    ----------
    image_pair : 2-tuple
        array of which to calculate overlaps metrics

    Returns
    -------
    coefficents : dict
        coeffiecients['dice'] : float
            Dice index

        coeffiecients['jaccard'] : float
            Jaccard index

        coeffiecients['cross_corr'] : float
            cross-correlation

        coeffiecients['coverage'] : float
            coverage index

    Examples
    --------
    >>> import numpy as np
    >>> a1 = np.array([0, 0, 0, 1, 1, 1])
    >>> a2 = np.array([0, 0, 1, 1, 0, 1])
    >>> tuple(calculate_overlap((a1, a2)).values())
    (0.6666666666666666, 0.5, 0.33333333333333326, 0.6666666666666666)
    >>> tuple(calculate_overlap((a1, a1)).values())
    (1.0, 1.0, 0.9999999999999998, 1.0)
    >>> tuple(calculate_overlap((a2, a2)).values())
    (1.0, 1.0, 0.9999999999999998, 1.0)
    '''  # noqa: E501  # pylint: disable=line-too-long
    if len(image_pair) != 2:
        raise IndexError('`calculate_overlap` requires 2 images, but '
                         f'{len(image_pair)} were provided')
    if len(image_pair[0]) != len(image_pair[1]):
        image_pair = _repeat_shorter(image_pair)
    image_pair = tuple(image.astype(bool) for image in image_pair)
    intersect = image_pair[0] * image_pair[1]
    vols = [np.sum(image) for image in image_pair]
    vol_intersect = np.sum(intersect)
    vol_sum = sum(vols)
    vol_union = vol_sum - vol_intersect
    coefficients = {
        'dice': 2 * vol_intersect / vol_sum,
        'jaccard': vol_intersect / vol_union,
        'cross_corr': np.corrcoef(image_pair)[0, 1],
        'coverage': vol_intersect / min(vols)
    }
    for name, coefficient in coefficients.items():
        if not 1 >= abs(coefficient) >= 0 and not np.isnan(coefficient):
            raise ValueError(f'Valid range for {name} is [0, 1] but value '
                             f'{coefficient} was calculated.')
    return coefficients


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


def _connect_xcp(wf, strat_pool, qc_file, original, final, t1w_bold,
                 brain_mask_key, output_key, pipe_num):
    # pylint: disable=invalid-name, too-many-arguments
    if (
        strat_pool.check_rpool('movement-parameters') and
        strat_pool.check_rpool(brain_mask_key)
    ):
        wf = _connect_motion(wf, strat_pool, qc_file, brain_mask_key, final,
                             pipe_num)
    else:
        qc_file.inputs.censor_indices = []
        for key in [*motion_params, 'movement_parameters',
                    'framewise_displacement_jenkinson', 'dvars_after',
                    'fdj_after']:
            setattr(qc_file.inputs, key, None)
    wf.connect([
        (original['anat'].node, qc_file, [
            (original['anat'].out, 'original_anat')]),
        (original['func'].node, qc_file, [
            (original['func'].out, 'original_func')]),
        (final['anat'].node, qc_file, [(final['anat'].out, 'final_anat')]),
        (final['func'].node, qc_file, [(final['func'].out, 'final_func')]),
        (t1w_bold.node, qc_file, [(t1w_bold.out, 'space_T1w_bold')])])
    outputs = {output_key: (qc_file, 'qc_file')}
    return wf, outputs


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


def generate_xcp_qc(space, desc, original_anat,
                    final_anat, original_func, final_func, space_T1w_bold,
                    movement_parameters, dvars, censor_indices,
                    framewise_displacement_jenkinson, dvars_after, fdj_after,
                    template=None):
    # pylint: disable=too-many-arguments, too-many-locals, invalid-name
    """Function to generate an RBC-style QC CSV

    Parameters
    ----------
    space : str
        'native' or 'template'

    desc : str
        description string

    original_anat : str
        path to original 'T1w' image

    final_anat : str
        path to 'desc-preproc_T1w' image

    original_func : str
        path to original 'bold' image

    final_bold : str
        path to 'desc-preproc_bold' image

    space_T1w_bold : str
        path to 'space-T1w_desc-mean_bold' image

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
        path to template

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
        'original_anat': nb.load(original_anat),
        'original_func': nb.load(original_func),
        'final_anat': nb.load(final_anat),
        'final_func': nb.load(final_func),
        'space-T1w_bold': nb.load(space_T1w_bold)
    }
    if template is not None:
        images['template'] = nb.load(template)

    # `sub` through `desc`
    from_bids = {
        **strings_from_bids(original_func),
        'space': space,
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

    if dvars:
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
    else:
        meanDV = _na_dict([f'{dv}DV{"Corr" + ts if dv == "motion" else ts}'
                           for dv in ['mean', 'motion']
                           for ts in ['Init', 'Final']])
        power_params = {'meanFD': 'n/a'}
        rms_params = _na_dict(['relMeansRMSMotion', 'relMaxRMSMotion'])

    # Overlap
    overlap_images = {variable: image.get_fdata().ravel() for
                      variable, image in images.items() if
                      variable in ['space-T1w_bold',
                                   'original_anat', 'template']}
    overlap_params = {}
    (overlap_params['coregDice'], overlap_params['coregJaccard'],
     overlap_params['coregCrossCorr'], overlap_params['coregCoverage']
     ) = [[item] for item in calculate_overlap(
        (overlap_images['space-T1w_bold'], overlap_images['original_anat'])
    ).values()]
    if space == 'native':
        for key in ['normDice', 'normJaccard', 'normCrossCorr',
                    'normCoverage']:
            overlap_params[key] = ['N/A: native space']
    elif template is not None:
        (overlap_params['normDice'], overlap_params['normJaccard'],
         overlap_params['normCrossCorr'], overlap_params['normCoverage']
         ) = [[item] for item in calculate_overlap(
            (images['final_func'].get_fdata().ravel(),
             overlap_images['template'])).values()]
    else:
        overlap_params = _na_dict(['normDice', 'normJaccard', 'normCrossCorr',
                                   'normCoverage'])

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


def _na_dict(keys):
    return {key: 'n/a' for key in keys}


def _prep_qc_xcp(strat_pool, pipe_num, space):
    qc_file = pe.Node(Function(input_names=['subject', 'scan',
                                            'space', 'desc', 'template',
                                            'original_func', 'final_func',
                                            'original_anat', 'final_anat',
                                            'space_T1w_bold',
                                            'movement_parameters',
                                            'censor_indices', 'dvars',
                                            'framewise_displacement_jenkinson',
                                            'dvars_after', 'fdj_after'],
                               output_names=['qc_file'],
                               function=generate_xcp_qc,
                               as_module=True),
                      name=f'xcpqc-{space}_{pipe_num}')
    qc_file.inputs.desc = 'preproc'
    qc_file.inputs.space = space
    original = {}
    final = {}
    original['anat'] = strat_pool.node_data('T1w')
    original['func'] = strat_pool.node_data('bold')
    final['anat'] = strat_pool.node_data('desc-preproc_T1w')
    t1w_bold = strat_pool.node_data('space-T1w_desc-mean_bold')
    return qc_file, original, final, t1w_bold


def qc_xcp_native(wf, cfg, strat_pool, pipe_num, opt=None):
    # pylint: disable=invalid-name, unused-argument
    """
    {'name': 'qc_xcp_native',
     'config': ['pipeline_setup', 'output_directory', 'quality_control'],
     'switch': ['generate_xcpqc_files'],
     'option_key': 'None',
     'option_val': 'None',
     'inputs': [('bold', 'subject', 'scan', 'T1w', 'max-displacement', 'dvars',
                'censor-indices', 'desc-preproc_bold',
                'desc-preproc_T1w', 'space-T1w_desc-mean_bold',
                'space-bold_desc-brain_mask', 'movement-parameters',
                'framewise-displacement-jenkinson', 'rels-displacement',
                'coordinate-transformation')],
     'outputs': ['desc-xcp_quality']}
    """
    space = 'native'
    qc_file, original, final, t1w_bold = _prep_qc_xcp(strat_pool, pipe_num,
                                                      space)
    final['func'] = strat_pool.node_data('desc-preproc_bold')
    return _connect_xcp(wf, strat_pool, qc_file, original, final, t1w_bold,
                        'space-bold_desc-brain_mask', 'desc-xcp_quality',
                        pipe_num)


def qc_xcp_skullstripped(wf, cfg, strat_pool, pipe_num, opt=None):
    # pylint: disable=invalid-name, unused-argument
    r"""
    Same as ``qc_xcp_native`` except no motion inputs.
    Node Block:
    {'name': 'qc_xcp_skullstripped',
     'config': ['pipeline_setup', 'output_directory', 'quality_control'],
     'switch': ['generate_xcpqc_files'],
     'option_key': 'None',
     'option_val': 'None',
     'inputs': [('bold', 'subject', 'scan', 'T1w', 'desc-preproc_bold',
                'desc-preproc_T1w', 'space-T1w_desc-mean_bold',
                'space-bold_desc-brain_mask')],
     'outputs': ['desc-xcp_quality']}
    """
    # If strat has dvars, it should be captured by 'qc_xcp_native'
    if 'dvars' in strat_pool.get('desc-preproc_bold').get(
            'json', {}).get('Sources', {}):
        return wf, {}
    return qc_xcp_native(wf, cfg, strat_pool, pipe_num, opt)


def qc_xcp_template(wf, space, strat_pool, pipe_num, reference):
    """Handle either T1 or EPI template XCP-QC."""
    if not (
        strat_pool.check_rpool(f'space-{space}_desc-preproc_bold') and
        strat_pool.check_rpool(reference)
    ):
        return wf, {}
    qc_file, original, final, t1w_bold = _prep_qc_xcp(strat_pool, pipe_num,
                                                      space)
    final['func'] = strat_pool.node_data(f'space-{space}_desc-preproc_bold')
    template = strat_pool.node_data(reference)
    wf.connect(template.node, template.out, qc_file, 'template')
    return _connect_xcp(wf, strat_pool, qc_file, original, final, t1w_bold,
                        f'space-{space}_desc-bold_mask',
                        f'space-{space}_desc-xcp_quality', pipe_num)


# pylint: disable=invalid-name, unused-argument
def qc_xcp_EPItemplate(wf, cfg, strat_pool, pipe_num, opt=None):
    """
    {'name': 'qc_xcp_EPItemplate',
     'config': ['pipeline_setup', 'output_directory', 'quality_control'],
     'switch': ['generate_xcpqc_files'],
     'option_key': 'None',
     'option_val': 'None',
     'inputs': [('bold', 'subject', 'scan', 'T1w',
                'space-T1w_desc-mean_bold', 'desc-preproc_T1w',
                'space-EPItemplate_desc-bold_mask')],
     'outputs': ['space-EPItemplate_desc-xcp_quality']}
    """
    return qc_xcp_template(wf, 'EPItemplate', strat_pool, pipe_num,
                           'EPI-template')


# pylint: disable=invalid-name, unused-argument
def qc_xcp_T1template(wf, cfg, strat_pool, pipe_num, opt=None):
    """
    {'name': 'qc_xcp_T1template',
     'config': ['pipeline_setup', 'output_directory', 'quality_control'],
     'switch': ['generate_xcpqc_files'],
     'option_key': 'None',
     'option_val': 'None',
     'inputs': [('bold', 'subject', 'scan', 'T1w',
                'space-T1w_desc-mean_bold', 'desc-preproc_T1w',
                'space-template_desc-bold_mask')],
     'outputs': ['space-template_desc-xcp_quality']}
    """
    return qc_xcp_template(wf, 'template', strat_pool, pipe_num,
                           'T1w-brain-template-funcreg')


def _repeat_shorter(images):
    '''
    Parameters
    ----------
    images : 2-tuple

    Returns
    -------
    images : 2-tuple

    Examples
    --------
    >>> _repeat_shorter((np.array([1, 2, 3]), np.array([1])))
    (array([1, 2, 3]), array([1, 1, 1]))
    >>> _repeat_shorter((np.array([0, 0]), (np.array([1, 1, 2, 2, 3, 4]))))
    (array([0, 0, 0, 0, 0, 0]), array([1, 1, 2, 2, 3, 4]))
    '''
    lens = (len(images[0]), len(images[1]))
    if lens[1] > lens[0] and lens[1] % lens[0] == 0:
        return (np.tile(images[0], lens[1] // lens[0]), images[1])
    if lens[0] > lens[1] and lens[0] % lens[1] == 0:
        return (images[0], np.tile(images[1], lens[0] // lens[1]))
    raise ValueError('operands could not be broadcast together with shapes '
                     f'({lens[0]},) ({lens[1]},)')


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
