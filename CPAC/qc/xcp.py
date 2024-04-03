"""
.. seealso::

      `User Guide: Generate eXtensible Connectivity Pipeline-style quality control files <https://fcp-indi.github.io/docs/user/xcpqc>`_

Columns
-------
sub : str
    subject label :footcite:`BIDS21`
ses : str
    session label :footcite:`BIDS21`
task : str
    task label :footcite:`BIDS21`
run : int
    run index :footcite:`BIDS21`
desc : str
    description :footcite:`BIDS21`
regressors : str
    'Name' of regressors in the current fork
space : str
    space label :footcite:`BIDS21`
meanFD : float
    mean Jenkinson framewise displacement :footcite:`xcp_22,Jenk02` :func:`CPAC.generate_motion_statistics.calculate_FD_J` after preprocessing
relMeansRMSMotion : float
    "mean value of RMS motion" :footcite:`xcp_22,Ciri19`
relMaxRMSMotion : float
    "maximum vaue of RMS motion" :footcite:`xcp_22,Ciri19`
meanDVInit : float
    "mean DVARS" :footcite:`xcp_22,Ciri19`
meanDVFinal : float
    "mean DVARS" :footcite:`xcp_22,Ciri19`
nVolCensored : int
    "total number of volume(s) censored :footcite:`Ciri19`
nVolsRemoved : int
    number of volumes in derivative minus number of volumes in original
    functional scan
motionDVCorrInit : float
    "correlation of RMS and DVARS before regresion" :footcite:`Ciri19`
motionDVCorrFinal : float
    "correlation of RMS and DVARS after regresion" :footcite:`Ciri19`
coregDice : float
    "Coregsitration of Functional and T1w:[…] Dice index" :footcite:`xcp_22,Ciri19`
coregJaccard : float
    "Coregsitration of Functional and T1w:[…] Jaccard index" :footcite:`xcp_22,Ciri19`
coregCrossCorr : float
    "Coregsitration of Functional and T1w:[…] cross correlation" :footcite:`xcp_22,Ciri19`
coregCoverag : float
    "Coregsitration of Functional and T1w:[…] Coverage index" :footcite:`xcp_22,Ciri19`
normDice : float
    "Normalization of T1w/Functional to Template:[…] Dice index" :footcite:`xcp_22,Ciri19`
normJaccard : float
    "Normalization of T1w/Functional to Template:[…] Jaccard index" :footcite:`xcp_22,Ciri19`
normCrossCorr : float
    "Normalization of T1w/Functional to Template:[…] cross correlation" :footcite:`xcp_22,Ciri19`
normCoverage : float
    "Normalization of T1w/Functional to Template:[…] Coverage index" :footcite:`xcp_22,Ciri19`
"""  # noqa: E501  # pylint: disable=line-too-long
import os
import re
from io import BufferedReader

import nibabel as nb
import numpy as np
import pandas as pd
from bids.layout import parse_file_entities
from nipype.interfaces import afni, fsl

from CPAC.generate_motion_statistics.generate_motion_statistics import \
    DVARS_strip_t0, ImageTo1D
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.nodeblock import nodeblock
from CPAC.qc.qcmetrics import regisQ
from CPAC.utils.interfaces.function import Function

motion_params = ['dvars', 'framewise-displacement-jenkinson',
                 'desc-movementParametersUnfiltered_motion', 'desc-movementParameters_motion']


def _connect_motion(wf, nodes, strat_pool, qc_file, pipe_num):
    """
    Connect the motion metrics to the workflow.

    Parameters
    ----------
    wf : nipype.pipeline.engine.Workflow
        The workflow to connect the motion metrics to.

    nodes : dict
        Dictionary of nodes already collected from the strategy pool.

    strat_pool : CPAC.pipeline.engine.ResourcePool
        The current strategy pool.

    qc_file : nipype.pipeline.engine.Node
        A function node with the function ``generate_xcp_qc``.

    pipe_num : int

    Returns
    -------
    wf : nipype.pipeline.engine.Workflow
    """
    # pylint: disable=invalid-name, too-many-arguments
    try:
        nodes = {**nodes,
                 'censor-indices': strat_pool.node_data('censor-indices')}
        wf.connect(nodes['censor-indices'].node, nodes['censor-indices'].out,
                   qc_file, 'censor_indices')
    except LookupError:
        qc_file.inputs.censor_indices = []
    cal_DVARS = pe.Node(ImageTo1D(method='dvars'),
                        name=f'cal_DVARS_{pipe_num}',
                        mem_gb=0.4,
                        mem_x=(739971956005215 / 151115727451828646838272,
                               'in_file'),
                        throttle=True)
    cal_DVARS_strip = pe.Node(Function(input_names=['file_1D'],
                                       output_names=['out_file'],
                                       function=DVARS_strip_t0,
                                       as_module=True),
                              name=f'cal_DVARS_strip_{pipe_num}')
    wf.connect([
        (nodes['desc-preproc_bold'].node, cal_DVARS, [
            (nodes['desc-preproc_bold'].out, 'in_file')]),
        (nodes['space-bold_desc-brain_mask'].node, cal_DVARS, [
            (nodes['space-bold_desc-brain_mask'].out, 'mask')]),
        (cal_DVARS, cal_DVARS_strip, [('out_file', 'file_1D')]),
        (cal_DVARS_strip, qc_file, [('out_file', 'dvars_after')]),
        *[(nodes[node].node, qc_file, [
            (nodes[node].out, node.replace('-', '_'))
        ]) for node in motion_params if node in nodes]])
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


def generate_xcp_qc(sub, ses, task, run, desc, regressors, bold2t1w_mask,
                    t1w_mask, bold2template_mask, template_mask, original_func,
                    final_func, movement_parameters, dvars, censor_indices,
                    framewise_displacement_jenkinson, dvars_after, template):
    # pylint: disable=too-many-arguments, too-many-locals, invalid-name
    """Function to generate an RBC-style QC CSV

    Parameters
    ----------
    sub : str
        subject ID

    ses : str
        session ID

    task : str
        task ID

    run : str or int
        run ID

    desc : str
        description string

    regressors : str
        'Name' of regressors in fork

    original_func : str
        path to original 'bold' image

    final_bold : str
        path to 'space-template_desc-preproc_bold' image

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

    template : str
        path to registration template

    Returns
    -------
    str
        path to space-template_desc-xcp_quality TSV
    """
    columns = (
        'sub,ses,task,run,desc,regressors,space,meanFD,relMeansRMSMotion,'
        'relMaxRMSMotion,meanDVInit,meanDVFinal,nVolCensored,nVolsRemoved,'
        'motionDVCorrInit,motionDVCorrFinal,coregDice,coregJaccard,'
        'coregCrossCorr,coregCoverage,normDice,normJaccard,normCrossCorr,'
        'normCoverage'.split(',')
    )

    images = {
        'original_func': nb.load(original_func),
        'final_func': nb.load(final_func),
    }

    # `sub` through `space`
    from_bids = {
        'sub': sub, 'ses': ses, 'task': task, 'run': run, 'desc': desc,
        'regressors': regressors,
        'space': os.path.basename(template).split('.', 1)[0].split('_', 1)[0]
    }
    if from_bids['space'].startswith('tpl-'):
        from_bids['space'] = from_bids['space'][4:]

    # `nVolCensored` & `nVolsRemoved`
    n_vols_censored = len(
        censor_indices) if censor_indices is not None else 'unknown'
    shape_params = {'nVolCensored': n_vols_censored,
                    'nVolsRemoved': images['original_func'].shape[3] -
                    images['final_func'].shape[3]}

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
    meanDV['meanDVFinal'] = np.mean(np.loadtxt(dvars_after))
    try:
        meanDV['motionDVCorrFinal'] = dvcorr(dvars_after,
                                             framewise_displacement_jenkinson)
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


def get_bids_info(subject, scan, wf_name):
    """
    Function to gather BIDS information from a strat_pool

    Parameters
    ----------
    subject : str
        subject ID

    scan : str
        scan ID

    wf_name : str
        workflow name

    Returns
    -------
    subject : str
        subject ID

    session : str
        session ID

    task : str
        task ID

    run : str or int
        run ID

    Examples
    --------
    >>> subject, session, task, run = get_bids_info(
    ...     subject='DavidBowman', scan='rest_acq-1_run-1',
    ...     wf_name='cpac_DavidBowman_3')
    >>> subject
    'DavidBowman'
    >>> session
    '3'
    >>> task
    'rest'
    >>> run
    '1'
    >>> get_bids_info(subject='sub-colornest035', scan='rest_run-01',
    ...               wf_name='cpac_sub-colornest035_ses-1')
    ('colornest035', '1', 'rest', '01')
    """
    returns = ('subject', 'session', 'task', 'run')
    ses = wf_name.split('_')[-1]
    if not ses.startswith('ses-'):
        ses = f'ses-{ses}'
    if not subject.startswith('sub-'):
        subject = f'sub-{subject}'
    resource = '_'.join([subject, ses,
                         scan if 'task' in scan else f'task-{scan}'])
    entities = parse_file_entities(resource)
    returns = {key: entities.get(key) for key in returns}
    if any(value is None for value in returns.values()):
        entity_parts = resource.split('_')

        def get_entity_part(key):
            key = f'{key}-'
            matching_parts = [part for part in entity_parts if
                              part.startswith(key)]
            if matching_parts:
                return matching_parts[0].replace(key, '')
            return None

        for key, value in returns.items():
            if value is None:
                if key == 'task':
                    returns[key] = get_entity_part(key)
                else:
                    returns[key] = get_entity_part(key[:3])
    return tuple(str(returns.get(key)) for key in [
                 'subject', 'session', 'task', 'run'])


@nodeblock(
    name="qc_xcp",
    config=["pipeline_setup", "output_directory", "quality_control"],
    switch=["generate_xcpqc_files"],
    inputs=[
        (
            "subject",
            "scan",
            "bold",
            "desc-preproc_bold",
            "space-T1w_sbref",
            "space-T1w_desc-brain_mask",
            "max-displacement",
            "space-template_desc-preproc_bold",
            "space-bold_desc-brain_mask",
            ["T1w-brain-template-mask", "EPI-template-mask"],
            ["space-template_desc-bold_mask", "space-EPItemplate_desc-bold_mask"],
            "regressors",
            ["T1w-brain-template-funcreg", "EPI-brain-template-funcreg"],
            ["desc-movementParametersUnfiltered_motion", "desc-movementParameters_motion"],
            "dvars",
            "framewise-displacement-jenkinson",
        )
    ],
    outputs={
        "space-template_desc-xcp_quality": {"Template": "T1w-brain-template-mask"}
    },
)
def qc_xcp(wf, cfg, strat_pool, pipe_num, opt=None):
    # pylint: disable=invalid-name, unused-argument
    if cfg['nuisance_corrections', '2-nuisance_regression', 'run'
           ] and not strat_pool.check_rpool('regressors'):
        return wf, {}
    bids_info = pe.Node(Function(input_names=['subject', 'scan', 'wf_name'],
                                 output_names=['subject', 'session', 'task',
                                               'run'],
                                 imports=['from bids.layout import '
                                          'parse_file_entities'],
                                 function=get_bids_info, as_module=True),
                        name=f'bids_info_{pipe_num}')
    bids_info.inputs.wf_name = wf.name
    qc_file = pe.Node(Function(input_names=['sub', 'ses', 'task', 'run',
                                            'desc', 'bold2t1w_mask',
                                            't1w_mask', 'bold2template_mask',
                                            'template_mask', 'original_func',
                                            'final_func', 'template',
                                            'movement_parameters', 'dvars',
                                            'censor_indices', 'regressors',
                                            'framewise_displacement_jenkinson',
                                            'dvars_after'],
                               output_names=['qc_file'],
                               function=generate_xcp_qc,
                               as_module=True),
                      name=f'qcxcp_{pipe_num}')
    qc_file.inputs.desc = 'preproc'
    qc_file.inputs.regressors = strat_pool.node_data(
        'regressors').node.name.split('regressors_'
    )[-1][::-1].split('_', 1)[-1][::-1]
    bold_to_T1w_mask = pe.Node(interface=fsl.ImageMaths(),
                               name=f'binarize_bold_to_T1w_mask_{pipe_num}',
                               op_string='-bin ')
    nodes = {key: strat_pool.node_data(key) for key in [
        'bold', 'desc-preproc_bold', 'max-displacement',
        'scan', 'space-bold_desc-brain_mask', 'space-T1w_desc-brain_mask',
        'space-T1w_sbref', 'space-template_desc-preproc_bold',
        'subject', *motion_params] if strat_pool.check_rpool(key)}
    nodes['bold2template_mask'] = strat_pool.node_data([
        'space-template_desc-bold_mask', 'space-EPItemplate_desc-bold_mask'])
    nodes['template_mask'] = strat_pool.node_data(
        ['T1w-brain-template-mask', 'EPI-template-mask'])
    nodes['template'] = strat_pool.node_data(['T1w-brain-template-funcreg',
                                              'EPI-brain-template-funcreg'])
    resample_bold_mask_to_template = pe.Node(
        afni.Resample(), name=f'resample_bold_mask_to_anat_res_{pipe_num}',
        mem_gb=0, mem_x=(0.0115, 'in_file', 't'))
    resample_bold_mask_to_template.inputs.outputtype = 'NIFTI_GZ'
    wf = _connect_motion(wf, nodes, strat_pool, qc_file, pipe_num=pipe_num)
    wf.connect([
        (nodes['subject'].node, bids_info, [
            (nodes['subject'].out, 'subject')]),
        (nodes['scan'].node, bids_info, [(nodes['scan'].out, 'scan')]),
        (nodes['space-T1w_sbref'].node, bold_to_T1w_mask, [
            (nodes['space-T1w_sbref'].out, 'in_file')]),
        (nodes['space-T1w_desc-brain_mask'].node, qc_file, [
            (nodes['space-T1w_desc-brain_mask'].out, 't1w_mask')]),
        (bold_to_T1w_mask, qc_file, [('out_file', 'bold2t1w_mask')]),
        (nodes['template_mask'].node, qc_file, [
            (nodes['template_mask'].out, 'template_mask')]),
        (nodes['bold'].node, qc_file, [(nodes['bold'].out, 'original_func')]),
        (nodes['space-template_desc-preproc_bold'].node, qc_file, [
            (nodes['space-template_desc-preproc_bold'].out, 'final_func')]),
        (nodes['template'].node, qc_file, [
            (nodes['template'].out, 'template')]),
        (nodes['template_mask'].node, resample_bold_mask_to_template, [
             (nodes['template_mask'].out, 'master')]),
        (nodes['bold2template_mask'].node, resample_bold_mask_to_template,
            [(nodes['bold2template_mask'].out, 'in_file')]),
        (resample_bold_mask_to_template, qc_file, [
            ('out_file', 'bold2template_mask')]),
        (bids_info, qc_file, [
            ('subject', 'sub'),
            ('session', 'ses'),
            ('task', 'task'),
            ('run', 'run')])])

    return wf, {'space-template_desc-xcp_quality': (qc_file, 'qc_file')}
