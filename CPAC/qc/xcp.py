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
import os
import re
from io import BufferedReader

import nibabel as nb
import numpy as np
import pandas as pd
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.interfaces.function import Function


def generate_desc_qc(original_anat, final_anat, original_func, final_func,
                     n_vols_censored, space_T1w_bold):
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
    final_filename = final_func.split('/')[-1]
    bids_entities = final_filename.split('_')
    from_bids = dict(
        tuple(entity.split('-', 1)) if '-' in entity else
        ('suffix', entity) for entity in bids_entities
    )
    from_bids = {k: [from_bids[k]] for k in from_bids}
    if 'space' not in from_bids:
        from_bids['space'] = ['native']

    # `nVolCensored` & `nVolsRemoved`
    shape_params = {'nVolCensored': n_vols_censored,
                    'nVolsRemoved': images['final_func'].shape[3] -
                    images['original_func'].shape[3]}

    # `meanFD (Jenkinson)`
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
    power_params = {'meanFD': np.mean(np.loadtxt('_'.join([
        *final_func.split('_')[:-1],
        'framewise-displacement-jenkinson.1D'
    ])))}

    # `relMeansRMSMotion` & `relMaxRMSMotion`
    mot = np.genfromtxt('_'.join([
        *final_func.split('_')[:-1],
        'movement-parameters.1D'
    ])).T
    # Relative RMS of translation
    rms = np.sqrt(mot[3] ** 2 + mot[4] ** 2 + mot[5] ** 2)
    rms_params = {
        'relMeansRMSMotion': [np.mean(rms)],
        'relMaxRMSMotion': [np.max(rms)]
    }

    # `meanDVinit` & `meanDVFinal`
    # ?

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
        **overlap_params
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


def qc_xcp(wf, cfg, strat_pool, pipe_num, opt=None):
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

    qc_file = pe.Node(Function(input_names=['original_func', 'final_func',
                                            'original_anat', 'final_anat',
                                            'space_T1w_bold',
                                            'n_vols_censored'],
                               output_names=['qc_file'],
                               function=generate_desc_qc,
                               as_module=True),
                      name=f'xcpqc_{pipe_num}')

    try:
        censor_node, n_vols_censored = strat_pool.get_data('n_vols_censored')
        wf.connect(censor_node, n_vols_censored, qc_file, 'n_vols_censored')
    except LookupError:
        qc_file.inputs.n_vols_censored = 'unknown'

    wf.connect(original['anat']['node'], original['anat']['out'],
               qc_file, 'original_anat')
    wf.connect(original['func']['node'], original['func']['out'],
               qc_file, 'original_func')
    wf.connect(final['anat']['node'], final['anat']['out'],
               qc_file, 'final_anat')
    wf.connect(final['func']['node'], final['func']['out'],
               qc_file, 'final_func')
    wf.connect(t1w_bold['node'], t1w_bold['out'],
               qc_file, 'space_T1w_bold')

    outputs = {
        'xcpqc': (qc_file, 'qc_file'),
    }

    return (wf, outputs)
