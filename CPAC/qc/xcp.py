"""Generate eXtensible Connectivity Pipeline-style quality control files

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
    # TODO
relMaxRMSMotion : float
    # TODO
meanDVInit : float
    # TODO
meanDVFinal : float
    # TODO
nVolCensored : int
    # TODO
nVolsRemoved : int
    # TODO
motionDVCorrInit : float
    # TODO
motionDVCorrFinal : float
    # TODO
coregDice : float
    :cite:`cite-Ciri19` :cite:`cite-Penn19`
coregJaccard : float
    :cite:`cite-Ciri19` :cite:`cite-Penn19`
coregCrossCorr : float
    :cite:`cite-Ciri19` :cite:`cite-Penn19`
coregCoverag : float
    :cite:`cite-Ciri19` :cite:`cite-Penn19`
normDice : float
    :cite:`cite-Ciri19` :cite:`cite-Penn19`
normJaccard : float
    :cite:`cite-Ciri19` :cite:`cite-Penn19`
normCrossCorr : float
    :cite:`cite-Ciri19` :cite:`cite-Penn19`
normCoverage : float
    :cite:`cite-Ciri19` :cite:`cite-Penn19`

.. rubric References

.. bibliography:: /references/xcpqc_citation.bib
   :style: cpac_docs_style
   :cited:
   :keyprefix: cite-
"""  # noqa E501  # pylint: disable=line-too-long
import os
import re
from io import BufferedReader

import nibabel as nb
import numpy as np
import pandas as pd
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.interfaces.function import Function


# class OverlapInterface(IdentityInterface):
#     """XCP QC interface for overlap measures :cite:`cite-Ciri19` :cite:`cite-Penn19`."""  # noqa E501  # pylint: disable=line-too-long
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args,
#                          fields=['Dice', 'Jaccard', 'CrossCorr', 'Coverage'],
#                          **kwargs)


def generate_desc_qc(original, final):
    """Function to generate an RBC-style QC CSV

    Parameters
    ----------
    original : str
        path to original image

    final : str
        path to final image

    coreg : nipype.interfaces.base.specs.DynamicTraitedSpec
        Output of OverlapInterface (with attributes Dice, Jaccard,
        CrossCorr, and Coverage)

    norm : nipype.interfaces.base.specs.DynamicTraitedSpec
        Output of OverlapInterface (with attributes Dice, Jaccard,
        CrossCorr, and Coverage)

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
        'original': nb.load(original),
        'final': nb.load(final)
    }

    # `sub` through `space`
    final_filename = final.split('/')[-1]
    bids_entities = final_filename.split('_')
    from_bids = dict(
        tuple(entity.split('-', 1)) if '-' in entity else
        ('suffix', entity) for entity in bids_entities
    )
    from_bids = {k: [from_bids[k]] for k in from_bids}
    if 'space' not in from_bids:
        from_bids['space'] = ['native']

    # `nVolCensored` & `nVolsRemoved`
    if images['original'].shape == images['final'].shape:
        shape_diff = 0
    else:
        shape_diff = 'qc log not yet implemented'  # TODO
    shape_params = {shape_key: [shape_diff] for
                    shape_key in {'nVolCensored', 'nVolsRemoved'}}

    # `meanFD (Jenkinson)`
    if isinstance(final, BufferedReader):
        final = final.name
    qc_filepath = _generate_filename(final)

    desc_span = re.search(r'_desc-.*_', final)
    if desc_span:
        desc_span = desc_span.span()
        final = '_'.join([
            final[:desc_span[0]],
            final[desc_span[1]:]
        ])
    del desc_span
    power_params = {'meanFD': np.mean(np.loadtxt('_'.join([
        *final.split('_')[:-1],
        'framewise-displacement-jenkinson.1D'
    ])))}

    # `relMeansRMSMotion` & `relMaxRMSMotion`
    mot = np.genfromtxt('_'.join([
        *final.split('_')[:-1],
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
    images = {variable: image.get_fdata().ravel() for
              variable, image in images.items()}
    intersect = images['original'] * images['final']
    vols = {variable: np.sum(image) for variable, image in images.items()}
    vol_intersect = np.sum(intersect)
    vol_sum = sum(vols.values())
    vol_union = vol_sum - vol_intersect
    overlap_params = {
        'coregDice': 2 * vol_intersect / vol_sum,
        'coregJaccard': vol_intersect / vol_union,
        'coregCrossCorr': np.corrcoef(
            images['original'],
            images['final'])[0, 1],
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
     'inputs': ['bold', 'desc-preproc_bold'],
     'outputs': ['xcp-qc']}
    """
    original = {}
    final = {}
    original['node'], original['out'] = strat_pool.get_data('bold')
    final['node'], final['out'] = strat_pool.get_data('desc-preproc_bold')

    qc_file = pe.Node(Function(input_names=['original', 'final'],
                               output_names=['qc_file'],
                               function=generate_desc_qc,
                               as_module=True),
                      name=f'xcpqc_{pipe_num}')

    wf.connect(original['node'], original['out'],
               qc_file, 'original')
    wf.connect(final['node'], final['out'],
               qc_file, 'final')

    outputs = {
        'xcp-qc': (qc_file, 'qc_file'),
    }

    return (wf, outputs)
