"""Generate eXtensible Connectivity Pipeline-style quality control files"""
import os
import re
from io import BufferedReader

import nibabel as nb
import numpy as np
import pandas as pd
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.interfaces.function import Function


def generate_desc_qc(original, final):
    """Function to generate an RBC-style QC CSV

    Parameters
    ----------
    original : str
        path to original image

    final : str
        path to final image

    Returns
    -------
    str
        path to XCP QC TSV
    """
    columns = (
        'sub,ses,task,run,desc,space,meanFD,relMeansRMSMotion,'
        'relMaxRMSMotion,meanDVInit,meanDVFinal,nVolCensored,nVolsRemoved,'
        'motionDVCorrInit,motionDVCorrFinal,coregDice,coregJaccard,'
        'coregCrossCorr,coregCoverag,normDice,normJaccard,normCrossCorr,'
        'normCoverage'.split(',')
    )

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
    if nb.load(original).shape == nb.load(final).shape:
        shape_diff = 0
    else:
        shape_diff = 'qc log not yet implemented'  # TODO
    shape_params = {shape_key: [shape_diff] for
                    shape_key in {'nVolCensored', 'nVolsRemoved'}}

    # `meanFD`
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
    power_params = {'meanFD': pd.read_csv(
        '_'.join([
            *final.split('_')[:-1],
            'power-params.txt'
        ])
    ).to_dict(orient='list').get('MeanFD_Power')}

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

    qc_dict = {
        **from_bids,
        **power_params,
        **rms_params,
        **shape_params,
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
     'switch': ['generate_xcp_qc_files'],
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
