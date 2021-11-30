"""Util to set nuisance output datatypes"""
from nipype.interfaces.afni.utils import Calc
from CPAC.pipeline import nipype_pipeline_engine as pe

def report_as_float(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    Node block to set the final datatype of 'desc-cleaned_bold' to
    float

    Node Block:
    {"name": "report_as_float",
     "config": "None",
     "switch": "None",
     "option_key": "None",
     "option_val": "None",
     "inputs": ["desc-cleaned_bold"],
     "outputs": ["desc-cleaned_bold"]}
    '''
    set_float = pe.Node(Calc(), name=f'report_as_float_{pipe_num}')

    set_float.inputs.expr = 'a+0'
    set_float.inputs.args = '-datum float'
    set_float.inputs.overwrite = True

    cleaned_bold = {}
    cleaned_bold['node'], cleaned_bold['out'] = strat_pool.get_data(
        'desc-cleaned_bold')

    wf.connect([
        (cleaned_bold['node'], set_float, [
            (cleaned_bold['out'], 'in_file_a'),
            (cleaned_bold['out'], 'out_file')
        ])
    ])

    outputs = {'desc-cleaned_bold': (set_float, 'out_file')}

    return (wf, outputs)
