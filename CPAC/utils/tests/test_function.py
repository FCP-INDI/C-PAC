from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.interfaces.function import Function
from CPAC.utils.datasource import get_rest


rest_dict = {
    'rest_acq-1_run-1': {
        'scan': 's3://fcp-indi/sub-1019436_ses-1_task-rest_bold.nii.gz'
    }
}

resource = 'scan'
scan = 'rest_acq-1_run-1'


def test_function():

    f = pe.Node(Function(input_names=['scan',
                                      'rest_dict',
                                      'resource'],
                         output_names=['file_path'],
                         function=get_rest,
                         as_module=True),
                name='get_rest')

    f.inputs.set(
        resource=resource,
        rest_dict=rest_dict,
        scan=scan
    )

    results = f.run()
    assert rest_dict['rest_acq-1_run-1']['scan'] == results.outputs.file_path


def test_function_str():

    f = pe.Node(Function(input_names=['scan',
                                      'rest_dict',
                                      'resource'],
                         output_names=['file_path'],
                         function=get_rest),
                name='get_rest')

    f.inputs.set(
        resource=resource,
        rest_dict=rest_dict,
        scan=scan
    )

    results = f.run()
    assert rest_dict['rest_acq-1_run-1']['scan'] == results.outputs.file_path
