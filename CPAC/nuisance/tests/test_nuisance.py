from CPAC import nuisance
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.afni as afni
import pytest


def test_hello_world():

    print("hello world!")
    assert 0 == 0


def test_find_offending_time_points_fd_dvars_extend():

    find_censors = pe.Node(util.Function(input_names=['thresh_metric', 'out_file_path', 'fd_file_path', 'dvars_file_path',
                                                      'fd_threshold', 'dvars_threshold',
                                                      'number_of_previous_trs_to_remove',
                                                      'number_of_subsequent_trs_to_remove'],
                                         output_names=['out_file'],
                                         function=nuisance.find_offending_time_points),
                           name="find_censors_fd_dvars_extend")

    find_censors.inputs.thresh_metric = "FD+DVARS"
    find_censors.inputs.out_file_path = "censors_extend.1D"
    find_censors.inputs.fd_file_path = '/home/ccraddock/nuisance_test/fd.1D'
    find_censors.inputs.dvars_file_path = '/home/ccraddock/nuisance_test/dvars.1D'
    find_censors.inputs.fd_threshold = "1SD"
    find_censors.inputs.dvars_threshold = "1SD"
    find_censors.inputs.number_of_previous_trs_to_remove = 1
    find_censors.inputs.number_of_subsequent_trs_to_remove = 2
    find_censors.base_dir = '/home/ccraddock/nuisance_test/working_dir'

    res = find_censors.run()

    assert 0 == 0


def test_find_offending_time_points_fd():

    find_censors = pe.Node(util.Function(input_names=['thresh_metric', 'out_file_path', 'fd_file_path', 'dvars_file_path',
                                                      'fd_threshold', 'dvars_threshold',
                                                      'number_of_previous_trs_to_remove',
                                                      'number_of_subsequent_trs_to_remove'],
                                         output_names=['out_file'],
                                         function=nuisance.find_offending_time_points),
                           name="find_censors_fd")

    find_censors.inputs.thresh_metric = "FD"
    find_censors.inputs.out_file_path = "censors.1D"
    find_censors.inputs.fd_file_path = '/home/ccraddock/nuisance_test/fd.1D'
    find_censors.inputs.dvars_file_path = '/home/ccraddock/nuisance_test/dvars.1D'
    find_censors.inputs.fd_threshold = "1SD"
    find_censors.inputs.dvars_threshold = "1SD"
    find_censors.inputs.number_of_previous_trs_to_remove = 0
    find_censors.inputs.number_of_subsequent_trs_to_remove = 0
    find_censors.base_dir = '/home/ccraddock/nuisance_test/working_dir'

    res = find_censors.run()

    assert 0 == 0


def test_find_offending_time_points_dvars():
    find_censors = pe.Node(util.Function(input_names=['thresh_metric', 'out_file_path', 'fd_file_path', 'dvars_file_path',
                                                      'fd_threshold', 'dvars_threshold',
                                                      'number_of_previous_trs_to_remove',
                                                      'number_of_subsequent_trs_to_remove'],
                                         output_names=['out_file'],
                                         function=nuisance.find_offending_time_points),
                           name="find_censors_dvars")

    find_censors.inputs.thresh_metric = "DVARS"
    find_censors.inputs.out_file_path = "censors.1D"
    find_censors.inputs.fd_file_path = '/home/ccraddock/nuisance_test/fd.1D'
    find_censors.inputs.dvars_file_path = '/home/ccraddock/nuisance_test/dvars.1D'
    find_censors.inputs.fd_threshold = "1SD"
    find_censors.inputs.dvars_threshold = "1SD"
    find_censors.inputs.number_of_previous_trs_to_remove = 0
    find_censors.inputs.number_of_subsequent_trs_to_remove = 0
    find_censors.base_dir = '/home/ccraddock/nuisance_test/working_dir'

    res = find_censors.run()

    assert 0 == 0


def test_find_offending_time_points_fd_dvars():

    find_censors = pe.Node(util.Function(input_names=['thresh_metric', 'out_file_path', 'fd_file_path', 'dvars_file_path',
                                                      'fd_threshold', 'dvars_threshold',
                                                      'number_of_previous_trs_to_remove',
                                                      'number_of_subsequent_trs_to_remove'],
                                         output_names=['out_file'],
                                         function=nuisance.find_offending_time_points),
                           name="find_censors_fd_dvars")

    find_censors.inputs.thresh_metric = "FD+DVARS"
    find_censors.inputs.out_file_path = "censors.1D"
    find_censors.inputs.fd_file_path = '/home/ccraddock/nuisance_test/fd.1D'
    find_censors.inputs.dvars_file_path = '/home/ccraddock/nuisance_test/dvars.1D'
    find_censors.inputs.fd_threshold = "1SD"
    find_censors.inputs.dvars_threshold = "1SD"
    find_censors.inputs.number_of_previous_trs_to_remove = 0
    find_censors.inputs.number_of_subsequent_trs_to_remove = 0
    find_censors.base_dir = '/home/ccraddock/nuisance_test/working_dir'

    res = find_censors.run()

    assert 0 == 0


@pytest.mark.skip(reason="too slow")
def test_localstat():
    """
        3dLocalstat -prefix __WMeLOCAL_r${r} -nbhd 'SPHERE('${r}')' \
                -stat mean -mask  __mask_WMe${view} \
                -use_nonmask ${fn_epi}
    """
    radius = 2

    create_anaticor_regressor = pe.Node(interface=nuisance.Localstat(), name='create_anaticor_regressor')
    create_anaticor_regressor.inputs.in_file = '/home/ccraddock/nuisance_test/data_2mm.nii.gz'
    create_anaticor_regressor.inputs.out_file = 'local_WMe.nii.gz'
    create_anaticor_regressor.inputs.mask = '/home/ccraddock/nuisance_test/wm_mask_2mm.nii.gz'
    create_anaticor_regressor.inputs.neighborhood = "'SPHERE({0})'".format(radius)
    create_anaticor_regressor.inputs.statistic = "mean"
    create_anaticor_regressor.inputs.use_nonmask = True
    create_anaticor_regressor.interface.num_threads = 8
    create_anaticor_regressor.base_dir = '/home/ccraddock/nuisance_test/working_dir'

    res = create_anaticor_regressor.run()

    assert 0 == 0


def test_erode():

    erode = afni.Calc()
    erode.inputs.args = '-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
    erode.inputs.in_file_a = '/home/ccraddock/nuisance_test/working_dir/nuisance/' \
                             'gm_anat_to_2mm_flirt_applyxfm/gm_mask_flirt.nii.gz'
    erode.inputs.expr = 'a*(1-amongst(0,b,c,d,e,f,g))'
    erode.inputs.outputtype = 'NIFTI_GZ'
    erode.inputs.out_file = 'gm_mask_2mm_eroded.nii.gz'
    print("Erosion command line: {0}".format(erode.cmdline))

    assert "a+i" in erode.cmdline


def test_tproject():

    nuisance_regression = pe.Node(interface=nuisance.Tproject(), name='nuisance_regression')
    nuisance_regression.inputs.in_file = '/home/ccraddock/nuisance_test/functional.nii.gz'
    nuisance_regression.inputs.out_file = 'residuals.nii.gz'
    nuisance_regression.inputs.censor_file = '/home/ccraddock/nuisance_test/censors.1D'
    nuisance_regression.inputs.censor_idx = [51, 52, 53]
    nuisance_regression.inputs.censor_mode = 'NTRP'
    # nuisance_regression.inputs.catenation_file = '/home/ccraddock/nuisance_test/catenation.1D'
    nuisance_regression.inputs.noblock = True
    nuisance_regression.inputs.orthogonalize_file = '/home/ccraddock/nuisance_test/nuisance.1D'
    nuisance_regression.inputs.orthogonalize_polynomial = 1
    nuisance_regression.inputs.orthogonalize_dataset = '/home/ccraddock/nuisance_test/regressors.nii.gz'
    # nuisance_regression.inputs.stopband = [0.01, 0.05]
    nuisance_regression.inputs.bandpass = [0.005, 0.1]
    nuisance_regression.inputs.tr = 3.0
    nuisance_regression.inputs.mask = '/home/ccraddock/nuisance_test/func_mask.nii.gz'
    nuisance_regression.inputs.automask = True
    nuisance_regression.inputs.blur = 6
    nuisance_regression.inputs.normalize = True
    nuisance_regression.inputs.verb = True
    nuisance_regression.interface.num_threads = 8

    nuisance_regression.base_dir = '/home/ccraddock/nuisance_test/working_dir'

    res = nuisance_regression.run()

    assert 0 == 0


def test_mask_summarize_time_course():

    outfile = nuisance.mask_summarize_time_course("/home/ccraddock/nuisance_test/functional.nii.gz",
                                                  "/home/ccraddock/nuisance_test/func_mask.nii.gz",
                                                  "/home/ccraddock/nuisance_test/extracted_mean.nii.gz",
                                                  method="DetrendNormMean")

    outfile = nuisance.mask_summarize_time_course("/home/ccraddock/nuisance_test/functional.nii.gz",
                                                  "/home/ccraddock/nuisance_test/func_mask.nii.gz",
                                                  "/home/ccraddock/nuisance_test/extracted_mean.nii.gz",
                                                  method="PCA", num_pcs=5)

    outfile = nuisance.mask_summarize_time_course("/home/ccraddock/nuisance_test/functional.nii.gz",
                                         "/home/ccraddock/nuisance_test/func_mask.nii.gz",
                                         "/home/ccraddock/nuisance_test/extracted_pca_multilabels.tsv",
                                         method="PCA", num_pcs=5, mask_label=[[1, 2, 3]])
    assert 0 == 0

#@pytest.mark.skip(reason="too slow")
def test_nuisance_workflow_type1():

    """
    test_selector = {'Anaticor' : None | {radius = <radius in mm>},
        'aCompCor' : None | {num_pcs = <number of components to retain>,
                            tissues = 'WM' | 'CSF' | 'WM+CSF',
                            include_delayed = True | False,
                            include_squared = True | False,
                            include_delayed_squared = True | False},
        'WhiteMatter' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                       num_pcs = <number of components to retain>,
                       include_delayed = True | False,
                       include_squared = True | False,
                       include_delayed_squared = True | False},
        'Ventricles' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                       num_pcs = <number of components to retain>,
                       include_delayed = True | False,
                       include_squared = True | False,
                       include_delayed_squared = True | False},
        'GreyMatter' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                       num_pcs = <number of components to retain>,
                       include_delayed = True | False,
                       include_squared = True | False,
                       include_delayed_squared = True | False},
        'GlobalSignal' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                           num_pcs = <number of components to retain>,
                           include_delayed = True | False,
                           include_squared = True | False,
                           include_delayed_squared = True | False},
        'Motion' : None | {include_delayed = True | False,
                           include_squared = True | False,
                           include_delayed_squared = True | False},
        'Censor' : None | { thresh_metric = 'RMSD','DVARS', or 'RMSD+DVARS',
                            threshold = <threshold to be applied to metric, if using
                              RMSD+DVARS, this should be a tuple (RMSD thresh, DVARS thresh)>,
                            number_of_previous_trs_to_remove = True | False,
                            number_of_subsequent_trs_to_remove = True | False,
                            method = 'Kill', 'Zero', 'Interpolate', 'SpikeRegression'},
        'PolyOrt' : None | { degree = <polynomial degree up to which will be removed, e.g. 2 means
                                       constant + linear + quadratic, practically that is probably,
                                       the most that will be need esp. if band pass filtering>},
        'Bandpass' : None | { bottom_frequency = <frequency in hertz of the highpass part of the pass
                                                  band, frequencies below this will be removed>,
                              top_frequency = <frequency in hertz of the lowpass part of the pass
                                               band, frequencies above this will be removed>},
        }

    """

    selector_test1 = {'Anaticor': None,
                      'aCompCor': {'num_pcs': 5,
                                   'tissues': 'WM',
                                   'include_delayed': False,
                                   'include_squared': False,
                                   'include_delayed_squared': False},
                      'WhiteMatter': None,
                      'Ventricles': {'summary_method': 'DetrendNormMean',
                                     'num_pcs': None,
                                     'include_delayed': False,
                                     'include_squared': False,
                                     'include_delayed_squared': False},
                      'GreyMatter': None,
                      'GlobalSignal': None,
                      'Motion': {'include_delayed': True,
                                 'include_squared': True,
                                 'include_delayed_squared': True},
                      'Censor': {'thresh_metric': 'FD',
                                 'fd_threshold': '1.5 SD',
                                 'number_of_previous_trs_to_remove': 0,
                                 'number_of_subsequent_trs_to_remove': 0,
                                 'censor_method': 'SpikeRegression'},
                      'PolyOrt': {'degree': 2},
                      'Bandpass': None
                      }

    nuisance_regression_workflow = nuisance.create_nuisance_workflow(use_ants=True, selector=selector_test1)

    nuisance_regression_workflow.inputs.inputspec.lat_ventricles_mask_file_path = \
        '/home/ccraddock/nuisance_test/MNI152_T1_2mm_VentricleMask.nii.gz'

    nuisance_regression_workflow.inputs.inputspec.fd_file_path = '/home/ccraddock/nuisance_test/fd.1D'
    nuisance_regression_workflow.inputs.inputspec.dvars_file_path = '/home/ccraddock/nuisance_test/dvars.1d'
    nuisance_regression_workflow.inputs.inputspec.functional_file_path = '/home/ccraddock/nuisance_test/functional.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.wm_mask_file_path = '/home/ccraddock/nuisance_test/wm_mask.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.csf_mask_file_path = '/home/ccraddock/nuisance_test/csf_mask.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.gm_mask_file_path = '/home/ccraddock/nuisance_test/gm_mask.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.brain_mask_file_path = '/home/ccraddock/nuisance_test/func_mask.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.anat_to_mni_initial_xfm_file_path = \
        '/home/ccraddock/nuisance_test/anat_to_mni_initial_xfm.mat'
    nuisance_regression_workflow.inputs.inputspec.anat_to_mni_rigid_xfm_file_path = \
        '/home/ccraddock/nuisance_test/anat_to_mni_rigid_xfm.mat'
    nuisance_regression_workflow.inputs.inputspec.anat_to_mni_affine_xfm_file_path = \
        '/home/ccraddock/nuisance_test/anat_to_mni_affine_xfm.mat'
    nuisance_regression_workflow.inputs.inputspec.func_to_anat_linear_xfm_file_path =\
        '/home/ccraddock/nuisance_test/func_to_anat_linear_xfm.mat'

    nuisance_regression_workflow.inputs.inputspec.motion_parameters_file_path = \
        '/home/ccraddock/nuisance_test/motion_parameters.1D'

    nuisance_regression_workflow.inputs.inputspec.selector = selector_test1

    nuisance_regression_workflow.inputs.inputspec.brain_template_file_path = \
        '/home/ccraddock/nuisance_test/MNI152_T1_2mm_brain.nii.gz'

    nuisance_regression_workflow.base_dir = '/home/ccraddock/nuisance_test/working_dir'

    retval = nuisance_regression_workflow.run()

    assert 0 == 0

if __name__ == "__main__":
    test_nuisance_workflow_type1()
