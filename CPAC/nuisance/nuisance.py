import re
import os
import numpy as np
import nibabel as nb
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import CPAC

from nipype.interfaces import fsl
from nipype.interfaces import ants
from nipype.interfaces import c3
from nipype.interfaces import afni
from nipype.interfaces.afni import utils as afni_utils
from scipy.fftpack import fft, ifft
from CPAC import utils
from CPAC.utils.interfaces.function import Function
from CPAC.utils.interfaces.masktool import MaskTool
from CPAC.utils.interfaces.pc import PC
from CPAC.nuisance.utils import (
    find_offending_time_points,
    generate_summarize_tissue_mask,
    temporal_variance_mask)
from CPAC.nuisance.utils.compcor import (
    calc_compcor_components,
    cosine_filter,
    TR_string_to_float)
from CPAC.utils.datasource import check_for_s3
from .bandpass import (bandpass_voxels, afni_1dBandpass)


def gather_nuisance(functional_file_path,
                    selector,
                    grey_matter_summary_file_path=None,
                    white_matter_summary_file_path=None,
                    csf_summary_file_path=None,
                    acompcor_file_path=None,
                    tcompcor_file_path=None,
                    global_summary_file_path=None,
                    motion_parameters_file_path=None,
                    custom_file_paths=None,
                    censor_file_path=None):
    """
    Gathers the various nuisance regressors together into a single tab separated values file that is an appropriate for
    input into 3dTproject

    :param functional_file_path: path to file that the regressors are being calculated for, is used to calculate
           the length of the regressors for error checking and in particular for calculating spike regressors
    :param output_file_path: path to output TSV that will contain the various nuisance regressors as columns
    :param grey_matter_summary_file_path: path to TSV that includes summary of grey matter time courses, e.g. output of
        mask_summarize_time_course
    :param white_matter_summary_file_path: path to TSV that includes summary of white matter time courses, e.g. output
        of mask_summarize_time_course
    :param csf_summary_file_path: path to TSV that includes summary of csf time courses, e.g. output
        of mask_summarize_time_course
    :param acompcor_file_path: path to TSV that includes acompcor time courses, e.g. output
        of mask_summarize_time_course
    :param tcompcor_file_path: path to TSV that includes tcompcor time courses, e.g. output
        of mask_summarize_time_course
    :param global_summary_file_path: path to TSV that includes summary of global time courses, e.g. output
        of mask_summarize_time_course
    :param motion_parameters_file_path: path to TSV that includes motion parameters
    :param custom_file_paths: path to CSV/TSV files to use as regressors
    :param censor_file_path: path to TSV with a single column with 1's for indices that should be retained and 0's
              for indices that should be censored
    :return:
    """

    # Basic checks for the functional image
    if not functional_file_path or \
        (not functional_file_path.endswith(".nii") and
         not functional_file_path.endswith(".nii.gz")):

        raise ValueError("Invalid value for input_file ({}). Should be a nifti "
                         "file and should exist".format(functional_file_path))

    try:
        functional_image = nb.load(functional_file_path)
    except:
        raise ValueError("Invalid value for input_file ({}). Should be a nifti "
                         "file and should exist".format(functional_file_path))

    if len(functional_image.shape) < 4 or functional_image.shape[3] < 2:
        raise ValueError("Invalid input_file ({}). Expected 4D file."
                         .format(functional_file_path))

    regressor_length = functional_image.shape[3]

    selector = selector.selector

    if not isinstance(selector, dict):
        raise ValueError("Invalid type for selectors {0}, expecting dict"
                         .format(type(selector)))

    regressor_files = {
        'aCompCor': acompcor_file_path,
        'tCompCor': tcompcor_file_path,
        'GlobalSignal': global_summary_file_path,
        'GreyMatter': grey_matter_summary_file_path,
        'WhiteMatter': white_matter_summary_file_path,
        'CerebrospinalFluid': csf_summary_file_path,
        'Motion': motion_parameters_file_path,
    }

    regressors_order = [
        'Motion',
        'GlobalSignal',
        'aCompCor',
        'tCompCor',
        'CerebrospinalFluid',
        'WhiteMatter',
        'GreyMatter',
    ]

    motion_labels = ['RotY', 'RotX', 'RotZ', 'Y', 'X', 'Z']

    # Compile regressors into a matrix
    column_names = []
    nuisance_regressors = []

    for regressor_type in regressors_order:

        if regressor_type not in selector:
            continue

        regressor_file = regressor_files[regressor_type]

        regressor_selector = selector.get(regressor_type) or {}

        if 'summary' in regressor_selector:
            if type(regressor_selector['summary']) is str:
                regressor_selector['summary'] = {
                    'method': regressor_selector['summary'],
                }

        if not regressor_file or not os.path.isfile(regressor_file):
            raise ValueError("Regressor type {0} specified in selectors "
                             "but the corresponding file was not found!"
                             .format(regressor_type))

        try:
            regressors = np.loadtxt(regressor_file)
        except:
            print("Could not read regressor {0} from {1}."
                  .format(regressor_type, regressor_file))
            raise

        if regressors.shape[0] != regressor_length:
            raise ValueError("Number of time points in {0} ({1}) is "
                             "inconsistent with length of functional "
                             "file {2} ({3})"
                             .format(regressor_file,
                                     regressors.shape[0],
                                     functional_file_path,
                                     regressor_length))

        if regressor_type == "Motion":
            num_regressors = 6
        elif not regressor_selector.get('summary', {}).get('components'):
            num_regressors = 1
        else:
            num_regressors = regressor_selector['summary']['components']

        if len(regressors.shape) == 1:
            regressors = np.expand_dims(regressors, axis=1)

        regressors = regressors[:, 0:num_regressors]

        if regressors.shape[1] != num_regressors:
            raise ValueError("Expecting {0} regressors for {1}, but "
                             "found {2} in file {3}."
                             .format(num_regressors,
                                     regressor_type,
                                     regressors.shape[1],
                                     regressor_file))

        # Add in the regressors, making sure to also add in the column name
        for regressor_index in range(regressors.shape[1]):
            if regressor_type == "Motion":
                regressor_name = motion_labels[regressor_index]
            else:
                summary_method = regressor_selector['summary']
                if type(summary_method) is dict:
                    summary_method = summary_method['method']

                regressor_name = "{0}{1}{2}".format(regressor_type,
                                                    summary_method,
                                                    regressor_index)

            column_names.append(regressor_name)
            nuisance_regressors.append(regressors[:, regressor_index])

            if regressor_selector.get('include_delayed', False):
                column_names.append("{0}Delay".format(regressor_name))
                nuisance_regressors.append(
                    np.append([0.0], regressors[0:-1, regressor_index])
                )

            if regressor_selector.get('include_backdiff', False):
                column_names.append("{0}BackDiff".format(regressor_name))
                nuisance_regressors.append(
                    np.append([0.0], np.diff(regressors[:, regressor_index], n=1))
                )

            if regressor_selector.get('include_squared', False):
                column_names.append("{0}Sq".format(regressor_name))
                nuisance_regressors.append(
                    np.square(regressors[:, regressor_index])
                )

            if regressor_selector.get('include_delayed_squared', False):
                column_names.append("{0}DelaySq".format(regressor_name))
                nuisance_regressors.append(
                    np.square(
                        np.append([0.0], regressors[0:-1, regressor_index])
                    )
                )

            if regressor_selector.get('include_backdiff_squared', False):
                column_names.append("{0}BackDiffSq".format(regressor_name))
                nuisance_regressors.append(
                    np.square(
                        np.append([0.0], np.diff(regressors[:, regressor_index], n=1))
                    )
                )

    # Add custom regressors
    if custom_file_paths:
        for custom_file_path in custom_file_paths:

            try:
                custom_regressor = np.loadtxt(custom_file_path)
            except:
                raise ValueError("Could not read regressor {0} from {1}."
                                .format('Custom', custom_file_path))

            if (len(custom_regressor.shape) > 1 and custom_regressor.shape[1] > 1):
                raise ValueError(
                    "Invalid format for censor file {0}, should be a single "
                    "column containing 1s for volumes to keep and 0s for volumes "
                    "to censor.".format(custom_file_path)
                )

            column_names.append(custom_file_path)
            nuisance_regressors.append(custom_regressor)

    # Add spike regressors
    if selector.get('Censor', {}).get('method') == 'SpikeRegression':

        selector = selector['Censor']

        regressor_file = censor_file_path

        if not regressor_file:
            # ↓ This section is gross and temporary ↓
            num_thresh = len(selector['thresholds'])
            plural_s = '' if num_thresh == 1 else 's'
            thresh_list = [
                thresh.get('value') for thresh in selector['thresholds']
            ]
            print(f"{selector['method']} Censor "
                  "specified with "
                  f"{'no ' if num_thresh == 0 else ''}threshold"
                  f"{plural_s} {str(thresh_list)} in selectors but "
                  f" threshold was not reached.")
            # ↑ This section is gross and temporary ↑
            # All good to pass through if nothing to censor
            censor_volumes = np.ones((regressor_length,), dtype=int)
        else:
            try:
                censor_volumes = np.loadtxt(regressor_file)
            except:
                raise ValueError("Could not read regressor {0} from {1}."
                                 .format(regressor_type, regressor_file))

        if (len(censor_volumes.shape) > 1 and censor_volumes.shape[1] > 1) or \
           not np.all(np.isin(censor_volumes, [0, 1])):

            raise ValueError(
                "Invalid format for censor file {0}, should be a single "
                "column containing 1s for volumes to keep and 0s for volumes "
                "to censor.".format(regressor_file)
            )

        censor_volumes = censor_volumes.flatten()
        censor_indices = np.where(censor_volumes == 0)[0]

        out_of_range_censors = censor_indices >= regressor_length
        if np.any(out_of_range_censors):
            raise ValueError(
                "Censor volumes {0} are out of range"
                "on censor file {1}, calculated "
                "regressor length is {2}".format(
                    censor_indices[out_of_range_censors],
                    regressor_file,
                    regressor_length
                )
            )

        if len(censor_indices) > 0:

            # if number_of_previous_trs_to_censor and number_of_subsequent_trs_to_censor
            # are not set, assume they should be zero
            previous_trs_to_censor = \
                selector.get('number_of_previous_trs_to_censor', 0)

            subsequent_trs_to_censor = \
                selector.get('number_of_subsequent_trs_to_censor', 0)

            spike_regressors = np.zeros(regressor_length)

            for censor_index in censor_indices:

                censor_begin_index = censor_index - previous_trs_to_censor
                if censor_begin_index < 0:
                    censor_begin_index = 0

                censor_end_index = censor_index + subsequent_trs_to_censor
                if censor_end_index >= regressor_length:
                    censor_end_index = regressor_length - 1

                spike_regressors[censor_begin_index:censor_end_index + 1] = 1

            for censor_index in np.where(spike_regressors == 1)[0]:

                column_names.append("SpikeRegression{0}".format(censor_index))
                spike_regressor_index = np.zeros(regressor_length)
                spike_regressor_index[censor_index] = 1
                nuisance_regressors.append(spike_regressor_index.flatten())

    if len(nuisance_regressors) == 0:
        return None

    # Compile columns into regressor file
    output_file_path = os.path.join(os.getcwd(), "nuisance_regressors.1D")

    with open(output_file_path, "w") as ofd:

        # write out the header information
        ofd.write("# C-PAC {0}\n".format(CPAC.__version__))
        ofd.write("# Nuisance regressors:\n")
        ofd.write("# " + "\t".join(column_names) + "\n")

        nuisance_regressors = np.array(nuisance_regressors)
        np.savetxt(ofd, nuisance_regressors.T, fmt='%.18f', delimiter='\t')

    return output_file_path


def create_regressor_workflow(nuisance_selectors,
                              use_ants,
                              ventricle_mask_exist,
                              name='nuisance_regressors'):
    """
    Workflow for the removal of various signals considered to be noise from resting state
    fMRI data.  The residual signals for linear regression denoising is performed in a single
    model.  Therefore the residual time-series will be orthogonal to all signals.

    Parameters
    ----------
    :param nuisance_selectors: dictionary describing nuisance regression to be performed
    :param use_ants: flag indicating whether FNIRT or ANTS is used
    :param name: Name of the workflow, defaults to 'nuisance'
    :return: nuisance : nipype.pipeline.engine.Workflow
        Nuisance workflow.

    Notes
    -----

    Workflow Inputs
    ---------------
    Workflow Inputs::

        inputspec.functional_file_path : string (nifti file)
            Path to realigned and motion corrected functional image (nifti) file.

        inputspec.functional_brain_mask_file_path : string (nifti file)
            Whole brain mask corresponding to the functional data.

        inputspec.anatomical_file_path : string (nifti file)
            Corresponding preprocessed anatomical.
        inputspec.wm_mask_file_path : string (nifti file)
            Corresponding white matter mask.
        inputspec.csf_mask_file_path : string (nifti file)
            Corresponding cerebral spinal fluid mask.
        inputspec.gm_mask_file_path : string (nifti file)
            Corresponding grey matter mask.
        inputspec.lat_ventricles_mask_file_path : string (nifti file)
            Mask of lateral ventricles calculated from the Harvard Oxford Atlas.

        inputspec.mni_to_anat_linear_xfm_file_path: string (nifti file)
            FLIRT Linear MNI to Anat transform
        inputspec.anat_to_mni_initial_xfm_file_path: string (nifti file)
            ANTS initial transform from anat to MNI
        inputspec.anat_to_mni_rigid_xfm_file_path: string (nifti file)
            ANTS rigid (6 parameter, no scaling) transform from anat to MNI
        inputspec.anat_to_mni_affine_xfm_file_path: string (nifti file)
            ANTS affine (13 parameter, scales and shears) transform from anat to MNI

        inputspec.func_to_anat_linear_xfm_file_path: string (nifti file)
            FLIRT Linear Transform between functional and anatomical spaces

        inputspec.motion_parameter_file_path : string (text file)
            Corresponding rigid-body motion parameters. Matrix in the file should be of shape
            (`T`, `R`), `T` time points and `R` motion parameters.
        inputspec.fd_j_file_path : string (text file)
            Framewise displacement calculated from the volume alignment.
        inputspec.fd_p_file_path : string (text file)
            Framewise displacement calculated from the motion parameters.
        inputspec.dvars_file_path : string (text file)
            DVARS calculated from the functional data.

        inputspec.selector : Dictionary containing configuration parameters for nuisance regression.
            To not run a type of nuisance regression, it may be ommited from the dictionary.

            selector = {
                aCompCor: {
                    summary: {
                        filter: 'cosine', Principal components are estimated after using a discrete cosine filter with 128s cut-off,
                            Leave filter field blank, if selected aCompcor method is 'DetrendPC'
                        method: 'DetrendPC', aCompCor will extract the principal components from
                            detrended tissues signal,
                        components: number of components to retain,
                    },
                    tissues: list of tissues to extract regressors.
                        Valid values are: 'WhiteMatter', 'CerebrospinalFluid',
                    extraction_resolution: None | floating point value indicating isotropic
                        resolution (ex. 2 for 2mm x 2mm x 2mm that data should be extracted at,
                        the corresponding tissue mask will be resampled to this resolution. The
                        functional data will also be resampled to this resolution, and the
                        extraction will occur at this new resolution. The goal is to avoid
                        contamination from undesired tissue components when extracting nuisance
                        regressors,
                    erode_mask: True | False, whether or not the mask should be eroded to
                        further avoid a mask overlapping with a different tissue class,
                    include_delayed: True | False, whether or not to include a one-frame delay regressor,
                        default to False,
                    include_squared: True | False, whether or not to include a squared regressor,
                        default to False,
                    include_delayed_squared: True | False, whether or not to include a squared one-frame
                        delay regressor, default to False,
                    include_backdiff: True | False, whether or not to include a one-lag difference,
                        default to False,
                    include_backdiff_squared: True | False, whether or not to include a squared one-frame
                        delay regressor, default to False,
                },
                tCompCor: {
                    summary: {
                        filter: 'cosine', Principal components are estimated after using a discrete cosine filter with 128s cut-off,
                            Leave filter field blank, if selected tCompcor method is 'DetrendPC'
                        method: 'DetrendPC', tCompCor will extract the principal components from
                            detrended tissues signal,
                        components: number of components to retain,
                    },
                    threshold:
                        floating point number = cutoff as raw variance value,
                        floating point number followed by SD (ex. 1.5SD) = mean + a multiple of the SD,
                        floating point number followed by PCT (ex. 2PCT) = percentile from the top (ex is top 2%),
                    by_slice: True | False, whether or not the threshold criterion should be applied
                        by slice or across the entire volume, makes most sense for thresholds
                        using SD or PCT,
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                    include_backdiff: True | False (same as for aCompCor),
                    include_backdiff_squared: True | False (same as for aCompCor),
                },
                WhiteMatter: {
                    summary: {
                        method: 'PC', 'DetrendPC', 'Mean', 'NormMean' or 'DetrendNormMean',
                        components: number of components to retain, if PC,
                    },
                    extraction_resolution: None | floating point value (same as for aCompCor),
                    erode_mask: True | False (same as for aCompCor),
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                    include_backdiff: True | False (same as for aCompCor),
                    include_backdiff_squared: True | False (same as for aCompCor),
                },
                CerebrospinalFluid: {
                    summary: {
                        method: 'PC', 'DetrendPC', 'Mean', 'NormMean' or 'DetrendNormMean',
                        components: number of components to retain, if PC,
                    },
                    extraction_resolution: None | floating point value (same as for aCompCor),
                    erode_mask: True | False (same as for aCompCor),
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                    include_backdiff: True | False (same as for aCompCor),
                    include_backdiff_squared: True | False (same as for aCompCor),
                },
                GreyMatter: {
                    summary: {
                        method: 'PC', 'DetrendPC', 'Mean', 'NormMean' or 'DetrendNormMean',
                        components: number of components to retain, if PC,
                    },
                    extraction_resolution: None | floating point value (same as for aCompCor),
                    erode_mask: True | False (same as for aCompCor),
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                    include_backdiff: True | False (same as for aCompCor),
                    include_backdiff_squared: True | False (same as for aCompCor),
                },
                GlobalSignal: {
                    summary: {
                        method: 'PC', 'DetrendPC', 'Mean', 'NormMean' or 'DetrendNormMean',
                        components: number of components to retain, if PC,
                    },
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                    include_backdiff: True | False (same as for aCompCor),
                    include_backdiff_squared: True | False (same as for aCompCor),
                },
                Motion: None | {
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                    include_backdiff: True | False (same as for aCompCor),
                    include_backdiff_squared: True | False (same as for aCompCor),
                },
                Censor: {
                    method: 'Kill', 'Zero', 'Interpolate', 'SpikeRegression',
                    thresholds: list of dictionary, {
                        type: 'FD_J', 'FD_P', 'DVARS',
                        value: threshold value to be applied to metric
                    },
                    number_of_previous_trs_to_censor: integer, number of previous
                        TRs to censor (remove or regress, if spike regression)
                    number_of_subsequent_trs_to_censor: integer, number of
                        subsequent TRs to censor (remove or regress, if spike
                        regression)
                },
                PolyOrt: {
                    degree: integer, polynomial degree up to which will be removed,
                        e.g. 2 means constant + linear + quadratic, practically
                        that is probably, the most that will be need especially
                        if band pass filtering
                },
                Bandpass: {
                    bottom_frequency: floating point value, frequency in hertz of
                        the highpass part of the pass band, frequencies below this
                        will be removed,
                    top_frequency: floating point value, frequency in hertz of the
                        lowpass part of the pass band, frequencies above this
                        will be removed
                },
                Custom: [
                    {
                        file: file containing the regressors. It can be a CSV file,
                            with one regressor per column, or a Nifti image, with
                            one regressor per voxel.
                        convolve: perform the convolution operation of the given
                            regressor with the timeseries.
                    }
                ]
            }

    Workflow Outputs::

        outputspec.residual_file_path : string (nifti file)
            Path of residual file in nifti format
        outputspec.regressors_file_path : string (TSV file)
            Path of TSV file of regressors used. Column name indicates the regressors included .

    Nuisance Procedure:

    1. Compute nuisance regressors based on input selections.
    2. Calculate residuals with respect to these nuisance regressors in a
       single model for every voxel.

    High Level Workflow Graph:

    .. exec::
        from CPAC.nuisance import create_regressor_workflow
        wf = create_regressor_workflow({
            'PolyOrt': {'degree': 2},
            'tCompCor': {'summary': {'method': 'PC', 'components': 5}, 'threshold': '1.5SD', 'by_slice': True},
            'aCompCor': {'summary': {'method': 'PC', 'components': 5}, 'tissues': ['WhiteMatter', 'CerebrospinalFluid'], 'extraction_resolution': 2},
            'WhiteMatter': {'summary': {'method': 'PC', 'components': 5}, 'extraction_resolution': 2},
            'CerebrospinalFluid': {'summary': {'method': 'PC', 'components': 5}, 'extraction_resolution': 2, 'erode_mask': True},
            'GreyMatter': {'summary': {'method': 'PC', 'components': 5}, 'extraction_resolution': 2, 'erode_mask': True},
            'GlobalSignal': {'summary': 'Mean', 'include_delayed': True, 'include_squared': True, 'include_delayed_squared': True},
            'Motion': {'include_delayed': True, 'include_squared': True, 'include_delayed_squared': True},
            'Censor': {'method': 'Interpolate', 'thresholds': [{'type': 'FD_J', 'value': 0.5}, {'type': 'DVARS', 'value': 0.7}]}
        }, use_ants=False)

        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/nuisance.dot'
        )

    .. image:: ../../images/generated/nuisance.png
       :width: 1000

    Detailed Workflow Graph:

    .. image:: ../../images/generated/nuisance_detailed.png
       :width: 1000

    """

    nuisance_wf = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=[
        'selector',
        'functional_file_path',
        'anatomical_file_path',
        'anatomical_eroded_brain_mask_file_path',
        'gm_mask_file_path',
        'wm_mask_file_path',
        'csf_mask_file_path',
        'lat_ventricles_mask_file_path',
        'functional_brain_mask_file_path',
        'func_to_anat_linear_xfm_file_path',
        'anat_to_func_linear_xfm_file_path',
        'mni_to_anat_linear_xfm_file_path',
        'anat_to_mni_initial_xfm_file_path',
        'anat_to_mni_rigid_xfm_file_path',
        'anat_to_mni_affine_xfm_file_path',
        'motion_parameters_file_path',
        'fd_j_file_path',
        'fd_p_file_path',
        'dvars_file_path',
        'creds_path',
        'dl_dir',
        'tr',
    ]), name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['regressors_file_path']),
                         name='outputspec')

    functional_mean = pe.Node(interface=afni_utils.TStat(),
                        name='functional_mean')

    functional_mean.inputs.options = '-mean'
    functional_mean.inputs.outputtype = 'NIFTI_GZ'

    nuisance_wf.connect(inputspec, 'functional_file_path',
                        functional_mean, 'in_file')

    # Resources to create regressors
    pipeline_resource_pool = {
        "Anatomical": (inputspec, 'anatomical_file_path'),
        "AnatomicalErodedMask": (inputspec, 'anatomical_eroded_brain_mask_file_path'),
        "Functional": (inputspec, 'functional_file_path'),
        "Functional_mean" : (functional_mean, 'out_file'),
        "GlobalSignal": (inputspec, 'functional_brain_mask_file_path'),
        "WhiteMatter": (inputspec, 'wm_mask_file_path'),
        "CerebrospinalFluid": (inputspec, 'csf_mask_file_path'),
        "GreyMatter": (inputspec, 'gm_mask_file_path'),
        "Ventricles": (inputspec, 'lat_ventricles_mask_file_path'),

        "Transformations": {
            "func_to_anat_linear_xfm": (inputspec, "func_to_anat_linear_xfm_file_path"),
            "anat_to_func_linear_xfm": (inputspec, "anat_to_func_linear_xfm_file_path"),
            "mni_to_anat_linear_xfm": (inputspec, "mni_to_anat_linear_xfm_file_path"),
            "anat_to_mni_initial_xfm": (inputspec, "anat_to_mni_initial_xfm_file_path"),
            "anat_to_mni_rigid_xfm": (inputspec, "anat_to_mni_rigid_xfm_file_path"),
            "anat_to_mni_affine_xfm": (inputspec, "anat_to_mni_affine_xfm_file_path"),
        },

        "DVARS": (inputspec, 'dvars_file_path'),
        "FD_J": (inputspec, 'framewise_displacement_j_file_path'),
        "FD_P": (inputspec, 'framewise_displacement_p_file_path'),
        "Motion": (inputspec, 'motion_parameters_file_path'),
    }

    # Regressor map to simplify construction of the needed regressors
    regressors = {
        'GreyMatter': ['grey_matter_summary_file_path', (), 'ort'],
        'WhiteMatter': ['white_matter_summary_file_path', (), 'ort'],
        'CerebrospinalFluid': ['csf_summary_file_path', (), 'ort'],
        'aCompCor': ['acompcor_file_path', (), 'ort'],
        'tCompCor': ['tcompcor_file_path', (), 'ort'],
        'GlobalSignal': ['global_summary_file_path', (), 'ort'],
        'Custom': ['custom_file_paths', (), 'ort'],
        'VoxelCustom': ['custom_file_paths', (), 'dsort'],
        'DVARS': ['dvars_file_path', (), 'ort'],
        'FD_J': ['framewise_displacement_j_file_path', (), 'ort'],
        'FD_P': ['framewise_displacement_p_file_path', (), 'ort'],
        'Motion': ['motion_parameters_file_path', (), 'ort']
    }

    motion = ['DVARS', 'FD_J', 'FD_P', 'Motion']
    derived = ['tCompCor', 'aCompCor']
    tissues = ['GreyMatter', 'WhiteMatter', 'CerebrospinalFluid']

    for regressor_type, regressor_resource in regressors.items():

        if regressor_type not in nuisance_selectors:
            continue

        regressor_selector = nuisance_selectors[regressor_type]

        if regressor_type == 'Custom':

            custom_ort_check_s3_nodes = []
            custom_dsort_check_s3_nodes = []
            custom_dsort_convolve_nodes = []

            for custom_regressor in sorted(regressor_selector, key=lambda c: c['file']):
                custom_regressor_file = custom_regressor['file']

                custom_check_s3_node = pe.Node(Function(
                    input_names=[
                        'file_path',
                        'creds_path',
                        'dl_dir',
                        'img_type'
                    ],
                    output_names=[
                        'local_path'
                    ],
                    function=check_for_s3,
                    as_module=True),
                    name='custom_check_for_s3_%s' % name)

                custom_check_s3_node.inputs.set(
                    file_path=custom_regressor_file,
                    img_type='func'
                )

                if custom_regressor_file.endswith('.nii.gz') or \
                    custom_regressor_file.endswith('.nii'):

                    if custom_regressor.get('convolve'):
                        custom_dsort_convolve_nodes += [custom_check_s3_node]
                    else:
                        custom_dsort_check_s3_nodes += [custom_check_s3_node]

                else:
                    custom_ort_check_s3_nodes += [custom_check_s3_node]

            if len(custom_ort_check_s3_nodes) > 0:
                custom_ort_merge = pe.Node(
                    util.Merge(len(custom_ort_check_s3_nodes)),
                    name='custom_ort_merge'
                )

                for i, custom_check_s3_node in enumerate(custom_ort_check_s3_nodes):
                    nuisance_wf.connect(
                        custom_check_s3_node, 'local_path',
                        custom_ort_merge, "in{}".format(i + 1)
                    )

                pipeline_resource_pool['custom_ort_file_paths'] = \
                    (custom_ort_merge, 'out')

                regressors['Custom'][1] = \
                    pipeline_resource_pool['custom_ort_file_paths']

            if len(custom_dsort_convolve_nodes) > 0:
                custom_dsort_convolve_merge = pe.Node(
                    util.Merge(len(custom_dsort_convolve_nodes)),
                    name='custom_dsort_convolve_merge'
                )

                for i, custom_check_s3_node in enumerate(custom_dsort_convolve_nodes):
                    nuisance_wf.connect(
                        custom_check_s3_node, 'local_path',
                        custom_dsort_convolve_merge, "in{}".format(i + 1)
                    )

            if len(custom_dsort_check_s3_nodes) > 0:

                images_to_merge = len(custom_dsort_check_s3_nodes)
                if len(custom_dsort_convolve_nodes) > 0:
                    images_to_merge += 1

                custom_dsort_merge = pe.Node(
                    util.Merge(images_to_merge),
                    name='custom_dsort_merge'
                )

                for i, custom_check_s3_node in enumerate(custom_dsort_check_s3_nodes):
                    nuisance_wf.connect(
                        custom_check_s3_node, 'local_path',
                        custom_dsort_merge, "in{}".format(i + 1)
                    )

                if len(custom_dsort_convolve_nodes) > 0:
                    nuisance_wf.connect(
                        custom_dsort_convolve_merge, 'out',
                        custom_dsort_merge, "in{}".format(i + 1)
                    )

                pipeline_resource_pool['custom_dsort_file_paths'] = \
                    (custom_dsort_merge, 'out')

                regressors['VoxelCustom'][1] = \
                    pipeline_resource_pool['custom_dsort_file_paths']

            continue

        if regressor_type in motion:
            regressor_resource[1] = \
                pipeline_resource_pool[regressor_type]
            continue

        # Set summary method for tCompCor and aCompCor
        if regressor_type in derived:

            if 'summary' not in regressor_selector:
                regressor_selector['summary'] = {}

            if type(regressor_selector['summary']) is not dict:
                raise ValueError("Regressor {0} requires PC summary method, "
                                 "but {1} specified"
                                 .format(regressor_type,
                                         regressor_selector['summary']))

            if not regressor_selector['summary'].get('components'):
                regressor_selector['summary']['components'] = 1


        # If regressor is not present, build up the regressor
        if not regressor_resource[1]:

            # We don't have the regressor, look for it in the resource pool,
            # build a corresponding key, this is seperated in to a mask key
            # and an extraction key, which when concatenated provide the
            # resource key for the regressor
            regressor_descriptor = {'tissue': regressor_type}

            if regressor_type == 'aCompCor':
                if not regressor_selector.get('tissues'):
                    raise ValueError("Tissue type required for aCompCor, "
                                     "but none specified")

                regressor_descriptor = {
                    'tissue': regressor_selector['tissues']
                }


            if regressor_type == 'tCompCor':
                if not regressor_selector.get('threshold'):
                    raise ValueError("Threshold required for tCompCor, "
                                     "but none specified.")

                regressor_descriptor = {
                    'tissue': 'FunctionalVariance-{}'
                              .format(regressor_selector['threshold'])
                }

                if regressor_selector.get('by_slice'):
                    regressor_descriptor['tissue'] += '-BySlice'
                else:
                    regressor_selector['by_slice'] = False

                if regressor_selector.get('erode_mask_mm'):
                    erosion_mm = regressor_selector['erode_mask_mm']
                else:
                    erosion_mm = False

                if regressor_selector.get('degree'):
                    degree = regressor_selector['degree']
                else:
                    degree = 1

                temporal_wf = temporal_variance_mask(regressor_selector['threshold'],
                                                     by_slice=regressor_selector['by_slice'],
                                                     erosion=erosion_mm,
                                                     degree=degree)

                nuisance_wf.connect(*(pipeline_resource_pool['Functional'] + (temporal_wf, 'inputspec.functional_file_path')))

                if erosion_mm: # TODO: in func/anat space
                    # transform eroded anat brain mask to functional space
                    # convert_xfm
                    anat_to_func_linear_xfm = pe.Node(interface=fsl.ConvertXFM(), name='anat_to_func_linear_xfm')
                    anat_to_func_linear_xfm.inputs.invert_xfm = True
                    nuisance_wf.connect(*(pipeline_resource_pool['Transformations']['func_to_anat_linear_xfm'] + (anat_to_func_linear_xfm, 'in_file')))

                    # flirt
                    anat_to_func_mask = pe.Node(interface=fsl.FLIRT(), name='Functional_eroded_mask')
                    anat_to_func_mask.inputs.output_type = 'NIFTI_GZ'
                    anat_to_func_mask.inputs.apply_xfm = True
                    anat_to_func_mask.inputs.interp = 'nearestneighbour'
                    nuisance_wf.connect(anat_to_func_linear_xfm, 'out_file', anat_to_func_mask, 'in_matrix_file')
                    nuisance_wf.connect(*(pipeline_resource_pool['AnatomicalErodedMask'] + (anat_to_func_mask, 'in_file')))
                    nuisance_wf.connect(*(pipeline_resource_pool['GlobalSignal'] + (anat_to_func_mask, 'reference')))

                    # connect workflow
                    nuisance_wf.connect(anat_to_func_mask, 'out_file', temporal_wf, 'inputspec.mask_file_path')
                else:
                    nuisance_wf.connect(*(pipeline_resource_pool['GlobalSignal'] + (temporal_wf, 'inputspec.mask_file_path')))

                pipeline_resource_pool[regressor_descriptor['tissue']] = \
                    (temporal_wf, 'outputspec.mask')

            if type(regressor_selector['summary']) is not dict:
                regressor_selector['summary'] = {
                    "filter": regressor_selector['summary'],
                    "method": regressor_selector['summary']
                }

            # Add selector into regressor description

            if regressor_selector.get('extraction_resolution'):
                regressor_descriptor['resolution'] = \
                    str(regressor_selector['extraction_resolution']) + "mm"

            elif regressor_type in tissues:
                regressor_selector['extraction_resolution'] = "Functional"
                regressor_descriptor['resolution'] = "Functional"

            if regressor_selector.get('erode_mask'):
                regressor_descriptor['erosion'] = 'Eroded'

            if not regressor_selector.get('summary'):
                raise ValueError("Summary method required for {0}, "
                                 "but none specified".format(regressor_type))

            regressor_descriptor['extraction'] = \
                regressor_selector['summary']['method']

            if regressor_descriptor['extraction'] in ['DetrendPC', 'PC']:
                if not regressor_selector['summary'].get('components'):
                    raise ValueError("Summary method PC requires components, "
                                     "but received none.")

                regressor_descriptor['extraction'] += \
                    '_{0}'.format(regressor_selector['summary']['components'])

            if type(regressor_descriptor['tissue']) is not list:
                regressor_descriptor['tissue'] = \
                    [regressor_descriptor['tissue']]


            if regressor_selector.get('extraction_resolution') and \
                    regressor_selector["extraction_resolution"] != "Functional":

                functional_at_resolution_key = "Functional_{0}mm".format(
                    regressor_selector["extraction_resolution"]
                )

                anatomical_at_resolution_key = "Anatomical_{0}mm".format(
                    regressor_selector["extraction_resolution"]
                )

                if anatomical_at_resolution_key not in pipeline_resource_pool:

                    anat_resample = pe.Node(
                        interface=fsl.FLIRT(),
                        name='{}_flirt'
                             .format(anatomical_at_resolution_key)
                    )
                    anat_resample.inputs.apply_isoxfm = regressor_selector["extraction_resolution"]

                    nuisance_wf.connect(*(
                        pipeline_resource_pool['Anatomical'] +
                        (anat_resample, 'in_file')
                    ))

                    nuisance_wf.connect(*(
                        pipeline_resource_pool['Anatomical'] +
                        (anat_resample, 'reference')
                    ))

                    pipeline_resource_pool[anatomical_at_resolution_key] = \
                        (anat_resample, 'out_file')

                if functional_at_resolution_key not in pipeline_resource_pool:

                    func_resample = pe.Node(
                        interface=fsl.FLIRT(),
                        name='{}_flirt'
                             .format(functional_at_resolution_key)
                    )
                    func_resample.inputs.apply_xfm = True

                    nuisance_wf.connect(*(
                        pipeline_resource_pool['Transformations']['func_to_anat_linear_xfm'] +
                        (func_resample, 'in_matrix_file')
                    ))

                    nuisance_wf.connect(*(
                        pipeline_resource_pool['Functional'] +
                        (func_resample, 'in_file')
                    ))

                    nuisance_wf.connect(*(
                        pipeline_resource_pool[anatomical_at_resolution_key] +
                        (func_resample, 'reference')
                    ))

                    pipeline_resource_pool[functional_at_resolution_key] = \
                        (func_resample, 'out_file')

            # Create merger to summarize the functional timeseries
            regressor_mask_file_resource_keys = []
            for tissue in regressor_descriptor['tissue']:

                # Ignore non tissue masks
                if tissue not in tissues and \
                    not tissue.startswith('FunctionalVariance'):
                    regressor_mask_file_resource_keys += [tissue]
                    continue

                tissue_regressor_descriptor = regressor_descriptor.copy()
                tissue_regressor_descriptor['tissue'] = tissue

                # Generate resource masks
                (pipeline_resource_pool,
                 regressor_mask_file_resource_key) = \
                    generate_summarize_tissue_mask(
                        nuisance_wf,
                        pipeline_resource_pool,
                        tissue_regressor_descriptor,
                        regressor_selector,
                        use_ants=use_ants,
                        ventricle_mask_exist=ventricle_mask_exist
                    )

                regressor_mask_file_resource_keys += \
                    [regressor_mask_file_resource_key]

            # Keep tissus ordered, to avoid duplicates
            regressor_mask_file_resource_keys = \
                list(sorted(regressor_mask_file_resource_keys))

            # Create key for the final regressors
            regressor_file_resource_key = "_".join([
                "-".join(regressor_descriptor[key])
                if type(regressor_descriptor[key]) == list
                else regressor_descriptor[key]

                for key in ['tissue', 'resolution', 'erosion', 'extraction']
                if key in regressor_descriptor
            ])

            if regressor_file_resource_key not in pipeline_resource_pool:

                # Merge mask paths to extract voxel timeseries
                merge_masks_paths = pe.Node(
                    util.Merge(len(regressor_mask_file_resource_keys)),
                    name='{}_merge_masks'.format(regressor_type)
                )
                for i, regressor_mask_file_resource_key in \
                        enumerate(regressor_mask_file_resource_keys):

                    node, node_output = \
                        pipeline_resource_pool[regressor_mask_file_resource_key]

                    nuisance_wf.connect(
                        node, node_output,
                        merge_masks_paths, "in{}".format(i + 1)
                    )

                union_masks_paths = pe.Node(
                    MaskTool(outputtype='NIFTI_GZ'),
                    name='{}_union_masks'.format(regressor_type)
                )

                nuisance_wf.connect(
                    merge_masks_paths, 'out',
                    union_masks_paths, 'in_files'
                )

                functional_key = 'Functional'
                if regressor_selector.get('extraction_resolution') and \
                        regressor_selector["extraction_resolution"] != "Functional":

                    functional_key = 'Functional_{}mm'.format(
                        regressor_selector['extraction_resolution']
                    )

                summary_filter = regressor_selector['summary'].get('filter', '')
                summary_filter_input = pipeline_resource_pool[functional_key]

                summary_method = regressor_selector['summary']['method']
                summary_method_input = pipeline_resource_pool[functional_key]

                if 'DetrendPC' in summary_method:

                    compcor_imports = ['import os',
                                       'import scipy.signal as signal',
                                       'import nibabel as nb',
                                       'import numpy as np',
                                       'from CPAC.utils import safe_shape']

                    compcor_node = pe.Node(Function(input_names=['data_filename',
                                                                 'num_components',
                                                                 'mask_filename'],
                                                    output_names=[
                                                        'compcor_file'],
                                                    function=calc_compcor_components,
                                                    imports=compcor_imports),
                                           name='{}_DetrendPC'.format(regressor_type), mem_gb=2.0)

                    compcor_node.inputs.num_components = regressor_selector['summary']['components']

                    nuisance_wf.connect(
                        summary_method_input[0], summary_method_input[1],
                        compcor_node, 'data_filename'
                    )

                    nuisance_wf.connect(
                        union_masks_paths, 'out_file',
                        compcor_node, 'mask_filename'
                    )

                    summary_method_input = (compcor_node, 'compcor_file')

                else:
                    if 'cosine' in summary_filter:
                        cosfilter_imports = ['import os',
                                             'import numpy as np',
                                             'import nibabel as nb',
                                             'from nipype import logging']

                        cosfilter_node = pe.Node(util.Function(input_names=['input_image_path',
                                                                            'timestep'],
                                                               output_names=[
                                                                   'cosfiltered_img'],
                                                               function=cosine_filter,
                                                               imports=cosfilter_imports),
                                                 name='{}_cosine_filter'.format(regressor_type))
                        nuisance_wf.connect(
                            summary_filter_input[0], summary_filter_input[1],
                            cosfilter_node, 'input_image_path'
                        )
                        tr_string2float_node = pe.Node(util.Function(input_names=['tr'],
                                                                     output_names=[
                                                                         'tr_float'],
                                                                     function=TR_string_to_float),
                                                       name='{}_tr_string2float'.format(regressor_type))

                        nuisance_wf.connect(
                            inputspec, 'tr',
                            tr_string2float_node, 'tr'
                        )

                        nuisance_wf.connect(
                            tr_string2float_node, 'tr_float',
                            cosfilter_node, 'timestep'
                        )

                        summary_method_input = (
                            cosfilter_node, 'cosfiltered_img')

                    if 'Detrend' in summary_method:

                        detrend_node = pe.Node(
                            afni.Detrend(args='-polort 1', outputtype='NIFTI'),
                            name='{}_detrend'.format(regressor_type)
                        )

                        nuisance_wf.connect(
                            summary_method_input[0], summary_method_input[1],
                            detrend_node, 'in_file'
                        )

                        summary_method_input = (detrend_node, 'out_file')

                    if 'Norm' in summary_method:

                        l2norm_node = pe.Node(
                            afni.TStat(args='-l2norm', outputtype='NIFTI'),
                            name='{}_l2norm'.format(regressor_type)
                        )
                        nuisance_wf.connect(
                            summary_method_input[0], summary_method_input[1],
                            l2norm_node, 'in_file'
                        )
                        nuisance_wf.connect(
                            union_masks_paths, 'out_file',
                            l2norm_node, 'mask'
                        )

                        norm_node = pe.Node(
                            afni.Calc(expr='a/b', outputtype='NIFTI'),
                            name='{}_norm'.format(regressor_type)
                        )
                        nuisance_wf.connect(
                            summary_method_input[0], summary_method_input[1],
                            norm_node, 'in_file_a'
                        )
                        nuisance_wf.connect(
                            l2norm_node, 'out_file',
                            norm_node, 'in_file_b'
                        )

                        summary_method_input = (norm_node, 'out_file')

                    if 'Mean' in summary_method:

                        mean_node = pe.Node(
                            afni.ROIStats(quiet=False, args='-1Dformat'),
                            name='{}_mean'.format(regressor_type)
                        )
                        nuisance_wf.connect(
                            summary_method_input[0], summary_method_input[1],
                            mean_node, 'in_file'
                        )

                        nuisance_wf.connect(
                            union_masks_paths, 'out_file',
                            mean_node, 'mask_file'
                        )

                        summary_method_input = (mean_node, 'out_file')

                    if 'PC' in summary_method:

                        std_node = pe.Node(
                            afni.TStat(args='-nzstdev', outputtype='NIFTI'),
                            name='{}_std'.format(regressor_type)
                        )
                        nuisance_wf.connect(
                            summary_method_input[0], summary_method_input[1],
                            std_node, 'in_file'
                        )
                        nuisance_wf.connect(
                            union_masks_paths, 'out_file',
                            std_node, 'mask'
                        )

                        standardized_node = pe.Node(
                            afni.Calc(expr='a/b', outputtype='NIFTI'),
                            name='{}_standardized'.format(regressor_type)
                        )
                        nuisance_wf.connect(
                            summary_method_input[0], summary_method_input[1],
                            standardized_node, 'in_file_a'
                        )
                        nuisance_wf.connect(
                            std_node, 'out_file',
                            standardized_node, 'in_file_b'
                        )

                        pc_node = pe.Node(
                            PC(args='-vmean -nscale', pcs=regressor_selector['summary']['components'], outputtype='NIFTI_GZ'),
                            name='{}_pc'.format(regressor_type)
                        )

                        nuisance_wf.connect(
                            standardized_node, 'out_file',
                            pc_node, 'in_file'
                        )
                        nuisance_wf.connect(
                            union_masks_paths, 'out_file',
                            pc_node, 'mask'
                        )

                        summary_method_input = (pc_node, 'pcs_file')

                pipeline_resource_pool[regressor_file_resource_key] = \
                    summary_method_input

                # Add it to internal resource pool
                regressor_resource[1] = \
                    pipeline_resource_pool[regressor_file_resource_key]

    # Build regressors and combine them into a single file
    build_nuisance_regressors = pe.Node(Function(
        input_names=['functional_file_path',
                     'selector',
                     'grey_matter_summary_file_path',
                     'white_matter_summary_file_path',
                     'csf_summary_file_path',
                     'acompcor_file_path',
                     'tcompcor_file_path',
                     'global_summary_file_path',
                     'motion_parameters_file_path',
                     'custom_file_paths',
                     'censor_file_path'],
        output_names=['out_file'],
        function=gather_nuisance,
        as_module=True
    ), name="build_nuisance_regressors")

    nuisance_wf.connect(
        inputspec, 'functional_file_path',
        build_nuisance_regressors, 'functional_file_path'
    )

    nuisance_wf.connect(
        inputspec, 'selector',
        build_nuisance_regressors, 'selector'
    )

    # Check for any regressors to combine into files
    has_nuisance_regressors = any(
        regressor_node
        for regressor_key, (regressor_arg, regressor_node, regressor_target)
        in regressors.items()
        if regressor_target == 'ort'
    )

    if has_nuisance_regressors:
        for regressor_key, (
            regressor_arg, regressor_node, regressor_target
        ) in regressors.items():
            if regressor_target != 'ort':
                continue

            if regressor_key in nuisance_selectors:
                nuisance_wf.connect(
                    regressor_node[0], regressor_node[1],
                    build_nuisance_regressors, regressor_arg
                )

    # Check for any regressors to combine into files
    has_voxel_nuisance_regressors = any(
        regressor_node
        for regressor_key, (
            regressor_arg, regressor_node, regressor_target
        ) in regressors.items()
        if regressor_target == 'dsort'
    )

    if has_voxel_nuisance_regressors:

        voxel_nuisance_regressors = [
            (regressor_key, (regressor_arg, regressor_node, regressor_target))
            for regressor_key, (
                regressor_arg, regressor_node, regressor_target
            ) in regressors.items()
            if regressor_target == 'dsort'
        ]

        voxel_nuisance_regressors_merge = pe.Node(
            util.Merge(len(voxel_nuisance_regressors)),
            name='voxel_nuisance_regressors_merge'
        )

        for i, (regressor_key, (
            regressor_arg, regressor_node, regressor_target)
        ) in enumerate(voxel_nuisance_regressors):

            if regressor_target != 'dsort':
                continue

            node, node_output = regressor_node

            nuisance_wf.connect(
                node, node_output,
                voxel_nuisance_regressors_merge, "in{}".format(i + 1)
            )

    nuisance_wf.connect(build_nuisance_regressors, 'out_file',
                        outputspec, 'regressors_file_path')

    return nuisance_wf


def create_nuisance_regression_workflow(nuisance_selectors,
                                        name='nuisance_regression'):

    inputspec = pe.Node(util.IdentityInterface(fields=[
        'functional_file_path',
        'functional_brain_mask_file_path',
        'regressor_file',
        'fd_j_file_path',
        'fd_p_file_path',
        'dvars_file_path'
    ]), name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['residual_file_path']),
                         name='outputspec')

    nuisance_wf = pe.Workflow(name=name)

    if nuisance_selectors.get('Censor'):

        censor_methods = ['Kill', 'Zero', 'Interpolate', 'SpikeRegression']

        censor_selector = nuisance_selectors.get('Censor')
        if censor_selector.get('method') not in censor_methods:
            raise ValueError("Improper censoring method specified ({0}), "
                             "should be one of {1}."
                             .format(censor_selector.get('method'),
                                     censor_methods))

        find_censors = pe.Node(Function(
            input_names=['fd_j_file_path',
                         'fd_j_threshold',
                         'fd_p_file_path',
                         'fd_p_threshold',
                         'dvars_file_path',
                         'dvars_threshold',
                         'number_of_previous_trs_to_censor',
                         'number_of_subsequent_trs_to_censor'],
            output_names=['out_file'],
            function=find_offending_time_points,
            as_module=True
        ), name="find_offending_time_points")

        if not censor_selector.get('thresholds'):
            raise ValueError(
                'Censoring requested, but thresh_metric not provided.'
            )

        for threshold in censor_selector['thresholds']:

            if 'type' not in threshold or threshold['type'] not in ['DVARS', 'FD_J', 'FD_P']:
                raise ValueError(
                    'Censoring requested, but with invalid threshold type.'
                )

            if 'value' not in threshold:
                raise ValueError(
                    'Censoring requested, but threshold not provided.'
                )

            if threshold['type'] == 'FD_J':
                find_censors.inputs.fd_j_threshold = threshold['value']
                nuisance_wf.connect(inputspec, "fd_j_file_path",
                                    find_censors, "fd_j_file_path")

            if threshold['type'] == 'FD_P':
                find_censors.inputs.fd_p_threshold = threshold['value']
                nuisance_wf.connect(inputspec, "fd_p_file_path",
                                    find_censors, "fd_p_file_path")

            if threshold['type'] == 'DVARS':
                find_censors.inputs.dvars_threshold = threshold['value']
                nuisance_wf.connect(inputspec, "dvars_file_path",
                                    find_censors, "dvars_file_path")

        if censor_selector.get('number_of_previous_trs_to_censor') and \
                censor_selector['method'] != 'SpikeRegression':

            find_censors.inputs.number_of_previous_trs_to_censor = \
                censor_selector['number_of_previous_trs_to_censor']

        else:
            find_censors.inputs.number_of_previous_trs_to_censor = 0

        if censor_selector.get('number_of_subsequent_trs_to_censor') and \
                censor_selector['method'] != 'SpikeRegression':

            find_censors.inputs.number_of_subsequent_trs_to_censor = \
                censor_selector['number_of_subsequent_trs_to_censor']

        else:
            find_censors.inputs.number_of_subsequent_trs_to_censor = 0

    # Use 3dTproject to perform nuisance variable regression
    nuisance_regression = pe.Node(interface=afni.TProject(),
                                  name='nuisance_regression')

    nuisance_regression.inputs.out_file = 'residuals.nii.gz'
    nuisance_regression.inputs.outputtype = 'NIFTI_GZ'
    nuisance_regression.inputs.norm = False

    if nuisance_selectors.get('Censor'):
        if nuisance_selectors['Censor']['method'] == 'SpikeRegression':
            nuisance_wf.connect(find_censors, 'out_file',
                                nuisance_regression, 'censor')
        else:
            if nuisance_selectors['Censor']['method'] == 'Interpolate':
                nuisance_regression.inputs.cenmode = 'NTRP'
            else:
                nuisance_regression.inputs.cenmode = \
                    nuisance_selectors['Censor']['method'].upper()

            nuisance_wf.connect(find_censors, 'out_file',
                                nuisance_regression, 'censor')

    if nuisance_selectors.get('PolyOrt'):
        if not nuisance_selectors['PolyOrt'].get('degree'):
            raise ValueError("Polynomial orthogonalization requested, "
                             "but degree not provided.")

        nuisance_regression.inputs.polort = \
            nuisance_selectors['PolyOrt']['degree']

    else:
        nuisance_regression.inputs.polort = 0

    nuisance_wf.connect(inputspec, 'functional_file_path',
                        nuisance_regression, 'in_file')

    nuisance_wf.connect(inputspec, 'functional_brain_mask_file_path',
                        nuisance_regression, 'mask')

    if nuisance_selectors.get('Custom'):
        if nuisance_selectors['Custom'][0].get('file'):
            if nuisance_selectors['Custom'][0]['file'].endswith('.nii') or \
                    nuisance_selectors['Custom'][0]['file'].endswith('.nii.gz'):
                nuisance_wf.connect(inputspec, 'regressor_file',
                                    nuisance_regression, 'dsort')
            else:
                nuisance_wf.connect(inputspec, 'regressor_file',
                                    nuisance_regression, 'ort')
        else:
            nuisance_wf.connect(inputspec, 'regressor_file',
                                nuisance_regression, 'ort')
    else:
        # there's no regressor file generated if only Bandpass in nuisance_selectors
        if not ('Bandpass' in nuisance_selectors and len(nuisance_selectors.selector.keys()) == 1):
            nuisance_wf.connect(inputspec, 'regressor_file',
                                nuisance_regression, 'ort')

    nuisance_wf.connect(nuisance_regression, 'out_file',
                        outputspec, 'residual_file_path')

    return nuisance_wf


def filtering_bold_and_regressors(nuisance_selectors,
                                  name='filtering_bold_and_regressors'):

    inputspec = pe.Node(util.IdentityInterface(fields=[
        'functional_file_path',
        'regressors_file_path',
        'functional_brain_mask_file_path',
        'tr'
    ]), name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['residual_file_path',
                                                        'residual_regressor']),
                         name='outputspec')

    filtering_wf = pe.Workflow(name=name)
    bandpass_selector = nuisance_selectors.get('Bandpass')

    if bandpass_selector.get('method'):
        bandpass_method = bandpass_selector.get('method')
    else:
        bandpass_method = 'default'

    if bandpass_method == 'default':

        frequency_filter = pe.Node(
                    Function(input_names=['realigned_file',
                                        'regressor_file',
                                        'bandpass_freqs',
                                        'sample_period'],
                            output_names=['bandpassed_file',
                                        'regressor_file'],
                            function=bandpass_voxels,
                            as_module=True),
                    name='frequency_filter'
                )

        frequency_filter.inputs.bandpass_freqs = [
                    bandpass_selector.get('bottom_frequency'),
                    bandpass_selector.get('top_frequency')
                ]

        filtering_wf.connect(inputspec, 'functional_file_path',
                            frequency_filter, 'realigned_file')

        filtering_wf.connect(inputspec, 'regressors_file_path',
                            frequency_filter, 'regressor_file')

        filtering_wf.connect(frequency_filter, 'bandpassed_file',
                            outputspec, 'residual_file_path')

        filtering_wf.connect(frequency_filter, 'regressor_file',
                            outputspec, 'residual_regressor')

    elif bandpass_method == 'AFNI':

        bandpass_ts = pe.Node(interface=afni.Bandpass(),
                                    name='bandpass_ts')

        bandpass_ts.inputs.highpass = bandpass_selector.get('bottom_frequency')
        bandpass_ts.inputs.lowpass = bandpass_selector.get('top_frequency')
        bandpass_ts.inputs.outputtype = 'NIFTI_GZ'

        tr_string2float_node = pe.Node(util.Function(input_names=['tr'],
                                                     output_names=['tr_float'],
                                                     function=TR_string_to_float),
                                        name='tr_string2float')

        filtering_wf.connect(inputspec, 'tr',
                            tr_string2float_node, 'tr')

        filtering_wf.connect(tr_string2float_node, 'tr_float',
                            bandpass_ts, 'tr')

        filtering_wf.connect(inputspec, 'functional_file_path',
                            bandpass_ts, 'in_file')

        filtering_wf.connect(inputspec, 'functional_brain_mask_file_path',
                            bandpass_ts, 'mask')

        filtering_wf.connect(bandpass_ts, 'out_file',
                            outputspec, 'residual_file_path')

        bandpass_regressor = pe.Node(Function(input_names=['in_file',
                                                           'highpass',
                                                           'lowpass',
                                                           'tr'],
                                              output_names=['out_file'],
                                              function=afni_1dBandpass),
                                     name='bandpass_regressor')

        bandpass_regressor.inputs.highpass = bandpass_selector.get('bottom_frequency')
        bandpass_regressor.inputs.lowpass = bandpass_selector.get('top_frequency')

        filtering_wf.connect(inputspec, 'regressors_file_path',
                            bandpass_regressor, 'in_file')

        filtering_wf.connect(tr_string2float_node, 'tr_float',
                            bandpass_regressor, 'tr')

        filtering_wf.connect(bandpass_regressor, 'out_file',
                            outputspec, 'residual_regressor')

    return filtering_wf
