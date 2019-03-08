import re
import os
import numpy as np
import nibabel as nb

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants
from nipype.interfaces import afni

import CPAC
import CPAC.utils as utils
from CPAC.utils.function import Function
from CPAC.nuisance import (
    find_offending_time_points,
    create_temporal_variance_mask,
    summarize_timeseries,
    generate_summarize_tissue_mask,
)

from scipy.fftpack import fft, ifft

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants


def gather_nuisance(functional_file_path,
                    selector,
                    grey_matter_summary_file_path=None,
                    white_matter_summary_file_path=None,
                    csf_summary_file_path=None,
                    acompcor_file_path=None,
                    tcompcor_file_path=None,
                    global_summary_file_path=None,
                    motion_parameters_file_path=None,
                    dvars_file_path=None,
                    framewise_displacement_file_path=None,
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
    :param dvars_file_path: path to TSV that includes DVARS as a single column
    :param framewise_displacement_file_path: path to TSV that includes framewise displacement as a single column
    :param censor_file_path: path to TSV with a single column with 1's for indices that should be retained and 0's
              for indices that should be censored
    :return: 
    """

    # Basic checks for the functional image
    if not functional_file_path or \
        not os.path.isfile(functional_file_path) or \
        (not functional_file_path.endswith(".nii") and
         not functional_file_path.endswith(".nii.gz")):

        raise ValueError("Invalid value for input_file ({}). Should be a nifti "
                         "file and should exist".format(functional_file_path))

    functional_image = nb.load(functional_file_path)

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
                    'method': regressor_selector['summary']
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
                column_names.append("{0}Back".format(regressor_name))
                nuisance_regressors.append(
                    np.append([0.0], regressors[0:-1, regressor_index])
                )

            if regressor_selector.get('include_squared', False):
                column_names.append("{0}Sq".format(regressor_name))
                nuisance_regressors.append(
                    np.square(regressors[:, regressor_index])
                )

            if regressor_selector.get('include_delayed_squared', False):
                column_names.append("{0}BackSq".format(regressor_name))
                nuisance_regressors.append(
                    np.square(
                        np.append([0.0], regressors[0:-1, regressor_index])
                    )
                )

    # Add spike regressors
    if selector.get('Censor', {}).get('method') == 'SpikeRegression':

        selector = selector['Censor']

        regressor_file = censor_file_path

        if not regressor_file:
            raise ValueError("Regressor type Censor specified in selectors but "
                             "the corresponding file was not found!")

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

            for censor_index in censor_indices:

                censor_begin_index = censor_index - previous_trs_to_censor
                if censor_begin_index < 0:
                    censor_begin_index = 0

                censor_end_index = censor_index + subsequent_trs_to_censor
                if censor_end_index >= regressor_length:
                    censor_end_index = regressor_length - 1

                spike_regressor = np.zeros((regressor_length, 1))
                spike_regressor[censor_begin_index:censor_end_index + 1, 0] = 1

                column_names.append("SpikeRegression{0}".format(censor_index))
                nuisance_regressors.append(spike_regressor.flatten())

    # Compile columns into regressor file
    output_file_path = os.path.join(os.getcwd(), "nuisance_regressors.1D")

    with open(output_file_path, "w") as ofd:

        # write out the header information
        ofd.write("# CPAC {0}\n".format(CPAC.__version__))
        ofd.write("# Nuisance regressors:\n")
        ofd.write("# " + "\t".join(column_names) + "\n")

        nuisance_regressors = np.array(nuisance_regressors)
        np.savetxt(ofd, nuisance_regressors.T, fmt='%.18f', delimiter='\t')

    return output_file_path


def create_nuisance_workflow(nuisance_selectors,
                             use_ants,
                             name='nuisance'):
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
        inputspec.fd_file_path : string (text file)
            Framewise displacement calculated from the motion parameters.
        inputspec.dvars_file_path : string (text file)
            DVARS calculated from the functional data.

        inputspec.selector : Dictionary containing configuration parameters for nuisance regression.
            To not run a type of nuisance regression, it may be ommited from the dictionary.
            selector = {
                aCompCor: {
                    symmary: {
                        method: 'DetrendPC', aCompCor will always extract the principal components from
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
                },
                tCompCor: {
                    symmary: {
                        method: 'PC', tCompCor will always extract the principal components from
                            BOLD signal,
                        components: number of components to retain,
                    },
                    threshold:
                        floating point number = cutoff as raw variance value,
                        floating point number followed by SD (ex. 1.5SD) = mean + a multiple of the SD,
                        floating point number followed by PCT (ex. 2PCT) = percentile from the top (ex is top 2%),
                    by_slice: True | False, whether or not the threshold criterion should be applied
                        by slice or across the entire volume, makes most sense for thresholds
                        using SD or PCT,
                    include_delayed: True | False,
                    include_squared: True | False,
                    include_delayed_squared: True | False,
                },
                WhiteMatter: {
                    symmary: {
                        method: 'PC', 'DetrendPC', 'Mean', 'NormMean' or 'DetrendNormMean',
                        components: number of components to retain, if PC,
                    },
                    extraction_resolution: None | floating point value (same as for aCompCor),
                    erode_mask: True | False (same as for aCompCor),
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                },
                CerebrospinalFluid: {
                    symmary: {
                        method: 'PC', 'DetrendPC', 'Mean', 'NormMean' or 'DetrendNormMean',
                        components: number of components to retain, if PC,
                    },
                    extraction_resolution: None | floating point value (same as for aCompCor),
                    erode_mask: True | False (same as for aCompCor),
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                },
                GreyMatter: {
                    symmary: {
                        method: 'PC', 'DetrendPC', 'Mean', 'NormMean' or 'DetrendNormMean',
                        components: number of components to retain, if PC,
                    },
                    extraction_resolution: None | floating point value (same as for aCompCor),
                    erode_mask: True | False (same as for aCompCor),
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                },
                GlobalSignal: {
                    symmary: {
                        method: 'PC', 'DetrendPC', 'Mean', 'NormMean' or 'DetrendNormMean',
                        components: number of components to retain, if PC,
                    },
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                },
                Motion: None | { 
                    include_delayed: True | False (same as for aCompCor),
                    include_squared: True | False (same as for aCompCor),
                    include_delayed_squared: True | False (same as for aCompCor),
                },
                Censor: {
                    method: 'Kill', 'Zero', 'Interpolate', 'SpikeRegression',
                    thresholds: list of dictionary, {
                        type: 'FD', 'DVARS',
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
                }
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

    Workflow Graph:

    .. image:: ../images/nuisance.dot.png
        :width: 500

    Detailed Workflow Graph:

    .. image:: ../images/nuisance_detailed.dot.png
        :width: 500    

    """

    nuisance_wf = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=[
        'selector',
        'functional_file_path',

        'anatomical_file_path',
        'gm_mask_file_path',
        'wm_mask_file_path',
        'csf_mask_file_path',
        'lat_ventricles_mask_file_path',

        'functional_brain_mask_file_path',

        'func_to_anat_linear_xfm_file_path',
        'mni_to_anat_linear_xfm_file_path',
        'anat_to_mni_initial_xfm_file_path',
        'anat_to_mni_rigid_xfm_file_path',
        'anat_to_mni_affine_xfm_file_path',

        'motion_parameters_file_path',
        'dvars_file_path',
        'fd_file_path',
    ]), name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['residual_file_path',
                                                        'no_bandpass_residual_file_path',
                                                        'regressors_file_path']),
                         name='outputspec')

    # Resources to create regressors
    pipeline_resource_pool = {
        "Anatomical": (inputspec, 'anatomical_file_path'),
        "Functional": (inputspec, 'functional_file_path'),
        "GlobalSignal": (inputspec, 'functional_brain_mask_file_path'),
        "WhiteMatter": (inputspec, 'wm_mask_file_path'),
        "CerebrospinalFluid": (inputspec, 'csf_mask_file_path'),
        "GreyMatter": (inputspec, 'gm_mask_file_path'),
        "Ventricles": (inputspec, 'lat_ventricles_mask_file_path'),

        "Transformations": {
            "func_to_anat_linear_xfm": (inputspec, "func_to_anat_linear_xfm_file_path"),
            "mni_to_anat_linear_xfm": (inputspec, "mni_to_anat_linear_xfm_file_path"),
            "anat_to_mni_initial_xfm": (inputspec, "anat_to_mni_initial_xfm_file_path"),
            "anat_to_mni_rigid_xfm": (inputspec, "anat_to_mni_rigid_xfm_file_path"),
            "anat_to_mni_affine_xfm": (inputspec, "anat_to_mni_affine_xfm_file_path"),
        }
    }

    # Regressor map to simplify construction of the needed regressors
    regressors = {
        'GreyMatter': ['grey_matter_summary_file_path', ()],
        'WhiteMatter': ['white_matter_summary_file_path', ()],
        'CerebrospinalFluid': ['csf_summary_file_path', ()],
        'aCompCor': ['acompcor_file_path', ()],
        'tCompCor': ['tcompcor_file_path', ()],
        'GlobalSignal': ['global_summary_file_path', ()],
        'DVARS': ['dvars_file_path', (inputspec, 'dvars_file_path')],
        'FD': ['framewise_displacement_file_path', (inputspec, 'framewise_displacement_file_path')],
        'Motion': ['motion_parameters_file_path', (inputspec, 'motion_parameters_file_path')]
    }

    derived = ['tCompCor', 'aCompCor']
    tissues = ['GreyMatter', 'WhiteMatter', 'CerebrospinalFluid']

    for regressor_type, regressor_resource in regressors.items():

        if regressor_type not in nuisance_selectors:
            continue

        regressor_selector = nuisance_selectors[regressor_type]

        # Set summary method for tCompCor and aCompCor
        if regressor_type in derived:

            if 'summary' not in regressor_selector:
                regressor_selector['summary'] = {}

            if type(regressor_selector['summary']) is not dict:
                raise ValueError("Regressor {0} requires PC summary method, "
                                 "but {1} specified"
                                 .format(regressor_type,
                                         regressor_selector['summary']))

            regressor_selector['summary']['method'] = \
                'DetrendPC' if regressor_type == 'aCompCor' else 'PC'

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

            if type(regressor_selector['summary']) is dict:
                regressor_descriptor['extraction'] = \
                    regressor_selector['summary']['method']
            else:
                regressor_descriptor['extraction'] = \
                    regressor_selector['summary']

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
                        use_ants=use_ants
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

                # Retrieve summary from voxels at provided mask
                summarize_timeseries_node = pe.Node(
                    Function(
                        input_names=[
                            'functional_path',
                            'masks_path',
                            'summary'
                        ],
                        output_names=['components_path'],
                        function=summarize_timeseries,
                        as_module=True,
                    ),
                    name='{}_summarization'.format(regressor_type)
                )

                summarize_timeseries_node.inputs.summary = \
                    regressor_selector['summary']

                # Merge mask paths to extract voxel timeseries
                merge_masks_paths = pe.Node(
                    util.Merge(len(regressor_mask_file_resource_keys)),
                    name='{}_marge_masks'.format(regressor_type)
                )
                for i, regressor_mask_file_resource_key in \
                        enumerate(regressor_mask_file_resource_keys):

                    node, node_output = \
                        pipeline_resource_pool[regressor_mask_file_resource_key]

                    nuisance_wf.connect(
                        node, node_output,
                        merge_masks_paths, "in{}".format(i + 1)
                    )

                nuisance_wf.connect(
                    merge_masks_paths, 'out',
                    summarize_timeseries_node, 'masks_path'
                )

                functional_key = 'Functional'
                if regressor_selector.get('extraction_resolution') and \
                        regressor_selector["extraction_resolution"] != "Functional":

                    functional_key = 'Functional_{}mm'.format(
                        regressor_selector['extraction_resolution']
                    )

                nuisance_wf.connect(*(
                    pipeline_resource_pool[functional_key] +
                    (summarize_timeseries_node, 'functional_path')
                ))

                pipeline_resource_pool[regressor_file_resource_key] = \
                    (summarize_timeseries_node, 'components_path')

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
                     'dvars_file_path',
                     'framewise_displacement_file_path',
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
        regressor_resource[1]
        for regressor_key, regressor_resource
        in regressors.items()
    )

    if has_nuisance_regressors:
        for regressor_key, (regressor_arg, regressor_node) in regressors.items():
            if regressor_key in nuisance_selectors:
                nuisance_wf.connect(
                    regressor_node[0], regressor_node[1],
                    build_nuisance_regressors, regressor_arg
                )

    if nuisance_selectors.get('Censor'):

        censor_methods = ['Kill', 'Zero', 'Interpolate', 'SpikeRegression']

        censor_selector = nuisance_selectors.get('Censor')
        if censor_selector.get('method') not in censor_methods:
            raise ValueError("Improper censoring method specified ({0}), "
                             "should be one of {1}."
                             .format(censor_selector.get('method'),
                                     censor_methods))

        find_censors = pe.Node(Function(
            input_names=['fd_file_path',
                         'fd_threshold',
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

            if 'type' not in threshold or threshold['type'] not in ['DVARS', 'FD']:
                raise ValueError(
                    'Censoring requested, but with invalid threshold type.'
                )

            if 'value' not in threshold:
                raise ValueError(
                    'Censoring requested, but threshold not provided.'
                )

            if threshold['type'] == 'FD':
                find_censors.inputs.fd_threshold = threshold['value']
                nuisance_wf.connect(inputspec, "fd_file_path",
                                    find_censors, "fd_file_path")

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
                                build_nuisance_regressors, 'censor_file_path')
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

    no_bandpass_nuisance_regression = None
    if nuisance_selectors.get('Bandpass'):

        no_bandpass_nuisance_regression = nuisance_regression.clone('no_bandpass_nuisance_regression')

        bandpass_selector = nuisance_selectors['Bandpass']
        bottom_frequency = bandpass_selector.get('bottom_frequency', 0.0)
        top_frequency = bandpass_selector.get('top_frequency', 9999.9)

        no_bandpass_nuisance_regression.inputs.bandpass = (float(bottom_frequency),
                                                        float(top_frequency))

    nuisance_wf.connect([
        (inputspec, nuisance_regression, [
            ('functional_file_path', 'in_file'),
            ('functional_brain_mask_file_path', 'mask'),
        ]),
    ])

    if has_nuisance_regressors:
        nuisance_wf.connect(build_nuisance_regressors, 'out_file',
                            nuisance_regression, 'ort')

    nuisance_wf.connect(nuisance_regression, 'out_file',
                        outputspec, 'residual_file_path')

    if nuisance_selectors.get('Bandpass'):
        nuisance_wf.connect([
            (inputspec, no_bandpass_nuisance_regression, [
                ('functional_file_path', 'in_file'),
                ('functional_brain_mask_file_path', 'mask'),
            ]),
        ])

        nuisance_wf.connect(build_nuisance_regressors, 'out_file',
                            no_bandpass_nuisance_regression, 'ort')

        nuisance_wf.connect(no_bandpass_nuisance_regression, 'out_file',
                            outputspec, 'no_bandpass_residual_file_path')
    else:
        nuisance_wf.connect(nuisance_regression, 'out_file',
                            outputspec, 'no_bandpass_residual_file_path')


    nuisance_wf.connect(build_nuisance_regressors, 'out_file',
                        outputspec, 'regressors_file_path')

    return nuisance_wf
