# Copyright (C) 2012-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
import re
import os
from CPAC.pipeline.nodeblock import nodeblock
import numpy as np
import nibabel as nb
# pylint: disable=wrong-import-order
from nipype.pipeline.engine.workflows import Workflow
from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.utility as util
import CPAC

from nipype import logging
from nipype.interfaces import fsl
from nipype.interfaces import ants
from nipype.interfaces import c3
from nipype.interfaces import afni
from nipype.interfaces.afni import utils as afni_utils
from scipy.fftpack import fft, ifft
from CPAC.pipeline.engine import ResourcePool
from CPAC.utils.configuration import Configuration
from CPAC.utils.interfaces.function import Function
from CPAC.utils.interfaces.masktool import MaskTool
from CPAC.utils.interfaces.pc import PC
from CPAC.utils.typing import LITERAL, TUPLE

from CPAC.registration.registration import warp_timeseries_to_T1template, \
    warp_timeseries_to_EPItemplate, apply_transform
from CPAC.aroma.aroma import create_aroma

from CPAC.nuisance.utils import (
    find_offending_time_points,
    generate_summarize_tissue_mask,
    temporal_variance_mask)
from CPAC.nuisance.utils.compcor import (
    calc_compcor_components,
    cosine_filter,
    TR_string_to_float)

from CPAC.seg_preproc.utils import erosion, mask_erosion

from CPAC.utils.datasource import check_for_s3
from CPAC.utils.utils import check_prov_for_regtool
from .bandpass import (bandpass_voxels, afni_1dBandpass)
logger = logging.getLogger('nipype.workflow')


def choose_nuisance_blocks(cfg, rpool, generate_only=False):
    '''
    Function to handle selecting appropriate blocks based on
    existing config and resource pool

    Parameters
    ----------
    cfg : CPAC.utils.configuration.Configuration

    generate_only : boolean
        generate but don't run

    Returns
    -------
    nuisance : list
    '''
    nuisance = []
    to_template_cfg = cfg.registration_workflows['functional_registration'][
        'func_registration_to_template']
    apply_transform_using = to_template_cfg['apply_transform']['using']
    input_interface = {
        'default': ('desc-preproc_bold', ['desc-preproc_bold', 'bold']),
        'abcd': ('desc-preproc_bold', 'bold'),
        'single_step_resampling_from_stc': ("desc-preproc_bold",
                                            "desc-stc_bold")
    }.get(apply_transform_using)
    if input_interface is not None:
        if 'T1_template' in to_template_cfg['target_template']['using']:
            nuisance.append((nuisance_regressors_generation_T1w,
                             input_interface))
        if 'EPI_template' in to_template_cfg['target_template']['using']:
            nuisance.append((nuisance_regressors_generation_EPItemplate,
                             input_interface))

        if not generate_only and cfg['nuisance_corrections',
                                     '2-nuisance_regression',
                                     'space'] == 'native':
            nuisance.append((nuisance_regression_native, input_interface))

    return nuisance


def erode_mask(name, segmentmap=True):

    wf = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['mask',
                                                       'erode_mm',
                                                       'erode_prop',
                                                       'brain_mask',
                                                       'mask_erode_mm']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['eroded_mask']),
                         name='outputspec')

    def form_mask_erosion_prop(erosion_prop):
        if not isinstance(erosion_prop, (int, float)):
            erosion_prop = 0
        return erosion_prop ** 3

    ero_imports = ['import scipy.ndimage as nd', 'import numpy as np',
                   'import nibabel as nb', 'import os',
                   'from CPAC.seg_preproc.utils import _erode']

    eroded_mask = pe.Node(util.Function(
        input_names=['roi_mask', 'skullstrip_mask', 'mask_erosion_mm',
                     'mask_erosion_prop'],
        output_names=['output_roi_mask', 'eroded_skullstrip_mask'],
        function=mask_erosion,
        imports=ero_imports),
                          name='erode_skullstrip_mask',
                          mem_gb=2.3,
                          mem_x=(4664065662093477 / 2417851639229258349412352,
                                 'roi_mask'))

    wf.connect(inputspec, 'brain_mask', eroded_mask, 'skullstrip_mask')
    wf.connect(inputspec, 'mask', eroded_mask, 'roi_mask')

    wf.connect(inputspec, ('erode_prop', form_mask_erosion_prop), eroded_mask,
               'mask_erosion_prop')
    wf.connect(inputspec, 'mask_erode_mm', eroded_mask, 'mask_erosion_mm')

    if not segmentmap:
        wf.connect(eroded_mask, 'output_roi_mask', outputspec, 'eroded_mask')
    if segmentmap:
        erosion_segmentmap = pe.Node(util.Function(input_names=['roi_mask',
                                                                'erosion_mm',
                                                                'erosion_prop'
                                                                ],
                                                   output_names=[
                                                       'eroded_roi_mask'],
                                                   function=erosion,
                                                   imports=ero_imports),
                                     name='erode_mask')

        wf.connect(eroded_mask, 'output_roi_mask', erosion_segmentmap, 'roi_mask')

        wf.connect(inputspec, 'erode_prop', erosion_segmentmap, 'erosion_prop')
        wf.connect(inputspec, 'erode_mm', erosion_segmentmap, 'erosion_mm')

        wf.connect(erosion_segmentmap, 'eroded_roi_mask',
                   outputspec, 'eroded_mask')

    return wf


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
    Gathers the various nuisance regressors together into a single tab-
    separated values file that is an appropriate for input into
    3dTproject

    :param functional_file_path: path to file that the regressors are
        being calculated for, is used to calculate the length of the
        regressors for error checking and in particular for calculating
        spike regressors
    :param output_file_path: path to output TSV that will contain the
        various nuisance regressors as columns
    :param grey_matter_summary_file_path: path to TSV that includes
        summary of grey matter time courses, e.g. output of
        mask_summarize_time_course
    :param white_matter_summary_file_path: path to TSV that includes
        summary of white matter time courses, e.g. output of
        mask_summarize_time_course
    :param csf_summary_file_path: path to TSV that includes summary of
        csf time courses, e.g. output of mask_summarize_time_course
    :param acompcor_file_path: path to TSV that includes acompcor time
        courses, e.g. output of mask_summarize_time_course
    :param tcompcor_file_path: path to TSV that includes tcompcor time
        courses, e.g. output of mask_summarize_time_course
    :param global_summary_file_path: path to TSV that includes summary
        of global time courses, e.g. output of mask_summarize_time_course
    :param motion_parameters_file_path: path to TSV that includes
        motion parameters
    :param custom_file_paths: path to CSV/TSV files to use as regressors
    :param censor_file_path: path to TSV with a single column with '1's
        for indices that should be retained and '0's for indices that
        should be censored
    :return: out_file (str), censor_indices (list)
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

    #selector = selector.selector

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

    censor_indices = []
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

    return output_file_path, censor_indices


def create_regressor_workflow(nuisance_selectors,
                              use_ants,
                              ventricle_mask_exist,
                              csf_mask_exist,
                              all_bold=False,
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
        outputspec.censor_indices : list
            Indices of censored volumes

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
        'anat_to_mni_linear_xfm_file_path',
        'motion_parameters_file_path',
        'fd_j_file_path',
        'fd_p_file_path',
        'dvars_file_path',
        'creds_path',
        'dl_dir',
        'tr',
    ]), name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(
        fields=['regressors_file_path', 'censor_indices']), name='outputspec')

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
            "anat_to_mni_linear_xfm": (inputspec, "anat_to_mni_linear_xfm_file_path")
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

            for file_num, custom_regressor in enumerate(sorted(
                regressor_selector, key=lambda c: c['file']
            )):
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
                    name=f'custom_check_for_s3_{name}_{file_num}')

                custom_check_s3_node.inputs.set(
                    file_path=custom_regressor_file,
                    img_type='func'
                )

                if (
                    custom_regressor_file.endswith('.nii.gz') or
                    custom_regressor_file.endswith('.nii')
                ):

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
                            .format(anatomical_at_resolution_key),
                        mem_gb=3.63,
                        mem_x=(3767129957844731 / 1208925819614629174706176,
                            'in_file')
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
                            .format(functional_at_resolution_key),
                        mem_gb=0.521,
                        mem_x=(4394984950818853 / 302231454903657293676544,
                            'in_file')
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
                        csf_mask_exist,
                        use_ants=use_ants,
                        ventricle_mask_exist=ventricle_mask_exist,
                        all_bold=all_bold
                    )

                regressor_mask_file_resource_keys += \
                    [regressor_mask_file_resource_key]

            # Keep tissues ordered, to avoid duplicates
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
                    name='{}_union_masks'.format(regressor_type),
                    mem_gb=2.1,
                    mem_x=(1708448960473801 / 1208925819614629174706176,
                        'in_files')
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
                                        name='{}_DetrendPC'.format(regressor_type),
                                        mem_gb=0.4,
                                        mem_x=(3811976743057169 /
                                                151115727451828646838272,
                                                'data_filename'))

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

                        cosfilter_node = pe.Node(
                            util.Function(input_names=['input_image_path',
                                                    'timestep'],
                                        output_names=['cosfiltered_img'],
                                        function=cosine_filter,
                                        imports=cosfilter_imports),
                            name=f'{regressor_type}_cosine_filter',
                            mem_gb=8.0,
                            throttle=True)
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
                            name='{}_norm'.format(regressor_type),
                            mem_gb=1.7,
                            mem_x=(1233286593342025 /
                                151115727451828646838272,
                                'in_file_a')
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
                            name='{}_mean'.format(regressor_type),
                            mem_gb=5.0
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
        output_names=['out_file', 'censor_indices'],
        function=gather_nuisance,
        as_module=True
    ), name="build_nuisance_regressors")

    nuisance_wf.connect(
        inputspec, 'functional_file_path',
        build_nuisance_regressors, 'functional_file_path'
    )

    build_nuisance_regressors.inputs.selector = nuisance_selectors

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

    nuisance_wf.connect([(build_nuisance_regressors, outputspec, [
        ('out_file', 'regressors_file_path'),
        ('censor_indices', 'censor_indices')])])

    return nuisance_wf

def create_nuisance_regression_workflow(nuisance_selectors,
                                        name='nuisance_regression'):

    inputspec = pe.Node(util.IdentityInterface(fields=[
        'selector',
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
                                  name='nuisance_regression',
                                  mem_gb=1.716,
                                  mem_x=(6278549929741219 /
                                         604462909807314587353088, 'in_file'))

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
        if not ('Bandpass' in nuisance_selectors and len(nuisance_selectors.keys()) == 1):
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
        'nuisance_selectors',
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
                    name='frequency_filter',
                    mem_gb=0.5,
                    mem_x=(3811976743057169 / 151115727451828646838272,
                           'realigned_file')
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


@nodeblock(
    name="ICA_AROMA_FSLreg",
    config=["nuisance_corrections", "1-ICA-AROMA"],
    switch=["run"],
    inputs=[
        "desc-preproc_bold",
        "from-bold_to-T1w_mode-image_desc-linear_xfm",
        "from-T1w_to-template_mode-image_xfm",
    ],
    outputs=["desc-preproc_bold", "desc-cleaned_bold"],
)
def ICA_AROMA_FSLreg(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance('from-T1w_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    if reg_tool != 'fsl':
        return (wf, None)

    aroma_preproc = create_aroma(tr=None, wf_name=f'create_aroma_{pipe_num}')

    aroma_preproc.inputs.params.denoise_type = \
        cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type']

    node, out = strat_pool.get_data("desc-preproc_bold")
    wf.connect(node, out, aroma_preproc, 'inputspec.denoise_file')

    node, out = strat_pool.get_data(
        "from-bold_to-T1w_mode-image_desc-linear_xfm")
    wf.connect(node, out, aroma_preproc, 'inputspec.mat_file')

    node, out = strat_pool.get_data("from-T1w_to-template_mode-image_xfm")
    wf.connect(node, out, aroma_preproc, 'inputspec.fnirt_warp_file')

    if cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type'] == 'nonaggr':
        node, out = (aroma_preproc, 'outputspec.nonaggr_denoised_file')
    elif cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type'] == 'aggr':
        node, out = (aroma_preproc, 'outputspec.aggr_denoised_file')

    outputs = {
        'desc-preproc_bold': (node, out),
        'desc-cleaned_bold': (node, out)
    }

    return (wf, outputs)


@nodeblock(
    name="ICA_AROMA_ANTsreg",
    config=["nuisance_corrections", "1-ICA-AROMA"],
    switch=["run"],
    inputs=[
        (
            "desc-preproc_bold",
            "sbref",
            "from-bold_to-template_mode-image_xfm",
            "from-template_to-bold_mode-image_xfm",
        ),
        "T1w-brain-template-funcreg",
    ],
    outputs=["desc-preproc_bold", "desc-cleaned_bold"],
)
def ICA_AROMA_ANTsreg(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-bold_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    if reg_tool != 'ants':
        return (wf, None)

    num_cpus = cfg.pipeline_setup['system_config'][
        'max_cores_per_participant']

    num_ants_cores = cfg.pipeline_setup['system_config']['num_ants_threads']

    aroma_preproc = create_aroma(tr=None, wf_name=f'create_aroma_{pipe_num}')
    aroma_preproc.inputs.params.denoise_type = \
        cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type']

    wf, outputs = warp_timeseries_to_T1template(wf, cfg, strat_pool, pipe_num)
    for key, val in outputs.items():
        node, out = val

    wf.connect(node, out, aroma_preproc, 'inputspec.denoise_file')

    apply_xfm = apply_transform(f'ICA-AROMA_ANTs_template_to_bold_{pipe_num}',
                                reg_tool=reg_tool, time_series=True,
                                num_cpus=num_cpus,
                                num_ants_cores=num_ants_cores)
    apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'ANTs_pipelines']['interpolation']

    if cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type'] == 'nonaggr':
        node, out = (aroma_preproc, 'outputspec.nonaggr_denoised_file')
    elif cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type'] == 'aggr':
        node, out = (aroma_preproc, 'outputspec.aggr_denoised_file')

    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("sbref")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data('from-template_to-bold_mode-image_xfm')
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        'desc-preproc_bold': (apply_xfm, 'outputspec.output_image'),
        'desc-cleaned_bold': (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="ICA_AROMA_FSLEPIreg",
    switch=[
        ["nuisance_corrections", "1-ICA-AROMA", "run"],
        [
            "registration_workflows",
            "functional_registration",
            "EPI_registration",
            "run",
        ],
    ],
    inputs=[
        ["desc-brain_bold", "desc-motion_bold", "desc-preproc_bold", "bold"],
        "from-bold_to-EPItemplate_mode-image_xfm",
    ],
    outputs=["desc-preproc_bold", "desc-cleaned_bold"],
)
def ICA_AROMA_FSLEPIreg(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance('from-bold_to-EPItemplate_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    if reg_tool != 'fsl':
        return (wf, None)

    aroma_preproc = create_aroma(tr=None, wf_name=f'create_aroma_{pipe_num}')

    aroma_preproc.inputs.params.denoise_type = \
        cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type']

    node, out = strat_pool.get_data(["desc-brain_bold", "desc-motion_bold",
                                     "desc-preproc_bold", "bold"])
    wf.connect(node, out, aroma_preproc, 'inputspec.denoise_file')

    node, out = strat_pool.get_data("from-bold_to-EPItemplate_mode-image_xfm")
    wf.connect(node, out, aroma_preproc, 'inputspec.fnirt_warp_file')

    if cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type'] == 'nonaggr':
        node, out = (aroma_preproc, 'outputspec.nonaggr_denoised_file')
    elif cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type'] == 'aggr':
        node, out = (aroma_preproc, 'outputspec.aggr_denoised_file')

    outputs = {
        'desc-preproc_bold': (node, out),
        'desc-cleaned_bold': (node, out)
    }

    return (wf, outputs)


@nodeblock(
    name="ICA_AROMA_ANTsEPIreg",
    switch=[
        ["nuisance_corrections", "1-ICA-AROMA", "run"],
        [
            "registration_workflows",
            "functional_registration",
            "EPI_registration",
            "run",
        ],
    ],
    inputs=[
        (
            "desc-preproc_bold",
            "sbref",
            "from-bold_to-EPItemplate_mode-image_xfm",
            "from-EPItemplate_to-bold_mode-image_xfm",
        ),
        "EPI-template",
    ],
    outputs=["desc-preproc_bold", "desc-cleaned_bold"],
)
def ICA_AROMA_ANTsEPIreg(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-bold_to-EPItemplate_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    if reg_tool != 'ants':
        return (wf, None)

    num_cpus = cfg.pipeline_setup['system_config'][
        'max_cores_per_participant']

    num_ants_cores = cfg.pipeline_setup['system_config']['num_ants_threads']

    aroma_preproc = create_aroma(tr=None, wf_name=f'create_aroma_{pipe_num}')
    aroma_preproc.inputs.params.denoise_type = \
        cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type']

    wf, outputs = warp_timeseries_to_EPItemplate(wf, cfg, strat_pool,
                                                 pipe_num)
    for key, val in outputs.items():
        node, out = val

    wf.connect(node, out, aroma_preproc, 'inputspec.denoise_file')

    apply_xfm = apply_transform(f'ICA-AROMA_ANTs_EPItemplate_to_bold_{pipe_num}',
                                reg_tool=reg_tool, time_series=True,
                                num_cpus=num_cpus,
                                num_ants_cores=num_ants_cores)
    apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'ANTs_pipelines']['interpolation']

    if cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type'] == 'nonaggr':
        node, out = (aroma_preproc, 'outputspec.nonaggr_denoised_file')
    elif cfg.nuisance_corrections['1-ICA-AROMA']['denoising_type'] == 'aggr':
        node, out = (aroma_preproc, 'outputspec.aggr_denoised_file')

    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("sbref")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data('from-EPItemplate_to-bold_mode-image_xfm')
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        'desc-preproc_bold': (apply_xfm, 'outputspec.output_image'),
        'desc-cleaned_bold': (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="erode_mask_T1w",
    switch=[
        ["nuisance_corrections", "2-nuisance_regression", "create_regressors"],
        [
            "nuisance_corrections",
            "2-nuisance_regression",
            "regressor_masks",
            "erode_anatomical_brain_mask",
            "run",
        ],
    ],
    inputs=[
        ("space-T1w_desc-brain_mask", ["label-CSF_desc-preproc_mask", "label-CSF_mask"])
    ],
    outputs=["space-T1w_desc-eroded_mask"],
)
def erode_mask_T1w(wf, cfg, strat_pool, pipe_num, opt=None):

    erode = erode_mask(f'erode_T1w_mask_{pipe_num}', segmentmap=False)
    erode.inputs.inputspec.mask_erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks'][
        'erode_anatomical_brain_mask']['brain_mask_erosion_mm']
    erode.inputs.inputspec.erode_prop = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks'][
        'erode_anatomical_brain_mask']['brain_mask_erosion_prop']

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, erode, 'inputspec.brain_mask')

    node, out = strat_pool.get_data(['label-CSF_desc-preproc_mask',
                                     'label-CSF_mask'])
    wf.connect(node, out, erode, 'inputspec.mask')

    outputs = {
        'space-T1w_desc-eroded_mask': (erode, 'outputspec.eroded_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="erode_mask_CSF",
    switch=[
        ["nuisance_corrections", "2-nuisance_regression", "create_regressors"],
        [
            "nuisance_corrections",
            "2-nuisance_regression",
            "regressor_masks",
            "erode_csf",
            "run",
        ],
    ],
    inputs=[
        (["label-CSF_desc-preproc_mask", "label-CSF_mask"], "space-T1w_desc-brain_mask")
    ],
    outputs=["label-CSF_desc-eroded_mask"],
)
def erode_mask_CSF(wf, cfg, strat_pool, pipe_num, opt=None):

    erode = erode_mask(f'erode_CSF_mask_{pipe_num}')
    erode.inputs.inputspec.erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_csf'][
        'csf_erosion_mm']
    erode.inputs.inputspec.erode_prop = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_csf'][
        'csf_erosion_prop']

    erode.inputs.inputspec.mask_erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_csf'][
        'csf_mask_erosion_mm']

    node, out = strat_pool.get_data(['label-CSF_desc-preproc_mask',
                                     'label-CSF_mask'])
    wf.connect(node, out, erode, 'inputspec.mask')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, erode, 'inputspec.brain_mask')

    outputs = {
        'label-CSF_desc-eroded_mask': (erode, 'outputspec.eroded_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="erode_mask_GM",
    switch=[
        ["nuisance_corrections", "2-nuisance_regression", "create_regressors"],
        [
            "nuisance_corrections",
            "2-nuisance_regression",
            "regressor_masks",
            "erode_gm",
            "run",
        ],
    ],
    inputs=[["label-GM_desc-preproc", "label-GM_mask"]],
    outputs=["label-GM_desc-eroded_mask"],
)
def erode_mask_GM(wf, cfg, strat_pool, pipe_num, opt=None):

    erode = erode_mask(f'erode_GM_mask_{pipe_num}')
    erode.inputs.inputspec.erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_gm'][
        'gm_erosion_mm']
    erode.inputs.inputspec.erode_prop = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_gm'][
        'gm_erosion_prop']

    erode.inputs.inputspec.mask_erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_gm'][
        'gm_mask_erosion_mm']

    node, out = strat_pool.get_data(['label-GM_desc-preproc_mask',
                                     'label-GM_mask'])
    wf.connect(node, out, erode, 'inputspec.mask')

    outputs = {
        'label-GM_desc-eroded_mask': (erode, 'outputspec.eroded_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="erode_mask_WM",
    switch=[
        ["nuisance_corrections", "2-nuisance_regression", "create_regressors"],
        [
            "nuisance_corrections",
            "2-nuisance_regression",
            "regressor_masks",
            "erode_wm",
            "run",
        ],
    ],
    inputs=[
        (["label-WM_desc-preproc_mask", "label-WM_mask"], "space-T1w_desc-brain_mask")
    ],
    outputs=["label-WM_desc-eroded_mask"],
)
def erode_mask_WM(wf, cfg, strat_pool, pipe_num, opt=None):

    erode = erode_mask(f'erode_WM_mask_{pipe_num}')
    erode.inputs.inputspec.erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_wm'][
        'wm_erosion_mm']
    erode.inputs.inputspec.erode_prop = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_wm'][
        'wm_erosion_prop']

    erode.inputs.inputspec.mask_erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_wm'][
        'wm_mask_erosion_mm']

    node, out = strat_pool.get_data(['label-WM_desc-preproc_mask',
                                     'label-WM_mask'])
    wf.connect(node, out, erode, 'inputspec.mask')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, erode, 'inputspec.brain_mask')

    outputs = {
        'label-WM_desc-eroded_mask': (erode, 'outputspec.eroded_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="nuisance_regressors_generation_EPItemplate",
    config=["nuisance_corrections", "2-nuisance_regression"],
    switch=["create_regressors"],
    option_key="Regressors",
    option_val="USER-DEFINED",
    inputs=[
        (
            "desc-preproc_bold",
            "desc-brain_bold",
            "space-bold_desc-brain_mask",
            "desc-movementParameters_motion",
            "framewise-displacement-jenkinson",
            "framewise-displacement-power",
            "dvars",
            ["space-bold_desc-eroded_mask", "space-bold_desc-brain_mask"],
            [
                "space-bold_label-CSF_desc-eroded_mask",
                "space-bold_label-CSF_desc-preproc_mask",
                "space-bold_label-CSF_mask",
            ],
            [
                "space-bold_label-WM_desc-eroded_mask",
                "space-bold_label-WM_desc-preproc_mask",
                "space-bold_label-WM_mask",
            ],
            [
                "space-bold_label-GM_desc-eroded_mask",
                "space-bold_label-GM_desc-preproc_mask",
                "space-bold_label-GM_mask",
            ],
            "from-EPItemplate_to-bold_mode-image_desc-linear_xfm",
            "from-bold_to-EPItemplate_mode-image_desc-linear_xfm",
        ),
        "lateral-ventricles-mask",
        "TR",
    ],
    outputs=["desc-confounds_timeseries", "censor-indices"],
)

def nuisance_regressors_generation_EPItemplate(wf, cfg, strat_pool, pipe_num,
                                               opt=None):
    return nuisance_regressors_generation(wf, cfg, strat_pool, pipe_num, opt,
                                          'bold')

@nodeblock(
    name="nuisance_regressors_generation_T1w",
    config=["nuisance_corrections", "2-nuisance_regression"],
    switch=["create_regressors"],
    option_key="Regressors",
    option_val="USER-DEFINED",
    inputs=[
        (
            "desc-preproc_bold",
            "space-bold_desc-brain_mask",
            "from-bold_to-T1w_mode-image_desc-linear_xfm",
            "desc-movementParameters_motion",
            "framewise-displacement-jenkinson",
            "framewise-displacement-power",
            "dvars",
            "desc-brain_T1w",
            ["space-T1w_desc-eroded_mask", "space-T1w_desc-brain_mask"],
            [
                "label-CSF_desc-eroded_mask",
                "label-CSF_desc-preproc_mask",
                "label-CSF_mask",
            ],
            [
                "label-WM_desc-eroded_mask",
                "label-WM_desc-preproc_mask",
                "label-WM_mask",
            ],
            [
                "label-GM_desc-eroded_mask",
                "label-GM_desc-preproc_mask",
                "label-GM_mask",
            ],
            "from-template_to-T1w_mode-image_desc-linear_xfm",
            "from-T1w_to-template_mode-image_desc-linear_xfm",
        ),
        "lateral-ventricles-mask",
        "TR",
    ],
    outputs=["desc-confounds_timeseries", "censor-indices"],
)
def nuisance_regressors_generation_T1w(wf, cfg, strat_pool, pipe_num, opt=None
                                       ):
    return nuisance_regressors_generation(wf, cfg, strat_pool, pipe_num, opt,
                                          'T1w')


def nuisance_regressors_generation(wf: Workflow, cfg: Configuration,
                                   strat_pool: ResourcePool,
                                   pipe_num: int, opt: dict,
                                   space: LITERAL['T1w', 'bold']
                                   ) -> TUPLE[Workflow, dict]:
    '''
    Parameters
    ----------
    wf : ~nipype.pipeline.engine.workflows.Workflow

    cfg : ~CPAC.utils.configuration.Configuration

    strat_pool : ~CPAC.pipeline.engine.ResourcePool

    pipe_num : int

    opt : dict

    space : str
        T1w or bold

    Returns
    -------
    wf : nipype.pipeline.engine.workflows.Workflow

    outputs : dict
    '''
    prefixes = [f'space-{space}_'] * 2
    reg_tool = None
    if space == 'T1w':
        prefixes[0] = ''
        if strat_pool.check_rpool(
            'from-template_to-T1w_mode-image_desc-linear_xfm'):
            xfm_prov = strat_pool.get_cpac_provenance(
                'from-template_to-T1w_mode-image_desc-linear_xfm')
            reg_tool = check_prov_for_regtool(xfm_prov)
    elif space == 'bold':
        xfm_prov = strat_pool.get_cpac_provenance(
            'from-EPItemplate_to-bold_mode-image_desc-linear_xfm')
        reg_tool = check_prov_for_regtool(xfm_prov)
    if reg_tool is not None:
        use_ants = reg_tool == 'ants'
    if cfg.switch_is_on(['functional_preproc',
                         'motion_estimates_and_correction',
                         'motion_estimate_filter', 'run']):
        wf_name = (f'nuisance_regressors_{opt["Name"]}_filt-'
                    f'{strat_pool.filter_name(cfg)}_{pipe_num}')
    else:
        wf_name = f'nuisance_regressors_{opt["Name"]}_{pipe_num}'

    ventricle = strat_pool.check_rpool('lateral-ventricles-mask')
    csf_mask = strat_pool.check_rpool([f'{prefixes[0]}label-CSF_'
                                       'desc-eroded_mask',
                                       f'{prefixes[0]}label-CSF_'
                                       'desc-preproc_mask',
                                       f'{prefixes[0]}label-CSF_mask'])

    regressors = create_regressor_workflow(opt, use_ants,
                                           ventricle_mask_exist=ventricle,
                                           all_bold=space == 'bold',
                                           csf_mask_exist=csf_mask,
                                           name=wf_name)

    node, out = strat_pool.get_data("desc-preproc_bold")
    wf.connect(node, out, regressors, 'inputspec.functional_file_path')

    node, out = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out,
               regressors, 'inputspec.functional_brain_mask_file_path')

    if strat_pool.check_rpool(f'desc-brain_{space}'):
        node, out = strat_pool.get_data(f'desc-brain_{space}')
        wf.connect(node, out, regressors, 'inputspec.anatomical_file_path')

    if strat_pool.check_rpool([f'{prefixes[1]}desc-eroded_mask',
                               f'{prefixes[1]}desc-brain_mask']):
        node, out = strat_pool.get_data([f'{prefixes[1]}desc-eroded_mask',
                                         f'{prefixes[1]}desc-brain_mask'])
        wf.connect(node, out, regressors,
                   'inputspec.anatomical_eroded_brain_mask_file_path')
    else:
        logger.warning('No %s-space brain mask found in resource pool.', space)

    if strat_pool.check_rpool([f'{prefixes[0]}label-CSF_desc-eroded_mask',
                               f'{prefixes[0]}label-CSF_desc-preproc_mask',
                               f'{prefixes[0]}label-CSF_mask']):
        node, out = strat_pool.get_data([f'{prefixes[0]}label-CSF_'
                                         'desc-eroded_mask',
                                         f'{prefixes[0]}label-CSF_'
                                         'desc-preproc_mask',
                                         f'{prefixes[0]}label-CSF_mask'])
        wf.connect(node, out, regressors, 'inputspec.csf_mask_file_path')
    else:
        logger.warning('No %s-space CSF mask found in resource pool.', space)

    if strat_pool.check_rpool([f'{prefixes[0]}label-WM_desc-eroded_mask',
                               f'{prefixes[0]}label-WM_desc-preproc_mask',
                               f'{prefixes[0]}label-WM_mask']):
        node, out = strat_pool.get_data([f'{prefixes[0]}label-WM_'
                                         'desc-eroded_mask',
                                         f'{prefixes[0]}label-WM_'
                                         'desc-preproc_mask',
                                         f'{prefixes[0]}label-WM_mask'])
        wf.connect(node, out, regressors, 'inputspec.wm_mask_file_path')
    else:
        logger.warning('No %s-space WM mask found in resource pool.', space)

    if strat_pool.check_rpool([f'{prefixes[0]}label-GM_desc-eroded_mask',
                               f'{prefixes[0]}label-GM_desc-preproc_mask',
                               f'{prefixes[0]}label-GM_mask']):
        node, out = strat_pool.get_data([f'{prefixes[0]}label-GM_'
                                         'desc-eroded_mask',
                                         f'{prefixes[0]}label-GM_'
                                         'desc-preproc_mask',
                                         f'{prefixes[0]}label-GM_mask'])
        wf.connect(node, out, regressors, 'inputspec.gm_mask_file_path')
    else:
        logger.warning('No %s-space GM mask found in resource pool.', space)

    if ventricle:
        node, out = strat_pool.get_data('lateral-ventricles-mask')
        wf.connect(node, out,
                   regressors, 'inputspec.lat_ventricles_mask_file_path')

    if space == 'T1w':
        if strat_pool.check_rpool(
                'from-bold_to-T1w_mode-image_desc-linear_xfm'):
            node, out = strat_pool.get_data(
                'from-bold_to-T1w_mode-image_desc-linear_xfm')
            wf.connect(node, out,
                       regressors,
                       'inputspec.func_to_anat_linear_xfm_file_path')

            # invert func2anat matrix to get anat2func_linear_xfm
            anat2func_linear_xfm = pe.Node(interface=fsl.ConvertXFM(),
                                           name=f'anat_to_func_linear_xfm_'
                                                f'{opt["Name"]}_{pipe_num}')
            anat2func_linear_xfm.inputs.invert_xfm = True
            wf.connect(node, out, anat2func_linear_xfm, 'in_file')

            wf.connect(anat2func_linear_xfm, 'out_file',
                       regressors,
                       'inputspec.anat_to_func_linear_xfm_file_path')

        if strat_pool.check_rpool(
                'from-template_to-T1w_mode-image_desc-linear_xfm'):
            node, out = strat_pool.get_data(
                'from-template_to-T1w_mode-image_desc-linear_xfm')
            wf.connect(node, out,
                       regressors,
                       'inputspec.mni_to_anat_linear_xfm_file_path')

        if strat_pool.check_rpool(
            'from-T1w_to-template_mode-image_desc-linear_xfm'):
            node, out = strat_pool.get_data(
                'from-T1w_to-template_mode-image_desc-linear_xfm')
            wf.connect(node, out,
                       regressors,
                       'inputspec.anat_to_mni_linear_xfm_file_path')

    elif space == 'bold':
        if strat_pool.check_rpool(
                'from-EPItemplate_to-bold_mode-image_desc-linear_xfm'):
            node, out = strat_pool.get_data(
                'from-EPItemplate_to-bold_mode-image_desc-linear_xfm')
            wf.connect(node, out,
                       regressors,
                       'inputspec.mni_to_anat_linear_xfm_file_path')
            wf.connect(node, out,
                       regressors,
                       'inputspec.anat_to_func_linear_xfm_file_path')

        if strat_pool.check_rpool(
                'from-bold_to-EPItemplate_mode-image_desc-linear_xfm'):
            node, out = strat_pool.get_data(
                'from-bold_to-EPItemplate_mode-image_desc-linear_xfm')
            wf.connect(node, out,
                       regressors,
                       'inputspec.anat_to_mni_linear_xfm_file_path')
            wf.connect(node, out,
                       regressors,
                       'inputspec.func_to_anat_linear_xfm_file_path')

    if strat_pool.check_rpool('desc-movementParameters_motion'):
        node, out = strat_pool.get_data('desc-movementParameters_motion')
        wf.connect(node, out,
                   regressors, 'inputspec.motion_parameters_file_path')

    if strat_pool.check_rpool('framewise-displacement-jenkinson'):
        node, out = strat_pool.get_data('framewise-displacement-jenkinson')
        wf.connect(node, out, regressors, 'inputspec.fd_j_file_path')

    if strat_pool.check_rpool('framewise-displacement-power'):
        node, out = strat_pool.get_data('framewise-displacement-power')
        wf.connect(node, out, regressors, 'inputspec.fd_p_file_path')

    if strat_pool.check_rpool('dvars'):
        node, out = strat_pool.get_data('dvars')
        wf.connect(node, out, regressors, 'inputspec.dvars_file_path')

    node, out = strat_pool.get_data('TR')
    wf.connect(node, out, regressors, 'inputspec.tr')

    outputs = {
        'desc-confounds_timeseries': (regressors, 'outputspec.regressors_file_path'),
        'censor-indices': (regressors, 'outputspec.censor_indices')
    }

    return (wf, outputs)


def nuisance_regression(wf, cfg, strat_pool, pipe_num, opt, space, res=None):
    '''Nuisance regression in native (BOLD) or template space

    Parameters
    ----------
    wf, cfg, strat_pool, pipe_num, opt
        pass through from Node Block

    space : str
        native or template

    Returns
    -------
    wf : nipype.pipeline.engine.workflows.Workflow

    outputs : dict
    '''
    opt = strat_pool.regressor_dct(cfg)
    bandpass = 'Bandpass' in opt
    bandpass_before = bandpass and cfg['nuisance_corrections',
                                       '2-nuisance_regression',
                                       'bandpass_filtering_order'] == 'Before'

    name_suff = (f'space-{space}_reg-{opt["Name"]}_{pipe_num}'
                 if res is None else
                 f'space-{space}_res-{res}_reg-{opt["Name"]}_{pipe_num}')
    nuis_name = f'nuisance_regression_{name_suff}'

    nuis = create_nuisance_regression_workflow(opt, name=nuis_name)
    if bandpass_before:
        nofilter_nuis = nuis.clone(name=f'{nuis.name}-noFilter')

    desc_keys = ('desc-preproc_bold', 'desc-cleaned_bold',
                 'desc-denoisedNofilt_bold')
    if space != 'native':
        new_label = f'space-{space}'
        if res:
            new_label = f'{new_label}_res-{res}'
        desc_keys = tuple(f'{new_label}_{key}' for key in desc_keys)

    if space == 'template':
        # sometimes mm dimensions match but the voxel dimensions don't
        # so here we align the mask to the resampled data before applying
        match_grid = pe.Node(afni.Resample(),
                             name='align_template_mask_to_template_data_'
                                  f'{name_suff}')
        match_grid.inputs.outputtype = 'NIFTI_GZ'
        match_grid.inputs.resample_mode = 'Cu'
        node, out = strat_pool.get_data('FSL-AFNI-brain-mask')
        wf.connect(node, out, match_grid, 'in_file')
        node, out = strat_pool.get_data(desc_keys[0])
        wf.connect(node, out, match_grid, 'master')
        wf.connect(match_grid, 'out_file',
                   nuis, 'inputspec.functional_brain_mask_file_path')
        if bandpass_before:
            wf.connect(match_grid, 'out_file',
                       nofilter_nuis,
                       'inputspec.functional_brain_mask_file_path')
    else:
        node, out = strat_pool.get_data('space-bold_desc-brain_mask')
        wf.connect(node, out,
                   nuis, 'inputspec.functional_brain_mask_file_path')
        if bandpass_before:
            wf.connect(node, out,
                       nofilter_nuis,
                       'inputspec.functional_brain_mask_file_path')

    node, out = strat_pool.get_data(['desc-confounds_timeseries', 'parsed_regressors'])
    wf.connect(node, out, nuis, 'inputspec.regressor_file')
    if bandpass_before:
        wf.connect(node, out, nofilter_nuis, 'inputspec.regressor_file')

    if strat_pool.check_rpool('framewise-displacement-jenkinson'):
        node, out = strat_pool.get_data('framewise-displacement-jenkinson')
        wf.connect(node, out, nuis, 'inputspec.fd_j_file_path')
        if bandpass_before:
            wf.connect(node, out, nofilter_nuis, 'inputspec.fd_j_file_path')

    if strat_pool.check_rpool('framewise-displacement-power'):
        node, out = strat_pool.get_data('framewise-displacement-power')
        wf.connect(node, out, nuis, 'inputspec.fd_p_file_path')
        if bandpass_before:
            wf.connect(node, out, nofilter_nuis, 'inputspec.fd_p_file_path')

    if strat_pool.check_rpool('dvars'):
        node, out = strat_pool.get_data('dvars')
        wf.connect(node, out, nuis, 'inputspec.dvars_file_path')
        if bandpass_before:
            wf.connect(node, out, nofilter_nuis, 'inputspec.dvars_file_path')

    if bandpass:
        filt = filtering_bold_and_regressors(opt, name=f'filtering_bold_and_'
                                             f'regressors_{name_suff}')
        filt.inputs.inputspec.nuisance_selectors = opt

        node, out = strat_pool.get_data(['desc-confounds_timeseries', 'parsed_regressors'])
        wf.connect(node, out, filt, 'inputspec.regressors_file_path')

        if space == 'template':
            wf.connect(match_grid, 'out_file',
                       filt, 'inputspec.functional_brain_mask_file_path')
        else:
            node, out = strat_pool.get_data('space-bold_desc-brain_mask')
            wf.connect(node, out,
                       filt, 'inputspec.functional_brain_mask_file_path')

        node, out = strat_pool.get_data('TR')
        wf.connect(node, out, filt, 'inputspec.tr')

        if cfg['nuisance_corrections', '2-nuisance_regression',
               'bandpass_filtering_order'] == 'After':

            node, out = strat_pool.get_data(desc_keys[0])
            wf.connect(node, out, nuis, 'inputspec.functional_file_path')

            wf.connect(nuis, 'outputspec.residual_file_path',
                       filt, 'inputspec.functional_file_path')

            outputs = {
                desc_keys[0]: (filt, 'outputspec.residual_file_path'),
                desc_keys[1]: (filt, 'outputspec.residual_file_path'),
                desc_keys[2]: (nuis, 'outputspec.residual_file_path'),
                'desc-confounds_timeseries': (filt, 'outputspec.residual_regressor')
            }

        elif bandpass_before:

            node, out = strat_pool.get_data(desc_keys[0])
            wf.connect(node, out, filt, 'inputspec.functional_file_path')
            wf.connect(node, out,
                       nofilter_nuis, 'inputspec.functional_file_path')

            wf.connect(filt, 'outputspec.residual_file_path',
                       nuis, 'inputspec.functional_file_path')


            outputs = {
                desc_keys[0]: (nuis, 'outputspec.residual_file_path'),
                desc_keys[1]: (nuis, 'outputspec.residual_file_path'),
                desc_keys[2]: (nofilter_nuis, 'outputspec.residual_file_path'),
                'desc-confounds_timeseries': (filt, 'outputspec.residual_regressor')}

    else:
        node, out = strat_pool.get_data(desc_keys[0])
        wf.connect(node, out, nuis, 'inputspec.functional_file_path')

        outputs = {desc_key: (nuis, 'outputspec.residual_file_path') for
                   desc_key in desc_keys}

    return (wf, outputs)

@nodeblock(
    name="ingress_regressors",
    switch=[["nuisance_corrections", "2-nuisance_regression", "run"],
            ["nuisance_corrections", "2-nuisance_regression", "ingress_regressors", "run"]],
    inputs=["pipeline-ingress_desc-confounds_timeseries"],
    outputs=["parsed_regressors"]
)
def ingress_regressors(wf, cfg, strat_pool, pipe_num, opt=None):

    regressors_list = cfg.nuisance_corrections['2-nuisance_regression']['ingress_regressors'][
        'Regressors']['Columns']
    
    # Will need to generalize the name
    node, out = strat_pool.get_data('pipeline-ingress_desc-confounds_timeseries')
    if not regressors_list:
        logger.warning("\n[!] Ingress regressors is on, but no regressors provided. "  
                                            "The whole regressors file will be applied, but it may be" 
                                            "too large for the timeseries data!")
        outputs = {
            'parsed_regressors': (node, out)
        }
    else:
        ingress_imports = ['import numpy as np',
                   'import numpy as np', 'import os',
                   'import CPAC', 'from nipype import logging',
                   'logger = logging.getLogger("nipype.workflow")']
        ingress_regressors = pe.Node(Function(
                input_names=['regressors_file',
                            'regressors_list'],
                output_names=['parsed_regressors'],
                function=parse_regressors,
                imports=ingress_imports
            ), name="parse_regressors_file")

        wf.connect(node, out, ingress_regressors, 'regressors_file')
        ingress_regressors.inputs.regressors_list = regressors_list

        outputs = {
            'parsed_regressors': (ingress_regressors, 'parsed_regressors')
        }
    
    return wf, outputs

def parse_regressors(regressors_file, regressors_list):

    """
    
    Parses regressors file from outdir ingress.
    
    Parameters
    ----------
    confounds / regressors file : string
        Path of regressors / confounds file.
    regressors list : list, can be empty
        List containing names of regressors to select

        
    Returns
    -------
    parsed_regressors: dataframe
        Regressors 
    
    """
    import pandas as pd

    with open(regressors_file, 'r') as f:
        full_file = pd.read_table(regressors_file)
        parsed_regressors = pd.DataFrame()
        header = []
        for regressor in regressors_list:
            # Look through first 3 rows in case the header is nonstandard
            if regressor in full_file.iloc[:3]:
                header.append(regressor)
                parsed_regressors[regressor] = full_file.loc[:,regressor]
            else:
                logger.warning(f"\n[!] Regressor {regressor} not found in {regressors_file}")
    if parsed_regressors.empty:
        raise Exception(f"Regressors not found in {regressors_file}")

    regressors_path = os.path.join(os.getcwd(), "parsed_regressors.1D")
    parsed_regressors = parsed_regressors.to_numpy()
    check_vals = np.any(np.isnan(parsed_regressors))
    if check_vals:
        raise Exception('\n[!] This regressors file contains "N/A" values.\n' 
                            '[!] Please choose a different dataset or ' 
                                        'remove regressors with those values.')
    with open(regressors_path, "w") as ofd:

        # write out the header information
        ofd.write("# C-PAC {0}\n".format(CPAC.__version__))
        ofd.write("# Ingressed nuisance regressors:\n")
        np.savetxt(ofd, parsed_regressors, fmt='%.18f', delimiter='\t')

    return regressors_path

@nodeblock(
    name="nuisance_regression_native",
    config=["nuisance_corrections", "2-nuisance_regression"],
    switch=["run"],
    option_key="space",
    option_val="native",
    inputs=[
        (
            "desc-preproc_bold",
            ["desc-confounds_timeseries", "parsed_regressors"],
            "space-bold_desc-brain_mask",
            "framewise-displacement-jenkinson",
            "framewise-displacement-power",
            "dvars",
        ),
        "TR",
    ],
    outputs={
        "desc-preproc_bold": {
            "Description": "Preprocessed BOLD image that was nuisance-regressed in native space"
        },
        "desc-cleaned_bold": {
            "Description": "Preprocessed BOLD image that was nuisance-regressed in native space"
        },
        "desc-denoisedNofilt_bold": {
            "Description": "Preprocessed BOLD image that was nuisance-regressed in native space, but without frequency filtering."
        },
        "desc-confounds_timeseries": {"Description": "Regressors that were applied in native space"},
    },
)
def nuisance_regression_native(wf, cfg, strat_pool, pipe_num, opt=None):

    return nuisance_regression(wf, cfg, strat_pool, pipe_num, opt, 'native')


@nodeblock(
    name="nuisance_regression_template",
    config=["nuisance_corrections", "2-nuisance_regression"],
    switch=["run"],
    option_key="space",
    option_val="template",
    inputs=[
        (
            "desc-stc_bold",
            "space-template_desc-preproc_bold",
            "space-template_res-derivative_desc-preproc_bold",
            "desc-movementParameters_motion",
            ["desc-confounds_timeseries", "parsed_regressors"],
            "FSL-AFNI-brain-mask",
            "framewise-displacement-jenkinson",
            "framewise-displacement-power",
            "dvars",
        ),
        "TR",
    ],
    outputs={
        "space-template_desc-preproc_bold": {
            "Description": "Preprocessed BOLD image that was nuisance-regressed in template space"
        },
        "space-template_res-derivative_desc-preproc_bold": {
            "Description": "Preprocessed BOLD image that was nuisance-regressed in template space"
        },
        "space-template_desc-cleaned_bold": {
            "Description": "Preprocessed BOLD image that was nuisance-regressed in template space"
        },
        "space-template_res-derivative_desc-cleaned_bold": {
            "Description": "Preprocessed BOLD image that was nuisance-regressed in template space"
        },
        "space-template_desc-denoisedNofilt_bold": {
            "Description": "Preprocessed BOLD image that was nuisance-regressed in template space, but without frequency filtering."
        },
        "space-template_res-derivative_desc-denoisedNofilt_bold": {
            "Description": "Preprocessed BOLD image that was nuisance-regressed in template space, but without frequency filtering."
        },
        "desc-confounds_timeseries": {"Description": "Regressors that were applied in template space"},
    },
)
def nuisance_regression_template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Apply nuisance regression to template-space image'''
    wf, outputs = nuisance_regression(wf, cfg, strat_pool, pipe_num, opt,
                                      'template')
    if strat_pool.check_rpool(
        'space-template_res-derivative_desc-preproc_bold'
    ):
        wf, res_outputs = nuisance_regression(wf, cfg, strat_pool, pipe_num,
                                              opt, 'template', 'derivative')
        outputs.update(res_outputs)
    return (wf, outputs)


@nodeblock(
    name="erode_mask_bold",
    switch=[
        [
            "nuisance_corrections",
            "2-nuisance_regression",
            "regressor_masks",
            "erode_anatomical_brain_mask",
            "run",
        ],
        [
            "registration_workflows",
            "functional_registration",
            "EPI_registration",
            "run",
        ],
    ],
    inputs=[
        (
            "space-bold_desc-brain_mask",
            ["space-bold_label-CSF_desc-preproc_mask", "space-bold_label-CSF_mask"],
        )
    ],
    outputs=["space-bold_desc-eroded_mask"],
)
def erode_mask_bold(wf, cfg, strat_pool, pipe_num, opt=None):

    erode = erode_mask(f'erode_T1w_mask_{pipe_num}', segmentmap=False)
    erode.inputs.inputspec.mask_erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks'][
        'erode_anatomical_brain_mask']['brain_mask_erosion_mm']
    erode.inputs.inputspec.erode_prop = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks'][
        'erode_anatomical_brain_mask']['brain_mask_erosion_prop']

    node, out = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out, erode, 'inputspec.brain_mask')

    node, out = strat_pool.get_data(['space-bold_label-CSF_desc-preproc_mask',
                                     'space-bold_label-CSF_mask'])
    wf.connect(node, out, erode, 'inputspec.mask')

    outputs = {
        'space-bold_desc-eroded_mask': (erode, 'outputspec.eroded_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="erode_mask_boldCSF",
    switch=[
        [
            "nuisance_corrections",
            "2-nuisance_regression",
            "regressor_masks",
            "erode_csf",
            "run",
        ],
        [
            "registration_workflows",
            "functional_registration",
            "EPI_registration",
            "run",
        ],
    ],
    inputs=[
        (
            ["space-bold_label-CSF_desc-preproc_mask", "space-bold_label-CSF_mask"],
            "space-bold_desc-brain_mask",
        )
    ],
    outputs=["space-bold_label-CSF_desc-eroded_mask"],
)
def erode_mask_boldCSF(wf, cfg, strat_pool, pipe_num, opt=None):

    erode = erode_mask(f'erode_CSF_mask_{pipe_num}')
    erode.inputs.inputspec.erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_csf'][
        'csf_erosion_mm']
    erode.inputs.inputspec.erode_prop = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_csf'][
        'csf_erosion_prop']

    erode.inputs.inputspec.mask_erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_csf'][
        'csf_mask_erosion_mm']

    node, out = strat_pool.get_data(['space-bold_label-CSF_desc-preproc_mask',
                                     'space-bold_label-CSF_mask'])
    wf.connect(node, out, erode, 'inputspec.mask')

    node, out = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out, erode, 'inputspec.brain_mask')

    outputs = {
        'space-bold_label-CSF_desc-eroded_mask': (erode, 'outputspec.eroded_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="erode_mask_boldGM",
    switch=[
        [
            "nuisance_corrections",
            "2-nuisance_regression",
            "regressor_masks",
            "erode_gm",
            "run",
        ],
        [
            "registration_workflows",
            "functional_registration",
            "EPI_registration",
            "run",
        ],
    ],
    inputs=[["space-bold_label-GM_desc-preproc", "space-bold_label-GM_mask"]],
    outputs=["space-bold_label-GM_desc-eroded_mask"],
)
def erode_mask_boldGM(wf, cfg, strat_pool, pipe_num, opt=None):

    erode = erode_mask(f'erode_GM_mask_{pipe_num}')
    erode.inputs.inputspec.erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_gm'][
        'gm_erosion_mm']
    erode.inputs.inputspec.erode_prop = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_gm'][
        'gm_erosion_prop']

    erode.inputs.inputspec.mask_erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_gm'][
        'gm_mask_erosion_mm']

    node, out = strat_pool.get_data(['space-bold_label-GM_desc-preproc_mask',
                                     'space-bold_label-GM_mask'])
    wf.connect(node, out, erode, 'inputspec.mask')

    outputs = {
        'space-bold_label-GM_desc-eroded_mask': (erode, 'outputspec.eroded_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="erode_mask_boldWM",
    switch=[
        [
            "nuisance_corrections",
            "2-nuisance_regression",
            "regressor_masks",
            "erode_wm",
            "run",
        ],
        [
            "registration_workflows",
            "functional_registration",
            "EPI_registration",
            "run",
        ],
    ],
    inputs=[
        (
            ["space-bold_label-WM_desc-preproc_mask", "space-bold_label-WM_mask"],
            "space-bold_desc-brain_mask",
        )
    ],
    outputs=["space-bold_label-WM_desc-eroded_mask"],
)
def erode_mask_boldWM(wf, cfg, strat_pool, pipe_num, opt=None):

    erode = erode_mask(f'erode_WM_mask_{pipe_num}')
    erode.inputs.inputspec.erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_wm'][
        'wm_erosion_mm']
    erode.inputs.inputspec.erode_prop = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_wm'][
        'wm_erosion_prop']

    erode.inputs.inputspec.mask_erode_mm = cfg.nuisance_corrections[
        '2-nuisance_regression']['regressor_masks']['erode_wm'][
        'wm_mask_erosion_mm']

    node, out = strat_pool.get_data(['space-bold_label-WM_desc-preproc_mask',
                                     'space-bold_label-WM_mask'])
    wf.connect(node, out, erode, 'inputspec.mask')

    node, out = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out, erode, 'inputspec.brain_mask')

    outputs = {
        'space-bold_label-WM_desc-eroded_mask': (erode,
                                                 'outputspec.eroded_mask')
    }

    return (wf, outputs)
