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
import os
import sys
import collections.abc
import fnmatch
import gzip
import json
import numbers
import pickle
import numpy as np
import yaml

from copy import deepcopy
from itertools import repeat
from voluptuous.error import Invalid

CONFIGS_DIR = os.path.abspath(os.path.join(
    __file__, *repeat(os.path.pardir, 2), 'resources/configs/'))
with open(os.path.join(CONFIGS_DIR, '1.7-1.8-nesting-mappings.yml'), 'r',
          encoding='utf-8') as _f:
    NESTED_CONFIG_MAPPING = yaml.safe_load(_f)
with open(os.path.join(CONFIGS_DIR, '1.7-1.8-deprecations.yml'), 'r',
          encoding='utf-8') as _f:
    NESTED_CONFIG_DEPRECATIONS = yaml.safe_load(_f)
YAML_BOOLS = {True: ('on', 't', 'true', 'y', 'yes'),
              False: ('f', 'false', 'n', 'no', 'off')}


def get_last_prov_entry(prov):
    while not isinstance(prov[-1], str):
        prov = prov[-1]
    return prov[-1]


def check_prov_for_regtool(prov):
    last_entry = get_last_prov_entry(prov)
    last_node = last_entry.split(':')[1]
    if 'ants' in last_node.lower():
        return 'ants'
    elif 'fsl' in last_node.lower():
        return 'fsl'
    else:
        # go further back in case we're checking against a concatenated
        # downstream xfm like bold-to-template (and prov is the provenance of
        # that downstream xfm)
        if 'from-T1w_to-template_mode-image_xfm:' in str(prov):
            splitprov = str(prov).split(
                'from-T1w_to-template_mode-image_xfm:')
            node_name = splitprov[1].split("']")[0]
            if 'ANTs' in node_name:
                return 'ants'
            elif 'FSL' in node_name:
                return 'fsl'
        elif 'from-bold_to-template_mode-image_xfm:' in str(prov):
            splitprov = str(prov).split(
                'from-bold_to-template_mode-image_xfm:')
            node_name = splitprov[1].split("']")[0]
            if 'ANTs' in node_name:
                return 'ants'
            elif 'FSL' in node_name:
                return 'fsl'
            else:
                return None
        elif 'from-T1w_to-symtemplate_mode-image_xfm:' in str(prov):
            splitprov = str(prov).split(
                'from-T1w_to-symtemplate_mode-image_xfm:')
            node_name = splitprov[1].split("']")[0]
            if 'ANTs' in node_name:
                return 'ants'
            elif 'FSL' in node_name:
                return 'fsl'
        elif 'from-bold_to-symtemplate_mode-image_xfm:' in str(prov):
            splitprov = str(prov).split(
                'from-bold_to-symtemplate_mode-image_xfm:')
            node_name = splitprov[1].split("']")[0]
            if 'ANTs' in node_name:
                return 'ants'
            elif 'FSL' in node_name:
                return 'fsl'
            else:
                return None
        else:
            return None


def check_prov_for_motion_tool(prov):
    last_entry = get_last_prov_entry(prov)
    last_node = last_entry.split(':')[1]
    if '3dvolreg' in last_node.lower():
        return '3dvolreg'
    elif 'mcflirt' in last_node.lower():
        return 'mcflirt'
    else:
        # check entire prov
        if '3dvolreg' in str(prov):
            return '3dvolreg'
        elif 'mcflirt' in str(prov):
            return 'mcflirt'
        else:
            return None





def get_flag(in_flag):
    return in_flag


def get_flag_wf(wf_name='get_flag'):
    from CPAC.pipeline import nipype_pipeline_engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['in_flag']),
                         name='inputspec')

    get_flag = pe.Node(util.Function(input_names=['in_flag'],
                                     function=get_flag),
                       name='get_flag')

    wf.connect(input_node, 'in_flag', get_flag, 'in_flag')


def read_json(json_file):
    try:
        with open(json_file, 'r') as f:
            json_dct = json.load(f)
    except json.decoder.JSONDecodeError as err:
        raise Exception(f'\n\n{err}\n\nJSON file: {json_file}\n')
    return json_dct


def create_id_string(cfg, unique_id, resource, scan_id=None,
                     template_desc=None, atlas_id=None, fwhm=None,
                     subdir=None, extension=None):
    """Create the unique key-value identifier string for BIDS-Derivatives
    compliant file names.

    This is used in the file renaming performed during the Datasink
    connections.

    Example
    -------
    >>> from CPAC.utils.configuration import Configuration
    >>> create_id_string(Configuration(), 'sub-1_ses-1',
    ...                  'res-derivative_desc-Mean-1_timeseries',
    ...                  scan_id='rest', atlas_id='Yeo_desc-7')
    'sub-1_ses-1_task-rest_atlas-Yeo7_res-3mm_desc-Mean1_timeseries'
    """
    import re
    from CPAC.utils.bids_utils import combine_multiple_entity_instances, \
                                      res_in_filename
    if atlas_id:
        if '_desc-' in atlas_id:
            atlas, desc = atlas_id.split('_desc-')
            if not re.match(r'.*[0-9]$', atlas) and re.match(r'[a-z].*', desc):
                atlas_id = f'{atlas}{desc[0].upper()}{desc[1:]}'
            else:
                atlas_id = atlas_id.replace('_desc-', '')
        resource = f'atlas-{atlas_id}_{resource}'

    part_id = unique_id.split('_')[0]
    ses_id = unique_id.split('_')[1]
    if 'sub-' not in part_id:
        part_id = f'sub-{part_id}'
    if 'ses-' not in ses_id:
        ses_id = f'ses-{ses_id}'
    if scan_id:
        out_filename = f'{part_id}_{ses_id}_task-{scan_id}_{resource}'
    else:
        out_filename = f'{part_id}_{ses_id}_{resource}'

    template_tag = template_desc.split(' -')[0] if template_desc else '*'
    for prefix in ['space-', 'from-', 'to-']:
        for bidstag in out_filename.split('_'):
            if prefix in bidstag and 'template' in bidstag:
                out_filename = out_filename.replace(
                    bidstag, f'{prefix}{template_tag}')

    if fwhm:
        for tag in resource.split('_'):
            if 'desc-' in tag and '-sm' in tag:
                newtag = tag.replace('-sm', f'-sm{fwhm}')
                out_filename = out_filename.replace(tag, newtag)
                break
        else:
            raise Exception('\n[!] FWHM provided but no desc-sm?\n')
    
    
    if extension is not None: 
       out_filename = out_filename + "." + str(extension)
       

    # drop space- entities from from native-space filenames
    if subdir == 'anat':
        out_filename = out_filename.replace('_space-T1w_', '_')
    if subdir == 'func':
        out_filename = out_filename.replace('_space-bold_', '_')
    return combine_multiple_entity_instances(
        res_in_filename(cfg, out_filename))


def write_output_json(json_data, filename, indent=3, basedir=None):
    if not basedir:
        basedir = os.getcwd()
    if '.gii' in filename:
        filename = os.path.splitext(filename)[0]
        filename = f'{filename}.json'
    if '.json' not in filename:
        filename = f'{filename}.json'
    
    json_file = os.path.join(basedir, filename)
    json_data = json.dumps(json_data, indent=indent, sort_keys=True)
    with open(json_file, 'wt') as f:
        f.write(json_data)
    return json_file


def get_zscore(map_node=False, wf_name='z_score'):
    """
    Workflow to calculate z-scores

    Parameters
    ----------
    wf_name : string
        name of the workflow

    Returns
    -------
    wf : workflow object

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/network_centrality/z_score.py>`_

    Workflow Inputs::

        inputspec.input_file : string
            path to input functional derivative file for which z score has to be calculated
        inputspec.mask_file : string
            path to whole brain functional mask file required to calculate zscore

    Workflow Outputs::

        outputspec.z_score_img : string
             path to image containing Normalized Input Image Z scores across full brain.

    .. exec::
        from CPAC.utils import get_zscore
        wf = get_zscore('mean')
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/zscore.dot'
        )

    High Level Workflow Graph:

    .. image:: ../../images/generated/zscore.png
       :width: 500


    Detailed Workflow Graph:

    .. image:: ../../images/generated/zscore_detailed.png
       :width: 500

    Example
    -------
    >>> wf = get_zscore('example_input')
    >>> wf.inputs.inputspec.input_file = '/home/data/graph_working_dir/calculate_centrality/degree_centrality_binarize.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/graphs/GraphGeneration/new_mask_3m.nii.gz'
    >>> wf.run()  # doctest: +SKIP
    """  # noqa: E501  # pylint: disable=line-too-long
    # pylint: disable=import-outside-toplevel,redefined-outer-name,reimported
    from CPAC.pipeline import nipype_pipeline_engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl as fsl

    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['input_file',
                                                       'mask_file']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['z_score_img']),
                         name='outputspec')

    if map_node:
        mean = pe.MapNode(interface=fsl.ImageStats(),
                          name='mean',
                          iterfield=['in_file'])

        standard_deviation = pe.MapNode(interface=fsl.ImageStats(),
                                        name='standard_deviation',
                                        iterfield=['in_file'])

        op_string = pe.MapNode(util.Function(input_names=['mean', 'std_dev'],
                                             output_names=['op_string'],
                                             function=get_operand_string),
                               name='op_string',
                               iterfield=['mean', 'std_dev'])

        z_score = pe.MapNode(interface=fsl.MultiImageMaths(),
                             name='z_score',
                             iterfield=['in_file', 'op_string'])

    else:
        mean = pe.Node(interface=fsl.ImageStats(), name='mean')

        standard_deviation = pe.Node(interface=fsl.ImageStats(),
                                     name='standard_deviation')

        op_string = pe.Node(util.Function(input_names=['mean', 'std_dev'],
                                          output_names=['op_string'],
                                          function=get_operand_string),
                            name='op_string')

        z_score = pe.Node(interface=fsl.MultiImageMaths(), name='z_score')

    # calculate the mean
    mean.inputs.op_string = '-k %s -m'
    wflow.connect(inputNode, 'input_file', mean, 'in_file')
    wflow.connect(inputNode, 'mask_file', mean, 'mask_file')

    # calculate the standard deviation
    standard_deviation.inputs.op_string = '-k %s -s'
    wflow.connect(inputNode, 'input_file', standard_deviation, 'in_file')
    wflow.connect(inputNode, 'mask_file', standard_deviation, 'mask_file')

    # calculate the z-score
    wflow.connect(mean, 'out_stat', op_string, 'mean')
    wflow.connect(standard_deviation, 'out_stat', op_string, 'std_dev')

    # z_score.inputs.out_file = input_name + '_zstd.nii.gz'

    wflow.connect(op_string, 'op_string', z_score, 'op_string')
    wflow.connect(inputNode, 'input_file', z_score, 'in_file')
    wflow.connect(inputNode, 'mask_file', z_score, 'operand_files')

    wflow.connect(z_score, 'out_file', outputNode, 'z_score_img')

    return wflow


def get_fisher_zscore(input_name, map_node=False, wf_name='fisher_z_score'):
    """
    Runs the compute_fisher_z_score function as part of a one-node workflow.
    """

    from CPAC.pipeline import nipype_pipeline_engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl as fsl

    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['correlation_file',
                                                       'timeseries_one_d']),
                        name='inputspec')

    outputNode = pe.Node(
        util.IdentityInterface(fields=['fisher_z_score_img']),
        name='outputspec')

    if map_node:
        # node to separate out
        fisher_z_score = pe.MapNode(
            util.Function(input_names=['correlation_file',
                                       'timeseries_one_d',
                                       'input_name'],
                          output_names=['out_file'],
                          function=compute_fisher_z_score),
            name='fisher_z_score',
            iterfield=['correlation_file'])
    else:
        fisher_z_score = pe.Node(
            util.Function(input_names=['correlation_file',
                                       'timeseries_one_d',
                                       'input_name'],
                          output_names=['out_file'],
                          function=compute_fisher_z_score),
            name='fisher_z_score')

    fisher_z_score.inputs.input_name = input_name

    wflow.connect(inputNode, 'correlation_file',
                  fisher_z_score, 'correlation_file')
    wflow.connect(inputNode, 'timeseries_one_d',
                  fisher_z_score, 'timeseries_one_d')
    wflow.connect(fisher_z_score, 'out_file',
                  outputNode, 'fisher_z_score_img')

    return wflow


def compute_fisher_z_score(correlation_file, timeseries_one_d, input_name):
    """
    Computes the fisher z transform of the input correlation map
    If the correlation map contains data for multiple ROIs then
    the function returns z score for each ROI as a seperate nifti
    file


    Parameters
    ----------

    correlation_file : string
        Input correlations file


    Returns
    -------

    out_file : list (nifti files)
        list of z_scores for mask or ROI
    """

    import nibabel as nb
    import numpy as np
    import os

    if isinstance(timeseries_one_d, str):
        if '.1D' in timeseries_one_d or '.csv' in timeseries_one_d:
            timeseries_file = timeseries_one_d

    else:
        for timeseries in timeseries_one_d:
            if '.1D' in timeseries or '.csv' in timeseries:
                timeseries_file = timeseries

    # get the specific roi number
    filename = correlation_file.split("/")[-1]
    filename = filename.replace(".nii", "")
    if ".gz" in filename:
        filename = filename.replace(".gz", "")

    corr_img = nb.load(correlation_file)
    corr_data = corr_img.get_fdata()

    hdr = corr_img.header

    # calculate the Fisher r-to-z transformation
    corr_data = np.log((1 + corr_data) / (1 - corr_data)) / 2.0

    z_score_img = nb.Nifti1Image(corr_data, header=hdr,
                                 affine=corr_img.affine)

    out_file = os.path.join(os.getcwd(), filename + '_fisher_zstd.nii.gz')

    z_score_img.to_filename(out_file)

    return out_file


def get_operand_string(mean, std_dev):
    """
    Method to get operand string for Fsl Maths

    Parameters
    ----------
    mean : string
        path to img containing mean
    std_dev : string
        path to img containing standard deviation

    Returns
    ------
    op_string : string
        operand string
    """

    str1 = "-sub %f -div %f" % (float(mean), float(std_dev))
    op_string = str1 + " -mas %s"
    return op_string


def get_roi_num_list(timeseries_file, prefix=None):
    # extracts the ROI labels from the 3dROIstats output CSV file
    with open(timeseries_file, "r") as f:
        roi_file_lines = f.read().splitlines()

    roi_err = "\n\n[!] The output of 3dROIstats, used in extracting the " \
              "timeseries, is either empty, or not in the expected " \
              "format.\n\nROI output file: {0}\n\nIf there are no rows " \
              "in the output file, double-check your ROI/mask selection." \
              "\n\n".format(str(timeseries_file))

    for line in roi_file_lines:
        if "Mean_" in line:
            try:
                roi_list = line.split(",")
                # clear out any blank strings/non ROI labels in the list
                roi_list = [x for x in roi_list if "Mean" in x]
                # rename labels
                roi_list = [
                    x.replace("Mean", "ROI").replace(" ", "").replace("#", "")
                    for x in roi_list]
            except:
                raise Exception(roi_err)
            break
    else:
        raise Exception(roi_err)

    if prefix:
        temp_rois = []
        for roi in roi_list:
            roi = prefix + "_" + str(roi)
            temp_rois.append(roi)
        roi_list = temp_rois

    return roi_list


def safe_shape(*vol_data):
    """
    Checks if the volume (first three dimensions) of multiple ndarrays
    are the same shape.

    Parameters
    ----------
    vol_data0, vol_data1, ..., vol_datan : ndarray
        Volumes to check

    Returns
    -------
    same_volume : bool
        True only if all volumes have the same shape.
    """
    same_volume = True

    first_vol_shape = vol_data[0].shape[:3]
    for vol in vol_data[1:]:
        same_volume &= (first_vol_shape == vol.shape[:3])

    return same_volume


def extract_one_d(list_timeseries):
    if isinstance(list_timeseries, str):
        if '.1D' in list_timeseries or '.csv' in list_timeseries:
            return list_timeseries

    for timeseries in list_timeseries:
        if '.1D' in timeseries or '.csv' in timeseries:
            return timeseries

    raise Exception("Unable to retrieve roi timeseries 1D or csv"
                    " file. Files found:" + list_timeseries)


def extract_txt(list_timeseries):
    """
    Method to extract txt file containing
    roi timeseries required for dual regression
    """
    if isinstance(list_timeseries, str):
        if list_timeseries.endswith('.txt'):
            return list_timeseries

    out_file = None
    for timeseries in list_timeseries:
        if timeseries.endswith('.txt'):
            out_file = timeseries

    if not out_file:
        raise Exception("Unable to retrieve roi timeseries txt"
                        " file required for dual regression."
                        " Existing files are:%s" % (list_timeseries))

    return out_file


def zscore(data, axis):
    data = data.copy()
    data -= data.mean(axis=axis, keepdims=True)
    data /= data.std(axis=axis, keepdims=True)
    np.copyto(data, 0.0, where=np.isnan(data))
    return data


def correlation(matrix1, matrix2,
                match_rows=False, z_scored=False, symmetric=False):
    d1 = matrix1.shape[-1]
    d2 = matrix2.shape[-1]

    assert d1 == d2
    assert matrix1.ndim <= 2
    assert matrix2.ndim <= 2
    if match_rows:
        assert matrix1.shape == matrix2.shape

    var = np.sqrt(d1 * d2)

    if not z_scored:
        matrix1 = zscore(matrix1, matrix1.ndim - 1)
        matrix2 = zscore(matrix2, matrix2.ndim - 1)

    if match_rows:
        return np.einsum('...i,...i', matrix1, matrix2) / var

    if matrix1.ndim >= matrix2.ndim:
        r = np.dot(matrix1, matrix2.T) / var
    else:
        r = np.dot(matrix2, matrix1.T) / var

    r = np.clip(r, -1.0, 1.0)

    if symmetric:
        return (r + r.T) / 2

    return r


def check(params_dct, subject_id, scan_id, val_to_check, throw_exception):
    if val_to_check not in params_dct:
        if throw_exception:
            raise Exception("Missing Value for {0} for participant "
                            "{1}".format(val_to_check, subject_id))
        return None

    if isinstance(params_dct[val_to_check], dict):
        ret_val = params_dct[val_to_check][scan_id]
    else:
        ret_val = params_dct[val_to_check]

    if ret_val == 'None':
        if throw_exception:
            raise Exception("'None' Parameter Value for {0} for participant "
                            "{1}".format(val_to_check, subject_id))
        else:
            ret_val = None

    if ret_val == '' and throw_exception:
        raise Exception("Missing Value for {0} for participant "
                        "{1}".format(val_to_check, subject_id))

    return ret_val


def check_random_state(seed):
    """
    Turn seed into a np.random.RandomState instance
    Code from scikit-learn (https://github.com/scikit-learn/scikit-learn)

    Parameters
    ----------
    seed : None | int | instance of RandomState
        If seed is None, return the RandomState singleton used by np.random.
        If seed is an int, return a new RandomState instance seeded with seed.
        If seed is already a RandomState instance, return it.
        Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)


def try_fetch_parameter(scan_parameters, subject, scan, keys):
    scan_parameters = dict(
        (k.lower(), v)
        for k, v in scan_parameters.items()
    )

    for key in keys:

        key = key.lower()

        if key not in scan_parameters:
            continue

        if isinstance(scan_parameters[key], dict):
            value = scan_parameters[key][scan]
        else:
            value = scan_parameters[key]

        # Explicit none value
        if value == 'None':
            return None

        if value is not None:
            return value
    return None


def get_scan_params(subject_id, scan, pipeconfig_start_indx,
                    pipeconfig_stop_indx, data_config_scan_params=None):
    """
    Method to extract slice timing correction parameters
    and scan parameters.

    Parameters
    ----------
    subject_id : str
        subject id
    scan : str
        scan id
    pipeconfig_start_indx : int
        starting volume index as provided in the pipeline config yaml file
    pipeconfig_stop_indx : int
        ending volume index as provided in the pipeline config yaml file
    data_config_scan_params : str
        file path to scan parameter JSON file listed in data config yaml file

    Returns
    -------
    TR : a string
        TR value
    pattern : a string
        slice aquisition pattern string or file path
    ref_slice : an integer
        reference slice which is used to allign all other slices
    first_tr : an integer
        starting TR or starting volume index
    last_tr : an integer
        ending TR or ending volume index
    pe_direction : str
    effective_echo_spacing : float
    """
    
    import os
    import json
    import warnings

    check2 = lambda val: val if val == None or val == '' or \
                                isinstance(val, str) else int(val)

    # initialize vars to empty
    TR = ''
    TE = None
    pattern = ''
    ref_slice = ''
    first_tr = ''
    last_tr = ''
    unit = 's'
    pe_direction = ''
    effective_echo_spacing = None
    template = None

    if isinstance(pipeconfig_stop_indx, str):
        if "End" in pipeconfig_stop_indx or "end" in pipeconfig_stop_indx:
            pipeconfig_stop_indx = None
    if data_config_scan_params:
        if ".json" in data_config_scan_params:
            if not os.path.exists(data_config_scan_params):
                err = "\n[!] WARNING: Scan parameters JSON file listed in " \
                      "your data configuration file does not exist:\n{0}" \
                      "\n\n".format(data_config_scan_params)
                raise Exception(err)

            with open(data_config_scan_params, "r") as f:
                params_dct = json.load(f)

            # get details from the configuration
            # if this is a JSON file, the key values are the BIDS format
            # standard
            # TODO: better handling of errant key values!!!
            if "RepetitionTime" in params_dct.keys():
                TR = float(check(params_dct, subject_id, scan,
                                 "RepetitionTime", False))
            if "SliceTiming" in params_dct.keys():
                pattern = str(check(params_dct, subject_id, scan,
                                    "SliceTiming", False))
            elif "SliceAcquisitionOrder" in params_dct.keys():
                pattern = str(check(params_dct, subject_id, scan,
                                    "SliceAcquisitionOrder", False))
            if "PhaseEncodingDirection" in params_dct.keys():
                pe_direction = str(check(params_dct, subject_id, scan,
                                         "PhaseEncodingDirection", False))
            try:
                "EffectiveEchoSpacing" in params_dct.keys()
                effective_echo_spacing = float(
                        check(params_dct, subject_id, scan,
                            "EffectiveEchoSpacing", False))
            except TypeError:
                pass

        elif len(data_config_scan_params) > 0 and \
                isinstance(data_config_scan_params, dict):

            params_dct = data_config_scan_params

            # TODO: better handling of errant key values!!!
            # TODO: use schema validator to deal with it
            # get details from the configuration
            try: 
                TR = float(
                    try_fetch_parameter(
                        params_dct,
                        subject_id,
                        scan,
                        ['TR', 'RepetitionTime']
                    )
                )
            except TypeError:
                TR = None

            try: 
                template = str(
                    try_fetch_parameter(
                        params_dct,
                        subject_id,
                        scan,
                        ['Template', 'template']
                    )
                )
            except TypeError:
                template = None

            pattern = str(
                try_fetch_parameter(
                    params_dct,
                    subject_id,
                    scan,
                    ['acquisition', 'SliceTiming', 'SliceAcquisitionOrder']
                )
            )

            ref_slice = check(params_dct, subject_id, scan, 'reference',
                              False)
            if ref_slice:
                ref_slice = int(ref_slice)

            first_tr = check(params_dct, subject_id, scan, 'first_TR', False)
            if first_tr:
                first_tr = check2(first_tr)

            last_tr = check(params_dct, subject_id, scan, 'last_TR', False)
            if last_tr:
                last_tr = check2(last_tr)

            pe_direction = check(params_dct, subject_id, scan,
                                 'PhaseEncodingDirection', False)
            try:
                effective_echo_spacing = float(
                try_fetch_parameter(params_dct, subject_id, scan,
                                    ["EffectiveEchoSpacing"]))
            except TypeError:
                pass

        else:
            err = "\n\n[!] Could not read the format of the scan parameters " \
                  "information included in the data configuration file for " \
                  f"the participant {subject_id}.\n\n"
            raise Exception(err)
    if first_tr == '' or first_tr is None:
        first_tr = pipeconfig_start_indx

    if last_tr == '' or last_tr is None:
        last_tr = pipeconfig_stop_indx

    unit = 's'

    if 'None' in pattern or 'none' in pattern:
        pattern = None

    '''
    if not pattern:
        if pipeconfig_tpattern:
            if "Use NIFTI Header" in pipeconfig_tpattern:
                pattern = ''
            else:
                pattern = pipeconfig_tpattern
    '''

    # pattern can be one of a few keywords, a filename, or blank which
    # indicates that the images header information should be used
    tpattern_file = None

    valid_patterns = ['alt+z', 'altplus', 'alt+z2', 'alt-z', 'altminus',
                      'alt-z2', 'seq+z', 'seqplus', 'seq-z', 'seqminus']
    if pattern and pattern != '' and pattern not in valid_patterns:

        if isinstance(pattern, list) or \
                ("[" in pattern and "]" in pattern and "," in pattern):
            # if we got the slice timing as a list, from a BIDS-format scan
            # parameters JSON file

            if not isinstance(pattern, list):
                pattern = pattern.replace("[", "").replace("]", "").split(",")

            slice_timings = [float(x) for x in pattern]

            # write out a tpattern file for AFNI 3dTShift
            tpattern_file = os.path.join(os.getcwd(), "tpattern.txt")
            try:
                with open(tpattern_file, "wt") as f:
                    for time in slice_timings:
                        f.write("{0}\n".format(time).replace(" ", ""))
            except:
                err = "\n[!] Could not write the slice timing file meant as " \
                      "an input for AFNI 3dTshift (slice timing correction):" \
                      "\n{0}\n\n".format(tpattern_file)
                raise Exception(err)

        elif ".txt" in pattern and not os.path.exists(pattern):
            # if the user provided an acquisition pattern text file for
            # 3dTshift
            raise Exception("Invalid Pattern file path {0}, Please provide "
                            "the correct path".format(pattern))
        elif ".txt" in pattern:
            with open(pattern, "r") as f:
                lines = f.readlines()
            if len(lines) < 2:
                raise Exception('Invalid slice timing file format. The file '
                                'should contain only one value per row. Use '
                                'new line char as delimiter')
            tpattern_file = pattern
            slice_timings = [float(l.rstrip('\r\n')) for l in lines]
        else:
            # this only happens if there is a non-path string set in the data
            # config dictionary for acquisition pattern (like "alt+z"), except
            # the pattern is not listed in that list
            err = "\n[!] The slice timing acquisition pattern provided is " \
                  "not supported by AFNI 3dTshift:\n" \
                  "{0}\n".format(str(pattern))
            raise Exception(err)

        pattern = tpattern_file

        slice_timings.sort()
        max_slice_offset = slice_timings[-1]

        # checking if the unit of TR and slice timing match or not
        # if slice timing in ms convert TR to ms as well
        if TR and max_slice_offset > TR:
            warnings.warn("TR is in seconds and slice timings are in "
                          "milliseconds. Converting TR into milliseconds")
            TR = TR * 1000
            print("New TR value {0} ms".format(TR))
            unit = 'ms'

    else:
        # check to see, if TR is in milliseconds, convert it into seconds
        if TR and TR > 10:
            warnings.warn('TR is in milliseconds, Converting it into seconds')
            TR = TR / 1000.0
            print("New TR value {0} s".format(TR))
            unit = 's'

    # swap back in
    if TR:
        tr = "{0}{1}".format(str(TR), unit)
    else:
        tr = ""

    tpattern = pattern
    start_indx = first_tr
    stop_indx = last_tr

    return (tr if tr else None,
            tpattern if tpattern else None,
            template if template else None,
            ref_slice,
            start_indx,
            stop_indx,
            pe_direction,
            effective_echo_spacing)


def get_tr(tr):
    """
    Method to return TR in seconds
    """
    import re
    if tr:
        tr = re.search(r"\d+.\d+", str(tr)).group(0)
        tr = float(tr)
        if tr > 10:
            tr = tr / 1000.0
    else:
        tr = ""
    return tr


def check_tr(tr, in_file):
    # imageData would have to be the image data from the funcFlow workflow,
    # funcFlow outputspec.subject
    import nibabel as nib
    img = nib.load(in_file)

    # get header from image data, then extract TR information, TR is fourth
    # item in list returned by get_zooms()
    imageHeader = img.header
    imageZooms = imageHeader.get_zooms()
    header_tr = imageZooms[3]

    # If the TR information from header_tr (funcFlow) and convert_tr node
    # (TR from config file) do not match, prepare to update the TR information
    # from either convert_tr or header_tr using afni 3drefit, then append to
    # func_to_mni
    if header_tr != tr:
        if tr != None and tr != "":
            TR = tr
        else:
            TR = header_tr

        import warnings
        warnings.warn('Warning: The TR information does not match between '
                      'the config and subject list files.')

    return TR
        

def add_afni_prefix(tpattern):
    if ".txt" in tpattern:
        tpattern = "@{0}".format(tpattern)
    return tpattern


def write_to_log(workflow, log_dir, index, inputs, scan_id):
    """
    Method to write into log file the status of the workflow run.
    """

    import os
    import time
    import datetime

    from CPAC import __version__
    from nipype import logging

    iflogger = logging.getLogger('nipype.interface')

    version = __version__
    subject_id = os.path.basename(log_dir)

    if scan_id is None:
        scan_id = "scan_anat"

    strategy = ""
    ts = time.time()
    stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

    try:
        if workflow != 'DONE':
            wf_path = \
                os.path.dirname((os.getcwd()).split(workflow)[1]).strip("/")

            if wf_path and wf_path != "":
                if '/' in wf_path:
                    scan_id, strategy = wf_path.split('/', 1)
                    scan_id = scan_id.strip('_')
                    strategy = strategy.replace("/", "")
                else:
                    scan_id = wf_path.strip('_')

            file_path = os.path.join(log_dir, scan_id, workflow)

            try:
                os.makedirs(file_path)
            except Exception:
                iflogger.info(
                    "filepath already exist, filepath- {0}, "
                    "curr_dir - {1}".format(file_path, os.getcwd()))

        else:
            file_path = os.path.join(log_dir, scan_id)
    except Exception:
        print("ERROR in write log")
        raise

    try:
        os.makedirs(file_path)
    except Exception:
        iflogger.info(
            "filepath already exist, "
            "filepath: {0}, "
            "curr_dir: {1}".format(file_path, os.getcwd())
        )

    out_file = os.path.join(file_path, 'log_{0}.yml'.format(strategy))

    iflogger.info("CPAC custom log:")

    if isinstance(inputs, list):
        inputs = inputs[0]

    if os.path.exists(inputs):
        status_msg = "wf_status: DONE"
        iflogger.info(
            "version: {0}, "
            "timestamp: {1}, "
            "subject_id: {2}, "
            "scan_id: {3}, "
            "strategy: {4}, "
            "workflow: {5}, "
            "status: COMPLETED".format(
                str(version), str(stamp), subject_id,
                scan_id, strategy, workflow
            )
        )
    else:
        status_msg = "wf_status: ERROR"
        iflogger.info(
            "version: {0}, "
            "timestamp: {1}, "
            "subject_id: {2}, "
            "scan_id: {3}, "
            "strategy: {4}, "
            "workflow: {5}, "
            "status: ERROR".format(
                str(version), str(stamp), subject_id,
                scan_id, strategy, workflow
            )
        )

    with open(out_file, 'w') as f:
        f.write("version: {0}\n".format(str(version)))
        f.write("timestamp: {0}\n".format(str(stamp)))
        f.write("pipeline_index: {0}\n".format(index))
        f.write("subject_id: {0}\n".format(subject_id))
        f.write("scan_id: {0}\n".format(scan_id))
        f.write("strategy: {0}\n".format(strategy))
        f.write("workflow_name: {0}\n".format(workflow))
        f.write(status_msg)

    return out_file


def create_log(wf_name="log", scan_id=None):
    """
    Workflow to create log
    """

    from CPAC.pipeline import nipype_pipeline_engine as pe
    import nipype.interfaces.utility as util
    import CPAC.utils.interfaces.function as function

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['workflow',
                                                        'log_dir',
                                                        'index',
                                                        'inputs']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['out_file']),
                          name='outputspec')

    write_log = pe.Node(function.Function(input_names=['workflow',
                                                       'log_dir',
                                                       'index',
                                                       'inputs',
                                                       'scan_id'],
                                          output_names=['out_file'],
                                          function=write_to_log,
                                          as_module=True),
                        name='write_log')

    write_log.inputs.scan_id = scan_id

    wf.connect([
        (
            input_node, write_log, [
                ('workflow', 'workflow'),
                ('log_dir', 'log_dir'),
                ('index', 'index'),
                ('inputs', 'inputs')
            ]
        ),
        (
            write_log, output_node, [
                ('out_file', 'out_file')
            ]
        )
    ])

    return wf


def pick_wm(seg_prob_list):
    seg_prob_list.sort()
    return seg_prob_list[-1]


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


def extract_output_mean(in_file, output_name):
    '''
    function takes 'in_file', which should be an intermediary 1D file
    from individual-level analysis, containing the mean of the output across
    all voxels

    it then parses this value and writes it to a .csv file named
    output_means.csv located in the subject's output directory
    '''

    if os.path.exists(in_file):

        with open(in_file, 'r') as f:
            line = f.readline()

        line = line.split('[')[0].strip(' ')

        # get filename of input maskave 1D file
        filename = in_file.split("/")[-1]
        filename = filename[0:-3]

        split_fullpath = in_file.split("/")

        if "_mask_" in in_file and (
           "sca_roi" in in_file or "sca_tempreg" in in_file
        ):

            for dirname in split_fullpath:
                if "_mask_" in dirname:
                    maskname = dirname

            filename = split_fullpath[-1]

            if ".1D" in filename:
                filename = filename.replace(".1D", "")

            resource_name = "{0}_{1}_{2}".format(output_name, maskname,
                                                 filename)

        elif "_spatial_map_" in in_file and "dr_tempreg" in in_file:

            for dirname in split_fullpath:
                if "_spatial_map_" in dirname:
                    mapname = dirname

            filename = split_fullpath[-1]

            if ".1D" in filename:
                filename = filename.replace(".1D", "")

            resource_name = "{0}_{1}_{2}".format(output_name, mapname,
                                                 filename)

        elif "_mask_" in in_file and "centrality" in in_file:

            for dirname in split_fullpath:
                if "_mask_" in dirname:
                    maskname = dirname

            filename = split_fullpath[-1]

            if ".1D" in filename:
                filename = filename.replace(".1D", "")

            resource_name = "{0}_{1}_{2}".format(output_name, maskname,
                                                 filename)

        else:
            resource_name = output_name

        output_means_file = os.path.join(os.getcwd(),
                                         'mean_{0}.txt'.format(resource_name))

        with open(output_means_file, 'w') as f:
            f.write(line)

    return output_means_file


def create_output_mean_csv(subject_dir):
    '''
    this function finds all of the mean_{output}.txt files in the subject's
    output directory, collects the data and organizes them into one .csv
    file in the subject directory
    '''

    import os
    import csv

    output_vals = {}

    subID = subject_dir.split('/')[len(subject_dir.split('/')) - 1]
    means_dir = os.path.join(subject_dir, 'output_means')

    # extract the mean values
    for root, _, files in os.walk(means_dir):

        for filename in files:

            if 'mean_' in filename:

                output = filename.replace('mean_', '')
                output = output.replace('.txt', '')

                filepath = os.path.join(root, filename)

                if os.path.exists(filepath):
                    try:
                        mean_file = open(filepath, 'rU')
                        val = mean_file.readline()
                        val = val.strip('\n')
                    except:
                        print('\n\n[!] CPAC says: Could not open the output '
                              'mean text file.\n')
                        print('Path: ', filepath, '\n\n')
                        raise Exception

                else:
                    print('\n\n[!] CPAC says: Could not find the output mean '
                          'text file.\n')
                    print('Path not found: ', filepath, '\n\n')
                    raise Exception

                output_vals[output] = val

    # now take the extracted mean values and write them into the .csv file!
    csv_file_path = os.path.join(subject_dir, 'output_means_%s.csv' % subID)
    with open(csv_file_path, 'wt') as csv_file:

        output_items = list(output_vals.items())

        deriv_string = ','.join(v for v, _ in output_items)
        val_string = ','.join(v for _, v in output_items)

        csv_file.write(deriv_string + '\n')
        csv_file.write(val_string + '\n')


def check_command_path(path):
    import os
    return os.system("%s >/dev/null 2>&1" % path) != 32512


def check_system_deps(check_ants=False,
                      check_ica_aroma=False,
                      check_centrality_degree=False,
                      check_centrality_lfcd=False):
    '''
    Function to check system for neuroimaging tools AFNI, C3D, FSL,
    and (optionally) ANTs
    '''

    missing_install = []

    # Check AFNI
    if not check_command_path("3dcalc"):
        missing_install.append("AFNI")

    # Check FSL
    if not check_command_path("fslmaths"):
        missing_install.append("FSL")

    # Check ANTs/C3D
    if check_ants:
        if not check_command_path("c3d_affine_tool"):
            missing_install.append("C3D")
        if not check_command_path("antsRegistration"):
            missing_install.append("ANTS")

    if check_centrality_degree:
        if not check_command_path("3dDegreeCentrality"):
            missing_install.append("3dDegreeCentrality")
    if check_centrality_lfcd:
        if not check_command_path("3dLFCD"):
            missing_install.append("3dLFCD")

    # Check ICA-AROMA
    if check_ica_aroma:
        if not check_command_path("ICA_AROMA.py"):
            missing_install.append("ICA-AROMA")

    # If we're missing deps, raise Exception
    if len(missing_install) > 0:
        missing_string = ""
        for string in missing_install:
            missing_string = missing_string + string + "\n"
        err = "\n\n[!] CPAC says: It appears the following software " \
              "packages are not installed or configured properly:\n\n%s\n" \
              % missing_string
        raise Exception(err)


# Check pipeline config againts computer resources
def check_config_resources(c):
    # Import packages
    import psutil
    from multiprocessing import cpu_count

    # Init variables
    sys_virt_mem = psutil.virtual_memory()
    num_cores = cpu_count()

    # Check for pipeline memory for subject
    if c.pipeline_setup['system_config'][
        'maximum_memory_per_participant'] is None:
        # Get system memory and numSubsAtOnce
        sys_mem_gb = sys_virt_mem.total / (1024.0 ** 3)
        sub_mem_gb = sys_mem_gb / c.pipeline_setup['system_config'][
            'num_participants_at_once']
    else:
        sub_mem_gb = c.pipeline_setup['system_config'][
            'maximum_memory_per_participant']

    # If centrality is enabled, check to mem_sub >= mem_centrality
    if c.network_centrality['run']:
        if sub_mem_gb < c.network_centrality['memory_allocation']:
            err_msg = 'Memory allocated for subject: %d needs to be greater ' \
                      'than the memory allocated for centrality: %d. Fix ' \
                      'and try again.' % (c.pipeline_setup[
                          'system_config']['maximum_memory_per_participant'],
                                          c.network_centrality[
                          'memory_allocation'])
            raise Exception(err_msg)

    # Check for pipeline threads
    # Check if user specified cores
    if c.pipeline_setup['system_config']['max_cores_per_participant']:
        total_user_cores = c.pipeline_setup['system_config'][
                               'num_participants_at_once'] * \
                           c.pipeline_setup['system_config'][
                               'max_cores_per_participant']
        if total_user_cores > num_cores:
            raise SystemError('Configuration specifies more threads running '
                              'in parallel (%d) than number of threads '
                              'available (%d). Change this and try again' %
                              (total_user_cores, num_cores))
        num_cores_per_sub = c.pipeline_setup['system_config'][
            'max_cores_per_participant']
    else:
        num_cores_per_sub = num_cores / c.pipeline_setup['system_config'][
            'num_participants_at_once']

    # Now check ANTS
    if 'ANTS' in c.registration_workflows['anatomical_registration'][
        'registration']['using']:
        if c.pipeline_setup['system_config']['num_ants_threads'] is None:
            num_ants_cores = num_cores_per_sub
        elif c.pipeline_setup['system_config']['num_ants_threads'] > \
                c.pipeline_setup['system_config'][
                    'max_cores_per_participant']:
            err_msg = 'Number of threads for ANTS: %d is greater than the ' \
                      'number of threads per subject: %d. Change this and ' \
                      'try again.' % (
                      c.pipeline_setup['system_config']['num_ants_threads'],
                      c.pipeline_setup['system_config'][
                          'max_cores_per_participant'])
            raise Exception(err_msg)
        else:
            num_ants_cores = c.pipeline_setup['system_config'][
                'num_ants_threads']
    else:
        num_ants_cores = 1

    # Now check OMP
    if c.pipeline_setup['system_config']['num_OMP_threads'] is None:
        num_omp_cores = 1
    elif c.pipeline_setup['system_config']['num_OMP_threads'] > \
            c.pipeline_setup['system_config']['max_cores_per_participant']:
        err_msg = 'Number of threads for OMP: %d is greater than the ' \
                  'number of threads per subject: %d. Change this and ' \
                  'try again.' % (c.pipeline_setup['system_config'][
                                      'num_OMP_threads'],
                                  c.pipeline_setup['system_config'][
                                      'max_cores_per_participant'])
        raise Exception(err_msg)
    else:
        num_omp_cores = c.pipeline_setup['system_config']['num_OMP_threads']

    # Return memory and cores
    return sub_mem_gb, num_cores_per_sub, num_ants_cores, num_omp_cores


def _check_nested_types(d, keys):
    '''Helper function to check types for *_nested_value functions'''
    if not isinstance(d, dict):
        raise TypeError(f'Expected dict, got {type(d).__name__}: {str(d)}')
    if not isinstance(keys, list) and not isinstance(keys, tuple):
        raise TypeError(
            f'Expected list, got {type(keys).__name__}: {str(keys)}')


def delete_nested_value(d, keys):
    '''Helper function to delete nested values

    Parameters
    ---------
    d: dict
    keys: list or tuple

    Returns
    -------
    dict
        updated

    Examples
    --------
    >>> delete_nested_value(
    ...     {'nested': {'key1': 'value', 'key2': 'value'}},
    ...     ['nested', 'key1'])
    {'nested': {'key2': 'value'}}
    '''
    _check_nested_types(d, keys)
    if len(keys) == 1:
        del d[keys[0]]
        return d
    if not len(keys):  # pylint: disable=len-as-condition
        return d
    d[keys[0]] = delete_nested_value(d.get(keys[0], {}), keys[1:])
    return d


def ordereddict_to_dict(value):
    '''
    this function convert ordereddict into regular dict
    '''
    for k, v in value.items():
        if isinstance(v, dict):
            value[k] = ordereddict_to_dict(v)
    return dict(value)


def repickle(directory):
    """
    Function to check all of the pickles in a given directory, recursively, and
    convert any Python 2 pickles found to Python 3 pickles.

    Parameters
    ----------
    directory : str

    Returns
    -------
    None
    """
    for root, _, files in os.walk(directory, followlinks=True):
        for fn in files:
            p = os.path.join(root, fn)
            if fn.endswith(".pkl"):
                if _pickle2(p):
                    try:
                        with open(p, 'rb') as fp:
                            f = pickle.load(fp, encoding='latin1')
                        with open(p, 'wb') as fp:
                            pickle.dump(f, fp)
                        print(
                            f"Converted pickle {fn} from a Python 2 pickle to "
                            "a Python 3 pickle."
                        )
                    except Exception as e:
                        print(
                            f"Could not convert Python 2 pickle {p} "
                            f"because {e}\n"
                        )
                else:
                    print(f"Pickle {fn} is a Python 3 pickle.")
            elif fn.endswith(".pklz"):
                if _pickle2(p, True):
                    try:
                        with gzip.open(p, 'rb') as fp:
                            f = pickle.load(fp, encoding='latin1')
                        with gzip.open(p, 'wb') as fp:
                            pickle.dump(f, fp)
                        print(
                            f"Converted pickle {fn} from a Python 2 pickle to "
                            "a Python 3 pickle."
                        )
                    except Exception as e:
                        print(
                            f"Could not convert Python 2 pickle {p} "
                            f"because {e}\n"
                        )
                else:
                    print(f"Pickle {fn} is a Python 3 pickle.")


def _pickle2(p, z=False):
    """
    Helper function to check if a pickle is a Python 2 pickle. Also prints
    other exceptions raised by trying to load the file at p.

    Parameters
    ----------
    p : str
        path to pickle

    z : bool
        if pickle is gzipped

    Returns
    -------
    pickle2 : bool
        True if p is a Python 2 pickle
    """
    if z:
        with gzip.open(p, 'rb') as fp:
            try:
                pickle.load(fp)
            except UnicodeDecodeError:
                return True
            except Exception as e:
                print(
                    f"Pickle {p} may be a Python 3 pickle, but raised "
                    f"exception {e}"
                )
    else:
        with open(p, 'rb') as fp:
            try:
                pickle.load(fp)
            except UnicodeDecodeError:
                return True
            except Exception as e:
                print(
                    f"Pickle {p} may be a Python 3 pickle, but raised "
                    f"exception {e}"
                )
    return False


def _changes_1_8_0_to_1_8_1(config_dict):
    '''
    Examples
    --------
    Starting with 1.8.0
    >>> zero = {'anatomical_preproc': {
    ...     'non_local_means_filtering': True,
    ...     'n4_bias_field_correction': True
    ... }, 'functional_preproc': {
    ...     'motion_estimates_and_correction': {
    ...         'calculate_motion_first': False
    ...     }
    ... }, 'segmentation': {
    ...     'tissue_segmentation': {
    ...         'ANTs_Prior_Based': {
    ...             'CSF_label': 0,
    ...             'left_GM_label': 1,
    ...             'right_GM_label': 2,
    ...             'left_WM_label': 3,
    ...             'right_WM_label': 4}}}}
    >>> updated_apb = _changes_1_8_0_to_1_8_1(zero)[
    ...     'segmentation']['tissue_segmentation']['ANTs_Prior_Based']
    >>> updated_apb['CSF_label']
    [0]
    >>> updated_apb['GM_label']
    [1, 2]
    >>> updated_apb['WM_label']
    [3, 4]

    Starting with 1.8.1
    >>> one = {'anatomical_preproc': {
    ...     'non_local_means_filtering': True,
    ...     'n4_bias_field_correction': True
    ... }, 'functional_preproc': {
    ...     'motion_estimates_and_correction': {
    ...         'calculate_motion_first': False
    ...     }
    ... }, 'segmentation': {
    ...     'tissue_segmentation': {
    ...         'ANTs_Prior_Based': {
    ...             'CSF_label': [0],
    ...             'GM_label': [1, 2],
    ...             'WM_label': [3, 4]}}}}
    >>> updated_apb = _changes_1_8_0_to_1_8_1(one)[
    ...     'segmentation']['tissue_segmentation']['ANTs_Prior_Based']
    >>> updated_apb['CSF_label']
    [0]
    >>> updated_apb['GM_label']
    [1, 2]
    >>> updated_apb['WM_label']
    [3, 4]
    '''
    for key_sequence in {
        ('anatomical_preproc', 'non_local_means_filtering'),
        ('anatomical_preproc', 'n4_bias_field_correction')
    }:
        config_dict = _now_runswitch(config_dict, key_sequence)
    for combiners in {
        ((
            ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
             'CSF_label'),
        ), ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
            'CSF_label')),
        ((
            ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
             'left_GM_label'),
            ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
             'right_GM_label')
        ), ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
            'GM_label')),
        ((
            ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
             'left_WM_label'),
            ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
             'right_WM_label')
        ), ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
            'WM_label'))
    }:
        config_dict = _combine_labels(config_dict, *combiners)
    try:
        calculate_motion_first = lookup_nested_value(
            config_dict,
            ['functional_preproc', 'motion_estimates_and_correction',
                'calculate_motion_first']
        )
    except KeyError:
        calculate_motion_first = None
    if calculate_motion_first is not None:
        del config_dict['functional_preproc'][
            'motion_estimates_and_correction']['calculate_motion_first']
        config_dict = set_nested_value(config_dict, [
            'functional_preproc', 'motion_estimates_and_correction',
            'motion_estimates', 'calculate_motion_first'
        ], calculate_motion_first)

    return config_dict


def _combine_labels(config_dict, list_to_combine, new_key):
    '''
    Helper function to combine formerly separate keys into a
    combined key.

    Parameters
    ----------
    config_dict: dict

    key_sequence: iterable of lists or tuples

    new_key: list or tuple

    Returns
    -------
    updated_config_dict: dict
    '''
    new_value = []
    any_old_values = False
    for _to_combine in list_to_combine:
        try:
            old_value = lookup_nested_value(config_dict, _to_combine)
        except KeyError:
            old_value = None
        if old_value is not None:
            any_old_values = True
            if isinstance(old_value, (list, set, tuple)):
                for value in old_value:
                    new_value.append(value)
            else:
                new_value.append(old_value)
            config_dict = delete_nested_value(config_dict, _to_combine)
    if any_old_values:
        return set_nested_value(config_dict, new_key, new_value)
    return config_dict


def concat_list(in_list1=None, in_list2=None):
    """
    Parameters
    ----------
    in_list1 : list or str
        file path or a list of file paths

    in_list2 : list or str
        file path or a list of file paths

    Returns
    -------
    out_list : list
        a list of file paths
    """

    if in_list1 is not None:
        if not isinstance(in_list1, list):
            in_list1 = [in_list1]
    else:
        in_list1 = []

    if in_list2 is not None:
        if not isinstance(in_list2, list):
            in_list2 = [in_list2]
    else:
        in_list2 = []

    out_list = in_list1 + in_list2

    return out_list


def list_item_replace(l,  # noqa: E741  # pylint: disable=invalid-name
                      old, new):
    '''Function to replace an item in a list

    Parameters
    ----------
    l : list or string

    old : any
        item to replace

    new : any
        new item

    Returns
    -------
    l : list or string
        updated

    Examples
    --------
    >>> list_item_replace(['AFNI', 'FSL'], 'AFNI', '3dSkullStrip')
    ['3dSkullStrip', 'FSL']
    >>> list_item_replace(['AFNI', 'FSL'], 'FSL', 'BET')
    ['AFNI', 'BET']
    '''
    if isinstance(l, list) and old in l:
        l[l.index(old)] = new
    elif isinstance(l, str):
        l = l.replace(old, new)  # noqa: E741
    return l


def lookup_nested_value(d, keys):
    '''Helper method to look up nested values

    Parameters
    ---------
    d: dict
    keys: list or tuple

    Returns
    -------
    yaml: str or dict

    Examples
    --------
    >>> lookup_nested_value({'nested': {'True': True}}, ['nested', 'True'])
    True
    >>> lookup_nested_value({'nested': {'None': None}}, ['nested', 'None'])
    ''
    '''
    if not isinstance(d, dict):
        return d
    if len(keys) == 1:
        value = d[keys[0]]
        if value is None:
            return ''
        return value
    else:
        try:
            return lookup_nested_value(d[keys[0]], keys[1:])
        except KeyError as e:
            e.args = (keys,)
            raise


def _now_runswitch(config_dict, key_sequence):
    '''
    Helper function to convert a formerly forkable value to a
    runswitch.

    Parameters
    ----------
    config_dict: dict

    key_sequence: list or tuple

    Returns
    -------
    updated_config_dict: dict
    '''
    try:
        old_forkable = lookup_nested_value(config_dict, key_sequence)
    except KeyError:
        return config_dict
    if isinstance(old_forkable, (bool, list)):
        return set_nested_value(
            config_dict, key_sequence, {'run': old_forkable})
    return config_dict


def _remove_somethings(value, things_to_remove):
    '''Helper function to remove instances of any in a given set of
    values from a list.

    Parameters
    ----------
    value : list

    things_to_remove : set

    Returns
    -------
    list
    '''
    if isinstance(value, list):
        for thing in things_to_remove:
            while thing in value:
                value.remove(thing)
    return value


def remove_False(d, k):
    '''Function to remove "Off" and False from a list at a given nested key.

    Parameters
    ----------
    d : dict

    k : list

    Returns
    -------
    d : dict
       updated

    Examples
    --------
    >>> remove_False({'a': {'b': [1, False, 2, "Off", 3]}}, ['a', 'b'])
    {'a': {'b': [1, 2, 3]}}
    '''
    value = _remove_somethings(lookup_nested_value(d, k), {False, 'Off'})
    return set_nested_value(d, k, value)


def remove_None(d, k):
    '''Function to remove "None" and None from a list at a given nested key.

    Parameters
    ----------
    d : dict

    k : list

    Returns
    -------
    d : dict
       updated

    Examples
    --------
    >>> remove_None({'a': {'b': [1, None, 2, "None", 3]}}, ['a', 'b'])
    {'a': {'b': [1, 2, 3]}}
    '''
    value = _remove_somethings(lookup_nested_value(d, k), {None, 'None'})
    return set_nested_value(d, k, value)


def replace_in_strings(d, replacements=None):
    '''Helper function to recursively replace substrings.

    Parameters
    ----------
    d : any

    replacements : list of 2-tuples
        0 : str
            substring to replace
        1 : str
            replacement substring

    Returns
    -------
    d : any
       same as input, but updated

    Examples
    --------
    >>> replace_in_strings({'key': 'test${resolution_for_func_preproc}'})
    {'key': 'test${func_resolution}'}
    '''
    if replacements is None:
        replacements = [(r'${resolution_for_func_preproc}',
                         r'${func_resolution}')]
    if isinstance(d, dict):
        return {k: replace_in_strings(d[k], replacements) for k in d}
    if isinstance(d, list):
        return [replace_in_strings(i, replacements) for i in d]
    if isinstance(d, str):
        for replacement in replacements:
            d = d.replace(replacement[0], replacement[1])
    return d


def set_nested_value(d, keys, value):
    '''Helper method to set nested values

    Parameters
    ---------
    d: dict
    keys: list or tuple
    value: any

    Returns
    -------
    dict
        updated

    Examples
    --------
    >>> set_nested_value({}, ['nested', 'keys'], 'value')
    {'nested': {'keys': 'value'}}
    '''
    _check_nested_types(d, keys)
    if len(keys) == 1:
        d.update({keys[0]: value})
        return d
    if not len(keys):  # pylint: disable=len-as-condition
        return d
    new_d = {
        keys[0]: set_nested_value(d.get(keys[0], {}), keys[1:], value)
    }
    d = update_nested_dict(d, new_d)
    return d


def update_config_dict(old_dict):
    '''Function to convert an old config dict to a new config dict

    Parameters
    ----------
    old_dict : dict

    Returns
    -------
    new_dict : dict
        1.8 nested config dictionary

    old_dict : dict
        remaining undefined mappings

    combined_dict : dict
        1.8 nested config dictionary plus remaining undefined mappings

    Examples
    --------
    >>> a, b, c = update_config_dict({
    ...     'pipelineName': 'example-pipeline', '2': None})
    >>> a
    {'pipeline_setup': {'pipeline_name': 'example-pipeline'}}
    >>> b
    {'2': None}
    >>> c
    {'pipeline_setup': {'pipeline_name': 'example-pipeline'}, '2': None}
    '''
    def _append_to_list(current_value, new_value):
        '''Helper function to add new_value to the current_value list
        or create a list if one does not exist. Skips falsy elements
        in new_value

        Parameters
        ----------
        current_value : list

        new_value : list, bool, None, or str

        Returns
        -------
        list

        Examples
        --------
        >>> _append_to_list([1], [2])
        [1, 2]
        >>> _append_to_list([1, 2], [2])
        [1, 2]
        >>> _append_to_list(None, [2])
        [2]
        >>> _append_to_list([1], [1, 2])
        [1, 2]
        >>> _append_to_list([1], [None, 2])
        [1, 2]
        '''
        if not isinstance(current_value, list):
            if current_value is not None:
                current_value = [current_value]
            else:
                current_value = []
        else:
            current_value = [v for v in current_value if v is not None]
        if isinstance(new_value, list):
            for i in new_value:
                if i and i not in current_value and i != 'Off':
                    current_value.append(i)
        elif (
            new_value and new_value not in current_value and
            new_value != 'Off'
        ):
            current_value.append(new_value)
        return current_value

    def _bool_to_str(old_value, value_if_true):
        '''Helper function to convert a True or a list containing a
        True to a given string

        Parameters
        ----------
        old_value : list, bool, None, or str

        value_if_true : str

        Returns
        -------
        str or None

        Examples
        --------
        >>> _bool_to_str([0], 'test_str')
        >>> _bool_to_str([1], 'test_str')
        'test_str'
        >>> _bool_to_str(0, 'test_str')
        >>> _bool_to_str(1, 'test_str')
        'test_str'
        >>> _bool_to_str([True, False], 'test_str')
        'test_str'
        >>> _bool_to_str(None, 'test_str')
        >>> _bool_to_str([0, None, False], 'test_str')
        >>> _bool_to_str([0, None, False, 1], 'test_str')
        'test_str'
        '''
        if isinstance(old_value, list):
            if any(bool(i) for i in old_value):
                return value_if_true
        elif bool(old_value):
            return value_if_true
        return None

    def _get_old_values(old_dict, new_dict, key):
        '''Helper function to get old and current values of a special key
        being updated.

        Parameters
        ----------
        old_dict : dict

        new_dict : dict

        key : str

        Returns
        -------
        old_dict : dict

        new_dict : dict

        old_value : any

        current_value : any
        '''
        old_value = old_dict.pop(key)
        current_value = lookup_nested_value(
            new_dict, NESTED_CONFIG_MAPPING[key]
        )
        return old_dict, new_dict, old_value, current_value

    new_dict = {}
    for key in old_dict.copy():
        if key in NESTED_CONFIG_MAPPING:
            # handle special cases
            special_cases = {
                'acpc_run_preprocessing',
                'acpc_template_brain',
                'ANTs_prior_based_segmentation',
                'func_reg_input',
                'runRegisterFuncToTemplate',
                'runRegisterFuncToEPI',
                'fsl_linear_reg_only',
                'functional_registration',
                'template_for_resample',
                'fnirtConfig',
                'run_smoothing',
                'runZScoring',
                'run_longitudinal'
            }
            if key in special_cases:
                try:
                    (
                        old_dict, new_dict, old_value, current_value
                    ) = _get_old_values(old_dict, new_dict, key)
                except KeyError:
                    continue

                # longitudinal_template_generation.run
                if key == 'run_longitudinal':
                    if 'anat' in old_value or 'func' in old_value:
                        current_value = True
                    else:
                        current_value = False

                # anatomical_preproc.acpc_alignment.run_before_preproc
                if key == 'acpc_run_preprocessing':
                    current_value = True if old_value.lower(
                    ) == 'before' else False if old_value.lower(
                    ) == 'after' else None

                # anatomical_preproc.acpc_alignment.acpc_target
                if key == 'acpc_template_brain':
                    if current_value in {'None', None, ''}:
                        new_dict = set_nested_value(
                            new_dict,
                            ['anatomical_preproc', 'acpc_alignment',
                             'acpc_target'],
                            'whole-head'
                        )

                # segmentation.tissue_segmentation.using
                elif key == 'ANTs_prior_based_segmentation':
                    new_value = _bool_to_str(old_value, 'ANTs_Prior_Based')
                    if new_value == 'ANTs_Prior_Based':
                        new_dict = set_nested_value(
                            new_dict,
                            NESTED_CONFIG_MAPPING[key][:-1] +
                            [new_value, 'run'],
                            old_value
                        )

                # registration_workflows.functional_registration.
                # coregistration.func_input_prep.input
                elif key == 'func_reg_input':
                    new_value = _replace_in_value_list(old_value, (' ', '_'))
                    current_value = _replace_in_value_list(
                        current_value, (' ', '_'))

                # registration_workflows.functional_registration.
                # func_registration_to_template.target_template.using
                elif key in {
                    'runRegisterFuncToTemplate', 'runRegisterFuncToEPI'
                }:
                    current_value = _replace_in_value_list(
                        current_value, (' ', '_'))
                    if key == 'runRegisterFuncToTemplate':
                        current_value = [
                            v for v in current_value if v not in {
                                'Off', 'False', False
                            }
                        ]
                        new_value = []
                        new_dict = set_nested_value(
                            new_dict,
                            ['registration_workflows',
                             'functional_registration',
                             'func_registration_to_template', 'run'],
                            bool(current_value)
                        )
                    if key == 'runRegisterFuncToEPI':
                        new_value = _bool_to_str(old_value, 'EPI_template')

                # registration_workflows.anatomical_registration.registration.
                # using
                elif key == 'fsl_linear_reg_only':
                    new_value = _bool_to_str(old_value, 'FSL-linear')

                # registration_workflows.functional_registration.
                # func_registration_to_template.target_template.
                # EPI_template.EPI_template_for_resample
                elif key == 'template_for_resample':
                    new_dict = set_nested_value(
                        new_dict,
                        ['registration_workflows', 'functional_registration',
                         'func_registration_to_template', 'target_template',
                         'EPI_template', 'EPI_template_for_resample'],
                        current_value
                    )

                # registration_workflows.functional_registration.
                # EPI_registration.FSL-FNIRT.fnirt_config
                elif key == 'fnirtConfig':
                    current_value = old_value
                    new_dict = set_nested_value(
                        new_dict,
                        ['registration_workflows', 'functional_registration',
                         'EPI_registration', 'FSL-FNIRT', 'fnirt_config'],
                        current_value
                    )

                # post_processing.spatial_smoothing.output
                elif key == 'run_smoothing':
                    new_value = [_bool_to_str(old_value, 'smoothed')]
                    if any(not bool(value) for value in old_value):
                        new_value.append('nonsmoothed')
                    current_value = new_value

                # post_processing.z-scoring.output
                elif key == 'runZScoring':
                    new_value = [_bool_to_str(old_value, 'z-scored')]
                    if any(not bool(value) for value in old_value):
                        new_value.append('raw')
                    current_value = new_value

                # make sure list values are cast as lists
                if key not in {  # if key not in non-list-valued keys
                    'acpc_run_preprocessing', 'acpc_template_brain',
                    'functional_registration', 'template_for_resample',
                    'fnirtConfig', 'run_longitudinal'
                }:
                    current_value = _append_to_list(current_value, new_value)

            # update remaining keys
            else:
                current_value = old_dict.pop(key)

            if current_value == 'None':
                current_value = None

            new_dict = set_nested_value(
                new_dict, NESTED_CONFIG_MAPPING[key], current_value)
        elif key in NESTED_CONFIG_DEPRECATIONS:
            old_dict.pop(key)
    return new_dict, old_dict, update_nested_dict(new_dict.copy(), old_dict)


def update_nested_dict(d_base, d_update, fully_specified=False):
    """Update dictionary of varying depth.

    Parameters
    ----------
    d_base : dict
        original dictionary

    d_update : dict
        dictionary with updates

    fully_specified : bool
        if True, overwrite instead of update

    Returns
    -------
    d_base : dict
        original dictionary with updates

    Examples
    --------
    >>> d_base = {'pipeline_name': 'cpac-default-pipeline',
    ...     'output_directory': {'path': '/output',
    ...     'write_func_outputs': False,
    ...     'write_debugging_outputs': False,
    ...     'output_tree': 'default',
    ...     'quality_control': {
    ...         'generate_quality_control_images': True,
    ...         'generate_xcpqc_files': True}},
    ...     'working_directory': {'path': '/tmp', 'remove_working_dir': True},
    ...     'log_directory': {'run_logging': True, 'path': '/logs'},
    ...     'system_config': {'maximum_memory_per_participant': 1,
    ...     'max_cores_per_participant': 1,
    ...     'num_ants_threads': 4,
    ...     'num_participants_at_once': 1},
    ...     'Amazon-AWS': {'aws_output_bucket_credentials': None,
    ...                    's3_encryption': False}}
    >>> d_update = {'pipeline_name': 'cpac_fmriprep-options',
    ...     'system_config': {'num_ants_threads': 1},
    ...     'Amazon-AWS': {'s3_encryption': True}}
    >>> str(update_nested_dict(d_base, d_update)) == str({
    ...     'pipeline_name': 'cpac_fmriprep-options', 'output_directory': {
    ...         'path': '/output', 'write_func_outputs': False,
    ...         'write_debugging_outputs': False, 'output_tree': 'default',
    ...         'quality_control': {
    ...             'generate_quality_control_images': True,
    ...             'generate_xcpqc_files': True}
    ...     }, 'working_directory': {
    ...        'path': '/tmp', 'remove_working_dir': True
    ...     }, 'log_directory': {'run_logging': True, 'path': '/logs'},
    ...     'system_config': {
    ...     'maximum_memory_per_participant': 1,
    ...     'max_cores_per_participant': 1,
    ...     'num_ants_threads': 1, 'num_participants_at_once': 1
    ... }, 'Amazon-AWS': {
    ...     'aws_output_bucket_credentials': None, 's3_encryption': True}})
    True
    >>> tse_base = {'timeseries_extraction': {'run': True, 'tse_roi_paths': {
    ...     '/cpac_templates/CC400.nii.gz': 'Avg',
    ...     '/cpac_templates/aal_mask_pad.nii.gz': 'Avg'
    ... }, 'realignment': 'ROI_to_func'}}
    >>> str(update_nested_dict(tse_base, {})) == str({
    ...     'timeseries_extraction': {'run': True, 'tse_roi_paths': {
    ...         '/cpac_templates/CC400.nii.gz': 'Avg',
    ...         '/cpac_templates/aal_mask_pad.nii.gz': 'Avg'
    ... }, 'realignment': 'ROI_to_func'}})
    True
    >>> str(update_nested_dict(tse_base, {'timeseries_extraction': {
    ...     'tse_roi_paths': {'/cpac_templates/rois_3mm.nii.gz': 'Voxel'}
    ... }})) == str({'timeseries_extraction': {'run': True, 'tse_roi_paths': {
    ...     '/cpac_templates/rois_3mm.nii.gz': 'Voxel'
    ... }, 'realignment': 'ROI_to_func'}})
    True
    >>> str(update_nested_dict(tse_base, {'timeseries_extraction': {
    ...     'roi_paths_fully_specified': False,
    ...     'tse_roi_paths': {'/cpac_templates/rois_3mm.nii.gz': 'Voxel'}
    ... }})) == str({'timeseries_extraction': {'run': True, 'tse_roi_paths': {
    ...     '/cpac_templates/CC400.nii.gz': 'Avg',
    ...     '/cpac_templates/aal_mask_pad.nii.gz': 'Avg',
    ...     '/cpac_templates/rois_3mm.nii.gz': 'Voxel'
    ... }, 'realignment': 'ROI_to_func'}})
    True
    >>> str(update_nested_dict(tse_base, {'timeseries_extraction': {
    ...     'roi_paths_fully_specified': False,
    ...     'tse_roi_paths': {'/cpac_templates/aal_mask_pad.nii.gz': 'Voxel'}
    ... }})) == str({'timeseries_extraction': {'run': True,
    ...     'tse_roi_paths': {
    ...         '/cpac_templates/CC400.nii.gz': 'Avg',
    ...         '/cpac_templates/aal_mask_pad.nii.gz': 'Voxel'
    ... }, 'realignment': 'ROI_to_func'}})
    True
    >>> str(update_nested_dict(tse_base, {'timeseries_extraction': {
    ...     'tse_roi_paths': {'/cpac_templates/aal_mask_pad.nii.gz': 'Voxel'}
    ... }})) == str({'timeseries_extraction': {'run': True, 'tse_roi_paths': {
    ...     '/cpac_templates/aal_mask_pad.nii.gz': 'Voxel'
    ... }, 'realignment': 'ROI_to_func'}})
    True
    """  # noqa: E501  # pylint: disable=line-too-long

    # short-circuit if d_update has `*_roi_paths` and
    # `roi_paths_fully_specified` children
    if fully_specified:
        return d_update
    d_new = {} if d_base is None else deepcopy(d_base)
    for k, v in d_update.items():
        if k.endswith('_roi_paths'):
            fully_specified = d_update.get('roi_paths_fully_specified', True)
        else:
            fully_specified = False
        if k != 'roi_paths_fully_specified':
            if isinstance(v, collections.abc.Mapping):
                d_new[k] = update_nested_dict(d_new.get(k, {}), v,
                                              fully_specified)
            else:
                d_new[k] = v
    return d_new


def update_pipeline_values_1_8(d_old):
    '''Function to update pipeline config values that changed from
    C-PAC 1.7 to 1.8.

    Parameters
    ----------
    d_old : dict

    Returns
    -------
    d : dict
       updated

    Examples
    --------
    >>> update_pipeline_values_1_8({'segmentation': {'tissue_segmentation': {
    ...     'using': ['FSL-FAST Thresholding', 'Customized Thresholding']}}})
    {'segmentation': {'tissue_segmentation': {'using': ['FSL-FAST'], 'FSL-FAST': {'thresholding': {'use': 'Custom'}}}}}
    >>> update_pipeline_values_1_8({'segmentation': {'tissue_segmentation': {
    ...     'using': ['FSL-FAST Thresholding']}}})
    {'segmentation': {'tissue_segmentation': {'using': ['FSL-FAST'], 'FSL-FAST': {'thresholding': {'use': 'Auto'}}}}}
    '''  # noqa: E501  # pylint: disable=line-too-long
    from CPAC.pipeline.schema import valid_options \
        # pylint: disable=import-outside-toplevel

    d = replace_in_strings(d_old.copy())

    d = _replace_changed_values(
        d,
        ['anatomical_preproc', 'brain_extraction', 'using'],
        [('AFNI', '3dSkullStrip'), ('FSL', 'BET'), ('unet', 'UNet')]
    )

    d = _replace_changed_values(
        d,
        ['functional_preproc', 'func_masking', 'using'],
        [('3dAutoMask', 'AFNI'), ('BET', 'FSL')]
    )

    try:
        seg_use_threshold = lookup_nested_value(d, [
            'segmentation', 'tissue_segmentation', 'using'])
    except KeyError:
        seg_use_threshold = []

    if not isinstance(seg_use_threshold, list):
        seg_use_threshold = [seg_use_threshold]
    if 'FSL-FAST Thresholding' in seg_use_threshold:
        if 'using' in d['segmentation'].get(
            'tissue_segmentation', {}
        ):
            d['segmentation'][
                'tissue_segmentation'
            ]['using'].append('FSL-FAST')
        else:
            d = set_nested_value(d, [
                'segmentation', 'tissue_segmentation',
                'using'], ['FSL-FAST'])
        seg_use_threshold.remove('FSL-FAST Thresholding')
    if 'Customized Thresholding' in seg_use_threshold:
        seg_use_threshold.remove('Customized Thresholding')
        d = set_nested_value(d, [
            'segmentation', 'tissue_segmentation',
            'FSL-FAST', 'thresholding', 'use'], 'Custom')
    else:
        d = set_nested_value(d, [
            'segmentation', 'tissue_segmentation',
            'FSL-FAST', 'thresholding', 'use'], 'Auto')

    for centr in ['degree_centrality', 'eigenvector_centrality',
                  'local_functional_connectivity_density']:
        centr_keys = ['network_centrality', centr, 'weight_options']
        try:
            centr_value = lookup_nested_value(d, centr_keys)
            if any(isinstance(v, bool) for v in centr_value):
                for i in range(2):
                    if centr_value[i] is True:
                        centr_value[i] = valid_options['centrality'][
                            'weight_options'][i]
                while False in centr_value:
                    centr_value.remove(False)
                d = set_nested_value(d, centr_keys, centr_value)
        except KeyError:
            continue

    seg_template_key = [
        'segmentation', 'tissue_segmentation',
        'Template_Based', 'template_for_segmentation']
    try:
        seg_template = lookup_nested_value(d, seg_template_key)
        for replacement in [
            ('EPI_template', valid_options['segmentation']['template'][0]),
            ('T1_template', valid_options['segmentation']['template'][1])
        ]:
            seg_template = list_item_replace(seg_template, *replacement)
        while 'Off' in seg_template:
            seg_template.remove('Off')
        while False in seg_template:
            seg_template.remove(False)
        d = set_nested_value(d, seg_template_key, seg_template)
        d = remove_None(d, seg_template_key)
    except KeyError:
        pass

    distcor_key = ['functional_preproc', 'distortion_correction', 'using']
    try:
        lookup_nested_value(d, distcor_key)
        d = remove_None(d, distcor_key)
    except KeyError:
        pass

    if 'functional_registration' in d and isinstance(
        d['functional_registration'], dict
    ):
        if '1-coregistration' in d['functional_registration']:
            coreg = d['functional_registration'].pop('1-coregistration')
            d = set_nested_value(
                d, ['registration_workflows', 'functional_registration',
                    'coregistration'],
                coreg
            )
        if not(bool(d['functional_registration'])):
            d.pop('functional_registration')

    return update_values_from_list(d)


def update_values_from_list(d_old, last_exception=None):
    '''Function to convert 1-length lists of an expected type to
    single items of that type, or to convert singletons of an expected
    list of a type into lists thereof. Also handles some type
    conversions against the schema.

    Parameters
    ----------
    d_old : dict

    last_exception: Exception or None
        if the same exception recurs, raise it.

    Returns
    -------
    d : dict
       updated

    Examples
    --------
    >>> update_values_from_list({'pipeline_setup': {
    ...     'pipeline_name': ['one_string']}})
    {'pipeline_setup': {'pipeline_name': 'one_string'}}
    >>> update_values_from_list({'nuisance_corrections': {
    ...     '1-ICA-AROMA': {'run': [False]}}})
    {'nuisance_corrections': {'1-ICA-AROMA': {'run': [False]}}}
    '''
    from CPAC.pipeline.schema import schema

    d = d_old.copy()

    try:
        schema(d)
    except Invalid as e:
        if (
            last_exception and last_exception.path == e.path and
            last_exception.msg == e.msg
        ):
            raise e
        observed = lookup_nested_value(d, e.path)
        if observed == 'None':
            return update_values_from_list(
                set_nested_value(d, e.path, None), e)

        expected = e.msg.split('expected')[-1].strip(
        ) if 'expected' in e.msg else 'unknown'

        if (
            expected != 'bool' and isinstance(observed, list) and
            len(observed) == 1
        ):
            try:
                return update_values_from_list(
                    set_nested_value(d, e.path, observed[0]), e)
            except TypeError:
                raise e

        if expected == 'bool':
            if isinstance(observed, int):  # pylint: disable=no-else-return
                return update_values_from_list(
                    set_nested_value(d, e.path, bool(observed)), e)
            elif isinstance(observed, list):
                if len(observed) == 0:  # pylint: disable=no-else-return
                    return update_values_from_list(set_nested_value(
                        d, e.path, False), e)
                else:
                    # maintain a list if list expected
                    list_expected = (e.path[-1] == 0)
                    e_path = e.path[:-1] if list_expected else e.path
                    if len(observed) == 1:  # pylint: disable=no-else-return
                        if isinstance(observed[0], int):
                            value = bool(observed[0])
                        elif observed[0].lower() in YAML_BOOLS[True]:
                            value = True
                        elif observed[0].lower() in YAML_BOOLS[False]:
                            value = False
                        return update_values_from_list(set_nested_value(
                            d, e_path, [value] if list_expected else value), e)
                    else:
                        return update_values_from_list(set_nested_value(
                            d, e_path, [bool(value) for value in observed]), e)
            elif observed.lower() in YAML_BOOLS[True]:
                return update_values_from_list(
                    set_nested_value(d, e.path, True), e)
            elif observed.lower() in YAML_BOOLS[False]:
                return update_values_from_list(
                    set_nested_value(d, e.path, False), e)
            else:
                return update_values_from_list(
                    set_nested_value(d, e_path, observed[0]), e)

        elif expected == 'a list':
            return update_values_from_list(
                set_nested_value(d, e.path, [observed]), e)
        else:
            raise e
    return d


def _replace_changed_values(d, nested_key, replacement_list):
    '''Helper function to replace values changed from C-PAC 1.7 to C-PAC 1.8.

    Parameters
    ----------
    d : dict

    nested_key : list of strings

    replacement_list : list of tuples
        0 : any
            value to replace
        1 : any
            replacement value

    Returns
    -------
    d : dict

    Examples
    --------
    >>> d = {'test': {'this': ['function']}}
    >>> _replace_changed_values(d, ['test', 'this'], [('function', 'success')])
    {'test': {'this': ['success']}}
    '''
    try:
        current_value = lookup_nested_value(d, nested_key)
    except KeyError:
        return d
    if isinstance(current_value, list):
        current_value = _replace_in_value_list(
            current_value, replacement_list)
    else:
        for replacement in replacement_list:
            current_value = list_item_replace(current_value, *replacement)
    return set_nested_value(d, nested_key, current_value)


def _replace_in_value_list(current_value, replacement_tuple):
    '''Helper function to make character replacements in
    `current_value` and drop falsy values.

    Parameters
    ----------
    current_value : list

    replacement_tuple : tuple or list of tuples
        0 : str
            character value to replace
        1 : str
            replacement character value

    Returns
    -------
    current_value : list

    Examples
    --------
    >>> current_value = ['EPI Template', 'T1_Template', 'None']
    >>> _replace_in_value_list(current_value, (' ', '_'))
    ['EPI_Template', 'T1_Template']
    >>> current_value = ['AFNI', 'FSL', 'None', None, False]
    >>> _replace_in_value_list(current_value, [
    ...     ('AFNI', '3dSkullStrip'), ('FSL', 'BET')])
    ['3dSkullStrip', 'BET']
    '''
    if isinstance(replacement_tuple, list):
        for rt in replacement_tuple:
            current_value = _replace_in_value_list(current_value, rt)
        return current_value
    if not isinstance(current_value, list):
        current_value = [current_value]
    return [
        v.replace(*replacement_tuple) for v in current_value
        if bool(v) and v not in {'None', 'Off', ''}
    ]
