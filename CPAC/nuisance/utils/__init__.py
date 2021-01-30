import os
import re
from collections import OrderedDict

import nipype.interfaces.ants as ants
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe
from nipype.interfaces import afni

from CPAC.nuisance.utils.compcor import calc_compcor_components
from CPAC.nuisance.utils.crc import encode as crc_encode
from CPAC.utils.interfaces.function import Function
from CPAC.registration.utils import check_transforms, generate_inverse_transform_flags


def find_offending_time_points(fd_j_file_path=None, fd_p_file_path=None, dvars_file_path=None,
                               fd_j_threshold=None, fd_p_threshold=None, dvars_threshold=None,
                               number_of_previous_trs_to_censor=0,
                               number_of_subsequent_trs_to_censor=0):
    """
    Applies criterion in method to find time points whose FD or DVARS (or both)
    are above threshold.

    :param fd_j_file_path: path to TSV containing framewise displacement as a
        single column. If not specified, it will not be used.
    :param fd_p_file_path: path to TSV containing framewise displacement as a
        single column. If not specified, it will not be used.
    :param dvars_file_path: path to TSV containing DVARS as a single column.
        If not specified, it will not be used.
    :param fd_j_threshold: threshold to apply to framewise displacement (Jenkinson),
        it can be a value such as 0.2 or a floating point multiple of the
        standard deviation specified as, e.g. '1.5SD'.
    :param fd_p_threshold: threshold to apply to framewise displacement (Power),
        it can be a value such as 0.2 or a floating point multiple of the
        standard deviation specified as, e.g. '1.5SD'.
    :param dvars_threshold: threshold to apply to DVARS, can be a value such
        as 0.5 or a floating point multiple of the standard deviation specified
        as, e.g. '1.5SD'.
    :param number_of_previous_trs_to_censor: extent of censorship window before
        the censor.
    :param number_of_subsequent_trs_to_censor: extent of censorship window after
        the censor.

    :return: File path to TSV file containing the volumes to be censored.
    """
    import numpy as np
    import os
    import re

    offending_time_points = set()
    time_course_len = 0

    types = ['FDJ', 'FDP', 'DVARS']
    file_paths = [fd_j_file_path, fd_p_file_path, dvars_file_path]
    thresholds = [fd_j_threshold, fd_p_threshold, dvars_threshold]
    #types = ['FDP', 'DVARS']
    #file_paths = [fd_p_file_path, dvars_file_path]
    #thresholds = [fd_p_threshold, dvars_threshold]

    for type, file_path, threshold in zip(types, file_paths, thresholds):

        if not file_path:
            continue

        if not os.path.isfile(file_path):
            raise ValueError(
                "File {0} could not be found."
                .format(file_path)
            )

        if not threshold:
            raise ValueError("Method requires the specification of a threshold, none received")

        metric = np.loadtxt(file_path)
        if type == 'DVARS':
            metric = np.array([0.0] + metric.tolist())

        if not time_course_len:
            time_course_len = metric.shape[0]
        else:
            assert time_course_len == metric.shape[0], "Threshold metric files does not have same size."

        try:
            threshold_sd = \
                re.match(r"([0-9]*\.*[0-9]*)\s*SD", str(threshold))

            if threshold_sd:
                threshold_sd = float(threshold_sd.groups()[0])
                threshold = metric.mean() + \
                    threshold_sd * metric.std()
            else:
                threshold = float(threshold)
        except:
            raise ValueError("Could not translate threshold {0} into a "
                             "meaningful value".format(threshold))

        offending_time_points |= \
            set(np.where(metric > threshold)[0].tolist())

    extended_censors = []
    for censor in offending_time_points:
        extended_censors += list(range(
            (censor - number_of_previous_trs_to_censor),
            (censor + number_of_subsequent_trs_to_censor + 1)))

    extended_censors = [
        censor
        for censor in np.unique(extended_censors)
        if 0 <= censor < time_course_len
    ]

    censor_vector = np.ones((time_course_len, 1))
    censor_vector[extended_censors] = 0

    out_file_path = os.path.join(os.getcwd(), "censors.tsv")
    np.savetxt(
        out_file_path, censor_vector, fmt='%d', header='censor', comments=''
    )

    return out_file_path


def compute_threshold(in_file, mask, threshold):
    return threshold


def compute_pct_threshold(in_file, mask, threshold_pct):
    import nibabel as nb
    import numpy as np
    m = nb.load(mask).get_data().astype(bool)
    if not np.any(m):
        return 0.0
    d = nb.load(in_file).get_data()[m]
    return np.percentile(
        d,
        100.0 - threshold_pct
    )


def compute_sd_threshold(in_file, mask, threshold_sd):
    import nibabel as nb
    import numpy as np
    m = nb.load(mask).get_data().astype(bool)
    if not np.any(m):
        return 0.0
    d = nb.load(in_file).get_data()[m]
    return d.mean() + threshold_sd * d.std()


def temporal_variance_mask(threshold, by_slice=False, erosion=False, degree=1):

    threshold_method = "VAR"

    if isinstance(threshold, str):
        regex_match = {
            "SD": r"([0-9]+(\.[0-9]+)?)\s*SD",
            "PCT": r"([0-9]+(\.[0-9]+)?)\s*PCT",
        }

        for method, regex in regex_match.items():
            matched = re.match(regex, threshold)
            if matched:
                threshold_method = method
                threshold_value = matched.groups()[0]

    try:
        threshold_value = float(threshold_value)
    except:
        raise ValueError("Error converting threshold value {0} from {1} to a "
                        "floating point number. The threshold value can "
                        "contain SD or PCT for selecting a threshold based on "
                        "the variance distribution, otherwise it should be a "
                        "floating point number.".format(threshold_value,
                                                        threshold))

    if threshold_value < 0:
        raise ValueError("Threshold value should be positive, instead of {0}."
                        .format(threshold_value))

    if threshold_method is "PCT" and threshold_value >= 100.0:
        raise ValueError("Percentile should be less than 100, received {0}."
                        .format(threshold_value))

    threshold = threshold_value

    wf = pe.Workflow(name='tcompcor')

    input_node = pe.Node(util.IdentityInterface(fields=['functional_file_path', 'mask_file_path']), name='inputspec')
    output_node = pe.Node(util.IdentityInterface(fields=['mask']), name='outputspec')

    # C-PAC default performs linear regression while nipype performs quadratic regression
    detrend = pe.Node(afni.Detrend(args='-polort {0}'.format(degree), outputtype='NIFTI'), name='detrend')
    wf.connect(input_node, 'functional_file_path', detrend, 'in_file')

    std = pe.Node(afni.TStat(args='-nzstdev', outputtype='NIFTI'), name='std')
    wf.connect(input_node, 'mask_file_path', std, 'mask')
    wf.connect(detrend, 'out_file', std, 'in_file')

    var = pe.Node(afni.Calc(expr='a*a', outputtype='NIFTI'), name='var')
    wf.connect(std, 'out_file', var, 'in_file_a')

    if by_slice:
        slices = pe.Node(fsl.Slice(), name='slicer')
        wf.connect(var, 'out_file', slices, 'in_file')

        mask_slices = pe.Node(fsl.Slice(), name='mask_slicer')
        wf.connect(input_node, 'mask_file_path', mask_slices, 'in_file')

        mapper = pe.MapNode(util.IdentityInterface(fields=['out_file', 'mask_file']), name='slice_mapper', iterfield=['out_file', 'mask_file'])
        wf.connect(slices, 'out_files', mapper, 'out_file')
        wf.connect(mask_slices, 'out_files', mapper, 'mask_file')

    else:
        mapper_list = pe.Node(util.Merge(1), name='slice_mapper_list')
        wf.connect(var, 'out_file', mapper_list, 'in1')

        mask_mapper_list = pe.Node(util.Merge(1), name='slice_mask_mapper_list')
        wf.connect(input_node, 'mask_file_path', mask_mapper_list, 'in1')

        mapper = pe.Node(util.IdentityInterface(fields=['out_file', 'mask_file']), name='slice_mapper')
        wf.connect(mapper_list, 'out', mapper, 'out_file')
        wf.connect(mask_mapper_list, 'out', mapper, 'mask_file')

    if threshold_method is "PCT":
        threshold_node = pe.MapNode(Function(input_names=['in_file', 'mask', 'threshold_pct'],
                                             output_names=['threshold'],
                                             function=compute_pct_threshold, as_module=True),
                                    name='threshold_value', iterfield=['in_file', 'mask'])
        threshold_node.inputs.threshold_pct = threshold_value
        wf.connect(mapper, 'out_file', threshold_node, 'in_file')
        wf.connect(mapper, 'mask_file', threshold_node, 'mask')

    elif threshold_method is "SD":
        threshold_node = pe.MapNode(Function(input_names=['in_file', 'mask', 'threshold_sd'],
                                             output_names=['threshold'],
                                             function=compute_sd_threshold, as_module=True),
                                    name='threshold_value', iterfield=['in_file', 'mask'])
        threshold_node.inputs.threshold_sd = threshold_value
        wf.connect(mapper, 'out_file', threshold_node, 'in_file')
        wf.connect(mapper, 'mask_file', threshold_node, 'mask')

    else:
        threshold_node = pe.MapNode(Function(input_names=['in_file', 'mask', 'threshold'],
                                             output_names=['threshold'],
                                             function=compute_threshold, as_module=True),
                                    name='threshold_value', iterfield=['in_file', 'mask'])
        threshold_node.inputs.threshold = threshold_value
        wf.connect(mapper, 'out_file', threshold_node, 'in_file')
        wf.connect(mapper, 'mask_file', threshold_node, 'mask')

    threshold_mask = pe.MapNode(interface=fsl.maths.Threshold(), name='threshold', iterfield=['in_file', 'thresh'])
    threshold_mask.inputs.args = '-bin'
    wf.connect(mapper, 'out_file', threshold_mask, 'in_file')
    wf.connect(threshold_node, 'threshold', threshold_mask, 'thresh')

    merge_slice_masks = pe.Node(interface=fsl.Merge(), name='merge_slice_masks')
    merge_slice_masks.inputs.dimension = 'z'
    wf.connect(
        threshold_mask, 'out_file',
        merge_slice_masks, 'in_files'
    )

    wf.connect(merge_slice_masks, 'merged_file', output_node, 'mask')

    return wf


def generate_summarize_tissue_mask(nuisance_wf,
                                   pipeline_resource_pool,
                                   regressor_descriptor,
                                   regressor_selector,
                                   use_ants=True,
                                   ventricle_mask_exist=True):
    """
    Add tissue mask generation into pipeline according to the selector.

    :param nuisance_wf: Nuisance regressor workflow.
    :param pipeline_resource_pool: dictionary of available resources.
    :param regressor_descriptor: dictionary of steps to build, including keys:
        'tissue', 'resolution', 'erosion'
    :param regressor_selector: dictionary with the original selector

    :return: the full path of the 3D nifti file containing the mask created by
        this operation.
    """

    steps = [
        key
        for key in ['tissue', 'resolution', 'erosion']
        if key in regressor_descriptor
    ]

    full_mask_key = "_".join(
        regressor_descriptor[s]
        for s in steps
    )

    for step_i, step in enumerate(steps):

        mask_key = "_".join(
            regressor_descriptor[s]
            for s in steps[:step_i+1]
        )

        if mask_key in pipeline_resource_pool:
            continue

        node_mask_key = re.sub(r"[^\w]", "_", mask_key)

        prev_mask_key = "_".join(
            regressor_descriptor[s]
            for s in steps[:step_i]
        )

        if step == 'tissue':
            pass

        elif step == 'resolution':
            mask_to_epi = pe.Node(interface=fsl.FLIRT(),
                                  name='{}_flirt'
                                       .format(node_mask_key))

            mask_to_epi.inputs.interp = 'nearestneighbour'

            if regressor_selector['extraction_resolution'] == "Functional":
                # apply anat2func matrix
                mask_to_epi.inputs.apply_xfm = True
                mask_to_epi.inputs.output_type = 'NIFTI_GZ'
                nuisance_wf.connect(*(
                    pipeline_resource_pool['Functional_mean'] +
                    (mask_to_epi, 'reference')
                ))
                nuisance_wf.connect(*(
                    pipeline_resource_pool['Transformations']['anat_to_func_linear_xfm'] +
                    (mask_to_epi, 'in_matrix_file')
                ))

            else:
                resolution = regressor_selector['extraction_resolution']
                mask_to_epi.inputs.apply_isoxfm = resolution

                nuisance_wf.connect(*(
                    pipeline_resource_pool['Anatomical_{}mm'
                                           .format(resolution)] +
                    (mask_to_epi, 'reference')
                ))

            nuisance_wf.connect(*(
                pipeline_resource_pool[prev_mask_key] +
                (mask_to_epi, 'in_file')
            ))

            pipeline_resource_pool[mask_key] = \
                (mask_to_epi, 'out_file')

            if full_mask_key.startswith('CerebrospinalFluid'):
                pipeline_resource_pool = generate_summarize_tissue_mask_ventricles_masking(
                    nuisance_wf,
                    pipeline_resource_pool,
                    regressor_descriptor,
                    regressor_selector,
                    node_mask_key,
                    use_ants,
                    ventricle_mask_exist
                )

        elif step == 'erosion':

            erode_mask_node = pe.Node(
                afni.Calc(args='-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k', expr='a*(1-amongst(0,b,c,d,e,f,g))', outputtype='NIFTI_GZ'),
                name='{}'.format(node_mask_key)
            )

            nuisance_wf.connect(*(
                pipeline_resource_pool[prev_mask_key] +
                (erode_mask_node, 'in_file_a')
            ))

            pipeline_resource_pool[mask_key] = \
                (erode_mask_node, 'out_file')

    return pipeline_resource_pool, full_mask_key


def generate_summarize_tissue_mask_ventricles_masking(nuisance_wf,
                                                      pipeline_resource_pool,
                                                      regressor_descriptor,
                                                      regressor_selector,
                                                      mask_key,
                                                      use_ants=True,
                                                      ventricle_mask_exist=True):

    # Mask CSF with Ventricles
    if '{}_Unmasked'.format(mask_key) not in pipeline_resource_pool:

        # reduce CSF mask to the lateral ventricles
        mask_csf_with_lat_ven = pe.Node(interface=afni.Calc(outputtype='NIFTI_GZ'), name='{}_Ventricles'.format(mask_key))
        mask_csf_with_lat_ven.inputs.expr = 'a*b'
        mask_csf_with_lat_ven.inputs.out_file = 'csf_lat_ven_mask.nii.gz'

        if ventricle_mask_exist :
            ventricles_key = 'VentriclesToAnat'
            if 'resolution' in regressor_descriptor:
                ventricles_key += '_{}'.format(regressor_descriptor['resolution'])

            if ventricles_key not in pipeline_resource_pool:

                transforms = pipeline_resource_pool['Transformations']

                if use_ants is True:

                    # perform the transform using ANTS
                    collect_linear_transforms = pe.Node(util.Merge(3), name='{}_ants_transforms'.format(ventricles_key))

                    nuisance_wf.connect(*(transforms['anat_to_mni_linear_xfm'] + (collect_linear_transforms, 'in1')))

                    # generate inverse transform flags, which depends on the number of transforms
                    inverse_transform_flags = pe.Node(util.Function(input_names=['transform_list'],
                                                                    output_names=['inverse_transform_flags'],
                                                                    function=generate_inverse_transform_flags),
                                                                    name='{0}_inverse_transform_flags'.format(ventricles_key))
                    nuisance_wf.connect(collect_linear_transforms, 'out', inverse_transform_flags, 'transform_list')

                    lat_ven_mni_to_anat = pe.Node(interface=ants.ApplyTransforms(), name='{}_ants'.format(ventricles_key))
                    lat_ven_mni_to_anat.inputs.interpolation = 'NearestNeighbor'
                    lat_ven_mni_to_anat.inputs.dimension = 3

                    nuisance_wf.connect(inverse_transform_flags, 'inverse_transform_flags', lat_ven_mni_to_anat, 'invert_transform_flags')
                    nuisance_wf.connect(collect_linear_transforms, 'out', lat_ven_mni_to_anat, 'transforms')

                    nuisance_wf.connect(*(pipeline_resource_pool['Ventricles'] + (lat_ven_mni_to_anat, 'input_image')))
                    nuisance_wf.connect(*(pipeline_resource_pool[mask_key] + (lat_ven_mni_to_anat, 'reference_image')))

                    pipeline_resource_pool[ventricles_key] = (lat_ven_mni_to_anat, 'output_image')

                else:
                    # perform the transform using FLIRT
                    lat_ven_mni_to_anat = pe.Node(interface=fsl.FLIRT(), name='{}_flirt'.format(ventricles_key))
                    lat_ven_mni_to_anat.inputs.interp = 'nearestneighbour'

                    nuisance_wf.connect(*(transforms['mni_to_anat_linear_xfm'] + (lat_ven_mni_to_anat, 'in_matrix_file')))
                    nuisance_wf.connect(*(pipeline_resource_pool['Ventricles'] + (lat_ven_mni_to_anat, 'in_file')))
                    nuisance_wf.connect(*(pipeline_resource_pool[mask_key] + (lat_ven_mni_to_anat, 'reference')))

                    pipeline_resource_pool[ventricles_key] = (lat_ven_mni_to_anat, 'out_file')

            nuisance_wf.connect(*(pipeline_resource_pool[ventricles_key] + (mask_csf_with_lat_ven, 'in_file_a')))
            nuisance_wf.connect(*(pipeline_resource_pool[mask_key] + (mask_csf_with_lat_ven, 'in_file_b')))

            pipeline_resource_pool['{}_Unmasked'.format(mask_key)] = pipeline_resource_pool[mask_key]
            pipeline_resource_pool[mask_key] = (mask_csf_with_lat_ven, 'out_file')
        else :
            pipeline_resource_pool['{}_Unmasked'.format(mask_key)] = pipeline_resource_pool[mask_key]

        return pipeline_resource_pool


class NuisanceRegressor(object):

    def __init__(self, selector):
        self.selector = selector

        if 'Bandpass' in self.selector:
            s = self.selector['Bandpass']
            if type(s) is not dict or \
               (not s.get('bottom_frequency') and
                not s.get('top_frequency')):

                del self.selector['Bandpass']

    def get(self, key, default=None):
        return self.selector.get(key, default)

    def __contains__(self, key):
        return key in self.selector

    def __getitem__(self, key):
        return self.selector[key]

    @staticmethod
    def _derivative_params(selector):
        nr_repr = ''
        if not selector:
            return nr_repr
        if selector.get('include_squared'):
            nr_repr += 'S'
        if selector.get('include_delayed'):
            nr_repr += 'D'
        if selector.get('include_delayed_squared'):
            nr_repr += 'B'
        if selector.get('include_backdiff'):
            nr_repr += 'V'
        if selector.get('include_backdiff_squared'):
            nr_repr += 'C'
        return nr_repr

    @staticmethod
    def _summary_params(selector):
        summ = selector['summary']

        methods = {
            'PC': 'PC',
            'DetrendPC': 'DPC',
            'Mean': 'M',
            'NormMean': 'NM',
            'DetrendMean': 'DM',
            'DetrendNormMean': 'DNM',
        }

        if type(summ) == dict:
            method = summ['method']
            rep = methods[method]
            if method in ['DetrendPC', 'PC']:
                rep += "%d" % summ['components']
        else:
            rep = methods[summ]

        return rep

    @staticmethod
    def encode(selector):
        regs = OrderedDict([
            ('GreyMatter', 'GM'),
            ('WhiteMatter', 'WM'),
            ('CerebrospinalFluid', 'CSF'),
            ('tCompCor', 'tC'),
            ('aCompCor', 'aC'),
            ('GlobalSignal', 'G'),
            ('Motion', 'M'),
            ('Custom', 'T'),
            ('PolyOrt', 'P'),
            ('Bandpass', 'BP'),
            ('Censor', 'C')
        ])

        tissues = ['GreyMatter', 'WhiteMatter', 'CerebrospinalFluid']

        selectors_representations = []

        # tC-1.5PT-PC5S-SDB
        # aC-WC-2mmE-PC5-SDB

        # WM-2mmE-PC5-SDB
        # CSF-2mmE-M-SDB
        # GM-2mmE-DNM-SDB

        # G-PC5-SDB
        # M-SDB
        # C-S-FD1.5SD-D1.5SD
        # P-2
        # B-T0.01-B0.1

        for r in regs.keys():
            if r not in selector:
                continue

            s = selector[r]

            pieces = [regs[r]]

            if r in tissues:
                if s.get('extraction_resolution') and s['extraction_resolution'] != 'Functional':
                    res = "%.2gmm" % s['extraction_resolution']
                    if s.get('erode_mask'):
                        res += 'E'
                    pieces += [res]

                pieces += [NuisanceRegressor._summary_params(s)]
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == 'tCompCor':

                threshold = ""
                if s.get('by_slice'):
                    threshold += 'S'
                t = s.get('threshold')
                if t:
                    if type(t) != str:
                        t = "%.2f" % t
                    threshold += t
                if s.get('erode_mask'):
                    threshold += 'E'
                if s.get('degree'):
                    d = s.get('degree')
                    threshold += str(d)

                pieces += [threshold]
                pieces += [NuisanceRegressor._summary_params(s)]
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == 'aCompCor':
                if s.get('tissues'):
                    pieces += ["+".join([regs[t] for t in sorted(s['tissues'])])]

                if s.get('extraction_resolution'):
                    res = "%.2gmm" % s['extraction_resolution']
                    if s.get('erode_mask'):
                        res += 'E'
                    pieces += [res]

                pieces += [NuisanceRegressor._summary_params(s)]
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == 'Custom':
                for ss in s:
                    pieces += [
                        os.path.basename(ss['file'])[0:5] +
                        crc_encode(ss['file'])
                    ]

            elif r == 'GlobalSignal':
                pieces += [NuisanceRegressor._summary_params(s)]
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == 'Motion':
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == 'PolyOrt':
                pieces += ['%d' % s['degree']]

            elif r == 'Bandpass':
                if s.get('bottom_frequency'):
                    pieces += ['B%.2g' % s['bottom_frequency']]
                if s.get('top_frequency'):
                    pieces += ['T%.2g' % s['top_frequency']]

            elif r == 'Censor':
                censoring = {
                    'Kill': 'K',
                    'Zero': 'Z',
                    'Interpolate': 'I',
                    'SpikeRegression': 'S',
                }

                thresholds = {
                    'FD_J': 'FD-J',
                    'FD_P': 'FD-P',
                    'DVARS': 'DV',
                }

                pieces += [censoring[s['method']]]

                trs_range = ['0', '0']
                if s.get('number_of_previous_trs_to_censor'):
                    trs_range[0] = '%d' % s['number_of_previous_trs_to_censor']
                if s.get('number_of_subsequent_trs_to_censor'):
                    trs_range[1] = '%d' % s['number_of_subsequent_trs_to_censor']

                pieces += ['+'.join(trs_range)]

                threshs = sorted(s['thresholds'], reverse=True, key=lambda d: d['type'])
                for st in threshs:
                    thresh = thresholds[st['type']]
                    if type(st['value']) == str:
                        thresh += st['value']
                    else:
                        thresh += "%.2g" % st['value']

                    pieces += [thresh]

            selectors_representations += ['-'.join([_f for _f in pieces if _f])]

        return "_".join(selectors_representations)

    def __repr__(self):
        return NuisanceRegressor.encode(self.selector)
