
def combine_inputs_into_list(input1, input2, input3):
    outputs_list = [input1, input2, input3]
    return outputs_list


def seperate_warps_list(warp_list, selection):
    selected_warp = None
    for warp in warp_list:
        if selection == 'Warp':
            if '3Warp' in warp or '2Warp' in warp or '1Warp' in warp :
                selected_warp = warp
        else:
            if selection in warp:
                selected_warp = warp
    return selected_warp


def check_transforms(transform_list):
    transform_number = list(filter(None, transform_list))
    return[(transform_number[index]) for index in range(len(transform_number))], len(transform_number)


def generate_inverse_transform_flags(transform_list):
    inverse_transform_flags=[]
    for transform in transform_list:
        # check `blip_warp_inverse` file name and rename it
        if 'WARPINV' in transform:    
            inverse_transform_flags.append(False)
        if 'updated_affine' in transform:
            inverse_transform_flags.append(True)
        if 'Initial' in transform:
            inverse_transform_flags.append(True)
        if 'Rigid' in transform: 
            inverse_transform_flags.append(True)
        if 'Affine' in transform:
            inverse_transform_flags.append(True)
        if 'InverseWarp' in transform:
            inverse_transform_flags.append(False)
    return inverse_transform_flags


def hardcoded_reg(moving_brain, reference_brain, moving_skull,
                  reference_skull, reference_mask, moving_mask, ants_para, fixed_image_mask=None, interp=None):

    #TODO: expand transforms to cover all in ANTs para
    
    regcmd = ["antsRegistration"]
    for para_index in range(len(ants_para)):
        for para_type in ants_para[para_index]:
            if para_type == 'dimensionality':
                if ants_para[para_index][para_type] not in [2, 3, 4]:
                    err_msg = 'Dimensionality specified in ANTs parameters: %d, is not supported. ' \
                        'Change to 2, 3, or 4 and try again' % ants_para[para_index][para_type]
                    raise Exception(err_msg)
                else:
                    regcmd.append("--dimensionality")
                    regcmd.append(str(ants_para[para_index][para_type]))

            elif para_type == 'verbose':
                if ants_para[para_index][para_type] not in [0,1]:
                    err_msg = 'Verbose output option in ANTs parameters: %d, is not supported. ' \
                        'Change to 0 or 1 and try again' % ants_para[para_index][para_type]
                    raise Exception(err_msg)
                else:
                    regcmd.append("--verbose")
                    regcmd.append(str(ants_para[para_index][para_type]))

            elif para_type == 'float':
                if ants_para[para_index][para_type] not in [0,1]:
                    err_msg = 'Float option in ANTs parameters: %d, is not supported. ' \
                        'Change to 0 or 1 and try again' % ants_para[para_index][para_type]
                    raise Exception(err_msg)
                else:
                    regcmd.append("--float")
                    regcmd.append(str(ants_para[para_index][para_type]))

            elif para_type == 'collapse-output-transforms':
                if ants_para[para_index][para_type] not in [0, 1]:
                    err_msg = 'collapse-output-transforms specified in ANTs parameters: %d, is not supported. ' \
                        'Change to 0 or 1 and try again' % ants_para[para_index][para_type]
                    raise Exception(err_msg)
                else:
                    regcmd.append("--collapse-output-transforms")
                    regcmd.append(str(ants_para[para_index][para_type]))
            elif para_type == 'initial-moving-transform':
                if ants_para[para_index][para_type]['initializationFeature'] is None:
                    err_msg = 'Please specifiy initializationFeature of ANTs parameters in pipeline config. '
                    raise Exception(err_msg)
                else:
                    regcmd.append("--initial-moving-transform")
                    regcmd.append("[{0},{1},{2}]".format(
                        reference_brain, moving_brain, ants_para[para_index][para_type]['initializationFeature']))

            elif para_type == 'transforms':
                for trans_index in range(len(ants_para[para_index][para_type])):
                    for trans_type in ants_para[para_index][para_type][trans_index]:
                        regcmd.append("--transform")
                        if trans_type == 'Rigid' or trans_type == 'Affine':
                            if ants_para[para_index][para_type][trans_index][trans_type]['gradientStep'] is None:
                                err_msg = 'Please specifiy % s Gradient Step of ANTs parameters in pipeline config. ' % trans_type
                                raise Exception(err_msg)
                            else:
                                regcmd.append("{0}[{1}]".format(
                                    trans_type, ants_para[para_index][para_type][trans_index][trans_type]['gradientStep']))

                        if trans_type == 'SyN':
                            if ants_para[para_index][para_type][trans_index][trans_type]['gradientStep'] is None:
                                err_msg = 'Please specifiy % s Gradient Step of ANTs parameters in pipeline config. ' % trans_type
                                raise Exception(err_msg)
                            else:
                                SyN_para = []
                                SyN_para.append("{0}".format(
                                    ants_para[para_index][para_type][trans_index][trans_type]['gradientStep']))
                                if ants_para[para_index][para_type][trans_index][trans_type]['updateFieldVarianceInVoxelSpace'] is not None:
                                    SyN_para.append("{0}".format(
                                        ants_para[para_index][para_type][trans_index][trans_type]['updateFieldVarianceInVoxelSpace']))
                                if ants_para[para_index][para_type][trans_index][trans_type]['totalFieldVarianceInVoxelSpace'] is not None:
                                    SyN_para.append("{0}".format(
                                        ants_para[para_index][para_type][trans_index][trans_type]['totalFieldVarianceInVoxelSpace']))
                                SyN_para = ','.join([str(elem)
                                                    for elem in SyN_para])
                                regcmd.append("{0}[{1}]".format(
                                    trans_type, SyN_para))

                        if ants_para[para_index][para_type][trans_index][trans_type]['metric']['type'] == 'MI':
                            if ants_para[para_index][para_type][trans_index][trans_type]['metric']['metricWeight'] is None or ants_para[para_index][para_type][trans_index][trans_type]['metric']['numberOfBins'] is None:
                                err_msg = 'Please specifiy metricWeight and numberOfBins for metric MI of ANTs parameters in pipeline config.'
                                raise Exception(err_msg)
                            else:
                                MI_para = []
                                MI_para.append("{0},{1}".format(ants_para[para_index][para_type][trans_index][trans_type]['metric']
                                                                ['metricWeight'], ants_para[para_index][para_type][trans_index][trans_type]['metric']['numberOfBins']))
                                if 'samplingStrategy' in ants_para[para_index][para_type][trans_index][trans_type]['metric'] and ants_para[para_index][para_type][trans_index][trans_type]['metric']['samplingStrategy'] in ['None', 'Regular', 'Random']:
                                    MI_para.append("{0}".format(
                                        ants_para[para_index][para_type][trans_index][trans_type]['metric']['samplingStrategy']))
                                if 'samplingPercentage' in ants_para[para_index][para_type][trans_index][trans_type]['metric'] and ants_para[para_index][para_type][trans_index][trans_type]['metric']['samplingPercentage'] is not None:
                                    MI_para.append("{0}".format(
                                        ants_para[para_index][para_type][trans_index][trans_type]['metric']['samplingPercentage']))
                                MI_para = ','.join([str(elem) for elem in MI_para])
                                regcmd.append("--metric")
                                regcmd.append("MI[{0},{1},{2}]".format(
                                    reference_brain, moving_brain, MI_para))

                        if ants_para[para_index][para_type][trans_index][trans_type]['metric']['type'] == 'CC':
                            if ants_para[para_index][para_type][trans_index][trans_type]['metric']['metricWeight'] is None or ants_para[para_index][para_type][trans_index][trans_type]['metric']['radius'] is None:
                                err_msg = 'Please specifiy metricWeight and radius for metric CC of ANTs parameters in pipeline config.'
                                raise Exception(err_msg)
                            else:
                                CC_para = []
                                CC_para.append("{0},{1}".format(ants_para[para_index][para_type][trans_index][trans_type]['metric']
                                                                ['metricWeight'], ants_para[para_index][para_type][trans_index][trans_type]['metric']['radius']))
                                if 'samplingStrategy' in ants_para[para_index][para_type][trans_index][trans_type]['metric'] and ants_para[para_index][para_type][trans_index][trans_type]['metric']['samplingStrategy'] in ['None', 'Regular', 'Random']:
                                    CC_para.append("{0}".format(
                                        ants_para[para_index][para_type][trans_index][trans_type]['metric']['samplingStrategy']))
                                if 'samplingPercentage' in ants_para[para_index][para_type][trans_index][trans_type]['metric'] and ants_para[para_index][para_type][trans_index][trans_type]['metric']['samplingPercentage'] is not None:
                                    CC_para.append("{0}".format(
                                        ants_para[para_index][para_type][trans_index][trans_type]['metric']['samplingPercentage']))
                                CC_para = ','.join([str(elem) for elem in CC_para])
                                regcmd.append("--metric")
                                regcmd.append("CC[{0},{1},{2}]".format(
                                    reference_skull, moving_skull, CC_para))

                        if 'convergence' in ants_para[para_index][para_type][trans_index][trans_type]:
                            convergence_para = []
                            if ants_para[para_index][para_type][trans_index][trans_type]['convergence']['iteration'] is None:
                                err_msg = 'Please specifiy convergence iteration of ANTs parameters in pipeline config.'
                                raise Exception(err_msg)
                            else:
                                convergence_para.append("{0}".format(
                                    ants_para[para_index][para_type][trans_index][trans_type]['convergence']['iteration']))
                                if 'convergenceThreshold' in ants_para[para_index][para_type][trans_index][trans_type]['convergence'] and ants_para[para_index][para_type][trans_index][trans_type]['convergence']['convergenceThreshold'] is not None:
                                    convergence_para.append("{0}".format(
                                        ants_para[para_index][para_type][trans_index][trans_type]['convergence']['convergenceThreshold']))
                                if 'convergenceWindowSize' in ants_para[para_index][para_type][trans_index][trans_type]['convergence'] and ants_para[para_index][para_type][trans_index][trans_type]['convergence']['convergenceWindowSize'] is not None:
                                    convergence_para.append("{0}".format(
                                        ants_para[para_index][para_type][trans_index][trans_type]['convergence']['convergenceWindowSize']))
                                convergence_para = ','.join(
                                    [str(elem) for elem in convergence_para])
                                regcmd.append("--convergence")
                                regcmd.append("[{0}]".format(convergence_para))

                        if 'smoothing-sigmas' in ants_para[para_index][para_type][trans_index][trans_type] and ants_para[para_index][para_type][trans_index][trans_type]['smoothing-sigmas'] is not None:
                            regcmd.append("--smoothing-sigmas")
                            regcmd.append("{0}".format(
                                ants_para[para_index][para_type][trans_index][trans_type]['smoothing-sigmas']))

                        if 'shrink-factors' in ants_para[para_index][para_type][trans_index][trans_type] and ants_para[para_index][para_type][trans_index][trans_type]['shrink-factors'] is not None:
                            regcmd.append("--shrink-factors")
                            regcmd.append("{0}".format(
                                ants_para[para_index][para_type][trans_index][trans_type]['shrink-factors']))

                        if 'use-histogram-matching' in ants_para[para_index][para_type][trans_index][trans_type]:
                            if ants_para[para_index][para_type][trans_index][trans_type]['use-histogram-matching']:
                                regcmd.append("--use-histogram-matching")
                                regcmd.append("1")
                            else:
                                continue

                        if 'winsorize-image-intensities' in ants_para[para_index][para_type][trans_index][trans_type] and ants_para[para_index][para_type][trans_index][trans_type]['winsorize-image-intensities']['lowerQuantile'] is not None and ants_para[para_index][para_type][trans_index][trans_type]['winsorize-image-intensities']['upperQuantile'] is not None:
                            regcmd.append("--winsorize-image-intensities")
                            regcmd.append("[{0},{1}]".format(ants_para[para_index][para_type][trans_index][trans_type]['winsorize-image-intensities']
                                                            ['lowerQuantile'], ants_para[para_index][para_type][trans_index][trans_type]['winsorize-image-intensities']['upperQuantile']))
    
            elif para_type == 'masks':
                # lesion preproc has 
                if fixed_image_mask is not None:
                    regcmd.append("--masks")
                    regcmd.append(str(fixed_image_mask))
                else: 
                    if ants_para[para_index][para_type]['fixed_image_mask'] == False and ants_para[para_index][para_type]['moving_image_mask'] == True:
                        err_msg = 'Masks option in ANTs parameters: %d is not supported. ' \
                            'Please set `fixed_image_mask` as True. ' \
                                'Or set both `fixed_image_mask` and `moving_image_mask` as False' % ants_para[para_index][para_type]
                        raise Exception(err_msg)
                    elif ants_para[para_index][para_type]['fixed_image_mask'] == True and ants_para[para_index][para_type]['moving_image_mask'] == True:
                        regcmd.append("--masks")
                        regcmd.append('['+str(reference_mask)+','+str(moving_mask)+']')
                    elif ants_para[para_index][para_type]['fixed_image_mask'] == True and ants_para[para_index][para_type]['moving_image_mask'] == False:
                        regcmd.append("--masks")
                        regcmd.append('['+str(reference_mask)+']')
                    else:
                        continue

    if interp is not None: 
        regcmd.append("--interpolation")
        regcmd.append("{0}".format(interp))

    regcmd.append("--output")
    regcmd.append("[transform,transform_Warped.nii.gz]")

    # write out the actual command-line entry for testing/validation later
    command_file = os.path.join(os.getcwd(), 'command.txt')
    with open(command_file, 'wt') as f:
        f.write(' '.join(regcmd))

    try:
        retcode = subprocess.check_output(regcmd)
    except Exception as e:
        raise Exception('[!] ANTS registration did not complete successfully.'
                        '\n\nError details:\n{0}\n{1}\n'.format(e, e.output))

    warp_list = []
    warped_image = None

    files = [f for f in os.listdir('.') if os.path.isfile(f)]

    for f in files:
        if ("transform" in f) and ("Warped" not in f):
            warp_list.append(os.getcwd() + "/" + f)
        if "Warped" in f:
            warped_image = os.getcwd() + "/" + f

    if not warped_image:
        raise Exception("\n\n[!] No registration output file found. ANTS "
                        "registration may not have completed "
                        "successfully.\n\n")

    return warp_list, warped_image

def change_itk_transform_type(input_affine_file):
    """
    this function takes in the affine.txt produced by the c3d_affine_tool
    (which converted an FSL FLIRT affine.mat into the affine.txt)

    it then modifies the 'Transform Type' of this affine.txt so that it is
    compatible with the antsApplyTransforms tool and produces a new affine
    file titled 'updated_affine.txt'
    """

    new_file_lines = []

    with open(input_affine_file) as f:
        for line in f:
            if 'Transform:' in line:
                if 'MatrixOffsetTransformBase_double_3_3' in line:
                    transform_line = 'Transform: AffineTransform_double_3_3\n'
                    new_file_lines.append(transform_line)
            else:
                new_file_lines.append(line)

    updated_affine_file = os.path.join(os.getcwd(), 'updated_affine.txt')

    with open(updated_affine_file, 'wt') as f:
        for line in new_file_lines:
            f.write(line)

    return updated_affine_file


def run_ants_apply_warp(moving_image, reference, initial=None, rigid=None,
                        affine=None, nonlinear=None, func_to_anat=None,
                        anatomical_brain=None, dim=3, interp='Linear',
                        inverse=False):
    """Apply a transform using ANTs transforms."""

    import os
    import subprocess

    # if inverse:
    #     inverse = 1
    # else:
    #     inverse = 0

    if func_to_anat:
        # this assumes the func->anat affine transform is FSL-based and needs
        # to be converted to ITK format via c3d_affine_tool
        cmd = ['c3d_affine_tool', '-ref', anatomical_brain, '-src',
               moving_image, func_to_anat, '-fsl2ras', '-oitk', 'affine.txt']
        retcode = subprocess.check_output(cmd)
        func_to_anat = change_itk_transform_type(os.path.join(os.getcwd(),
                                                              'affine.txt'))

    out_image = os.path.join(os.getcwd(), moving_image[moving_image.rindex('/')+1:moving_image.rindex('.nii.gz')]+'_warp.nii.gz')

    cmd = ['antsApplyTransforms', '-d', str(dim), '-i', moving_image, '-r',
           reference, '-o', out_image, '-n', interp]

    if nonlinear:
        cmd.append('-t')
        if inverse:
            cmd.append('[{0}, {1}]'.format(os.path.abspath(nonlinear), '1'))
        else:
            cmd.append(os.path.abspath(nonlinear))

    if affine:
        cmd.append('-t')
        if inverse:
            cmd.append('[{0}, {1}]'.format(os.path.abspath(affine), '1'))
        else:
            cmd.append(os.path.abspath(affine))

    if rigid:
        cmd.append('-t')
        if inverse:
            cmd.append('[{0}, {1}]'.format(os.path.abspath(rigid), '1'))
        else:
            cmd.append(os.path.abspath(rigid))

    if initial:
        cmd.append('-t')
        if inverse:
            cmd.append('[{0}, {1}]'.format(os.path.abspath(initial), '1'))
        else:
            cmd.append(os.path.abspath(initial))

    if func_to_anat:
        cmd.append('-t')
        if inverse:
            cmd.append('[{0}, {1}]'.format(os.path.abspath(func_to_anat), '1'))
        else:
            cmd.append(os.path.abspath(func_to_anat))

    retcode = subprocess.check_output(cmd)

    return out_image


def cpac_ants_apply_nonlinear_inverse_warp(cpac_dir, moving_image, reference,
                                           dim=3, interp='Linear'):
    """Run antsApplyTransforms for inverse warping when given a C-PAC output
    directory."""

    import os

    cpac_dir = os.path.abspath(cpac_dir)

    for dir in os.listdir(cpac_dir):
        if 'ants_initial_xfm' in dir:
            pass

    #run_ants_apply_warp()
