def seperate_warps_list(warp_list, selection):

     return warp_list[selection]



def change_itk_transform_type(input_affine_file):

    '''
    this function takes in the affine.txt produced by the c3d_affine_tool
    (which converted an FSL FLIRT affine.mat into the affine.txt)

    it then modifies the 'Transform Type' of this affine.txt so that it is
    compatible with the antsApplyTransforms tool and produces a new affine
    file titled 'updated_affine.txt'
    '''

    import os

    new_file_lines = []

    with open(input_affine_file) as f:

        for line in f:

            if 'Transform:' in line:

                if 'MatrixOffsetTransformBase_double_3_3' in line:

                    transform_line = 'Transform: AffineTransform_double_3_3'
                    new_file_lines.append(transform_line)

            else:

                new_file_lines.append(line)


    updated_affine_file = os.path.join(os.getcwd(), 'updated_affine.txt')

    outfile = open(updated_affine_file, 'wt')

    for line in new_file_lines:

        print >>outfile, line.strip('\n')

    outfile.close()


    return updated_affine_file
