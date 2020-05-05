# -*- coding: utf-8 -*-
import os
# import ntpath
import numpy as np
import six
from multiprocessing.dummy import Pool as ThreadPool

from CPAC.utils.nifti_utils import nifti_image_input
import nibabel as nib
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
from nipype.interfaces.fsl import ConvertXFM


def read_ants_mat(ants_mat_file):
    if not os.path.exists(ants_mat_file):
        raise ValueError(str(ants_mat_file) + " does not exist.")

    with open(ants_mat_file) as f:
        for line in f:
            tmp = line.split(':')
            if tmp[0] == 'Parameters':
                oth_transform = np.reshape(
                    np.fromstring(tmp[1], float, sep=' '), (-1, 3))
            if tmp[0] == 'FixedParameters':
                translation = np.fromstring(tmp[1], float, sep=' ')
    return translation, oth_transform


def read_mat(input_mat):
    if isinstance(input_mat, np.ndarray):
        mat = input_mat
    elif isinstance(input_mat, six.string_types):
        if os.path.exists(input_mat):
            mat = np.loadtxt(input_mat)
        else:
            raise IOError("ERROR norm_transformation: " + input_mat +
                          " file does not exist")
    else:
        raise TypeError("ERROR norm_transformation: input_mat should be" +
                        " either a str or a numpy.ndarray matrix")

    if mat.shape != (4, 4):
        raise ValueError("ERROR norm_transformation: the matrix should be 4x4")

    # Translation vector
    translation = mat[0:3, 3]
    # 3x3 matrice of rotation, scaling and skewing
    oth_transform = mat[0:3, 0:3]

    return translation, oth_transform


def norm_transformations(translation, oth_transform):
    tr_norm = np.linalg.norm(translation)
    affine_norm = np.linalg.norm(oth_transform - np.identity(3), 'fro')
    return pow(tr_norm, 2) + pow(affine_norm, 2)


def norm_transformation(input_mat):
    """
    Calculate the squared norm of the translation + squared Frobenium norm
    of the difference between other affine transformations and the identity
    from an fsl FLIRT transformation matrix
    Parameters
    ----------
    input_mat: str or numpy.ndarray
        Either the path to text file matrix or a matrix already imported.

    Returns
    -------
        numpy.float64
            squared norm of the translation + squared Frobenius norm of the
            difference between other affine transformations and the identity
    """
    if isinstance(input_mat, np.ndarray):
        mat = input_mat
    elif isinstance(input_mat, six.string_types):
        if os.path.exists(input_mat):
            mat = np.loadtxt(input_mat)
        else:
            raise IOError("ERROR norm_transformation: " + input_mat +
                          " file does not exist")
    else:
        raise TypeError("ERROR norm_transformation: input_mat should be" +
                        " either a str (file_path) or a numpy.ndarray matrix")

    if mat.shape != (4, 4):
        raise ValueError("ERROR norm_transformation: the matrix should be 4x4")

    # Translation vector
    translation = mat[0:3, 3]
    # 3x3 matrice of rotation, scaling and skewing
    oth_affine_transform = mat[0:3, 0:3]
    tr_norm = np.linalg.norm(translation)
    affine_norm = np.linalg.norm(oth_affine_transform - np.identity(3), 'fro')
    return pow(tr_norm, 2) + pow(affine_norm, 2)


def template_convergence(mat_file, mat_type='matrix',
                         convergence_threshold=np.finfo(np.float64).eps):
    """
    Calculate the distance between transformation matrix with a matrix of no
    transformation
    Parameters
    ----------
    mat_file: str
        path to an fsl flirt matrix
    mat_type: str
        'matrix'(default), 'ITK'
        The type of matrix used to represent the transformations
    convergence_threshold: float
        (numpy.finfo(np.float64).eps (default)) threshold for the convergence
        The threshold is how different from no transformation is the
        transformation matrix.

    Returns
    -------

    """
    if mat_type == 'matrix':
        translation, oth_transform = read_mat(mat_file)
    elif mat_type == 'ITK':
        translation, oth_transform = read_ants_mat(mat_file)
    else:
        raise ValueError("ERROR template_convergence: this matrix type does " +
                         "not exist")
    distance = norm_transformations(translation, oth_transform)
    print("distance = " + str(abs(distance)) + ' | ')

    return abs(distance) <= convergence_threshold


def create_temporary_template(img_list, out_path, avg_method='median'):
    """
    Average all the 3D images of the list into one 3D image
    WARNING---the function assumes that all the images have the same header,
    the output image will have the same header as the first image of the list
    Parameters---
    ----------
    img_list: list of str
        list of images paths
    avg_method: str
        function names from numpy library such as 'median', 'mean', 'std' ...

    Returns
    -------
    tmp_template: Nifti1Image
    """
    if not img_list:
        raise ValueError('ERROR create_temporary_template: image list is empty')

    if len(img_list) == 1:
        return img_list[0]

    # ALIGN CENTERS

    avg_data = getattr(np, avg_method)(
        np.asarray([nifti_image_input(img).get_data() for img in img_list]), 0)

    nii = nib.Nifti1Image(avg_data, nifti_image_input(img_list[0]).affine)
    nib.save(nii, out_path)
    return out_path


def register_img_list(img_list, ref_img, dof=12,
                      interp='trilinear', cost='corratio', thread_pool=2):
    """
    Register a list of images to the reference image.
    Parameters
    ----------
    img_list: list of str
        list of images paths
    ref_img: str
        path to the reference image to which the images will be registered
    dof: integer (int of long)
        number of transform degrees of freedom (FLIRT) (12 by default)
    interp: str
        ('trilinear' (default) or 'nearestneighbour' or 'sinc' or 'spline')
        final interpolation method used in reslicing
    cost: str
        ('mutualinfo' or 'corratio' (default) or 'normcorr' or 'normmi' or
         'leastsq' or 'labeldiff' or 'bbr')
        cost function
    thread_pool: int or multiprocessing.dummy.Pool
        (default 2) number of threads. You can also provide a Pool so the
        node will be added to it to be run.

    Returns
    -------
    multiple_linear_reg: list of Node
        each Node 'node' has been run and
        node.inputs.out_file contains the path to the registered image
        node.inputs.out_matrix_file contains the path to the transformation
        matrix
    """
    if not img_list:
        raise ValueError('ERROR register_img_list: image list is empty')

    output_img_list = [os.path.join(os.getcwd(), os.path.basename(img))
                       for img in img_list]

    output_mat_list = [os.path.join(os.getcwd(),
                                    str(os.path.basename(img).split('.')[0])
                                    + '.mat')
                       for img in img_list]

    def flirt_node(in_img, output_img, output_mat):
        linear_reg = fsl.FLIRT()
        linear_reg.inputs.in_file = in_img
        linear_reg.inputs.out_file = output_img
        linear_reg.inputs.out_matrix_file = output_mat

        linear_reg.inputs.cost = cost
        linear_reg.inputs.dof = dof
        linear_reg.inputs.interp = interp
        linear_reg.inputs.reference = ref_img

        return linear_reg

    if isinstance(thread_pool, int):
        pool = ThreadPool(thread_pool)
    else:
        pool = thread_pool

    node_list = [flirt_node(img, out_img, out_mat)
                 for (img, out_img, out_mat) in zip(
                 img_list, output_img_list, output_mat_list)]
    pool.map(lambda node: node.run(), node_list)

    if isinstance(thread_pool, int):
        pool.close()
        pool.join()

    return node_list


def template_creation_flirt(img_list, init_reg=None, avg_method='median', dof=12,
                            interp='trilinear', cost='corratio', mat_type='matrix',
                            convergence_threshold=-1, thread_pool=2):
    """
    Parameters
    ----------
    img_list : list of str
        list of images paths
    init_reg : list of Node
        (default None so no initial registration performed)
        the output of the function register_img_list with another reference
        Reuter et al. 2012 (NeuroImage) section "Improved template estimation"
        doi:10.1016/j.neuroimage.2012.02.084 uses a ramdomly
        selected image from the input dataset
    avg_method : str
        function names from numpy library such as 'median', 'mean', 'std' ...
    dof : integer (int of long)
        number of transform degrees of freedom (FLIRT) (12 by default)
    interp : str
        ('trilinear' (default) or 'nearestneighbour' or 'sinc' or 'spline')
        final interpolation method used in reslicing
    cost : str
        ('mutualinfo' or 'corratio' (default) or 'normcorr' or 'normmi' or
         'leastsq' or 'labeldiff' or 'bbr')
        cost function
    mat_type : str
        'matrix'(default), 'ITK'
        The type of matrix used to represent the transformations
    convergence_threshold : float
        (numpy.finfo(np.float64).eps (default)) threshold for the convergence
        The threshold is how different from no transformation is the
        transformation matrix.
    thread_pool : int or multiprocessing.dummy.Pool
        (default 2) number of threads. You can also provide a Pool so the
        node will be added to it to be run.

    Returns
    -------
    template : str
        path to the final template

    """
    # DEBUG to skip the longitudinal template generation which takes a lot of time.
    # return 'CECI_EST_UN_TEST'

    if isinstance(thread_pool, int):
        pool = ThreadPool(thread_pool)
    else:
        pool = thread_pool

    if convergence_threshold == -1:
        convergence_threshold = np.finfo(np.float64).eps

    if not img_list:
        print('ERROR create_temporary_template: image list is empty')

    if len(img_list) == 1:
        print("img_list contains only 1 image, no template calculated")
        return img_list[0]

    final_warp_list = []
    final_warp_list_filenames = [os.path.join(
        os.getcwd(), str(os.path.basename(img).split('.')[0]) + '_anat_to_template.mat') for img in img_list]

    # I added this part because it is mentioned in the paper but I actually never used it
    # You could run a first register_img_list() with a selected image as starting point and
    # give the output to this function
    if init_reg is not None:
        if isinstance(init_reg, list):
            image_list = [node.inputs.out_file for node in init_reg]
            mat_list = [node.inputs.out_matrix_file for node in init_reg]
            final_warp_list = mat_list
            # test if every transformation matrix has reached the convergence
            convergence_list = [template_convergence(
                mat, mat_type, convergence_threshold) for mat in mat_list]
            converged = all(convergence_list)
        else:
            raise ValueError("init_reg must be a list of FLIRT nipype nodes files")
    else:
        image_list = img_list
        converged = False

    tmp_template = os.path.join(os.getcwd(), 'tmp_template.nii.gz')

    """ First is calculated an average image of the dataset to be the temporary template
    and the loop stops when this temporary template is close enough (with a transformation
    distance smaller than the threshold) to all the images of the precedent iteration.
    """
    while not converged:
        tmp_template = create_temporary_template(image_list,
                                                 out_path=tmp_template,
                                                 avg_method=avg_method)
        reg_list_node = register_img_list(image_list,
                                          ref_img=tmp_template,
                                          dof=dof,
                                          interp=interp,
                                          cost=cost)

        mat_list = [node.inputs.out_matrix_file for node in reg_list_node]
        if len(final_warp_list) == 0:
            final_warp_list = mat_list
        for index, mat in enumerate(mat_list):
            cmd = "convert_xfm -omat %s -inverse %s" % (final_warp_list_filenames[index], final_warp_list[index])
            os.system(cmd)
            final_warp_list[index] = final_warp_list_filenames[index]

        image_list = [node.inputs.out_file for node in reg_list_node]
        # test if every transformation matrix has reached the convergence
        convergence_list = [template_convergence(
            mat, mat_type, convergence_threshold) for mat in mat_list]
        converged = all(convergence_list)

    if isinstance(thread_pool, int):
        pool.close()
        pool.join()

    template = tmp_template
    return template, final_warp_list


def subject_specific_template(workflow_name='subject_specific_template',
                              method='flirt'):
    """
    Parameters
    ----------
    workflow_name
    method

    Returns
    -------
    """
    imports = [
        'import os',
        'import numpy as np',
        'from multiprocessing.dummy import Pool as ThreadPool',
        'from nipype.interfaces.fsl import ConvertXFM',
        'from CPAC.longitudinal_pipeline.longitudinal_preproc import ('
        '   create_temporary_template,'
        '   register_img_list,'
        '   template_convergence'
        ')'
    ]
    if method == 'flirt':
        template_gen_node = pe.Node(
            util.Function(
                input_names=[
                    'img_list',
                    'init_reg', 
                    'avg_method', 
                    'dof',
                    'interp', 
                    'cost',
                    'mat_type',
                    'convergence_threshold',
                    'thread_pool'],
                output_names=['template', 
                    'final_warp_list'],
                imports=imports,
                function=template_creation_flirt
            ),
            name=workflow_name
        )
    else:
        raise ValueError(str(method)
                         + 'this method has not yet been implemented')

    return template_gen_node

