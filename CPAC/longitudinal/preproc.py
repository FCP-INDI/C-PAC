# -*- coding: utf-8 -*-
# Copyright (C) 2020-2024  C-PAC Developers

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
"""Preprocessing for longitudinal pipelines."""

from collections import Counter
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing.pool import Pool
import os
from typing import Literal, Optional

import numpy as np
import nibabel as nib
from nipype.interfaces import fsl

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.interfaces import Function
from CPAC.utils.monitoring import IFLOGGER
from CPAC.utils.nifti_utils import nifti_image_input


def read_ants_mat(ants_mat_file):
    """Read a matrix, returning (translation) and (other transformations) matrices."""
    if not os.path.exists(ants_mat_file):
        raise ValueError(str(ants_mat_file) + " does not exist.")

    with open(ants_mat_file) as f:
        for line in f:
            tmp = line.split(":")
            if tmp[0] == "Parameters":
                oth_transform = np.reshape(
                    np.fromstring(tmp[1], float, sep=" "), (-1, 3)
                )
            if tmp[0] == "FixedParameters":
                translation = np.fromstring(tmp[1], float, sep=" ")
    return translation, oth_transform


def read_mat(input_mat):
    """Read a matrix, returning (translation) and (other transformations) matrices."""
    if isinstance(input_mat, np.ndarray):
        mat = input_mat
    elif isinstance(input_mat, str):
        if os.path.exists(input_mat):
            mat = np.loadtxt(input_mat)
        else:
            raise IOError(
                "ERROR norm_transformation: " + input_mat + " file does not exist"
            )
    else:
        raise TypeError(
            "ERROR norm_transformation: input_mat should be"
            + " either a str or a numpy.ndarray matrix"
        )

    if mat.shape != (4, 4):
        msg = "ERROR norm_transformation: the matrix should be 4x4"
        raise ValueError(msg)

    # Translation vector
    translation = mat[0:3, 3]
    # 3x3 matrice of rotation, scaling and skewing
    oth_transform = mat[0:3, 0:3]

    return translation, oth_transform


def norm_transformations(translation, oth_transform):
    """Calculate the sum of squares of norm translation and Frobenius norm."""
    tr_norm = np.linalg.norm(translation)
    affine_norm = np.linalg.norm(oth_transform - np.identity(3), "fro")
    return pow(tr_norm, 2) + pow(affine_norm, 2)


def norm_transformation(input_mat):
    """Calculate the sum of squares of norm translation and Frobenius norm.

    Calculate the squared norm of the translation + squared Frobenium norm
    of the difference between other affine transformations and the identity
    from an fsl FLIRT transformation matrix.

    Parameters
    ----------
    input_mat : str or ~numpy.ndarray
        Either the path to text file matrix or a matrix already imported.

    Returns
    -------
    ~numpy.float64
        squared norm of the translation + squared Frobenius norm of the
        difference between other affine transformations and the identity
    """
    if isinstance(input_mat, np.ndarray):
        mat = input_mat
    elif isinstance(input_mat, str):
        if os.path.exists(input_mat):
            mat = np.loadtxt(input_mat)
        else:
            msg = f"ERROR norm_transformation: {input_mat} file does not exist"
            raise IOError(msg)
    else:
        msg = (
            "ERROR norm_transformation: input_mat should be either a str"
            " (file_path) or a numpy.ndarray matrix"
        )
        raise TypeError(msg)

    if mat.shape != (4, 4):
        msg = "ERROR norm_transformation: the matrix should be 4x4"
        raise ValueError(msg)

    # Translation vector
    translation = mat[0:3, 3]
    # 3x3 matrice of rotation, scaling and skewing
    oth_affine_transform = mat[0:3, 0:3]
    tr_norm = np.linalg.norm(translation)
    affine_norm = np.linalg.norm(oth_affine_transform - np.identity(3), "fro")
    return pow(tr_norm, 2) + pow(affine_norm, 2)


def template_convergence(
    mat_file: str,
    mat_type: Literal["matrix", "ITK"] = "matrix",
    convergence_threshold: float | np.float64 = np.finfo(np.float64).eps,
) -> bool:
    """Check that the deistance between matrices is smaller than the threshold.

    Calculate the distance between transformation matrix with a matrix of no transformation.

    Parameters
    ----------
    mat_file
        path to an fsl flirt matrix
    mat_type
        The type of matrix used to represent the transformations
    convergence_threshold
        The threshold is how different from no transformation is the
        transformation matrix.
    """
    if mat_type == "matrix":
        translation, oth_transform = read_mat(mat_file)
    elif mat_type == "ITK":
        translation, oth_transform = read_ants_mat(mat_file)
    else:
        msg = f"template_convergence: matrix type {mat_type} does not exist"
        raise ValueError(msg)
    distance = norm_transformations(translation, oth_transform)
    IFLOGGER.info("distance = %s", abs(distance))

    return abs(distance) <= convergence_threshold


def create_temporary_template(
    input_brain_list,
    input_skull_list,
    output_brain_path,
    output_skull_path,
    avg_method="median",
):
    """Average all the 3D images of the list into one 3D image.

    Warnings
    --------
    The function assumes that all the images have the same header,
    the output image will have the same header as the first image of the list.

    Parameters
    ----------
    input_brain_list : list of str
        list of brain image paths
    input_skull_list : list of str
        list of skull image paths
    output_brain_path : ~nibabel.Nifti1Image
        temporary longitudinal brain template
    output_skull_path : ~nibabel.Nifti1Image
        temporary longitudinal skull template
    avg_method : str
        function names from numpy library such as 'median', 'mean', 'std' ...

    Returns
    -------
    output_brain_path : ~nibabel.Nifti1Image
        temporary longitudinal brain template
    output_skull_path : ~nibabel.Nifti1Image
        temporary longitudinal skull template
    """
    if not input_brain_list or not input_skull_list:
        msg = "ERROR create_temporary_template: image list is empty"
        raise ValueError(msg)

    if len(input_brain_list) == 1 and len(input_skull_list) == 1:
        return input_brain_list[0], input_skull_list[0]

    # ALIGN CENTERS
    avg_brain_data = getattr(np, avg_method)(
        np.asarray([nifti_image_input(img).get_fdata() for img in input_brain_list]), 0
    )

    avg_skull_data = getattr(np, avg_method)(
        np.asarray([nifti_image_input(img).get_fdata() for img in input_skull_list]), 0
    )

    nii_brain = nib.Nifti1Image(
        avg_brain_data, nifti_image_input(input_brain_list[0]).affine
    )
    nii_skull = nib.Nifti1Image(
        avg_skull_data, nifti_image_input(input_skull_list[0]).affine
    )

    nib.save(nii_brain, output_brain_path)
    nib.save(nii_skull, output_skull_path)

    return output_brain_path, output_skull_path


def register_img_list(
    input_brain_list,
    ref_img,
    dof=12,
    interp="trilinear",
    cost="corratio",
    thread_pool=2,
    unique_id_list=None,
):
    """
    Register a list of images to the reference image.

    Parameters
    ----------
    input_brain_list : list of str
        list of brain image paths
    ref_img : str
        path to the reference image to which the images will be registered
    dof : integer (int of long)
        number of transform degrees of freedom (FLIRT) (12 by default)
    interp : str
        ('trilinear' (default) or 'nearestneighbour' or 'sinc' or 'spline')
        final interpolation method used in reslicing
    cost : str
        ('mutualinfo' or 'corratio' (default) or 'normcorr' or 'normmi' or
         'leastsq' or 'labeldiff' or 'bbr')
        cost function
    thread_pool : int or multiprocessing.dummy.Pool
        (default 2) number of threads. You can also provide a Pool so the
        node will be added to it to be run.
    duplicated_basename : boolean
        whether there exists duplicated basename which may happen in non-BIDS dataset
    unique_id_list : list
        a list of unique IDs in data

    Returns
    -------
    node_list : list of Node
        each Node 'node' has been run and
        node.inputs.out_file contains the path to the registered image
        node.inputs.out_matrix_file contains the path to the transformation
        matrix
    """
    if not input_brain_list:
        msg = "ERROR register_img_list: image list is empty"
        raise ValueError(msg)

    basename_list = [
        str(os.path.basename(img).split(".")[0]) for img in input_brain_list
    ]
    counter = Counter(basename_list)
    duplicated_basename_list = [i for i, j in counter.items() if j > 1]

    if not duplicated_basename_list:
        output_img_list = [
            os.path.join(os.getcwd(), os.path.basename(img)) for img in input_brain_list
        ]

        output_mat_list = [
            os.path.join(os.getcwd(), str(os.path.basename(img).split(".")[0]) + ".mat")
            for img in input_brain_list
        ]
    else:
        output_img_list = [
            os.path.join(
                os.getcwd(),
                str(os.path.basename(img).split(".")[0])
                + "_"
                + unique_id_list[i]
                + ".nii.gz",
            )
            for i, img in enumerate(input_brain_list)
        ]

        output_mat_list = [
            os.path.join(
                os.getcwd(),
                str(os.path.basename(img).split(".")[0])
                + "_"
                + unique_id_list[i]
                + ".mat",
            )
            for i, img in enumerate(input_brain_list)
        ]

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

    node_list = [
        flirt_node(img, out_img, out_mat)
        for (img, out_img, out_mat) in zip(
            input_brain_list, output_img_list, output_mat_list
        )
    ]
    pool.map(lambda node: node.run(), node_list)

    if isinstance(thread_pool, int):
        pool.close()
        pool.join()

    return node_list


def template_creation_flirt(
    input_brain_list: list[str],
    input_skull_list: list[str],
    init_reg: Optional[list[pe.Node]] = None,
    avg_method: Literal["median", "mean", "std"] = "median",
    dof: Literal[12, 9, 7, 6] = 12,
    interp: Literal["trilinear", "nearestneighbour", "sinc", "spline"] = "trilinear",
    cost: Literal[
        "corratio", "mutualinfo", "normmi", "normcorr", "leastsq", "labeldiff", "bbr"
    ] = "corratio",
    mat_type: Literal["matrix", "ITK"] = "matrix",
    convergence_threshold: float | np.float64 = -1,
    max_iter: int = 5,
    thread_pool: int | Pool = 2,
    unique_id_list: Optional[list[str]] = None,
) -> tuple[str, str, list[str], list[str], list[str]]:
    """Create a temporary template from a list of images.

    Parameters
    ----------
    input_brain_list
        list of brain images paths
    input_skull_list
        list of skull images paths
    init_reg
        (default None so no initial registration performed)
        the output of the function register_img_list with another reference
        Reuter et al. 2012 (NeuroImage) section "Improved template estimation"
        doi:10.1016/j.neuroimage.2012.02.084 uses a ramdomly
        selected image from the input dataset
    avg_method
        function names from numpy library
    dof
        number of transform degrees of freedom (FLIRT)
    interp
        final interpolation method used in reslicing
    cost
        cost function
    mat_type
        The type of matrix used to represent the transformations
    convergence_threshold
        (numpy.finfo(np.float64).eps (default)) threshold for the convergence
        The threshold is how different from no transformation is the
        transformation matrix.
    max_iter
        Maximum number of iterations if transformation does not converge
    thread_pool : int or multiprocessing.dummy.Pool
        (default 2) number of threads. You can also provide a Pool so the
        node will be added to it to be run.
    unique_id_list : list of str
        list of unique IDs in data config

    Returns
    -------
    template : str
        path to the final template

    """
    # DEBUG to skip the longitudinal template generation which takes a lot of time.
    # return 'CECI_EST_UN_TEST'

    if not input_brain_list or not input_skull_list:
        msg = "ERROR create_temporary_template: image list is empty"
        raise ValueError(msg)

    warp_list = []

    # check if image basename_list are the same
    basename_list = [
        str(os.path.basename(img).split(".")[0]) for img in input_brain_list
    ]
    counter = Counter(basename_list)
    duplicated_basename_list = [i for i, j in counter.items() if j > 1]

    if (
        not duplicated_basename_list
    ):  # if duplicated_basename_list is empty, no duplicated basenames
        warp_list_filenames = [
            os.path.join(
                os.getcwd(),
                str(os.path.basename(img).split(".")[0]) + "_anat_to_template.mat",
            )
            for img in input_brain_list
        ]
    elif len(unique_id_list) == len(input_brain_list):
        warp_list_filenames = [
            os.path.join(
                os.getcwd(),
                str(os.path.basename(img).split(".")[0])
                + "_"
                + unique_id_list[i]
                + "_anat_to_template.mat",
            )
            for i, img in enumerate(input_brain_list)
        ]

    if isinstance(thread_pool, int):
        pool = ThreadPool(thread_pool)
    else:
        pool = thread_pool

    if convergence_threshold == -1:
        convergence_threshold = np.finfo(np.float64).eps

    if len(input_brain_list) == 1 or len(input_skull_list) == 1:
        IFLOGGER.warning(
            "input_brain_list or input_skull_list contains only 1 image, "
            "no need to calculate template"
        )
        warp_list.append(np.identity(4, dtype=float))  # return an identity matrix
        return (
            input_brain_list[0],
            input_skull_list[0],
            input_brain_list,
            input_skull_list,
            warp_list,
        )

    # Chris: I added this part because it is mentioned in the paper but I actually never used it
    # You could run a first register_img_list() with a selected image as starting point and
    # give the output to this function
    if init_reg is not None:
        if isinstance(init_reg, list):
            output_brain_list = [node.inputs.out_file for node in init_reg]
            mat_list = [node.inputs.out_matrix_file for node in init_reg]
            warp_list = mat_list
            # test if every transformation matrix has reached the convergence
            convergence_list = [
                template_convergence(mat, mat_type, convergence_threshold)
                for mat in mat_list
            ]
            converged = all(convergence_list)
        else:
            msg = "init_reg must be a list of FLIRT nipype nodes files"
            raise ValueError(msg)
    else:
        output_brain_list = input_brain_list
        output_skull_list = input_skull_list
        converged = False

    temporary_brain_template = os.path.join(
        os.getcwd(), "temporary_brain_template.nii.gz"
    )
    temporary_skull_template = os.path.join(
        os.getcwd(), "temporary_skull_template.nii.gz"
    )

    """ First is calculated an average image of the dataset to be the temporary template
    and the loop stops when this temporary template is close enough (with a transformation
    distance smaller than the threshold) to all the images of the precedent iteration.
    """
    iterator = 1
    iteration = 0
    if max_iter == -1:
        # make iteration < max_iter always True
        iterator = 0
        iteration = -2
    while not converged and iteration < max_iter:
        iteration += iterator
        temporary_brain_template, temporary_skull_template = create_temporary_template(
            input_brain_list=output_brain_list,
            input_skull_list=output_skull_list,
            output_brain_path=temporary_brain_template,
            output_skull_path=temporary_skull_template,
            avg_method=avg_method,
        )

        reg_list_node = register_img_list(
            input_brain_list=output_brain_list,
            ref_img=temporary_brain_template,
            dof=dof,
            interp=interp,
            cost=cost,
            unique_id_list=unique_id_list,
        )

        mat_list = [node.inputs.out_matrix_file for node in reg_list_node]

        # TODO clean code, refactor variables
        if len(warp_list) == 0:
            warp_list = mat_list

        for index, mat in enumerate(mat_list):
            cmd = (
                "flirt -in %s -ref %s -applyxfm -init %s -dof %s -interp %s -cost %s -out %s"
                % (
                    output_skull_list[index],
                    temporary_skull_template,
                    mat,
                    dof,
                    interp,
                    cost,
                    os.path.join(
                        os.getcwd(), os.path.basename(output_skull_list[index])
                    ),
                )
            )
            os.system(cmd)

            output_skull_list[index] = os.path.join(
                os.getcwd(), os.path.basename(output_skull_list[index])
            )

            # why inverse?
            cmd = "convert_xfm -omat %s -inverse %s" % (
                warp_list_filenames[index],
                warp_list[index],
            )
            os.system(cmd)

            warp_list[index] = warp_list_filenames[index]

        output_brain_list = [node.inputs.out_file for node in reg_list_node]

        # test if every transformation matrix has reached the convergence
        convergence_list = [
            template_convergence(mat, mat_type, convergence_threshold)
            for mat in mat_list
        ]
        converged = all(convergence_list)

    if isinstance(thread_pool, int):
        pool.close()
        pool.join()

    brain_template = temporary_brain_template
    skull_template = temporary_skull_template

    # register T1 to longitudinal template space
    reg_list_node = register_img_list(
        input_brain_list,
        ref_img=temporary_brain_template,
        dof=dof,
        interp=interp,
        cost=cost,
        unique_id_list=unique_id_list,
    )

    warp_list = [node.inputs.out_matrix_file for node in reg_list_node]

    return (
        brain_template,
        skull_template,
        output_brain_list,
        output_skull_list,
        warp_list,
    )


def subject_specific_template(
    workflow_name="subject_specific_template", method="flirt"
):
    """Create a subject specific template from a list of images.

    Parameters
    ----------
    workflow_name : str

    method : str

    Returns
    -------
    template_gen_node : ~nipype.pipeline.engine.Node
    """
    imports = [
        "import os",
        "import warnings",
        "import numpy as np",
        "from collections import Counter",
        "from multiprocessing.dummy import Pool as ThreadPool",
        "from nipype.interfaces.fsl import ConvertXFM",
        "from CPAC.longitudinal_pipeline.longitudinal_preproc import ("
        "   create_temporary_template,"
        "   register_img_list,"
        "   template_convergence"
        ")",
    ]
    if method == "flirt":
        template_gen_node = pe.Node(
            Function(
                input_names=[
                    "input_brain_list",
                    "input_skull_list",
                    "init_reg",
                    "avg_method",
                    "dof",
                    "interp",
                    "cost",
                    "mat_type",
                    "convergence_threshold",
                    "max_iter",
                    "thread_pool",
                    "unique_id_list",
                ],
                output_names=[
                    "brain_template",
                    "skull_template",
                    "output_brain_list",
                    "output_skull_list",
                    "warp_list",
                ],
                imports=imports,
                function=template_creation_flirt,
            ),
            name=workflow_name,
        )
    else:
        raise ValueError(str(method) + "this method has not yet been implemented")

    return template_gen_node
