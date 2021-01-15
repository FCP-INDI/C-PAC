import os
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from .utils import *
from CPAC.vmhc import *
from nipype.interfaces.afni import preprocess
from CPAC.registration import create_wf_calculate_ants_warp, \
                              output_func_to_standard
from CPAC.image_utils import spatial_smooth


def create_vmhc(workflow, num_strat, strat, pipeline_config_object,
        func_key='functional_nuisance_residuals', output_name='vmhc'):
    """
    Compute the map of brain functional homotopy, the high degree of synchrony
    in spontaneous activity between geometrically corresponding interhemispheric (i.e., homotopic) regions.



    Parameters
    ----------

    None

    Returns
    -------

    vmhc_workflow : Workflow

        Voxel Mirrored Homotopic Connectivity Analysis Workflow

    strat : Strategy

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/vmhc/vmhc.py>`_

    Workflow Inputs::

        inputspec.brain : string (existing nifti file)
            Anatomical image(without skull)

        inputspec.symmetric_brain : string (existing nifti file)
            MNI152_T1_2mm_symmetric_brain.nii.gz

        inputspec.rest_res_filt : string (existing nifti file)
            Band passed Image with nuisance signal regressed out(and optionally scrubbed). Recommended
            bandpass filter (0.001,0.1) )

        inputspec.reorient : string (existing nifti file)
            RPI oriented anatomical data

        inputspec.example_func2highres_mat : string (existing affine transformation .mat file)
            Specifies an affine transform that should be applied to the example_func before non linear warping

        inputspec.standard_for_func: string (existing nifti file)
            MNI152_T1_standard_resolution_brain.nii.gz

        inputspec.symmetric_skull : string (existing nifti file)
            MNI152_T1_2mm_symmetric.nii.gz

        inputspec.twomm_brain_mask_dil : string (existing nifti file)
            MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz

        inputspec.config_file_twomm_symmetric : string (existing .cnf file)
            T1_2_MNI152_2mm_symmetric.cnf

        inputspec.rest_mask : string (existing nifti file)
            A mask functional volume(derived by dilation from motion corrected functional volume)

        fwhm_input.fwhm : list (float)
            For spatial smoothing the Z-transformed correlations in MNI space.
            Generally the value of this parameter is 1.5 or 2 times the voxel size of the input Image.

        inputspec.mean_functional : string (existing nifti file)
            The mean functional image for use in the func-to-anat registration matrix conversion
            to ITK (ANTS) format, if the user selects to use ANTS.


    Workflow Outputs::

        outputspec.highres2symmstandard : string (nifti file)
            Linear registration of T1 image to symmetric standard image

        outputspec.highres2symmstandard_mat : string (affine transformation .mat file)
            An affine transformation .mat file from linear registration and used in non linear registration

        outputspec.highres2symmstandard_warp : string (nifti file)
            warp file from Non Linear registration of T1 to symmetrical standard brain

        outputspec.fnirt_highres2symmstandard : string (nifti file)
            Non Linear registration of T1 to symmetrical standard brain

        outputspec.highres2symmstandard_jac : string (nifti file)
            jacobian determinant image from Non Linear registration of T1 to symmetrical standard brain

        outputspec.rest_res_2symmstandard : string (nifti file)
            nonlinear registration (func to standard) image

        outputspec.VMHC_FWHM_img : string (nifti file)
            pearson correlation between res2standard and flipped res2standard

        outputspec.VMHC_Z_FWHM_img : string (nifti file)
            Fisher Z transform map

        outputspec.VMHC_Z_stat_FWHM_img : string (nifti file)
            Z statistic map

    Order of commands:

    - Perform linear registration of Anatomical brain in T1 space to symmetric standard space. For details
    see `flirt <http://www.fmrib.ox.ac.uk/fsl/flirt/index.html>`_::

        flirt
        -ref MNI152_T1_2mm_symmetric_brain.nii.gz
        -in mprage_brain.nii.gz
        -out highres2symmstandard.nii.gz
        -omat highres2symmstandard.mat
        -cost corratio
        -searchcost corratio
        -dof 12
        -interp trilinear

    - Perform nonlinear registration (higres to standard) to symmetric standard brain. For details
    see `fnirt <http://fsl.fmrib.ox.ac.uk/fsl/fnirt/>`_::

        fnirt
        --in=head.nii.gz
        --aff=highres2symmstandard.mat
        --cout=highres2symmstandard_warp.nii.gz
        --iout=fnirt_highres2symmstandard.nii.gz
        --jout=highres2symmstandard_jac.nii.gz
        --config=T1_2_MNI152_2mm_symmetric.cnf
        --ref=MNI152_T1_2mm_symmetric.nii.gz
        --refmask=MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz
        --warpres=10,10,10

    - Perform spatial smoothing on the input functional image(inputspec.rest_res_filt).  For details
    see `PrinciplesSmoothing <http://imaging.mrc-cbu.cam.ac.uk/imaging/PrinciplesSmoothing>`_
    `fslmaths <http://www.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro/index.htm>`_::

        fslmaths rest_res_filt.nii.gz
        -kernel gauss FWHM/ sqrt(8*ln(2))
        -fmean -mas rest_mask.nii.gz
        rest_res_filt_FWHM.nii.gz

    - Apply nonlinear registration (func to standard). For details see
    `applywarp <http://www.fmrib.ox.ac.uk/fsl/fnirt/warp_utils.html#applywarp>`_::

        applywarp
        --ref=MNI152_T1_2mm_symmetric.nii.gz
        --in=rest_res_filt_FWHM.nii.gz
        --out=rest_res_2symmstandard.nii.gz
        --warp=highres2symmstandard_warp.nii.gz
        --premat=example_func2highres.mat


    - Copy and L/R swap the output of applywarp command (rest_res_2symmstandard.nii.gz). For details
    see  `fslswapdim <http://fsl.fmrib.ox.ac.uk/fsl/fsl4.0/avwutils/index.html>`_::

        fslswapdim
        rest_res_2symmstandard.nii.gz
        -x y z
        tmp_LRflipped.nii.gz


    - Calculate pearson correlation between rest_res_2symmstandard.nii.gz and flipped
    rest_res_2symmstandard.nii.gz(tmp_LRflipped.nii.gz). For details see
    `3dTcorrelate <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTcorrelate.html>`_::

        3dTcorrelate
        -pearson
        -polort -1
        -prefix VMHC_FWHM.nii.gz
        rest_res_2symmstandard.nii.gz
        tmp_LRflipped.nii.gz

    .. exec::
        import nipype.pipeline.engine as pe
        from CPAC.utils.interfaces.function import Function
        from CPAC.utils.test_mocks import configuration_strategy_mock
        from CPAC.vmhc import create_vmhc
        pipeline_config_object, strat = configuration_strategy_mock(
            method="ANTS")
        pipeline_config_object.update('smoothing_method', ['AFNI'])
        pipeline_config_object.update('maxCoresPerParticipant', '1')
        strat.set_leaf_properties(
            pe.Node(Function([]), 'vhmc_leaf_node'), 'vhmc_outfile')
        wf = create_vmhc(
            pe.Workflow('vmhc'), 0, strat, pipeline_config_object)[0]
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/vmhc.dot'
        )

    Workflow:

    .. image:: ../../images/generated/vmhc.png
        :width: 500

    Workflow Detailed:

    .. image:: ../../images/generated/vmhc_detailed.png
        :width: 500


    References
    ----------

    .. [1] Zuo, X.-N., Kelly, C., Di Martino, A., Mennes, M., Margulies, D. S., Bangaru, S.,
           Grzadzinski, R., et al. (2010). Growing together and growing apart: regional and
           sex differences in the lifespan developmental trajectories of functional homotopy.
           The Journal of neuroscience : the official journal of the Society for Neuroscience,
           30(45), 15034-43. doi:10.1523/JNEUROSCI.2612-10.2010


    Examples
    --------

    >>> vmhc_w = create_vmhc()
    >>> vmhc_w.inputs.inputspec.symmetric_brain = 'MNI152_T1_2mm_symmetric_brain.nii.gz'
    >>> vmhc_w.inputs.inputspec.symmetric_skull = 'MNI152_T1_2mm_symmetric.nii.gz'
    >>> vmhc_w.inputs.inputspec.twomm_brain_mask_dil = 'MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz'
    >>> vmhc_w.inputs.inputspec.config_file_twomm = 'T1_2_MNI152_2mm_symmetric.cnf'
    >>> vmhc_w.inputs.inputspec.standard_for_func= 'MNI152_T1_2mm.nii.gz'
    >>> vmhc_w.inputs.fwhm_input.fwhm = [4.5, 6]
    >>> vmhc_w.get_node('fwhm_input').iterables = ('fwhm', [4.5, 6])
    >>> vmhc_w.inputs.inputspec.rest_res = os.path.abspath('/home/data/s1001/rest_res_filt.nii.gz')
    >>> vmhc_w.inputs.inputspec.reorient = os.path.abspath('/home/data/s1001/anat/mprage_RPI.nii.gz')
    >>> vmhc_w.inputs.inputspec.brain = os.path.abspath('/home/data/s1001/anat/mprage_brain.nii.gz')
    >>> vmhc_w.inputs.inputspec.example_func2highres_mat = os.path.abspath('/home/data/s1001/func2highres.mat')
    >>> vmhc_w.inputs.inputspec.rest_mask = os.path.abspath('/home/data/s1001/func/original/rest_mask.nii.gz')
    >>> vmhc_w.run() # doctest: +SKIP

    """

    nodes = strat.get_nodes_names()

    if not isinstance(func_key, str):
        raise ValueError('func_key should be a string, not a {0}'.format(type(func_key)))

    # we begin by smoothing the input file, which should be the current leaf node
    smooth_key = '{0}_smooth'.format(func_key)
    if smooth_key not in strat:
        spatial_smooth(workflow, 'leaf', 'functional_brain_mask', smooth_key,
            strat, num_strat, pipeline_config_object, input_image_type='func_4d')

    # next write it to symmetric MNI space
    func_symm_mni_key = 'func_preproc_symm_mni'
    if func_symm_mni_key not in strat:
        output_func_to_standard(workflow, smooth_key, 'template_skull_for_func_preproc',
            func_symm_mni_key, strat, num_strat, pipeline_config_object, input_image_type='func_4d')

    # write out a swapped version of the file
    # copy and L/R swap file
    copy_and_L_R_swap = pe.Node(interface=fsl.SwapDimensions(),
                      name='copy_and_L_R_swap_{0}'.format(num_strat))

    copy_and_L_R_swap.inputs.new_dims = ('-x', 'y', 'z')

    func_node, func_file = strat[func_symm_mni_key]
    workflow.connect(func_node, func_file,
                     copy_and_L_R_swap, 'in_file')

    # calculate correlation between original and swapped images
    pearson_correlation = pe.Node(interface=preprocess.TCorrelate(),
                      name='pearson_correlation_{0}'.format(num_strat))

    pearson_correlation.inputs.pearson = True
    pearson_correlation.inputs.polort = -1
    pearson_correlation.inputs.outputtype = 'NIFTI_GZ'

    workflow.connect(func_node, func_file,
                     pearson_correlation, 'xset')

    workflow.connect(copy_and_L_R_swap, 'out_file',
                     pearson_correlation, 'yset')

    # add the outputs to the resource pool
    strat.update_resource_pool({
        'vmhc_raw_score': (pearson_correlation, 'out_file')
    })

    strat.append_name(copy_and_L_R_swap.name)
    strat.append_name(pearson_correlation.name)

    return workflow, strat
