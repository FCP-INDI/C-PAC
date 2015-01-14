import sys
import os
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from utils import *
from CPAC.vmhc import *
from nipype.interfaces.afni import preprocess
from CPAC.registration import create_wf_calculate_ants_warp, \
                              create_wf_c3d_fsl_to_itk, \
                              create_wf_collect_transforms, \
                              create_wf_apply_ants_warp

def create_vmhc(use_ants, name='vmhc_workflow'):

    """
    Compute the map of brain functional homotopy, the high degree of synchrony in spontaneous activity between geometrically corresponding interhemispheric (i.e., homotopic) regions.



    Parameters
    ----------

    None

    Returns
    -------

    vmhc_workflow : workflow

        Voxel Mirrored Homotopic Connectivity Analysis Workflow



    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/vmhc/vmhc.py>`_ 

    Workflow Inputs::

        inputspec.brain : string (existing nifti file)
            Anatomical image(without skull)

        inputspec.symmetric_brain : string (existing nifti file)
            MNI152_T1_2mm_symmetric_brain.nii.gz
 
        inputspec.rest_res_filt : string (existing nifti file)
            Band passed Image with nuisance signal regressed out(and optionally scrubbed). Recommended bandpass filter (0.001,0.1) )

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

    - Perform linear registration of Anatomical brain in T1 space to symmetric standard space. For details see `flirt <http://www.fmrib.ox.ac.uk/fsl/flirt/index.html>`_::

        flirt
        -ref MNI152_T1_2mm_symmetric_brain.nii.gz
        -in mprage_brain.nii.gz
        -out highres2symmstandard.nii.gz
        -omat highres2symmstandard.mat
        -cost corratio
        -searchcost corratio
        -dof 12
        -interp trilinear    
        
    - Perform nonlinear registration (higres to standard) to symmetric standard brain. For details see `fnirt <http://fsl.fmrib.ox.ac.uk/fsl/fnirt/>`_::
    
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

    - Perform spatial smoothing on the input functional image(inputspec.rest_res_filt).  For details see `PrinciplesSmoothing <http://imaging.mrc-cbu.cam.ac.uk/imaging/PrinciplesSmoothing>`_ `fslmaths <http://www.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro/index.htm>`_::

        fslmaths rest_res_filt.nii.gz
        -kernel gauss FWHM/ sqrt(8-ln(2))
        -fmean -mas rest_mask.nii.gz
        rest_res_filt_FWHM.nii.gz
        
    - Apply nonlinear registration (func to standard). For details see  `applywarp <http://www.fmrib.ox.ac.uk/fsl/fnirt/warp_utils.html#applywarp>`_::
        
        applywarp
        --ref=MNI152_T1_2mm_symmetric.nii.gz
        --in=rest_res_filt_FWHM.nii.gz
        --out=rest_res_2symmstandard.nii.gz
        --warp=highres2symmstandard_warp.nii.gz
        --premat=example_func2highres.mat
        
        
    - Copy and L/R swap the output of applywarp command (rest_res_2symmstandard.nii.gz). For details see  `fslswapdim <http://fsl.fmrib.ox.ac.uk/fsl/fsl4.0/avwutils/index.html>`_::

        fslswapdim
        rest_res_2symmstandard.nii.gz
        -x y z
        tmp_LRflipped.nii.gz


    - Calculate pearson correlation between rest_res_2symmstandard.nii.gz and flipped rest_res_2symmstandard.nii.gz(tmp_LRflipped.nii.gz). For details see  `3dTcorrelate <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTcorrelate.html>`_::
        
        3dTcorrelate
        -pearson
        -polort -1
        -prefix VMHC_FWHM.nii.gz
        rest_res_2symmstandard.nii.gz
        tmp_LRflipped.nii.gz
    
    
    - Fisher Z Transform the correlation. For details see `3dcalc <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_::
        
        3dcalc
        -a VMHC_FWHM.nii.gz
        -expr 'log((a+1)/(1-a))/2'
        -prefix VMHC_FWHM_Z.nii.gz
    
        
    - Calculate the number of volumes(nvols) in flipped rest_res_2symmstandard.nii.gz(tmp_LRflipped.nii.gz) ::
        
        -Use Nibabel to do this
        
        
    - Compute the Z statistic map ::
        
        3dcalc
        -a VMHC_FWHM_Z.nii.gz
        -expr 'a*sqrt('${nvols}'-3)'
        -prefix VMHC_FWHM_Z_stat.nii.gz
    
    
    Workflow:
    
    .. image:: ../images/vmhc_graph.dot.png
        :width: 500 
    
    Workflow Detailed:
    
    .. image:: ../images/vmhc_detailed_graph.dot.png
        :width: 500 
    

    References
    ----------
    
    .. [1] Zuo, X.-N., Kelly, C., Di Martino, A., Mennes, M., Margulies, D. S., Bangaru, S., Grzadzinski, R., et al. (2010). Growing together and growing apart: regional and sex differences in the lifespan developmental trajectories of functional homotopy. The Journal of neuroscience : the official journal of the Society for Neuroscience, 30(45), 15034-43. doi:10.1523/JNEUROSCI.2612-10.2010


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
    >>> vmhc_w.inputs.inputspec.rest_res = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/rest_res_filt.nii.gz')
    >>> vmhc_w.inputs.inputspec.reorient = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/anat/mprage_RPI.nii.gz')
    >>> vmhc_w.inputs.inputspec.brain = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/anat/mprage_brain.nii.gz')
    >>> vmhc_w.inputs.inputspec.example_func2highres_mat = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/reg/example_func2highres.mat')
    >>> vmhc_w.inputs.inputspec.rest_mask = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/rest_mask.nii.gz')
    >>> vmhc_w.run() # doctest: +SKIP

    """

    vmhc = pe.Workflow(name=name)
    inputNode = pe.Node(util.IdentityInterface(fields=['rest_res',
                                                'example_func2highres_mat',
                                                'rest_mask',
                                                'standard_for_func',
                                                'mean_functional',
                                                'brain',
                                                'fnirt_nonlinear_warp',
                                                'ants_symm_initial_xfm',
                                                'ants_symm_rigid_xfm',
                                                'ants_symm_affine_xfm',
                                                'ants_symm_warp_field']),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=['rest_res_2symmstandard',
                                                'VMHC_FWHM_img',
                                                'VMHC_Z_FWHM_img',
                                                'VMHC_Z_stat_FWHM_img'
                                                ]),
                        name='outputspec')


    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')


    if use_ants == False:

        ## Apply nonlinear registration (func to standard)
        nonlinear_func_to_standard = pe.Node(interface=fsl.ApplyWarp(),
                          name='nonlinear_func_to_standard')

    elif use_ants == True:

        # ANTS warp image etc.
        fsl_to_itk_vmhc = create_wf_c3d_fsl_to_itk(0, name='fsl_to_itk_vmhc')

        collect_transforms_vmhc = create_wf_collect_transforms(0, name='collect_transforms_vmhc')

        apply_ants_xfm_vmhc = create_wf_apply_ants_warp(0,name='apply_ants_xfm_vmhc')

        # this has to be 3 instead of default 0 because it is a 4D file
        apply_ants_xfm_vmhc.inputs.inputspec.input_image_type = 3



    ## copy and L/R swap file
    copy_and_L_R_swap = pe.Node(interface=fsl.SwapDimensions(),
                      name='copy_and_L_R_swap')
    copy_and_L_R_swap.inputs.new_dims = ('-x', 'y', 'z')

    ## caculate vmhc
    pearson_correlation = pe.Node(interface=preprocess.TCorrelate(),
                      name='pearson_correlation')
    pearson_correlation.inputs.pearson = True
    pearson_correlation.inputs.polort = -1
    pearson_correlation.inputs.outputtype = 'NIFTI_GZ'

    z_trans = pe.Node(interface=preprocess.Calc(),
                         name='z_trans')
    z_trans.inputs.expr = 'log((1+a)/(1-a))/2'
    z_trans.inputs.outputtype = 'NIFTI_GZ'

    z_stat = pe.Node(interface=preprocess.Calc(),
                        name='z_stat')
    z_stat.inputs.outputtype = 'NIFTI_GZ'

    NVOLS = pe.Node(util.Function(input_names=['in_files'],
                                  output_names=['nvols'],
                    function=get_img_nvols),
                    name='NVOLS')

    generateEXP = pe.Node(util.Function(input_names=['nvols'],
                                        output_names=['expr'],
                          function=get_operand_expression),
                          name='generateEXP')


    smooth = pe.Node(interface=fsl.MultiImageMaths(),
                        name='smooth')


    if use_ants == False:

        vmhc.connect(inputNode, 'rest_res',
                     smooth, 'in_file')
        vmhc.connect(inputnode_fwhm, ('fwhm', set_gauss),
                     smooth, 'op_string')
        vmhc.connect(inputNode, 'rest_mask',
                     smooth, 'operand_files')
        vmhc.connect(smooth, 'out_file',
                     nonlinear_func_to_standard, 'in_file')
        vmhc.connect(inputNode, 'standard_for_func',
                     nonlinear_func_to_standard, 'ref_file')
        vmhc.connect(inputNode, 'fnirt_nonlinear_warp',
                     nonlinear_func_to_standard, 'field_file')
        ## func->anat matrix (bbreg)
        vmhc.connect(inputNode, 'example_func2highres_mat',
                     nonlinear_func_to_standard, 'premat')
        vmhc.connect(nonlinear_func_to_standard, 'out_file',
                     copy_and_L_R_swap, 'in_file')
        vmhc.connect(nonlinear_func_to_standard, 'out_file',
                     pearson_correlation, 'xset')

    elif use_ants == True:

        # connections for ANTS stuff

        # functional apply warp stuff

        vmhc.connect(inputNode, 'rest_res',
                     smooth, 'in_file')
        vmhc.connect(inputnode_fwhm, ('fwhm', set_gauss),
                     smooth, 'op_string')
        vmhc.connect(inputNode, 'rest_mask',
                     smooth, 'operand_files')

        vmhc.connect(smooth, 'out_file',
                     apply_ants_xfm_vmhc, 'inputspec.input_image')

        vmhc.connect(inputNode, 'ants_symm_initial_xfm',
                     collect_transforms_vmhc, 'inputspec.linear_initial')

        vmhc.connect(inputNode, 'ants_symm_rigid_xfm',
                     collect_transforms_vmhc, 'inputspec.linear_rigid')

        vmhc.connect(inputNode, 'ants_symm_affine_xfm',
                     collect_transforms_vmhc, 'inputspec.linear_affine')

        vmhc.connect(inputNode, 'ants_symm_warp_field',
                     collect_transforms_vmhc, 'inputspec.warp_file')

        ## func->anat matrix (bbreg)
        vmhc.connect(inputNode, 'example_func2highres_mat',
                     fsl_to_itk_vmhc, 'inputspec.affine_file')

        vmhc.connect(inputNode, 'brain', fsl_to_itk_vmhc,
                     'inputspec.reference_file')

        vmhc.connect(inputNode, 'mean_functional', fsl_to_itk_vmhc,
                     'inputspec.source_file')

        vmhc.connect(fsl_to_itk_vmhc, 'outputspec.itk_transform', 
                     collect_transforms_vmhc, 'inputspec.fsl_to_itk_affine')

        vmhc.connect(inputNode, 'standard_for_func',
                     apply_ants_xfm_vmhc, 'inputspec.reference_image')

        vmhc.connect(collect_transforms_vmhc, \
                     'outputspec.transformation_series', \
                     apply_ants_xfm_vmhc, 'inputspec.transforms')

        vmhc.connect(apply_ants_xfm_vmhc, 'outputspec.output_image',
                     copy_and_L_R_swap, 'in_file')

        vmhc.connect(apply_ants_xfm_vmhc, 'outputspec.output_image',
                     pearson_correlation, 'xset')


    vmhc.connect(copy_and_L_R_swap, 'out_file',
                 pearson_correlation, 'yset')
    vmhc.connect(pearson_correlation, 'out_file',
                 z_trans, 'in_file_a')
    vmhc.connect(copy_and_L_R_swap, 'out_file',
                 NVOLS, 'in_files')
    vmhc.connect(NVOLS, 'nvols',
                 generateEXP, 'nvols')
    vmhc.connect(z_trans, 'out_file',
                 z_stat, 'in_file_a')
    vmhc.connect(generateEXP, 'expr',
                 z_stat, 'expr')

    if use_ants == False:

        vmhc.connect(nonlinear_func_to_standard, 'out_file',
                     outputNode, 'rest_res_2symmstandard')

    elif use_ants == True:

        # ANTS warp outputs to outputnode

        vmhc.connect(apply_ants_xfm_vmhc, 'outputspec.output_image',
                     outputNode, 'rest_res_2symmstandard')


    vmhc.connect(pearson_correlation, 'out_file',
                 outputNode, 'VMHC_FWHM_img')
    vmhc.connect(z_trans, 'out_file',
                 outputNode, 'VMHC_Z_FWHM_img')
    vmhc.connect(z_stat, 'out_file',
                 outputNode, 'VMHC_Z_stat_FWHM_img')


    return vmhc
