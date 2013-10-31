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
from CPAC.interfaces.afni import preprocess
from CPAC.registration import create_ants_nonlinear_xfm, create_apply_ants_xfm

def create_vmhc(use_ants):

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

        inputspec.brain_symmetric : string (existing nifti file)
            MNI152_T1_2mm_brain_symmetric.nii.gz
 
        inputspec.rest_res_filt : string (existing nifti file)
            Band passed Image with nuisance signal regressed out(and optionally scrubbed). Recommended bandpass filter (0.001,0.1) )

        inputspec.reorient : string (existing nifti file)
            RPI oriented anatomical data

        inputspec.example_func2highres_mat : string (existing affine transformation .mat file)
            Specifies an affine transform that should be applied to the example_func before non linear warping

        inputspec.standard : string (existing nifti file)
            MNI152_T1_standard_resolution_brain.nii.gz

        inputspec.symm_standard : string (existing nifti file)
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
        -ref MNI152_T1_2mm_brain_symmetric.nii.gz
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
    >>> vmhc_w.inputs.inputspec.brain_symmetric = 'MNI152_T1_2mm_brain_symmetric.nii.gz'
    >>> vmhc_w.inputs.inputspec.symm_standard = 'MNI152_T1_2mm_symmetric.nii.gz'
    >>> vmhc_w.inputs.inputspec.twomm_brain_mask_dil = 'MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz'
    >>> vmhc_w.inputs.inputspec.config_file_twomm = 'T1_2_MNI152_2mm_symmetric.cnf'
    >>> vmhc_w.inputs.inputspec.standard = 'MNI152_T1_2mm.nii.gz'
    >>> vmhc_w.inputs.fwhm_input.fwhm = [4.5, 6]
    >>> vmhc_w.get_node('fwhm_input').iterables = ('fwhm', [4.5, 6])
    >>> vmhc_w.inputs.inputspec.rest_res = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/rest_res_filt.nii.gz')
    >>> vmhc_w.inputs.inputspec.reorient = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/anat/mprage_RPI.nii.gz')
    >>> vmhc_w.inputs.inputspec.brain = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/anat/mprage_brain.nii.gz')
    >>> vmhc_w.inputs.inputspec.example_func2highres_mat = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/reg/example_func2highres.mat')
    >>> vmhc_w.inputs.inputspec.rest_mask = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein/subjects/s1001/func/original/rest_mask.nii.gz')
    >>> vmhc_w.run() # doctest: +SKIP

    """

    vmhc = pe.Workflow(name='vmhc_workflow')
    inputNode = pe.Node(util.IdentityInterface(fields=['brain',
                                                'brain_symmetric',
                                                'rest_res',
                                                'reorient',
                                                'example_func2highres_mat',
                                                'symm_standard',
                                                'twomm_brain_mask_dil',
                                                'config_file_twomm',
                                                'rest_mask',
                                                'standard',
                                                'mean_functional']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['highres2symmstandard',
                                                'highres2symmstandard_mat',
                                                'highres2symmstandard_warp',
                                                'fnirt_highres2symmstandard',
                                                'highres2symmstandard_jac',
                                                'rest_res_2symmstandard',
                                                'VMHC_FWHM_img',
                                                'VMHC_Z_FWHM_img',
                                                'VMHC_Z_stat_FWHM_img'
                                                ]),
                        name='outputspec')


    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')


    if use_ants == False:

        ## Linear registration of T1 --> symmetric standard
        linear_T1_to_symmetric_standard = pe.Node(interface=fsl.FLIRT(),
                        name='linear_T1_to_symmetric_standard')
        linear_T1_to_symmetric_standard.inputs.cost = 'corratio'
        linear_T1_to_symmetric_standard.inputs.cost_func = 'corratio'
        linear_T1_to_symmetric_standard.inputs.dof = 12
        linear_T1_to_symmetric_standard.inputs.interp = 'trilinear'

        ## Perform nonlinear registration
        ##(higres to standard) to symmetric standard brain
        nonlinear_highres_to_symmetric_standard = pe.Node(interface=fsl.FNIRT(),
                      name='nonlinear_highres_to_symmetric_standard')
        nonlinear_highres_to_symmetric_standard.inputs.fieldcoeff_file = True
        nonlinear_highres_to_symmetric_standard.inputs.jacobian_file = True
        nonlinear_highres_to_symmetric_standard.inputs.warp_resolution = (10, 10, 10)

        # needs new inputs. needs input from resources for the field coeff of the template->symmetric.
        # and needs the field coeff of the anatomical-to-template registration

        ## Apply nonlinear registration (func to standard)
        nonlinear_func_to_standard = pe.Node(interface=fsl.ApplyWarp(),
                          name='nonlinear_func_to_standard')

    elif use_ants == True:

        # ANTS warp image etc.

        calculate_ants_xfm_vmhc = create_ants_nonlinear_xfm(name='calculate_ants_xfm_vmhc')

        apply_ants_xfm_vmhc = create_apply_ants_xfm(4,0,name='apply_ants_xfm_vmhc')


    ## copy and L/R swap file
    copy_and_L_R_swap = pe.Node(interface=fsl.SwapDimensions(),
                      name='copy_and_L_R_swap')
    copy_and_L_R_swap.inputs.new_dims = ('-x', 'y', 'z')

    ## caculate vmhc
    pearson_correlation = pe.Node(interface=preprocess.ThreedTcorrelate(),
                      name='pearson_correlation')
    pearson_correlation.inputs.pearson = True
    pearson_correlation.inputs.polort = -1
    pearson_correlation.inputs.outputtype = 'NIFTI_GZ'

    z_trans = pe.Node(interface=preprocess.Threedcalc(),
                         name='z_trans')
    z_trans.inputs.expr = '\'log((1+a)/(1-a))/2\''
    z_trans.inputs.outputtype = 'NIFTI_GZ'
    z_stat = pe.Node(interface=preprocess.Threedcalc(),
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

        vmhc.connect(inputNode, 'brain',
                     linear_T1_to_symmetric_standard, 'in_file')
        vmhc.connect(inputNode, 'brain_symmetric',
                     linear_T1_to_symmetric_standard, 'reference')
        vmhc.connect(inputNode, 'reorient',
                     nonlinear_highres_to_symmetric_standard, 'in_file')
        vmhc.connect(linear_T1_to_symmetric_standard, 'out_matrix_file',
                     nonlinear_highres_to_symmetric_standard, 'affine_file')
        vmhc.connect(inputNode, 'symm_standard',
                     nonlinear_highres_to_symmetric_standard, 'ref_file')
        vmhc.connect(inputNode, 'twomm_brain_mask_dil',
                     nonlinear_highres_to_symmetric_standard, 'refmask_file')
        vmhc.connect(inputNode, 'config_file_twomm',
                     nonlinear_highres_to_symmetric_standard, 'config_file')


        vmhc.connect(inputNode, 'rest_res',
                     smooth, 'in_file')
        vmhc.connect(inputnode_fwhm, ('fwhm', set_gauss),
                     smooth, 'op_string')
        vmhc.connect(inputNode, 'rest_mask',
                     smooth, 'operand_files')
        vmhc.connect(smooth, 'out_file',
                     nonlinear_func_to_standard, 'in_file')
        vmhc.connect(inputNode, 'symm_standard',
                     nonlinear_func_to_standard, 'ref_file')
        vmhc.connect(nonlinear_highres_to_symmetric_standard, 'fieldcoeff_file',
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

        # registration calculation stuff -- might go out the window
        vmhc.connect(inputNode, 'brain',
                     calculate_ants_xfm_vmhc, 'inputspec.anatomical_brain')
        vmhc.connect(inputNode, 'brain_symmetric',
                     calculate_ants_xfm_vmhc, 'inputspec.reference_brain')

        # functional apply warp stuff
        vmhc.connect(inputNode, 'rest_res',
                     smooth, 'in_file')
        vmhc.connect(inputnode_fwhm, ('fwhm', set_gauss),
                     smooth, 'op_string')
        vmhc.connect(inputNode, 'rest_mask',
                     smooth, 'operand_files')
        vmhc.connect(smooth, 'out_file',
                     apply_ants_xfm_vmhc, 'inputspec.in_file')
        vmhc.connect(inputNode, 'brain',
                     apply_ants_xfm_vmhc, 'inputspec.conversion_reference')
        vmhc.connect(inputNode, 'mean_functional',
                     apply_ants_xfm_vmhc, 'inputspec.conversion_source')
        vmhc.connect(inputNode, 'symm_standard',
                     apply_ants_xfm_vmhc, 'inputspec.warp_reference')
        vmhc.connect(calculate_ants_xfm_vmhc, 'outputspec.affine_transformation',
                     apply_ants_xfm_vmhc, 'inputspec.ants_affine')
        vmhc.connect(calculate_ants_xfm_vmhc, 'outputspec.warp_field',
                     apply_ants_xfm_vmhc, 'inputspec.nonlinear_field')
        ## func->anat matrix (bbreg)
        vmhc.connect(inputNode, 'example_func2highres_mat',
                     apply_ants_xfm_vmhc, 'inputspec.func_anat_affine')
        vmhc.connect(apply_ants_xfm_vmhc, 'outputspec.out_file',
                     copy_and_L_R_swap, 'in_file')
        vmhc.connect(apply_ants_xfm_vmhc, 'outputspec.out_file',
                     pearson_correlation, 'xset')



    vmhc.connect(copy_and_L_R_swap, 'out_file',
                 pearson_correlation, 'yset')
    vmhc.connect(pearson_correlation, 'out_file',
                 z_trans, 'infile_a')
    vmhc.connect(copy_and_L_R_swap, 'out_file',
                 NVOLS, 'in_files')
    vmhc.connect(NVOLS, 'nvols',
                 generateEXP, 'nvols')
    vmhc.connect(z_trans, 'out_file',
                 z_stat, 'infile_a')
    vmhc.connect(generateEXP, 'expr',
                 z_stat, 'expr')

    if use_ants == False:

        vmhc.connect(linear_T1_to_symmetric_standard, 'out_file',
                     outputNode, 'highres2symmstandard')
        vmhc.connect(linear_T1_to_symmetric_standard, 'out_matrix_file',
                     outputNode, 'highres2symmstandard_mat')
        vmhc.connect(nonlinear_highres_to_symmetric_standard, 'jacobian_file',
                     outputNode, 'highres2symmstandard_jac')
        vmhc.connect(nonlinear_highres_to_symmetric_standard, 'fieldcoeff_file',
                     outputNode, 'highres2symmstandard_warp')
        vmhc.connect(nonlinear_highres_to_symmetric_standard, 'warped_file',
                     outputNode, 'fnirt_highres2symmstandard')
        vmhc.connect(nonlinear_func_to_standard, 'out_file',
                     outputNode, 'rest_res_2symmstandard')

    elif use_ants == True:

        # ANTS warp outputs to outputnode

        vmhc.connect(calculate_ants_xfm_vmhc, 'outputspec.affine_transformation',
                     outputNode, 'highres2symmstandard_mat')
        vmhc.connect(calculate_ants_xfm_vmhc, 'outputspec.warp_field',
                     outputNode, 'highres2symmstandard_warp')
        vmhc.connect(calculate_ants_xfm_vmhc, 'outputspec.output_brain',
                     outputNode, 'fnirt_highres2symmstandard')
        vmhc.connect(apply_ants_xfm_vmhc, 'outputspec.out_file',
                     outputNode, 'rest_res_2symmstandard')


    vmhc.connect(pearson_correlation, 'out_file',
                 outputNode, 'VMHC_FWHM_img')
    vmhc.connect(z_trans, 'out_file',
                 outputNode, 'VMHC_Z_FWHM_img')
    vmhc.connect(z_stat, 'out_file',
                 outputNode, 'VMHC_Z_stat_FWHM_img')


    return vmhc
