
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
from CPAC.easy_thresh import easy_thresh



def create_group_analysis(ftest=False, wf_name='groupAnalysis'):
    """
    FSL `FEAT <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT>`_
    BASED Group Analysis

    Parameters
    ----------
    ftest : boolean, optional(default=False)
        Ftest help investigate several contrasts at the same time
        for example to see whether any of them (or any combination of them) is 
        significantly non-zero. Also, the F-test allows you to compare the 
        contribution of each contrast to the model and decide on significant 
        and non-significant ones
 
    wf_name : string 
        Workflow name
    
    Returns 
    -------
    grp_analysis : workflow object
        Group Analysis workflow object
    
    Notes
    -----
    `Source <https://github.com/openconnectome/C-PAC/blob/master/CPAC/group_analysis/group_analysis_preproc.py>`_
 
    Workflow Inputs::
        
        inputspec.mat_file : string (existing file)
           Mat file containing  matrix for design 
        
        inputspec.con_file : string (existing file)
           Contrast file containing contrast vectors 
        
        inputspec.grp_file : string (existing file)
           file containing matrix specifying the groups the covariance is split into
        
        inputspec.zmap_files : string (existing nifti file)
           derivative or the zmap file for which the group analysis is to be run
        
        inputspec.z_threshold : float
            Z Statistic threshold value for cluster thresholding. It is used to 
            determine what level of activation would be statistically significant. 
            Increasing this will result in higher estimates of required effect.
        
        inputspec.p_threshold : float
            Probability threshold for cluster thresholding.
            
        inputspec.fts_file : string (existing file)
           file containing matrix specifying f-contrasts
           
        inputspec.paramerters : string (tuple)
            tuple containing which MNI and FSLDIR path information
                      
    Workflow Outputs::
    
        outputspec.merged : string (nifti file)
            4D volume file after merging all the derivative 
            files from each specified subject.
            
        outputspec.zstats : list (nifti files)
            Z statistic image for each t contrast
            
        outputspec.zfstats : list (nifti files)
            Z statistic image for each f contrast
        
        outputspec.fstats : list (nifti files)
            F statistic for each contrast  
        
        outputspec.cluster_threshold : list (nifti files)
           the thresholded Z statistic image for each t contrast
        
        outputspec.cluster_index : list (nifti files)
            image of clusters for each t contrast; the values 
            in the clusters are the index numbers as used 
            in the cluster list.
        
        outputspec.cluster_localmax_txt : list (text files)
            local maxima text file for each t contrast, 
            defines the coordinates of maximum value in the cluster
        
        outputspec.overlay_threshold : list (nifti files)
            3D color rendered stats overlay image for t contrast
            After reloading this image, use the Statistics Color 
            Rendering GUI to reload the color look-up-table
        
        outputspec.overlay_rendered_image : list (nifti files)
           2D color rendered stats overlay picture for each t contrast
            
        outputspec.cluster_threshold_zf : list (nifti files)
           the thresholded Z statistic image for each f contrast
        
        outputspec.cluster_index_zf : list (nifti files)
            image of clusters for each f contrast; the values 
            in the clusters are the index numbers as used 
            in the cluster list.
            
        outputspec.cluster_localmax_txt_zf : list (text files)
            local maxima text file for each f contrast, 
            defines the coordinates of maximum value in the cluster
        
        outputspec.overlay_threshold_zf : list (nifti files)
            3D color rendered stats overlay image for f contrast
            After reloading this image, use the Statistics Color 
            Rendering GUI to reload the color look-up-table
        
        outputspec.overlay_rendered_image_zf : list (nifti files)
           2D color rendered stats overlay picture for each f contrast
    
    Order of commands:

    - Merge all the Z-map 3D images into 4D image file.  For details see `fslmerge <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Fslutils>`_::
    
        fslmerge -t sub01/sca/seed1/sca_Z_FWHM_merged.nii 
                    sub02/sca/seed1/sca_Z_FWHM.nii.gz ....  
                    merge.nii.gz
                    
        arguments 
            -t : concatenate images in time
            
    - Create mask specific for analysis. For details see `fslmaths <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Fslutils>`_::
    
        fslmaths merged.nii.gz 
                -abs -Tmin -bin mean_mask.nii.gz
        
        arguments 
             -Tmin  : min across time
             -abs   : absolute value
             -bin   : use (current image>0) to binarise
    
    - FSL FLAMEO to perform higher level analysis.  For details see `flameo <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT>`_::
        
        flameo --copefile = merged.nii.gz --covsplitfile = anova_with_meanFD.grp --designfile = anova_with_meanFD.mat 
               --fcontrastsfile = anova_with_meanFD.fts --ld=stats --maskfile = mean_mask.nii.gz --runmode=ols 
               --tcontrastsfile = anova_with_meanFD.con
           
        arguments
            --copefile        : cope regressor data file
            --designfile      : design matrix file
            --maskfile        : mask file
            --tcontrastsfile  : file containing an ASCII matrix specifying the t contrasts
            --fcontrastsfile  : file containing an ASCII matrix specifying the f contrasts
            --runmode         : Interference to perform (mixed effects - OLS)
            
    - Run FSL Easy thresh 
        
      Easy thresh is a simple script for carrying out cluster-based thresholding and colour activation overlaying::
        
        easythresh <raw_zstat> <brain_mask> <z_thresh> <prob_thresh> <background_image> <output_root> [--mm]
      
      A seperate workflow called easythresh is called to run easythresh steps.
      
    High Level Workflow Graph:
    
    .. image:: ../images/group_analysis.dot.png
       :width: 800
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/group_analysis_detailed.dot.png
       :width: 800

    Examples
    --------
    
    >>> from group_analysis_preproc import create_group_analysis
    >>> preproc = create_group_analysis()
    >>> preproc.inputs.inputspec.mat_file = '../group_models/anova_with_meanFD/anova_with_meanFD.mat'
    >>> preproc.inputs.inputspec.con_file = '../group_models/anova_with_meanFD/anova_with_meanFD.con'
    >>> preproc.inputs.inputspec.grp_file = '../group_models/anova_with_meanFD/anova_with_meanFD.grp'
    >>> preproc.inputs.inputspec.zmap_files = ['subjects/sub01/seeds_rest_Dickstein_DLPFC/sca_Z_FWHM.nii.gz', 
                                               'subjects/sub02/seeds_rest_Dickstein_DLPFC/sca_Z_FWHM.nii.gz']
    >>> preproc.inputs.inputspec.z_threshold = 2.3
    >>> preproc.inputs.inputspec.p_threshold = 0.05
    >>> preproc.inputs.inputspec.parameters = ('/usr/local/fsl/', 'MNI152')
    >>> preproc.run()  -- SKIP doctest
            
    """
    grp_analysis = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['merged_file',
                                                       'merge_mask',
                                                       'mat_file',
                                                       'con_file',
                                                       'grp_file',
                                                       'fts_file',
                                                       'z_threshold',
                                                       'p_threshold',
                                                       'parameters']),
                         name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['merged',
                                                        'zstats',
                                                        'zfstats',
                                                        'fstats',
                                                        'cluster_threshold',
                                                        'cluster_index',
                                                        'cluster_localmax_txt',
                                                        'overlay_threshold',
                                                        'rendered_image',
                                                        'cluster_localmax_txt_zf',
                                                        'cluster_threshold_zf',
                                                        'cluster_index_zf',
                                                        'overlay_threshold_zf',
                                                        'rendered_image_zf']),
                         name='outputspec')



    '''
    merge_to_4d = pe.Node(interface=fsl.Merge(),
                          name='merge_to_4d')
    merge_to_4d.inputs.dimension = 't'
    

    ### create analysis specific mask
    #-Tmin: min across time
    # -abs: absolute value
    #-bin: use (current image>0) to binarise
    merge_mask = pe.Node(interface=fsl.ImageMaths(),
                         name='merge_mask')
    merge_mask.inputs.op_string = '-abs -Tmin -bin'
    '''

    fsl_flameo = pe.Node(interface=fsl.FLAMEO(),
                         name='fsl_flameo')
    fsl_flameo.inputs.run_mode = 'ols'

    ### create analysis specific mask
    # fslmaths merged.nii.gz -abs -bin -Tmean -mul volume out.nii.gz
    #-Tmean: mean across time
    # create group_reg file
    # this file can provide an idea of how well the subjects 
    # in our analysis overlay with each other and the MNI brain.
    # e.g., maybe there is one subject with limited coverage.
    # not attached to sink currently 
    merge_mean_mask = pe.Node(interface=fsl.ImageMaths(),
                              name='merge_mean_mask')

    #function node to get the operation string for fslmaths command
    get_opstring = pe.Node(util.Function(input_names=['in_file'],
                                         output_names=['out_file'],
                                         function=get_operation),
                       name='get_opstring')

    #connections
    '''
    grp_analysis.connect(inputnode, 'zmap_files',
                         merge_to_4d, 'in_files')
    grp_analysis.connect(merge_to_4d, 'merged_file',
                         merge_mask, 'in_file')
    '''
    grp_analysis.connect(inputnode, 'merged_file',
                         fsl_flameo, 'cope_file')
    grp_analysis.connect(inputnode, 'merge_mask',
                         fsl_flameo, 'mask_file')
    grp_analysis.connect(inputnode, 'mat_file',
                         fsl_flameo, 'design_file')
    grp_analysis.connect(inputnode, 'con_file',
                         fsl_flameo, 't_con_file')
    grp_analysis.connect(inputnode, 'grp_file',
                         fsl_flameo, 'cov_split_file')

    if ftest:

        #calling easythresh for zfstats file
        grp_analysis.connect(inputnode, 'fts_file',
                             fsl_flameo, 'f_con_file')

        easy_thresh_zf = easy_thresh('easy_thresh_zf')

        grp_analysis.connect(fsl_flameo, 'zfstats',
                             easy_thresh_zf, 'inputspec.z_stats')
        grp_analysis.connect(inputnode, 'merge_mask',
                             easy_thresh_zf, 'inputspec.merge_mask')
        grp_analysis.connect(inputnode, 'z_threshold',
                             easy_thresh_zf, 'inputspec.z_threshold')
        grp_analysis.connect(inputnode, 'p_threshold',
                             easy_thresh_zf, 'inputspec.p_threshold')
        grp_analysis.connect(inputnode, 'parameters',
                             easy_thresh_zf, 'inputspec.parameters')
        grp_analysis.connect(easy_thresh_zf, 'outputspec.cluster_threshold',
                             outputnode, 'cluster_threshold_zf')
        grp_analysis.connect(easy_thresh_zf, 'outputspec.cluster_index',
                             outputnode, 'cluster_index_zf')
        grp_analysis.connect(easy_thresh_zf, 'outputspec.cluster_localmax_txt',
                             outputnode, 'cluster_localmax_txt_zf')
        grp_analysis.connect(easy_thresh_zf, 'outputspec.overlay_threshold',
                             outputnode, 'overlay_threshold_zf')
        grp_analysis.connect(easy_thresh_zf, 'outputspec.rendered_image',
                             outputnode, 'rendered_image_zf')

    #calling easythresh for zstats files
    easy_thresh_z = easy_thresh('easy_thresh_z')
    grp_analysis.connect(fsl_flameo, 'zstats',
                         easy_thresh_z, 'inputspec.z_stats')
    grp_analysis.connect(inputnode, 'merge_mask',
                         easy_thresh_z, 'inputspec.merge_mask')
    grp_analysis.connect(inputnode, 'z_threshold',
                         easy_thresh_z, 'inputspec.z_threshold')
    grp_analysis.connect(inputnode, 'p_threshold',
                         easy_thresh_z, 'inputspec.p_threshold')
    grp_analysis.connect(inputnode, 'parameters',
                         easy_thresh_z, 'inputspec.parameters')

    grp_analysis.connect(inputnode, 'merged_file',
                         get_opstring, 'in_file')
    grp_analysis.connect(inputnode, 'merged_file',
                         merge_mean_mask, 'in_file')
    grp_analysis.connect(get_opstring, 'out_file',
                         merge_mean_mask, 'op_string')

    grp_analysis.connect(fsl_flameo, 'zfstats',
                         outputnode, 'zfstats')
    grp_analysis.connect(fsl_flameo, 'fstats',
                         outputnode, 'fstats')
    grp_analysis.connect(inputnode, 'merged_file',
                         outputnode, 'merged')
    grp_analysis.connect(fsl_flameo, 'zstats',
                         outputnode, 'zstats')



    grp_analysis.connect(easy_thresh_z, 'outputspec.cluster_threshold',
                         outputnode, 'cluster_threshold')
    grp_analysis.connect(easy_thresh_z, 'outputspec.cluster_index',
                         outputnode, 'cluster_index')
    grp_analysis.connect(easy_thresh_z, 'outputspec.cluster_localmax_txt',
                         outputnode, 'cluster_localmax_txt')
    grp_analysis.connect(easy_thresh_z, 'outputspec.overlay_threshold',
                         outputnode, 'overlay_threshold')
    grp_analysis.connect(easy_thresh_z, 'outputspec.rendered_image',
                         outputnode, 'rendered_image')

    return grp_analysis



def get_operation(in_file):
    """
    Method to create operation string 
    for fslmaths
    
    Parameters
    ----------
    in_file : file
       input volume
    
    Returns
    -------
    op_string : string
        operation string for fslmaths
    
    Raises
    ------
    IOError
      If unable to load the input volume
    
    """
    try:
        from nibabel import load
        img = load(in_file)
        hdr = img.get_header()
        n_vol = int(hdr.get_data_shape()[3])
        op_string = '-abs -bin -Tmean -mul %d' % (n_vol)
        return op_string
    except:
        raise IOError("Unable to load the input nifti image")
