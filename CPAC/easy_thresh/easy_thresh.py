
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util


def easy_thresh(wf_name):
    """
    Workflow for carrying out cluster-based thresholding 
    and colour activation overlaying
    
    Parameters
    ----------
    wf_name : string 
        Workflow name
        
    Returns
    -------
    easy_thresh : object 
        Easy thresh workflow object
    
    Notes
    -----
    
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/easy_thresh/easy_thresh.py>`_
        
    Workflow Inputs::
    
        inputspec.z_stats : string (nifti file)
            z_score stats output for t or f contrast from flameo
        
        inputspec.merge_mask : string (nifti file)
            mask generated from 4D Merged derivative file
        
        inputspec.z_threshold : float
            Z Statistic threshold value for cluster thresholding. It is used to 
            determine what level of activation would be statistically significant. 
            Increasing this will result in higher estimates of required effect.
        
        inputspec.p_threshold : float
            Probability threshold for cluster thresholding.
        
        inputspec.paramerters : string (tuple)
            tuple containing which MNI and FSLDIR path information
            
    Workflow Outputs::
    
        outputspec.cluster_threshold : string (nifti files)
           the thresholded Z statistic image for each t contrast
        
        outputspec.cluster_index : string (nifti files)
            image of clusters for each t contrast; the values 
            in the clusters are the index numbers as used 
            in the cluster list.
        
        outputspec.overlay_threshold : string (nifti files)
            3D color rendered stats overlay image for t contrast
            After reloading this image, use the Statistics Color 
            Rendering GUI to reload the color look-up-table
        
        outputspec.overlay_rendered_image : string (nifti files)
           2D color rendered stats overlay picture for each t contrast
        
        outputspec.cluster_localmax_txt : string (text files)
            local maxima text file, defines the coordinates of maximum value
            in the cluster
    
    
    Order of commands:
    
    - Estimate smoothness of the image::
        
        smoothest --mask= merge_mask.nii.gz --zstat=.../flameo/stats/zstat1.nii.gz
        
        arguments
        --mask  :  brain mask volume
        --zstat :  filename of zstat/zfstat image
    
    - Create mask. For details see `fslmaths <http://www.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro/index.htm#fslutils>`_::
        
        fslmaths ../flameo/stats/zstat1.nii.gz 
                 -mas merge_mask.nii.gz 
                 zstat1_mask.nii.gz
        
        arguments
        -mas   : use (following image>0) to mask current image

    - Copy Geometry image dimensions, voxel dimensions, voxel dimensions units string, image orientation/origin or qform/sform info) from one image to another::
    
        fslcpgeom MNI152_T1_2mm_brain.nii.gz zstat1_mask.nii.gz
    
    - Cluster based thresholding. For details see `FEAT <http://www.fmrib.ox.ac.uk/fsl/feat5/detail.html#poststats>`_::
        
        cluster --dlh = 0.0023683100 
                --in = zstat1_mask.nii.gz 
                --oindex = zstat1_cluster_index.nii.gz 
                --olmax = zstat1_cluster_localmax.txt
                --othresh = zstat1_cluster_threshold.nii.gz 
                --pthresh = 0.0500000000 
                --thresh = 2.3000000000 
                --volume = 197071
                
        arguments 
        --in    :    filename of input volume
        --dlh   :    smoothness estimate = sqrt(det(Lambda))
        --oindex  :  filename for output of cluster index
        --othresh :  filename for output of thresholded image
        --olmax   :  filename for output of local maxima text file
        --volume  :  number of voxels in the mask
        --pthresh :  p-threshold for clusters
        --thresh  :  threshold for input volume
        
     Z statistic image is thresholded to show which voxels or clusters of voxels are activated at a particular significance level.
     A Z statistic threshold is used to define contiguous clusters. Then each cluster's estimated significance level (from GRF-theory) 
     is compared with the cluster probability threshold. Significant clusters are then used to mask the original Z statistic image.
    
    - Get the maximum intensity value of the output thresholded image. This used is while rendering the Z statistic image:: 
        
        fslstats zstat1_cluster_threshold.nii.gz -R
        
        arguments
        -R  : output <min intensity> <max intensity>

    - Rendering. For details see `FEAT <http://www.fmrib.ox.ac.uk/fsl/feat5/detail.html#poststats>`_::
         
        overlay 1 0 MNI152_T1_2mm_brain.nii.gz 
               -a zstat1_cluster_threshold.nii.gz 
               2.30 15.67 
               zstat1_cluster_threshold_overlay.nii.gz
               
        slicer zstat1_cluster_threshold_overlay.nii.gz 
               -L  -A 750 
               zstat1_cluster_threshold_overlay.png
    
      The Z statistic range selected for rendering is automatically calculated by default, 
      to run from red (minimum Z statistic after thresholding) to yellow (maximum Z statistic, here 
      maximum intensity).
      
    High Level Workflow Graph:
    
    .. exec::
        from CPAC.easy_thresh import easy_thresh
        wf = easy_thresh('easy_thresh_wf')
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/easy_thresh.dot'
        )

    .. image:: ../../images/generated/easy_thresh.png
       :width: 800
    
    
    Detailed Workflow Graph:
    
    .. image:: ../../images/generated/generated.easy_thresh_detailed.png
       :width: 800
               
    Examples
    --------
    
    >>> import easy_thresh
    >>> preproc = easy_thresh.easy_thresh("new_workflow")
    >>> preproc.inputs.inputspec.z_stats= 'flameo/stats/zstat1.nii.gz'
    >>> preproc.inputs.inputspec.merge_mask = 'merge_mask/alff_Z_fn2standard_merged_mask.nii.gz'
    >>> preproc.inputs.inputspec.z_threshold = 2.3
    >>> preproc.inputs.inputspec.p_threshold = 0.05
    >>> preproc.inputs.inputspec.parameters = ('/usr/local/fsl/', 'MNI152')
    >>> preporc.run()  -- SKIP doctest
    
    """

    easy_thresh = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['z_stats',
                                                       'merge_mask',
                                                       'z_threshold',
                                                       'p_threshold',
                                                       'parameters']),
                         name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['cluster_threshold',
                                                        'cluster_index',
                                                        'cluster_localmax_txt',
                                                        'overlay_threshold',
                                                        'rendered_image']),
                         name='outputspec')

    # fsl easythresh
    # estimate image smoothness
    smooth_estimate = pe.MapNode(interface=fsl.SmoothEstimate(),
                                    name='smooth_estimate',
                                    iterfield=['zstat_file'])

    # run clustering after fixing stats header for talspace
    zstat_mask = pe.MapNode(interface=fsl.MultiImageMaths(),
                                  name='zstat_mask',
                                  iterfield=['in_file'])
    # operations to perform
    #-mas use (following image>0) to mask current image
    zstat_mask.inputs.op_string = '-mas %s'

    # fslcpgeom
    #copy certain parts of the header information (image dimensions,
    #voxel dimensions, voxel dimensions units string, image orientation/origin
    #or qform/sform info) from one image to another
    geo_imports = ['import subprocess']
    copy_geometry = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'],
                                             output_names=['out_file'],
                                             function=copy_geom,
                                             imports=geo_imports),
                                             name='copy_geometry',
                                             iterfield=['infile_a', 'infile_b'])

    ##cluster-based thresholding
    #After carrying out the initial statistical test, the resulting
    #Z statistic image is then normally thresholded to show which voxels or
    #clusters of voxels are activated at a particular significance level.
    #A Z statistic threshold is used to define contiguous clusters.
    #Then each cluster's estimated significance level (from GRF-theory) is
    #compared with the cluster probability threshold. Significant clusters
    #are then used to mask the original Z statistic image for later production
    #of colour blobs.This method of thresholding is an alternative to
    #Voxel-based correction, and is normally more sensitive to activation.
#    cluster = pe.MapNode(interface=fsl.Cluster(),
#                            name='cluster',
#                            iterfield=['in_file', 'volume', 'dlh'])
#    #output of cluster index (in size order)
#    cluster.inputs.out_index_file = True
#    #thresholded image
#    cluster.inputs.out_threshold_file = True
#    #local maxima text file
#    #defines the cluster cordinates
#    cluster.inputs.out_localmax_txt_file = True

    cluster_imports = ['import os', 'import re', 'import subprocess as sb']
    cluster = pe.MapNode(util.Function(input_names=['in_file',
                                                    'volume',
                                                    'dlh',
                                                    'threshold',
                                                    'pthreshold',
                                                    'parameters'],
                                       output_names=['index_file',
                                                     'threshold_file',
                                                     'localmax_txt_file'],
                                       function=call_cluster,
                                       imports=cluster_imports),
                         name='cluster',
                         iterfield=['in_file', 'volume', 'dlh'])

    # max and minimum intensity values
    image_stats = pe.MapNode(interface=fsl.ImageStats(),
                             name='image_stats',
                             iterfield=['in_file'])
    image_stats.inputs.op_string = '-R'

    # create tuple of z_threshold and max intensity value of threshold file
    create_tuple = pe.MapNode(util.Function(input_names=['infile_a',
                                                         'infile_b'],
                                            output_names=['out_file'],
                                            function=get_tuple),
                              name='create_tuple',
                              iterfield=['infile_b'])

    # colour activation overlaying
    overlay = pe.MapNode(interface=fsl.Overlay(),
                         name='overlay',
                         iterfield=['stat_image', 'stat_thresh'])
    overlay.inputs.transparency = True
    overlay.inputs.auto_thresh_bg = True
    overlay.inputs.out_type = 'float'

    # colour rendering
    slicer = pe.MapNode(interface=fsl.Slicer(), name='slicer',
                        iterfield=['in_file'])
    # set max picture width
    slicer.inputs.image_width = 750
    # set output all axial slices into one picture
    slicer.inputs.all_axial = True

    # function mapnode to get the standard fsl brain image based on parameters
    # as FSLDIR,MNI and voxel size
    get_bg_imports = ['import os', 'from nibabel import load']
    get_backgroundimage = pe.MapNode(util.Function(input_names=['in_file',
                                                                'file_parameters'],
                                                   output_names=['out_file'],
                                                   function=get_standard_background_img,
                                                   imports=get_bg_imports),
                                     name='get_bckgrndimg1',
                                     iterfield=['in_file'])

    # function node to get the standard fsl brain image
    # outputs single file
    get_backgroundimage2 = pe.Node(util.Function(input_names=['in_file',
                                                              'file_parameters'],
                                                 output_names=['out_file'],
                                                 function=get_standard_background_img,
                                                 imports=get_bg_imports),
                                   name='get_backgrndimg2')

    # connections
    easy_thresh.connect(inputnode, 'z_stats', smooth_estimate, 'zstat_file')
    easy_thresh.connect(inputnode, 'merge_mask', smooth_estimate, 'mask_file')

    easy_thresh.connect(inputnode, 'z_stats', zstat_mask, 'in_file')
    easy_thresh.connect(inputnode, 'merge_mask', zstat_mask, 'operand_files')

    easy_thresh.connect(zstat_mask, 'out_file', get_backgroundimage,
                        'in_file')
    easy_thresh.connect(inputnode, 'parameters', get_backgroundimage,
                        'file_parameters')

    easy_thresh.connect(get_backgroundimage, 'out_file', copy_geometry,
                        'infile_a')
    easy_thresh.connect(zstat_mask, 'out_file', copy_geometry, 'infile_b')

    easy_thresh.connect(copy_geometry, 'out_file', cluster, 'in_file')
    easy_thresh.connect(inputnode, 'z_threshold', cluster, 'threshold')
    easy_thresh.connect(inputnode, 'p_threshold', cluster, 'pthreshold')
    easy_thresh.connect(smooth_estimate, 'volume', cluster, 'volume')
    easy_thresh.connect(smooth_estimate, 'dlh', cluster, 'dlh')
    easy_thresh.connect(inputnode, 'parameters', cluster, 'parameters')

    easy_thresh.connect(cluster, 'threshold_file', image_stats, 'in_file')

    easy_thresh.connect(image_stats, 'out_stat', create_tuple, 'infile_b')
    easy_thresh.connect(inputnode, 'z_threshold', create_tuple, 'infile_a')

    easy_thresh.connect(cluster, 'threshold_file', overlay, 'stat_image')
    easy_thresh.connect(create_tuple, 'out_file', overlay, 'stat_thresh')

    easy_thresh.connect(inputnode, 'merge_mask', get_backgroundimage2,
                        'in_file' )
    easy_thresh.connect(inputnode, 'parameters', get_backgroundimage2,
                        'file_parameters')

    easy_thresh.connect(get_backgroundimage2, 'out_file', overlay,
                        'background_image')

    easy_thresh.connect(overlay, 'out_file', slicer, 'in_file')

    easy_thresh.connect(cluster, 'threshold_file', outputnode,
                        'cluster_threshold')
    easy_thresh.connect(cluster, 'index_file', outputnode, 'cluster_index')
    easy_thresh.connect(cluster, 'localmax_txt_file', outputnode,
                        'cluster_localmax_txt')
    easy_thresh.connect(overlay, 'out_file', outputnode, 'overlay_threshold')
    easy_thresh.connect(slicer, 'out_file', outputnode, 'rendered_image')

    return easy_thresh


def call_cluster(in_file, volume, dlh, threshold, pthreshold, parameters):

    out_name = re.match('z(\w)*stat(\d)+', os.path.basename(in_file))

    filename, ext = os.path.splitext(os.path.basename(in_file))
    ext=  os.path.splitext(filename)[1] + ext
    filename = os.path.splitext(filename)[0]

    if out_name:
        out_name= out_name.group(0)
    else:
        out_name = filename

    FSLDIR = parameters[0]

    index_file = os.path.join(os.getcwd(), 'cluster_mask_' + out_name + ext)
    threshold_file = os.path.join(os.getcwd(), 'thresh_' + out_name + ext)
    localmax_txt_file = os.path.join(os.getcwd(), 'cluster_' + out_name + '.txt')

    cmd_path = os.path.join(FSLDIR, 'bin/cluster')

    f = open(localmax_txt_file,'wb')

    cmd = sb.Popen([ cmd_path,
                    '--dlh=' + str(dlh),
                    '--in=' + in_file,
                    '--oindex=' + index_file,
                    '--othresh=' + threshold_file,
                    '--pthresh=' + str(pthreshold),
                    '--thresh=' +  str(threshold),
                    '--volume=' + str(volume)],
                    stdout= f)

    stdout_value, stderr_value = cmd.communicate()
    f.close()

    return index_file, threshold_file, localmax_txt_file


def copy_geom(infile_a, infile_b):
    """
    Method to call fsl fslcpgeom command to copy 
    certain parts of the header information (image dimensions, 
    voxel dimensions, voxel dimensions units string, image 
    orientation/origin or qform/sform info) from one image to another
    
    Parameters
    -----------
    infile_a : nifti file
        input volume from which the geometry is copied from
    
    infile_b : nifti file
       input volume from which the geometry is copied to
       
    Returns
    -------    
    out_file : nifti file
        Input volume infile_b with modified geometry information
        in the header.
        
    Raises
    ------
    Exception 
        If fslcpgeom fails
    
    """

    try:
        out_file = infile_b
        cmd = ['fslcpgeom', infile_a, out_file]
        retcode = subprocess.check_output(cmd)
        return out_file
    except Exception:
        raise Exception("Error while using fslcpgeom to copy geometry")


def get_standard_background_img(in_file, file_parameters):
    """
    Method to get the standard brain image from FSL 
    standard data directory
    
    Parameters
    ----------
    in_file : nifti file
        Merged 4D Zmap volume
    file_parameters : tuple
       Value FSLDIR and MNI from config file
    
    Returns
    -------
    standard_path : string
        Standard FSL Image file path
    
    Raises
    ------
    Exception 
        If nibabel cannot load the input nifti volume
    
    """

    try:

        img = load(in_file)
        hdr = img.get_header()
        group_mm = int(hdr.get_zooms()[2])
        FSLDIR, MNI = file_parameters
        standard_path = \
            os.path.join(FSLDIR, 'data/standard/',
                         '{0}_T1_{1}mm_brain.nii.gz'.format(MNI, group_mm))
        return os.path.abspath(standard_path)

    except Exception:
        raise Exception("Error while loading background image")


def get_tuple(infile_a, infile_b):
    """
    Simple method to return tuple of z_threhsold
    maximum intensity values of Zstatistic image
    for input to the overlay.
    
    Parameters
    ----------
    z_theshold : float
        z threshold value
    intensity_stat : tuple of float values
        minimum and maximum intensity values
    
    Returns
    -------
    img_min_max : tuple (float)
        tuple of zthreshold and maximum intensity 
        value of z statistic image
    
    """
    out_file = (infile_a, infile_b[1])
    return out_file
