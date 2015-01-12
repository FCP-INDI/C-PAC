
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.interfaces.freesurfer as fs

from nipype import logging


def create_surface_registration(wf_name='surface_registration'):
    """
    Workflow to generate surface from anatomical data and register 
    the structural data to FreeSurfer anatomical and assign it 
    to a surface vertex. 
    
    Parameters
    ----------
    wf_name : string
        name of the workflow
    
    Returns 
    -------
    wflow : workflow object
        workflow object
        
    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/timeseries/timeseries_analysis.py>`_
        
    Workflow Inputs::
        
        inputspec.rest : string  (nifti file)
            path to input functional data
        inputspec.brain : string (nifti file)
            path to skull stripped anatomical image
        inputspec.recon_subjects : string
            path to subjects directory
        inputspec.subject_id : string
            subject id 
            
    Workflow Outputs::
    
        outputspec.reconall_subjects_dir : string
            freesurfer subjects directory
        outputspec.reconall_subjects_id : string
            subject id for which anatomical data is taken
        outputspec.out_reg_file : string
            path to bbregister output registration file
        outputspec.lh_surface_file : string (mgz file)
            path to left hemisphere cortical surface file 
        outputspec.rh_surface_file : string (mgz file)
            path to right hemisphere cortical surface file 
    
    Order of commands:
     
    - Generate surfaces and rois of structural data from T1 anatomical Image using FreeSurfer's reconall . For details see `ReconAll <https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all>`_::

        recon-all -all -subjid 0010001 -sd working_dir/SurfaceRegistration/anat_reconall
        
    - Register input volume to the FreeSurfer anatomical using FreeSurfer's bbregister. The input is the output of the recon-all command. For details see `BBRegister <http://surfer.nmr.mgh.harvard.edu/fswiki/bbregister>`_::
        
        bbregister --t2 --init -fsl --reg structural_bbreg_me.dat --mov structural.nii --s 0010001
        
    - Assign values from a volume to each surface vertex using FreeSurfer's mri_vol2surf . For details see `mri_vol2surf <http://surfer.nmr.mgh.harvard.edu/fswiki/mri_vol2surf>`_::
       
        For left hemisphere
        mri_vol2surf --mov structural.nii --reg structural_bbreg_me.dat --interp trilin --projfrac 0.5 --hemi lh --o surface_file.nii.gz 

        For right hemisphere        
        mri_vol2surf --mov structural.nii --reg structural_bbreg_me.dat --interp trilin --projfrac 0.5 --hemi rh --o surface_file.nii.gz
    
    
    High Level Workflow Graph:
    
    .. image:: ../images/surface_registration.dot.png
       :width: 1000
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/surface_registration_detailed.dot.png
       :width: 1000
       
    Example
    -------
    >>> import CPAC.timeseries.timeseries_analysis as t
    >>> wf = t.create_surface_registration()
    >>> wf.inputs.inputspec.rest = '/home/data/sub001/rest.nii.gz'
    >>> wf.inputs.inputspec.brain = '/home/data/sub001/anat.nii.gz'
    >>> wf.inputs.inputspec.recon_subjects = '/home/data/working_dir/SurfaceRegistration/anat_reconall'
    >>> wf.inputs.inputspec.subject_id = 'sub001'
    >>> wf.base_dir = './'
    >>> wf.run() 
        
    """
    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['recon_subjects',
                                                       'brain',
                                                       'subject_id',
                                                       'rest']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['reconall_subjects_dir',
                                                         'reconall_subjects_id',
                                                         'out_reg_file',
                                                         'lh_surface_file',
                                                         'rh_surface_file']),
                         name='outputspec')


    reconall = pe.Node(interface=fs.ReconAll(),
                       name="reconall")
    reconall.inputs.directive = 'all'


    wflow.connect(inputNode, 'brain',
                  reconall, 'T1_files')
    wflow.connect(inputNode, 'subject_id',
                  reconall, 'subject_id')
    wflow.connect(inputNode, 'recon_subjects',
                  reconall, 'subjects_dir')

    wflow.connect(reconall, 'subjects_dir',
                  outputNode, 'reconall_subjects_dir')
    wflow.connect(reconall, 'subject_id',
                  outputNode, 'reconall_subjects_id')


    bbregister = pe.Node(interface=fs.BBRegister(init='fsl',
                                                 contrast_type='t2',
                                                 registered_file=True,
                                                 out_fsl_file=True),
                        name='bbregister')

    wflow.connect(inputNode, 'rest',
                  bbregister, 'source_file')
    wflow.connect(reconall, 'subjects_dir',
                  bbregister, 'subjects_dir')
    wflow.connect(reconall, 'subject_id',
                  bbregister, 'subject_id')

    wflow.connect(bbregister, 'out_reg_file',
                  outputNode, 'out_reg_file')


    sample_to_surface_lh = pe.Node(interface=fs.SampleToSurface(hemi="lh"),
                                    name='sample_to_surface_lh')
    sample_to_surface_lh.inputs.no_reshape = True
    sample_to_surface_lh.inputs.interp_method = 'trilinear'
    sample_to_surface_lh.inputs.sampling_method = "point"
    sample_to_surface_lh.inputs.sampling_range = 0.5
    sample_to_surface_lh.inputs.sampling_units = "frac"

    wflow.connect(bbregister, 'out_reg_file',
                  sample_to_surface_lh, 'reg_file')
    wflow.connect(inputNode, 'rest',
                  sample_to_surface_lh, 'source_file')

    wflow.connect(sample_to_surface_lh, 'out_file',
               outputNode, 'lh_surface_file')


    sample_to_surface_rh = pe.Node(interface=fs.SampleToSurface(hemi="rh"),
                                    name='sample_to_surface_rh')
    sample_to_surface_rh.inputs.no_reshape = True
    sample_to_surface_rh.inputs.interp_method = 'trilinear'
    sample_to_surface_rh.inputs.sampling_method = "point"
    sample_to_surface_rh.inputs.sampling_range = 0.5
    sample_to_surface_rh.inputs.sampling_units = "frac"

    wflow.connect(bbregister, 'out_reg_file',
                  sample_to_surface_rh, 'reg_file')
    wflow.connect(inputNode, 'rest',
                  sample_to_surface_rh, 'source_file')

    wflow.connect(sample_to_surface_rh, 'out_file',
               outputNode, 'rh_surface_file')

    return wflow


def get_voxel_timeseries(wf_name='voxel_timeseries'):
    """
    Workflow to extract time series for each voxel
    in the data that is present in the input mask
    
    Parameters
    ----------
    wf_name : string
        name of the workflow
    
    Returns 
    -------
    wflow : workflow object
        workflow object
    
    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/timeseries/timeseries_analysis.py>`_
    
    Workflow Inputs::
        
        inputspec.rest : string  (nifti file)
            path to input functional data
        inputspec.output_type : string (list of boolean)
            list of boolean for csv and npz file formats
        input_mask.masks : string (nifti file)
            path to ROI mask
        
    Workflow Outputs::
    
        outputspec.mask_outputs: string (1D, csv and/or npz files)
            list of time series matrices stored in csv and/or
            npz files.By default it outputs mean of voxels 
            across each time point in a afni compatible 1D file.
    
        High Level Workflow Graph:
    
    Example
    -------
    >>> import CPAC.timeseries.timeseries_analysis as t
    >>> wf = t.get_voxel_timeseries()
    >>> wf.inputs.inputspec.rest = '/home/data/rest.nii.gz'
    >>> wf.inputs.input_mask.mask = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz'
    >>> wf.inputs.inputspec.output_type = [True,True]
    >>> wf.base_dir = './'
    >>> wf.run()
    
    """

    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['rest',
                                                       'output_type']),
                        name='inputspec')
    inputNode_mask = pe.Node(util.IdentityInterface(fields=['mask']),
                                name='input_mask')

    outputNode = pe.Node(util.IdentityInterface(fields=['mask_outputs']),
                        name='outputspec')

    timeseries_voxel = pe.Node(util.Function(input_names=['data_file',
                                                         'template',
                                                         'output_type'],
                                            output_names=['out_file'],
                                            function=gen_voxel_timeseries),
                              name='timeseries_voxel')

    wflow.connect(inputNode, 'rest',
                  timeseries_voxel, 'data_file')
    wflow.connect(inputNode, 'output_type',
                  timeseries_voxel, 'output_type')
    wflow.connect(inputNode_mask, 'mask',
                  timeseries_voxel, 'template')

    wflow.connect(timeseries_voxel, 'out_file',
                  outputNode, 'mask_outputs')

    return wflow


def get_roi_timeseries(wf_name='roi_timeseries'):

    """
    Workflow to extract timeseries for each node in the ROI mask.
    For each node, mean across all the timepoint is calculated and stored 
    in csv and npz format.
    
    Parameters
    ----------
    wf_name : string
        name of the workflow
    
    Returns 
    -------
    wflow : workflow object
        workflow object
    
    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/timeseries/timeseries_analysis.py>`_
    
    Workflow Inputs::
        
        inputspec.rest : string  (nifti file)
            path to input functional data
        inputspec.output_type : string (list of boolean)
            list of boolean for csv and npz file formats
        input_roi.roi : string (nifti file)
            path to ROI mask
        
    Workflow Outputs::
    
        outputspec.roi_outputs : string (list of files)
            Voxel time series stored in 1D (column wise timeseries for each node), 
            csv and/or npz files. By default it outputs timeseries in a 1D file.
            The 1D file is compatible with afni interfaces.
            
    Example
    -------
    >>> import CPAC.timeseries.timeseries_analysis as t
    >>> wf = t.get_roi_timeseries()
    >>> wf.inputs.inputspec.rest = '/home/data/rest.nii.gz'
    >>> wf.inputs.input_roi.roi = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-2mm.nii.gz'
    >>> wf.inputs.inputspec.output_type = [True,True]
    >>> wf.base_dir = './'
    >>> wf.run()
    
    """

    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['rest',
                                                       'output_type']),
                        name='inputspec')

    inputnode_roi = pe.Node(util.IdentityInterface(fields=['roi']),
                                name='input_roi')

    outputNode = pe.Node(util.IdentityInterface(fields=['roi_outputs']),
                        name='outputspec')

    timeseries_roi = pe.Node(util.Function(input_names=['data_file',
                                                        'template',
                                                        'output_type'],
                                                  output_names=['out_file'],
                                                  function=gen_roi_timeseries),
                                                  name='timeseries_roi')
    wflow.connect(inputNode, 'rest',
                  timeseries_roi, 'data_file')
    wflow.connect(inputNode, 'output_type',
                  timeseries_roi, 'output_type')
    wflow.connect(inputnode_roi, 'roi',
                  timeseries_roi, 'template')

    wflow.connect(timeseries_roi, 'out_file',
                  outputNode, 'roi_outputs')


    return wflow


def get_spatial_map_timeseries(wf_name='spatial_map_timeseries'):
    """
    Workflow to regress each provided spatial
    map to the subjects functional 4D file in order
    to return a timeseries for each of the maps

    Parameters
    ----------
    wf_name : string
        name of the workflow

    Returns
    -------
    wflow : workflow object
        workflow object

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/timeseries/timeseries_analysis.py>`_

    Workflow Inputs::

        inputspec.subject_rest : string  (nifti file)
            path to input functional data
        inputspec.subject_mask : string (nifti file)
            path to subject functional mask
        inputspec.spatial_map : string (nifti file)
            path to Spatial Maps
        inputspec.demean : Boolean
            control whether to demean model and data

    Workflow Outputs::

        outputspec.subject_timeseries: string (txt file)
            list of time series stored in a space separated
            txt file
            the columns are spatial maps, rows are timepoints


    Example
    -------
    >>> import CPAC.timeseries.timeseries_analysis as t
    >>> wf = t.get_spatial_map_timeseries()
    >>> wf.inputs.inputspec.subject_rest = '/home/data/rest.nii.gz'
    >>> wf.inputs.inputspec.subject_mask = '/home/data/rest_mask.nii.gz'
    >>> wf.inputs.inputspec.ICA_map = '/home/data/spatialmaps/spatial_map.nii.gz'
    >>> wf.inputs.inputspec.demean = True
    >>> wf.base_dir = './'
    >>> wf.run()

    """

    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface
                        (fields=['subject_rest',
                                 'subject_mask',
                                 'spatial_map',
                                 'demean']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface
                         (fields=['subject_timeseries']),
                          name='outputspec')

    spatialReg = pe.Node(interface=fsl.GLM(),
                         name='spatial_regression')

    spatialReg.inputs.out_file = 'spatial_map_timeseries.txt'

    wflow.connect(inputNode, 'subject_rest',
                spatialReg, 'in_file')
    wflow.connect(inputNode, 'subject_mask',
                spatialReg, 'mask')
    wflow.connect(inputNode, 'spatial_map',
                spatialReg, 'design')
    wflow.connect(inputNode, 'demean',
                spatialReg, 'demean')

    wflow.connect(spatialReg, 'out_file',
                  outputNode, 'subject_timeseries')

    return wflow


def get_vertices_timeseries(wf_name='vertices_timeseries'):

    """
    Workflow to get vertices time series from a FreeSurfer surface file 
    
    Parameters
    ----------
    wf_name : string
        name of the workflow
    
    Returns 
    -------
    wflow : workflow object
        workflow object
    
    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/timeseries/timeseries_analysis.py>`_
    
    Workflow Inputs::
            
        inputspec.lh_surface_file : string (nifti file)
            left hemishpere surface file
        inputspec.rh_surface_file : string (nifti file)
            right hemisphere surface file
        
    Workflow Outputs::
        
        outputspec.surface_outputs: string (csv and/or npz files)
            list of timeseries matrices stored in csv and/or
            npz files
        
    Example
    -------
    >>> import CPAC.timeseries.timeseries_analysis as t
    >>> wf = t.get_vertices_timeseries()
    >>> wf.inputs.inputspec.lh_surface_file = '/home/data/outputs/SurfaceRegistration/lh_surface_file.nii.gz'
    >>> wf.inputs.inputspec.rh_surface_file = '/home/data/outputs/SurfaceRegistration/rh_surface_file.nii.gz'
    >>> wf.base_dir = './'
    >>> wf.run()
    """

    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['lh_surface_file',
                                                       'rh_surface_file']),
                        name='inputspec')

    timeseries_surface = pe.Node(util.Function(input_names=['rh_surface_file',
                                                            'lh_surface_file'],
                                                output_names=['out_file'],
                                                function=gen_vertices_timeseries),
                                name='timeseries_surface')


    outputNode = pe.Node(util.IdentityInterface(fields=['surface_outputs']),
                        name='outputspec')

    wflow.connect(inputNode, 'rh_surface_file',
                  timeseries_surface, 'rh_surface_file')
    wflow.connect(inputNode, 'lh_surface_file',
                  timeseries_surface, 'lh_surface_file')

    wflow.connect(timeseries_surface, 'out_file',
                  outputNode, 'surface_outputs')


    return wflow
    
    
    
def get_normalized_moments(wf_name='normalized_moments'):

    """
    Workflow to calculate the normalized moments for skewedness calculations
    
    Parameters
    ----------
    wf_name : string
        name of the workflow
    
    Returns 
    -------
    wflow : workflow object
        workflow object
    
    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/timeseries/timeseries_analysis.py>`_
    
    Workflow Inputs::
            
        inputspec.spatial_timeseries : string (nifti file)
            spatial map timeseries
        
    Workflow Outputs::
        
        outputspec.moments: list
            list of moment values
        
    Example
    -------
    >>> import CPAC.timeseries.timeseries_analysis as t
    >>> wf = t.get_normalized_moments()
    >>> wf.inputs.inputspec.spatial_timeseries = '/home/data/outputs/SurfaceRegistration/lh_surface_file.nii.gz'
    >>> wf.base_dir = './'
    >>> wf.run()
    """

    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['spatial_timeseries']), name='inputspec')

    # calculate normalized moments
    # output of this node is a list, 'moments'
    norm_moments = pe.Node(util.CalculateNormalizedMoments(moment='3'), name='norm_moments')

    outputNode = pe.Node(util.IdentityInterface(fields=['moments_outputs']), name='outputspec')

    wflow.connect(inputNode, 'spatial_timeseries', norm_moments, 'timeseries_file')

    wflow.connect(norm_moments, 'moments', outputNode, 'moments_outputs')


    return wflow



def gen_roi_timeseries(data_file,
                       template,
                       output_type):
    """
    Method to extract mean of voxel across
    all timepoints for each node in roi mask

    Parameters
    ----------
    datafile : string
        path to input functional data
    template : string
        path to input roi mask in functional native space
    output_type : list
        list of two boolean values suggesting
        the output types - numpy npz file and csv
        format

    Returns
    -------
    out_list : list
        list of 1D file, txt file, csv file and/or npz file containing
        mean timeseries for each scan corresponding
        to each node in roi mask

    Raises
    ------
    Exception

    """
    import nibabel as nib
    import csv
    import numpy as np
    import os
    import shutil

    unit_data = nib.load(template).get_data()
    # Cast as rounded-up integer
    unit_data = np.int64(np.ceil(unit_data))
    datafile = nib.load(data_file)
    img_data = datafile.get_data()
    vol = img_data.shape[3]

    if unit_data.shape != img_data.shape[:3]:
        raise Exception('Invalid Shape Error.'\
                        'Please check the voxel dimensions.'\
                        'Data and roi should have'\
                        'same shape')

    nodes = np.unique(unit_data).tolist()
    sorted_list = []
    node_dict = {}
    out_list = []


    # extracting filename from input template
    tmp_file = os.path.splitext(
                    os.path.basename(template))[0]
    tmp_file = os.path.splitext(tmp_file)[0]
    oneD_file = os.path.abspath('roi_' + tmp_file + '.1D')
    txt_file = os.path.abspath('roi_' + tmp_file + '.txt')
    csv_file = os.path.abspath('roi_' + tmp_file + '.csv')
    numpy_file = os.path.abspath('roi_' + tmp_file + '.npz')
    
    nodes.sort()
    for n in nodes:
        if n > 0:
            node_array = img_data[unit_data == n]
            node_str = 'node_%s' % (n)
            avg = np.mean(node_array, axis=0)
            avg = np.round(avg, 6)
            list1 = [n] + avg.tolist()
            sorted_list.append(list1)
            node_dict[node_str] = avg.tolist()


    # writing to 1Dfile
    print "writing 1D file.."
    f = open(oneD_file, 'w')
    writer = csv.writer(f, delimiter='\t')

    value_list = []

    new_keys = sorted([int(float(key.split('node_')[1])) for key in node_dict.keys()])

    roi_number_list = [str(n) for n in new_keys]

    roi_number_str = []
    for number in roi_number_list:

        roi_number_str.append("#" + number)


    print "new keys: ", new_keys
    print "roi number str: ", roi_number_str
    for key in new_keys:
        value_list.append(node_dict['node_%s' % key])

    column_list = zip(*value_list)


    writer.writerow(roi_number_str)

    for column in column_list:
        writer.writerow(list(column))
    f.close()
    out_list.append(oneD_file)

    # copy the 1D contents to txt file
    shutil.copy(oneD_file, txt_file)
    out_list.append(txt_file)

    # if csv is required
    if output_type[0]:
        print "writing csv file.."
        f = open(csv_file, 'wt')
        writer = csv.writer(f, delimiter=',',
                                quoting=csv.QUOTE_MINIMAL)
        headers = ['node/volume'] + np.arange(vol).tolist()
        writer.writerow(headers)
        writer.writerows(sorted_list)
        f.close()
        out_list.append(csv_file)

    # if npz file is required
    if output_type[1]:
        print "writing npz file.."
        np.savez(numpy_file, roi_data=value_list, roi_numbers=roi_number_list)
        out_list.append(numpy_file)

    return out_list


def gen_voxel_timeseries(data_file,
                         template,
                         output_type):
    """
    Method to extract timeseries for each voxel
    in the data that is present in the input mask
    
    Parameters
    ----------
    datafile : string (nifti file)
        path to input functional data
    template : string (nifti file) 
        path to input mask in functional native space
    output_type :list
        list of two boolean values suggesting
        the output types - numpy npz file and csv 
        format
        
    Returns
    -------
    out_list : list of files
        Based on ouput_type options method returns a list containing 
        path to npz and csv file having timeseries of each voxel in 
        the data that is present in the input mask.The row header 
        corresponds to voxel's xyz cordinates and column headers corresponds 
        to the volume index in the csv. By default it outputs afni compatible 
        1D file with mean of timeseries of voxels across timepoints.
        
    Raises
    ------
    Exception
        
    """
    import nibabel as nib
    import numpy as np
    import csv
    import os

    unit = nib.load(template)
    unit_data = unit.get_data()
    datafile = nib.load(data_file)
    img_data = datafile.get_data()
    header_data = datafile.get_header()
    qform = header_data.get_qform()
    sorted_list = []
    vol_dict = {}
    out_list = []

    if unit_data.shape != img_data.shape[:3]:
        raise Exception('Invalid Shape Error.'\
                        'Please check the voxel dimensions.'\
                        'Data and mask should have same shape')



    tmp_file = os.path.splitext(
                  os.path.basename(template))[0]
    tmp_file = os.path.splitext(tmp_file)[0]
    oneD_file = os.path.abspath('mask_' + tmp_file + '.1D')
    f = open(oneD_file, 'wt')

    x, y, z = unit_data.shape

    node_array = img_data[unit_data != 0]
    node_array = node_array.T
    time_points = node_array.shape[0]
    for t in range(0, time_points):
        str = 'vol %s' % (t)
        vol_dict[str] = node_array[t]
        print >> f, np.round(np.mean(node_array[t]), 6)
        val = node_array[t].tolist()
        val.insert(0, t)
        sorted_list.append(val)

    f.close()
    out_list.append(oneD_file)

    if output_type[0]:
        csv_file = os.path.abspath('mask_' + tmp_file + '.csv')
        f = open(csv_file, 'wt')
        writer = csv.writer(f, delimiter=',',
                            quoting=csv.QUOTE_MINIMAL)
        one = np.array([1])
        headers = ['volume/xyz']
        cordinates = np.argwhere(unit_data != 0)
        for val in range(np.alen(cordinates)):
            ijk_mat = np.concatenate([cordinates[val], one])
            ijk_mat = ijk_mat.T
            product = np.dot(qform, ijk_mat)
            val = tuple(product.tolist()[0:3])
            headers.append(val)
        writer.writerow(headers)
        writer.writerows(sorted_list)
        f.close()
        out_list.append(csv_file)

    if output_type[1]:
        numpy_file = os.path.abspath('mask_' + tmp_file + '.npz')
        np.savez(numpy_file, **dict(vol_dict))
        out_list.append(numpy_file)

    return out_list

def gen_vertices_timeseries(rh_surface_file,
                        lh_surface_file):

    """
    Method to extract timeseries from vertices
    of a freesurfer surface file 
    
    Parameters
    ----------
    rh_surface_file : string (mgz/mgh file)
        left hemisphere FreeSurfer surface file
    lh_surface_file : string (mgz/mgh file)
        right hemisphere FreeSurfer surface file
        
    Returns
    -------
    out_list : string (list of file)
        list of vertices timeseries csv files  
    
    """

    import gradunwarp
    import numpy as np
    import os

    out_list = []
    rh_file = os.path.splitext(
                    os.path.basename(rh_surface_file))[0] + '_rh.csv'
    mghobj1 = gradunwarp.mgh.MGH()

    mghobj1.load(rh_surface_file)
    vol = mghobj1.vol
    (x, y) = vol.shape
#        print "rh shape", x, y

    np.savetxt(rh_file, vol, delimiter='\t')
    out_list.append(rh_file)

    lh_file = os.path.splitext(os.path.basename(lh_surface_file))[0] + '_lh.csv'
    mghobj2 = gradunwarp.mgh.MGH()

    mghobj2.load(lh_surface_file)
    vol = mghobj2.vol
    (x, y) = vol.shape
#        print "lh shape", x, y

    np.savetxt(lh_file,
               vol,
               delimiter=',')
    out_list.append(lh_file)

    return out_list
