import sys
from CPAC.interfaces.afni import preprocess
import os
import commands
import numpy as np
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from CPAC.sca.utils import *


def create_sca(name_sca='sca'):
    """
    Map of the correlations of the Region of Interest(Seed in native or MNI space) with the rest of brain voxels.
    The map is normalized to contain Z-scores, mapped in standard space and treated with spatial smoothing.

    Parameters
    ----------

    name_sca : a string
        Name of the SCA workflow

    Returns
    -------

    sca_workflow : workflow

        Seed Based Correlation Analysis Workflow



    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/sca/sca.py>`_

    Workflow Inputs::


        inputspec.rest_res_filt : string (existing nifti file)
            Band passed Image with Global Signal , white matter, csf and motion regression. Recommended bandpass filter (0.001,0.1) )

        inputspec.timeseries_one_d : string (existing nifti file)
            1D 3dTcorr1D compatible timeseries file. 1D file can be timeseries from a mask or from a parcellation containing ROIs




    Workflow Outputs::

        outputspec.correlation_file : string (nifti file)
            Correlations of the functional file and the input time series

        outputspec.Z_score : string (nifti file)
            Fisher Z transformed correlations of the seed


    SCA Workflow Procedure:

    1. Compute pearson correlation between input timeseries 1D file and input functional file
       Use 3dTcorr1D to compute that. Input timeseries can be a 1D file containing parcellation ROI's
       or a 3D mask

    2. Compute Fisher Z score of the correlation computed in step above. If a mask is provided then a
       a single Z score file is returned, otherwise z-scores for all ROIs are returned as a list of
       nifti files



    Workflow:

    .. image:: ../images/sca_graph.dot.png
        :width: 500

    Detailed Workflow:

    .. image:: ../images/sca_detailed_graph.dot.png
        :width: 500


    Examples
    --------

    >>> sca_w = create_sca("sca_wf")
    >>> sca_w.inputs.inputspec.rest_res_filt = '/home/data/subject/func/rest_bandpassed.nii.gz'
    >>> sca_w.inputs.inputspec.timeseries_one_d = '/home/data/subject/func/ts.1D'
    >>> sca_w.run() # doctest: +SKIP

    """

    sca = pe.Workflow(name=name_sca)
    inputNode = pe.Node(util.IdentityInterface(fields=['timeseries_one_d',
                                                'functional_file',
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'correlation_file',
                                                    'Z_score',
                                                    ]),
                        name='outputspec')



    # # 2. Compute voxel-wise correlation with Seed Timeseries
    corr = pe.Node(interface=preprocess.ThreedTcorrOneD(),
                      name='corr')
    corr.inputs.options = "-pearson"
    corr.inputs.outputtype = 'NIFTI_GZ'

    # # 3. Z-transform correlations
    z_score = pe.Node(util.Function(input_names=['correlation_file', 'timeseries_one_d'],
                               output_names=['out_file'],
                 function=compute_fisher_z_score), name='z_score')


    sca.connect(inputNode, 'timeseries_one_d',
                corr, 'y_one_d')
    sca.connect(inputNode, 'functional_file',
                corr, 'xset')
    sca.connect(corr, 'out_file',
                z_score, 'correlation_file')
    sca.connect(inputNode, 'timeseries_one_d',
                z_score, 'timeseries_one_d')
    sca.connect(corr, 'out_file',
                outputNode, 'correlation_file')
    sca.connect(z_score, 'out_file',
                outputNode, 'Z_score')

    return sca


def create_ca(wflow_name='ca', which= 'SR'):
    """
    Map of the correlations of the signal in the provided time series with all
    the voxels in the brain.

    Parameters
    ----------

    wflow_name : a string
        Name of the temporal regression workflow

    which: a string
        SR: Spatial Regression, RT: ROI Timeseries
        
    Returns
    -------

    wflow : workflow

        temporal multiple regression Workflow



    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/sca/sca.py>`_

    Workflow Inputs::

        inputspec.subject_rest : string (existing nifti file)
            Band passed Image with Global Signal , white matter, csf and motion regression. Recommended bandpass filter (0.001,0.1) )

        inputspec.subject_timeseries : string (existing txt file)
            text file containing the timeseries to be regressed on the subjects
            functional file 
            timeseries are organized by columns, timepoints by rows
        
        inputspec.subject_mask : string (existing nifti file)
            path to subject functional mask
            
        inputspec.demean : Boolean
            control whether to demean model and data
            
        inputspec.normalize : Boolean
            control whether to normalize the input timeseries to unit standard deviation



    Workflow Outputs::

        outputspec.component_map : string (nifti file)
            GLM parameter estimate image for each timeseries in the input file

        outputspec.component_map_z : string (nifti file)
            Normalized version of the GLM parameter estimates


    Temporal Regression Workflow Procedure:
    
    Enter all timeseries into a general linear model and regress these 
    timeseries to the subjects functional file to get spatial maps of voxels
    showing activation patterns related to those in the timeseries.

    References
    ----------
    `http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/DualRegression/UserGuide <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/DualRegression/UserGuide>`_

    Examples
    --------

    >>> ca_w = create_ca('temporal regression')
    >>> ca_w.inputs.inputspec.subject_rest = '/home/data/subject/func/rest_bandpassed.nii.gz'
    >>> ca_w.inputs.inputspec.subject_timeseries = '/home/data/subject/func/timeseries.txt'
    >>> ca_w.inputs.inputspec.subject_mask = '/home/data/spatialmaps/spatial_map.nii.gz'
    >>> ca_w.inputs.inputspec.demean = True
    >>> ca_w.inputs.inputspec.normalize = True
    >>> ca_w.run() # doctest: +SKIP

    """
    wflow = pe.Workflow(name=wflow_name)

    inputNode = pe.Node(util.IdentityInterface
                        (fields=['subject_rest',
                                 'subject_timeseries',
                                 'subject_mask',
                                 'demean',
                                 'normalize']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface
                         (fields=['component_map',
                                  'component_map_z',
                                  'component_map_z_stack']),
                          name='outputspec')
    
    check_timeseries = pe.Node(util.Function(input_names = ['in_file'],
                                             output_names = ['out_file'],
                                             function = check_ts),
                               name = 'check_timeseries')

    wflow.connect(inputNode, 'subject_timeseries',
                  check_timeseries, 'in_file')
    
    temporalReg = pe.Node(interface=fsl.FSLGLM(),
                          name='temporal_regression')
    
    temporalReg.inputs.output_file = 'component_map.nii.gz'
    temporalReg.inputs.out_z_name = 'component_map_z.nii.gz'

    wflow.connect(inputNode, 'subject_rest',
                  temporalReg, 'in_file')
    wflow.connect(check_timeseries, 'out_file', 
                  temporalReg, 'design_file')
    wflow.connect(inputNode, 'demean',
                  temporalReg, 'demean')
    wflow.connect(inputNode, 'normalize',
                  temporalReg, 'des_norm')
    wflow.connect(inputNode, 'subject_mask',
                  temporalReg, 'mask')


    wflow.connect(temporalReg, 'out_file',
                  outputNode, 'component_map')
    wflow.connect(temporalReg, 'out_z',
                  outputNode, 'component_map_z')

    
    split = pe.Node(interface = fsl.Split(),
                    name = 'split_volumes')
    split.inputs.dimension = 't'
    split.inputs.out_base_name = 'component_map_z_'
    
    wflow.connect(temporalReg, 'out_z',
                  split, 'in_file')

    if which == 'SR':
        wflow.connect(split, 'out_files',
                      outputNode, 'component_map_z_stack')
        
    elif which == 'RT':
        get_roi_order = pe.Node(util.Function(input_names=['maps',
                                                           'timeseries'],
                                              output_names=['labels',
                                                            'maps'],
                                              function= map_to_roi),
                                name='get_roi_order')
        
        wflow.connect(split, 'out_files',
                      get_roi_order, 'maps')
        
        wflow.connect(inputNode, 'subject_timeseries',
                      get_roi_order, 'timeseries')
        
        rename_maps = pe.MapNode(interface = util.Rename(), 
                                 name = 'rename_maps',
                                 iterfield = ['in_file',
                                              'format_string'])
        rename_maps.inputs.keep_ext = True
    
        wflow.connect(get_roi_order, 'labels',
                      rename_maps, 'format_string')
        wflow.connect(get_roi_order, 'maps',
                      rename_maps, 'in_file')
        
        wflow.connect(rename_maps, 'out_file',
                      outputNode, 'component_map_z_stack')
        
    
    return wflow

