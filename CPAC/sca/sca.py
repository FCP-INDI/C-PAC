import sys
from CPAC.interfaces.afni import preprocess
import os
import commands
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



    ## 2. Compute voxel-wise correlation with Seed Timeseries
    corr = pe.Node(interface=preprocess.ThreedTcorrOneD(),
                      name='corr')
    corr.inputs.options = "-pearson"
    corr.inputs.outputtype = 'NIFTI_GZ'

    ## 3. Z-transform correlations
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
