
def get_cent_zscore(wf_name = 'z_score'):
    
    """
    Workflow to calculate z-scores
    
    Parameters
    ----------
    wf_name : string
        name of the workflow
        
    Returns
    -------
    wf : workflow object
    
    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/network_centrality/z_score.py>`_
    
    
    Workflow Inputs::
        
        inputspec.input_file : string
            path to input functional derivative file for which z score has to be calculated
        inputspec.mask_file : string
            path to whole brain functional mask file required to calculate zscore
    
    Workflow Outputs::
        
        outputspec.z_score_img : string
             path to image containing Normalized Input Image Z scores across full brain.
    
    High Level Workflow Graph:
    
    .. image:: ../images/zscore.dot.png
       :width: 500
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/zscore_detailed.dot.png
       :width: 500
    
    Example
    -------
    >>> import get_zscore as z
    >>> wf = z.get_zscore()
    >>> wf.inputs.inputspec.input_file = '/home/data/graph_working_dir/calculate_centrality/degree_centrality_binarize.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/graphs/GraphGeneration/new_mask_3m.nii.gz'
    >>> wf.run()
    
    """
    
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl as fsl
    
    wflow = pe.Workflow(name = wf_name)
    
    inputNode = pe.Node(util.IdentityInterface(fields=['input_file',
                                                       'mask_file']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['z_score_img']),
                          name='outputspec')
    
    mean = pe.MapNode(interface=fsl.ImageStats(),
                      name='mean',
                      iterfield=['in_file'])
    mean.inputs.op_string = '-k %s -m'
    
    wflow.connect(inputNode, 'input_file',
                  mean, 'in_file')
    wflow.connect(inputNode, 'mask_file',
                  mean, 'mask_file')

    standard_deviation = pe.MapNode(interface=fsl.ImageStats(),
                                 name='standard_deviation',
                                 iterfield=['in_file'])
    standard_deviation.inputs.op_string = '-k %s -s'

    wflow.connect(inputNode, 'input_file',
                  standard_deviation, 'in_file')
    wflow.connect(inputNode, 'mask_file',
                  standard_deviation, 'mask_file')
    
    op_string = pe.MapNode(util.Function(input_names=['mean',
                                                   'std_dev'],
                                      output_names=['op_string'],
                           function=get_operand_string),
                           name='op_string',
                           iterfield=['mean',
                                      'std_dev'])
    
    wflow.connect(mean, 'out_stat',
                  op_string, 'mean')
    wflow.connect(standard_deviation, 'out_stat',
                 op_string, 'std_dev')
    
    z_score = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='z_score',
                        iterfield=['in_file',
                                   'op_string'])
    
    wflow.connect(op_string, 'op_string',
                  z_score, 'op_string')
    wflow.connect(inputNode, 'input_file',
                  z_score, 'in_file')
    wflow.connect(inputNode, 'mask_file',
                  z_score, 'operand_files')
    
    wflow.connect(z_score, 'out_file',
                  outputNode, 'z_score_img')
    
    return wflow
    
    
def get_operand_string(mean, std_dev):
    """
    Method to get operand string for Fsl Maths
    
    Parameters
    ----------
    mean : string
        path to img containing mean
    std_dev : string
        path to img containing standard deviation
    
    Returns
    ------
    op_string : string
        operand string
    """
    
    str1 = "-sub %f -div %f" % (float(mean), float(std_dev))
    op_string = str1 + " -mas %s"
    return op_string
