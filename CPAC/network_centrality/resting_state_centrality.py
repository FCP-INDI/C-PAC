
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from CPAC.network_centrality import *

def create_resting_state_graphs(generate_graph = False, wf_name = 'resting_state_graph'):
    """
    Workflow to calculate degree and eigenvector centrality
    measures for the resting state data.
    
    Parameters
    ----------
    generate_graph : boolean
        when true the workflow plots the adjacency matrix graph 
        and converts the adjacency matrix into compress sparse 
        matrix and stores it in a .mat file. By default its False
    wf_name : string
        name of the workflow
        
    Returns 
    -------
    wf : workflow object
        resting state graph workflow object
          
    Notes
    -----
    
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/network_centrality/resting_state_centrality.py>`_
    
    Workflow Inputs::
    
        inputspec.subject: string (nifti file)
            path to resting state input data for which centrality measure is to be calculated
            
        inputspec.template : string (existing nifti file)
            path to mask/parcellation unit 
        
        inputspec.threshold_option: string (int)
            threshold options:  0 for probability p_value, 1 for sparsity threshold, any other for threshold value
        
        inputspec.threshold: string (float)
            pvalue/sparsity_threshold/threshold value
           
        centrality_options.weight_options : string (list of boolean)
            list of two booleans for binarize and weighted options respectively
        
        centrality_options.method_options : string (list of boolean)
            list of two booleans for Degree and Eigenvector centrality method options respectively
        
    Workflow Outputs::
    
        outputspec.centrality_outputs : string (list of nifti files)
            path to list of centrality outputs for binarized or/and weighted and
            degree or/and eigen_vector 
        
        outputspec.threshold_matrix : string (numpy file)
            path to file containing thresholded correlation matrix
        
        outputspec.correlation_matrix : string (numpy file)
            path to file containing correlation matrix
        
        outputspec.graph_outputs : string (mat and png files)
            path to matlab compatible sparse adjacency matrix files 
            and adjacency graph images 
    
    Order of commands:
    
    - load the data and template, based on template type (parcellation unit ar mask)
      extract timeseries
    
    - Calculate the correlation matrix for the image data for each voxel in the mask or node
      in the parcellation unit
    
    - Based on threshold option (p_value or sparsity_threshold), calculate the threshold value
    
    - Threshold the correlation matrix
     
    - Based on weight options for edges in the network (binarize or weighted), calculate Degree 
      or Vector Based centrality measures
     
    
    High Level Workflow Graph:
    
    .. image:: ../images/resting_state_centrality.dot.png
       :width: 1000
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/resting_state_centrality_detailed.dot.png
       :width: 1000
    
    Examples
    --------
    
    >>> import resting_state_centrality as graph
    >>> wflow = graph.create_resting_state_graphs()
    >>> wflow.inputs.centrality_options.method_options=[True, True]
    >>> wflow.inputs.centrality_options.weight_options=[True, True]
    >>> wflow.inputs.inputspec.subject = '/home/work/data/rest_mc_MNI_TR_3mm.nii.gz'
    >>> wflow.inputs.inputspec.template = '/home/work/data/mask_3mm.nii.gz'
    >>> wflow.inputs.inputspec.threshold_option = 1
    >>> wflow.inputs.inputspec.threshold = 0.0744
    >>> wflow.base_dir = 'graph_working_directory'
    >>> wflow.run()
    
    """
    wf = pe.Workflow(name = wf_name)
    
    inputspec = pe.Node(util.IdentityInterface(fields=['subject',
                                                       'template',
                                                       'threshold',
                                                       'threshold_option']),
                        name='inputspec')
    
    outputspec = pe.Node(util.IdentityInterface(fields=['centrality_outputs',
                                                        'threshold_matrix',
                                                        'correlation_matrix',
                                                        'graph_outputs']),
                         name = 'outputspec')
    
    centrality_options = pe.Node(util.IdentityInterface(fields = ['weight_options',
                                                                  'method_options']),
                                 name = 'centrality_options')
    
   
    read_data = pe.Node(util.Function( input_names = ['datafile', 
                                                      'template'],
                                       output_names= ['timeseries_data', 
                                                      'affine', 
                                                      'mask_data',
                                                      'mask_type',
                                                      'scans'],
                                      function = load),
                        name = 'read_data')
    
    wf.connect(inputspec, 'subject', 
               read_data, 'datafile')
    wf.connect(inputspec, 'template', 
               read_data, 'template')
    
    
    find_correlation = pe.Node(util.Function(input_names = ['timeseries_data'],
                                        output_names = ['correlation_matrix'],
                                        function = calculate_correlation ),
                          name = 'find_correlation')
    
    wf.connect(read_data, 'timeseries_data', 
               find_correlation, 'timeseries_data')
    
    wf.connect(find_correlation, 'correlation_matrix',
               outputspec, 'correlation_matrix')
    
    threshold_correlation = pe.Node(util.Function(input_names = ['corr_matrix',  
                                                                 'option', 
                                                                 'threshold',
                                                                 'scans'],
                                                output_names = ['threshold_matrix'],
                                                function = threshold_rmatrix),
                                  name = 'threshold_correlation')
    
    wf.connect(read_data, 'scans', 
               threshold_correlation, 'scans')
    wf.connect(find_correlation, 'correlation_matrix',
               threshold_correlation, 'corr_matrix')
    wf.connect(inputspec, 'threshold', 
               threshold_correlation, 'threshold')
    wf.connect(inputspec, 'threshold_option', 
               threshold_correlation, 'option')
    
    wf.connect(threshold_correlation, 'threshold_matrix',
               outputspec, 'threshold_matrix')
    
    calculate_centrality = pe.Node(util.Function(input_names = ['threshold_matrix',
                                                                'correlation_matrix',
                                                                'template_data',
                                                                'affine',
                                                                'weight_options',
                                                                'method_options',
                                                                'template_type'],
                                                 output_names = ['out_list'],
                                                 function = get_centrality),
                                   name = 'calculate_centrality')
    
    wf.connect(centrality_options, 'method_options',
               calculate_centrality, 'method_options')
    wf.connect(centrality_options, 'weight_options',
               calculate_centrality, 'weight_options')
    wf.connect(find_correlation, 'correlation_matrix',
               calculate_centrality, 'correlation_matrix')
    wf.connect(threshold_correlation, 'threshold_matrix',
               calculate_centrality, 'threshold_matrix')
    wf.connect(read_data, 'mask_data',
                calculate_centrality, 'template_data')
    wf.connect(read_data, 'affine',
               calculate_centrality, 'affine')
    wf.connect(read_data, 'mask_type',
               calculate_centrality, 'template_type')
    
    wf.connect(calculate_centrality, 'out_list',
               outputspec, 'centrality_outputs')
    
    if generate_graph:
    
        generate_graph = pe.Node(util.Function(input_names = ['threshold_matrix',
                                                              'correlation_matrix',
                                                              'weight_options',
                                                              'template_data',
                                                              'template_type'],
                                                  output_names= ['out_list'],
                                                  function= generate_adjacency_graph),
                                    name = 'generate_graph')
        
        
        wf.connect(centrality_options, 'weight_options',
                   generate_graph, 'weight_options')
        wf.connect(find_correlation, 'correlation_matrix',
                   generate_graph, 'correlation_matrix')
        wf.connect(threshold_correlation, 'threshold_matrix',
                   generate_graph, 'threshold_matrix')
        wf.connect(read_data, 'mask_data',
                   generate_graph, 'template_data')
        wf.connect(read_data, 'mask_type',
                   generate_graph, 'template_type')
        
        wf.connect(generate_graph, 'out_list',
                   outputspec, 'graph_outputs')
    
    return wf




def load(datafile, template):
    
    """
    Method to read data from datafile and mask/parcellation unit
    and store the mask data, timeseries, affine matrix, mask type
    and scans. The output of this method is used by all other nodes.
    
    Parameters
    ----------
    datafile : string (nifti file)
        path to subject data file
    template : string (nifti file)
        path to mask/parcellation unit
        
    Returns
    -------
    timeseries_data: string (numpy npy file)
        path to file containing timeseries of the input data
    affine: string (numpy npy file)
        path to file containing affine matrix of teh input data
    mask_data: string (numpy npy file)
        path to file containing mask/parcellation unit matrix
    template_type: string 
        0 for mask, 1 for parcellation unit 
    scans: string (int)
        total no of scans in the input data
        
    Raises
    ------
    Exception
    """

    import os
    import nibabel as nib
    import numpy as np

    
    try:    
        if isinstance(datafile, list):
            img = nib.load(datafile[0])
        else:
            img = nib.load(datafile) 
        
        data = img.get_data()
        aff = img.get_affine()    
        mask = nib.load(template).get_data()   
        scans = data.shape[3]
        
    except:
        print "Error in loading images for graphs"
        raise
    
    
    
    if mask.shape != data.shape[:3]:
        raise Exception("Invalid Shape Error. mask and data file have"\
                        "different shape please check the voxel size of the two files")

    
    #check for parcellation
    nodes = np.unique(mask).tolist()
    nodes.sort()
    print "sorted nodes", nodes
    
    #extract timeseries
    if len(nodes)>2:
        flag=1
        for n in nodes:
            if n > 0:
                node_array = data[mask == n]
                avg = np.mean(node_array, axis =0)
                if flag:
                    timeseries = avg
                    flag=0
                else:
                    timeseries = np.vstack((timeseries, avg))
        #template_type is 1 for parcellation
        template_type = 1
    else:
        #template_type is 0 for mask
        template_type = 0
        mask = mask.astype('bool')
        timeseries = data[mask]
        
    #saving files
    img_name = os.path.splitext(os.path.splitext(os.path.basename(datafile))[0])[0] + '.npy'
    mask_name = os.path.splitext(os.path.splitext(os.path.basename(template))[0])[0]  + '.npy'
    
    mask_data = os.path.join(os.getcwd(), mask_name)
    timeseries_data = os.path.join(os.getcwd(),img_name)
    affine = os.path.join(os.getcwd(),'affine.npy')
        
    np.save(mask_data, mask)
    np.save(timeseries_data, timeseries)
    np.save(affine, aff)
    
    return timeseries_data, affine, mask_data, template_type, scans



def calculate_correlation(timeseries_data):
    """
    Method to calculate correlation between 
    each voxel or node of data present in the 
    template
    
    Parameters
    ----------
    timeseries_data : string (numpy matrix file)
        Path to file containing data matrix
    
    Returns
    -------
    corr_file : string (mat file)
        path to file containing the correlation matrix
    
    """
    
    import os
    import numpy as np
    from CPAC.network_centrality import load_mat
    
    timeseries = load_mat(timeseries_data)
    r_matrix = np.corrcoef(timeseries)
    cwd = os.getcwd()
    
    print "shape of correlation matrix", r_matrix.shape
    
    corr_mat_file = os.path.join(cwd, 'r_matrix.npy')
    np.save(corr_mat_file,r_matrix)
    
    return corr_mat_file



def threshold_rmatrix(corr_matrix, option, 
                      threshold, scans):
    
    """
    Method to threshold the correaltion matrix based on 
    any of the two threshold options- sparsity, probability
    or by simply providing correlation threshold. it is two 
    step process, first calculate the threshold and then apply
    the threshold to the correlation matrix.
    
    Parameters
    ----------
    corr_matrix : string (numpy npy file)
        patht o file containing correlation matrix
    option : string (int)
        list of threshold option: 0 for pvalue, 1 for sparsity, 
        any other for simply correlation threshold 
    threshold : string (float)
        pvalue/sparsity_threshold/correaltion_threshold
    scans : string (int)
        Total number of scans in input data
        
    Returns
    -------
    threshold_file : string (numpy npy file)
        file containing threshold correlation matrix
    
    Raises
    ------
    Exception
    """
    
    import numpy as np
    import os
    from CPAC.network_centrality import load_mat,\
                                        convert_pvalue_to_r,\
                                        convert_sparsity_to_r
    
    try:
        r_matrix = load_mat(corr_matrix)
       
        print "threshold_option -->", option
        
        try:
            if option == 0:
                r_value = convert_pvalue_to_r(scans, threshold)
            elif option == 1:
                r_value = convert_sparsity_to_r(r_matrix, threshold)
            else:
                r_value = threshold
        except:
            print "Exception in calculating thresold value"
            raise
        
        print "correlation threshold value -> ", r_value
        print "thresholding the correlation matrix...."
        
        threshold_matrix = r_matrix > r_value
        
        threshold_file = os.path.join(os.getcwd(), 'threshold_matrix.npy')
        np.save(threshold_file, threshold_matrix.astype(np.float))
    
    except Exception:
        print "Exception while thresholding correlation matrix"
        raise
    
    return threshold_file




def get_centrality(weight_options, method_options, 
                   threshold_matrix, correlation_matrix,
                   template_data, affine, template_type):
    
    """
    Method to calculate centrality measures
    for resting state data. It is two step process
    first calculate the centrality and then map the 
    centrality matrix to a nifti file.
    
    Parameters
    ----------
    weight_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    method_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    threshold_matrix : string (numpy npy file)
        path to file containing thresholded correlation matrix 
    correlation_matrix : string (numpy npy file)
        path to file containing correlation matrix
    template_data : string (numpy npy file)
        path to file containing mask or parcellation unit data
    affine : string (numpy npy file)
        path to file containing affine matrix of input data
    template_type : string 
      o for mask, 1 for parcellation_unit
    
    Returns
    -------
    out_list : string (list of files)
        nifti images for each input weight options and 
        centrality methods
    
    Raises
    ------
    Exception
    """
    from CPAC.network_centrality import map_centrality_matrix,\
                                        get_centrality_matrix
    
    out_list =[]
    
    try:
    
        if not method_options or not weight_options :
            raise Exception("Invalid options, method option"\
                            "and weight options are required") 
    
        if True not in method_options or True not in weight_options:
            raise Exception("Invalid method or weight options"\
                             "atleast one True value is required")
    
        def get_image(matrix, template_type):
            centrality_image = map_centrality_matrix(matrix, 
                                                    affine, 
                                                    template_data,
                                                    template_type)
            out_list.append(centrality_image)
        
        
        centrality_matrix = get_centrality_matrix(threshold_matrix, 
                                                  correlation_matrix, 
                                                  weight_options,
                                                  method_options)
        
        for mat in centrality_matrix:
                get_image(mat, template_type)
        
    except Exception:
        print "Exception in calculating centrality measures"
        raise
         
    return out_list

    
def generate_adjacency_graph(correlation_matrix, threshold_matrix, 
                             weight_options, template_data,template_type):
    """
    Method to store the adjacency matrix as a compress sparse matrix which
    can be loaded inot matlab. The method also create a png image of the 
    graph.
    
    Parameters
    ----------
    correlation_matrix : string (numpy matrix file)
        path to correlation matrix file
    threshold_matrix : string (numpy matrix file)
        path to thresholded correlation matrix file
    weight_options: boolean
        True for weighted and False for binarize
    template_data : string (numpy matrix file)
        path to file containing parcellation unit
    template_type : string
        0 for mask, 1 for parcellation unit
        
    Returns
    -------
    out_matrix : string (mat file)
        compressed sparse matrix
    out_img : string (png file)
        path to graph image file
            
    """
    from pylab import imsave
    import os
    from scipy.sparse import lil_matrix, csc_matrix
    from scipy.io import savemat
    from CPAC.network_centrality import load_mat
    
    out_list =[]
    
    #if int(template_type)==1:
    
    thresh_matrix = load_mat(threshold_matrix)
    corr_matrix = load_mat(correlation_matrix)
    
    if isinstance(template_data, list):
        mask_name = os.path.splitext(os.path.basename(template_data[0]))[0]
    else:
        mask_name = os.path.splitext(os.path.basename(template_data))[0]
    
    
    def save(filename, key, matrix):
        import os
        
        out_matrix = os.path.join(os.getcwd(),filename + ".mat")
        out_img = os.path.join(os.getcwd(), filename + ".png")
        savemat(out_matrix, {key: matrix})
        print out_matrix
        print out_img
        out_list.append(out_matrix)
        imsave(out_img, matrix.todense())
        out_list.append(out_img)
    
    try:
        
        if weight_options[0]:
            spedgemat = lil_matrix(thresh_matrix)
            spcscmat = csc_matrix(spedgemat)
            del spedgemat
            filename = mask_name + "_adjacency_matrix"
            save(filename, 'unit_graph', spcscmat)
        
        elif weight_options[1]:
            matrix = thresh_matrix * corr_matrix
            spedgemat = lil_matrix (matrix)
            spcscmat = csc_matrix (spedgemat)
            del spedgemat
            filename = mask_name +"_weighted_adjacency_matrix"
            save(filename,'unit_graph', spcscmat)
            
    except Exception:
        print "Not Enough Memory available to generate sparse matrix"
        raise
    #else:
    #    print "No voxel based sparse matrix generation, matrix size is too huge"
        
    
    return out_list
    
    
    
    
    
    
    