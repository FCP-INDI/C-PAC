
# Convert probability threshold value to correlation threshold
def convert_pvalue_to_r(datafile, p_value, two_tailed=False):
    '''
    Method to calculate correlation threshold from p_value

    Parameters
    ----------
    datafile : string
        filepath to dataset to extract number of time pts from
    p_value : float
        significance threshold p-value
    two_tailed : boolean (optional); default=False
        flag to indicate whether to calculate the two-tailed t-test
        threshold for the returned correlation value

    Returns
    -------
    r_value : float
        correlation threshold value 
    '''

    import nibabel as nb
    import numpy as np
    import scipy.stats

    # Get two-tailed distribution
    if two_tailed:
        p_value = p_value / 2

    # Load in data and number of time pts
    img = nb.load(datafile).get_data()
    t_pts = img.shape[-1]

    # N-2 degrees of freedom with Pearson correlation (two sample means)
    deg_freedom = t_pts-2

    # Inverse Survival Function (Inverse of SF)
    # Note: survival function (SF) is also known as the complementary
    # cumulative distribution function (CCDF): F_(x) = p = P(X > x) = 1 - F(x)
    # The inverse will yield: x = F_^-1(p) = F_^-1(P(X > x))
    # where x is a value under the distribution of the random variable X
    # such that the probability of getting greater than x, is p
    t_value = scipy.stats.t.isf(p_value, deg_freedom)
    r_value = np.sqrt(t_value ** 2 / (deg_freedom + t_value ** 2))

    # Return correlation coefficient
    return r_value


# Function to actually do the list merging
def merge_lists(deg_list=[], eig_list=[], lfcd_list=[]):
    merged_list = []
    merged_list.extend(deg_list)
    merged_list.extend(eig_list)
    merged_list.extend(lfcd_list)

    return merged_list


# Separate sub-briks of niftis and save
def sep_nifti_subbriks(nifti_file, out_names):
    '''
    '''
    import os
    import nibabel as nib

    output_niftis = []

    nii_img = nib.load(nifti_file)
    nii_arr = nii_img.get_data()
    nii_affine = nii_img.get_affine()
    nii_dims = nii_arr.shape

    if nii_dims[-1] != len(out_names):
        err_msg = 'out_names must have same number of elements as '\
                  'nifti sub-briks'
        raise Exception(err_msg)

    for brik, out_name in enumerate(out_names):
        brik_arr = nii_arr[:, :, :, 0, brik]
        out_file = os.path.join(os.getcwd(), out_name+'.nii.gz')
        out_img = nib.Nifti1Image(brik_arr, nii_affine)
        out_img.to_filename(out_file)
        output_niftis.append(out_file)

    return output_niftis


# Check centrality parameters
def check_centrality_params(method_option, threshold_option, threshold):
    '''
    Function to check the centrality parameters
    '''

    # Check method option
    if type(method_option) is int:
        if method_option == 0:
            method_option = 'degree'
        elif method_option == 1:
            method_option = 'eigenvector'
        elif method_option == 2:
            method_option = 'lfcd'
        else:
            err_msg = 'Method option: %d not supported' % method_option
            raise Exception(err_msg)
    elif type(method_option) is not str:
        err_msg = 'Method option must be a string, but type: %s provided' \
                  % str(type(method_option))

    # Check threshold option
    if type(threshold_option) is list:
        threshold_option = threshold_option[0]
    if type(threshold_option) is int:
        if threshold_option == 0:
            threshold_option = 'significance'
        elif threshold_option == 1:
            threshold_option = 'sparsity'
        elif threshold_option == 2:
            threshold_option = 'correlation'
        else:
            err_msg = 'Threshold option: %s not supported for network centrality '\
                      'measure: %s; fix this in the pipeline config'\
                      % (str(threshold_option), str(method_option))
            raise Exception(err_msg)
    elif type(threshold_option) is not str:
        err_msg = 'Threshold option must be a string, but type: %s provided' \
                  % str(type(threshold_option))

    # Init lists of acceptable strings
    acceptable_methods = ['degree', 'eigenvector', 'lfcd']
    acceptable_thresholds = ['significance', 'sparsity', 'correlation']

    # Format input strings
    method_option = method_option.lower().replace('centrality', '').rstrip(' ')
    threshold_option = threshold_option.lower().replace('threshold', '').rstrip(' ')

    # Check for strings properly formatted
    if method_option not in acceptable_methods:
        err_msg = 'Method option: %s not supported' % method_option
        raise Exception(err_msg)

    # Check for strings properly formatted
    if threshold_option not in acceptable_thresholds:
        err_msg = 'Threshold option: %s not supported for network centrality '\
                  'measure: %s; fix this in the pipeline config'\
                  % (str(threshold_option), str(method_option))
        raise Exception(err_msg)

    # If it's significance/sparsity thresholding, check for (0,1]
    if threshold_option == 'significance' or threshold_option == 'sparsity':
        if threshold <= 0 or threshold > 1:
            err_msg = 'Threshold value must be a positive number greater than '\
                      '0 and less than or equal to 1.\nCurrently it is set '\
                      'at %f' % threshold
            raise Exception(err_msg)
    # If it's correlation, check for [-1,1]
    elif threshold_option == 'correlation':
        if threshold < -1 or threshold > 1:
            err_msg = 'Threshold value must be greater than or equal to -1 and '\
                      'less than or equal to 1.\n Current it is set at %f'\
                      % threshold
            raise Exception(err_msg)
    else:
        err_msg = 'Threshold option: %s not supported for network centrality '\
                  'measure: %s; fix this in the pipeline config'\
                  % (str(threshold_option), str(method_option))
        raise Exception(err_msg)
    # 
    if method_option == 'lfcd' and threshold_option == 'sparsity':
        err_msg = 'lFCD must use significance or correlation-type '\
                  'thresholding. Check the pipline configuration has '\
                  'this setting'
        raise Exception(err_msg)

    # Return valid method and threshold options
    return method_option, threshold_option
