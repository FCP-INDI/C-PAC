

def safe_shape(*vol_data):
    """
    Checks if the volume (first three dimensions) of multiple ndarrays
    are the same shape.
    
    Parameters
    ----------
    vol_data0, vol_data1, ..., vol_datan : ndarray
        Volumes to check
    
    Returns
    -------
    same_volume : bool
        True only if all volumes have the same shape.
    """
    same_volume = True
    
    first_vol_shape = vol_data[0].shape[:3]
    for vol in vol_data[1:]:
        same_volume &= (first_vol_shape == vol.shape[:3]) 
        
    return same_volume


def extract_one_d(list_timeseries):


    for timeseries in list_timeseries:

        if '1D' in timeseries:


            return timeseries

        else:

            print "Error : ROI/Voxel TimeSeries 1D file not found"

            return None

def set_gauss(fwhm):

    op_string = ""

    fwhm = float(fwhm)

    sigma = float(fwhm / 2.3548)

    op = "-kernel gauss %f -fmean -mas " % (sigma) + "%s"
    op_string = op

    return op_string


def modify_model_files(model_files, sublist):

    def read_model_file(file):

        dict1 ={}
        f = open(file, 'r')
        for line in open(file, 'r'):
            if line.strip() !='':
                if 'NumWaves' in line: 
                    dict1['/NumWaves']= line.split()[1]
                elif 'NumPoints' in line:
                    dict1['/NumPoints'] = line.split()[1]
                elif 'PPheights' in line:
                    dict1['/PPHeights'] = line.split()[1:]
                elif 'Matrix' in line:
                    dict1['/Matrix'] = []
                else:
                    dict1.get('/Matrix').append(line)
        return dict1

