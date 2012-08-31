

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



def get_group_analysis_inputs(c, sublist):
    
    modelist = []

    try:
        
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

        models = os.listdir(c.modelsDirectory)

        for model in models:
            model = model.rstrip('\r\n')
            modelist.append(model)

        print 'checking for models ', modelist

        if len(modelist) < 1:
            raise EmptyListError("Model List is Empty")

        c.dervTemplateList = [sublist] + ['derivative']
        

#       c.FTest True if f_test need to perform
#       F-tests enables user to investigate several contrasts at the same time. 
#       Also, it enables to compare the contribution of each contrast to the 
#       model and decide on significant and non-significant ones.       
        
        template_dict = dict(mat=os.path.join(c.modelsDirectory, c.mat),
                             con=os.path.join(c.modelsDirectory, c.con),
                             grp=os.path.join(c.modelsDirectory, c.grp),
                             derv=dervTemplate)
                    
        template_args = dict(mat=[c.matTemplateList],
                             con=[c.conTemplateList],
                             grp=[c.grpTemplateList],
                             derv=[c.dervTemplateList])
                        
        if c.fTest:
            template_dict['fts'] = os.path.join(c.modelsDirectory, c.fts)
            template_args['fts'] = [c.grpTemplateList]
            
        
        print 'checking for template_dict', template_dict

        print 'checking for template_args', template_args

        return  modelist, template_dict, template_args

    except EmptyListError as e:
        print e.msg

    else:
        print "Error while preparing templates for group analysis"
        raise

