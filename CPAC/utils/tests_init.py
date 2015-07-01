# CPAC/utils/tests_init.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module contains functions that assist in initializing CPAC
tests resources
'''

# Return tests data config file
def return_data_config():
    '''
    Function to read in the template data config from the
    CPAC_RESOURCE_DIR and populate it with actual filepaths

    Parameters
    ----------
    None

    Returns
    -------
    data_config : string
        filepath to the newly written data config yaml file
    '''

    # Import packages
    import os

    # Init variables
    resource_dir = return_resource_dir()
    templates_dir = return_resource_subfolder('templates')
    config_dir = return_resource_subfolder('config')

    # Create template data config and store in config_dir
    data_config_template = os.path.join(templates_dir,
                                        'data_config_template.yml')
    data_config = os.path.join(config_dir, 'data_config_test.yml')

    # Open the files
    tmp_f = open(data_config_template, 'r')
    res_f = open(data_config, 'w')

    # Replace 'RESOURCE_DIR' string with actual directory
    for line in tmp_f:
        res_f.write(line.replace('RESOURCE_DIR', resource_dir))
    tmp_f.close()
    res_f.close()

    # Return filepath
    return data_config


# Return tests pipeline config file
def return_pipeline_config():
    '''
    Function to return the file path of the pipeline config file from
    the CPAC_RESOURCE_DIR

    Parameters
    ----------
    None

    Returns
    -------
    pipeline_config : string
        filepath to the pipeline config yaml file
    '''

    # Import packages
    import os

    # Init variables
    config_dir = return_resource_subfolder('config')

    # Get pipeline config
    pipeline_config = os.path.join(config_dir, 'pipeline_config_test.yml')

    # Return filepath
    return pipeline_config


# Return tests subject list
def return_subject_list():
    '''
    Function to return the file path of the subject list file from
    the CPAC_RESOURCE_DIR

    Parameters
    ----------
    None

    Returns
    -------
    subject_list : string
        filepath to the subject list yaml file
    '''

    # Import packages
    import os

    # Init variables
    config_dir = return_resource_subfolder('config')

    # Get sublist
    subject_list = os.path.join(config_dir, 'CPAC_subject_list_test.yml')

    # Return filepath
    return subject_list


# Return test strategies obj file
def return_strats_obj():
    '''
    Function to return the file path of the strategies obj file from
    the CPAC_RESOURCE_DIR

    Parameters
    ----------
    None

    Returns
    -------
    strats_obj : string
        filepath to the strategies obj file
    '''

    # Import packages
    import os

    # Init variables
    settings_dir = return_resource_subfolder('settings')

    # Get strategies obj
    strats_obj = os.path.join(settings_dir, 'strategies_test.obj')

    # Return filepath
    return strats_obj


# Look for CPAC_RESOURCE_DIR to be in environment
def return_resource_dir():
    '''
    Function to return the filepath of the CPAC_RESOURCE_DIR; note the
    CPAC_RESOURCE_DIR environment variable must be set

    Parameters
    ----------
    None

    Returns
    -------
    resource_dir : string
        the file path on disk where the cpac resources folder is
    '''
    # Import packages
    import os

    # Init variables
    resource_dir = os.getenv('CPAC_RESOURCE_DIR')

    # Check if set
    if not resource_dir:
        err_msg = 'CPAC_RESOURCE_DIR environment variable not set!\n' \
                  'Set this to the directory of the cpac_resources folder'
        raise Exception(err_msg)
    else:
        return resource_dir


# Return any subfolder of the resource directory
def return_resource_subfolder(subfolder):
    '''
    Funnction to return subfolders of the CPAC_RESOURCE_DIR

    Parameters
    ----------
    subfolder : string
        subfolder name to return path of

    Returns
    -------
    resource_subfolder : string
        filepath to the resource subfolder
    '''

    # Import packages
    import os

    # Init variables
    resource_dir = return_resource_dir()

    # Check if its a sub-subfolder
    if subfolder == 'config' or subfolder == 'creds' or \
       subfolder == 'templates':
        resource_subfolder = os.path.join(resource_dir, 'settings', subfolder)
    else:
        resource_subfolder = os.path.join(resource_dir, subfolder)

    # Return subfolder
    return resource_subfolder


# Get the AWS credentials
def return_aws_creds():
    '''
    Function to return the AWS credentials files located in the
    CPAC_RESOURCE_DIR

    Parameters
    ----------
    None

    Returns
    -------
    aws_creds : string
        filepath to the AWS credentials with access key id and secret
        access key
    db_creds: string
        filepath to the AWS database credentials with the username
        and password along with other details to connect to database
    '''

    # Import packages
    import os

    # Init variables
    creds_dir = return_resource_subfolder('creds')

    # Get credentials path
    aws_creds = os.path.join(creds_dir, 'aws_creds.csv')
    db_creds = os.path.join(creds_dir, 'db_creds.csv')

    # Return the aws creds and databse creds
    return aws_creds, db_creds


# Get the default test bucket name
def return_bucket_name():
    '''
    Function to return the default S3 bucket name used in test suite

    Parameters
    ----------
    None

    Returns
    -------
    bucket_name : string
        default S3 bucket name for testing
    '''

    # Set default bucket name
    bucket_name = 'fcp-indi'

    # Return bucket name
    return bucket_name
