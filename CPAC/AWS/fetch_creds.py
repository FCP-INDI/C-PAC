# CPAC/AWS/fetch_creds.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module contains functions which return sensitive information from 
a csv file, with regards to connection to AWS services.
'''

# Function to return AWS secure environment variables
def return_aws_keys(creds_path):
    '''
    Method to return AWS access key id and secret access key using
    credentials found in a local file.

    Parameters
    ----------
    creds_path : string (filepath)
        path to the csv file downloaded from AWS; can either be root
        or user credentials

    Returns
    -------
    aws_access_key_id : string
        string of the AWS access key ID
    aws_secret_access_key : string
        string of the AWS secret access key
    '''

    # Init variables
    with open(creds_path, 'r') as creds_in:
        # Grab csv rows
        row1 = creds_in.readline()
        row2 = creds_in.readline()

    # Are they root or user keys
    if 'User Name' in row1:
        # And split out for keys
        aws_access_key_id = row2.split(',')[1]
        aws_secret_access_key = row2.split(',')[2]
    elif 'AWSAccessKeyId' in row1:
        # And split out for keys
        aws_access_key_id = row1.split('=')[1]
        aws_secret_access_key = row2.split('=')[1]
    else:
        err_msg = 'Credentials file not recognized, check file is correct'
        raise Exception(err_msg)

    # Strip any carriage return/line feeds
    aws_access_key_id = aws_access_key_id.replace('\r', '').replace('\n', '')
    aws_secret_access_key = aws_secret_access_key.replace('\r', '').replace('\n', '')

    # Return keys
    return aws_access_key_id, aws_secret_access_key


# Function to return an AWS S3 bucket
def return_bucket(creds_path, bucket_name):
    '''
    Method to a return a bucket object which can be used to interact
    with an AWS S3 bucket using credentials found in a local file.

    Parameters
    ----------
    creds_path : string (filepath)
        path to the csv file with 'Access Key Id' as the header and the
        corresponding ASCII text for the key underneath; same with the
        'Secret Access Key' string and ASCII text
    bucket_name : string
        string corresponding to the name of the bucket on S3

    Returns
    -------
    bucket : boto.s3.bucket.Bucket
        a boto s3 Bucket object which is used to interact with files
        in an S3 bucket on AWS
    '''

    # Import packages
    try:
        import boto3
        import botocore
    except ImportError as exc:
        err_msg = 'Boto3 package is not installed - install boto3 and '\
                  'try again.'
        raise Exception(err_msg)

    # Try and get AWS credentials if a creds_path is specified
    if creds_path:
        try:
            aws_access_key_id, aws_secret_access_key = \
                return_aws_keys(creds_path)
        except Exception as exc:
            err_msg = 'There was a problem extracting the AWS credentials '\
                      'from the credentials file provided: %s. Error:\n%s'\
                      % (creds_path, exc)
            raise Exception(err_msg)
        # Init connection
        print 'Connecting to S3 bucket: %s with credentials from %s ...'\
              % (bucket_name, creds_path)
        # Use individual session for each instance of DataSink
        # Better when datasinks are being used in multi-threading, see:
        # http://boto3.readthedocs.org/en/latest/guide/resources.html#multithreading
        session = boto3.session.Session(aws_access_key_id=aws_access_key_id,
                                        aws_secret_access_key=aws_secret_access_key)
        s3_resource = session.resource('s3', use_ssl=True)

    # Otherwise, connect anonymously
    else:
        print 'Connecting to AWS: %s anonymously...' % bucket_name
        session = boto3.session.Session()
        s3_resource = session.resource('s3', use_ssl=True)
        s3_resource.meta.client.meta.events.register('choose-signer.s3.*',
                                                     botocore.handlers.disable_signing)

    # Explicitly declare a secure SSL connection for bucket object
    bucket = s3_resource.Bucket(bucket_name)

    # And try fetch the bucket with the name argument
    try:
        s3_resource.meta.client.head_bucket(Bucket=bucket_name)
    except botocore.exceptions.ClientError as exc:
        error_code = int(exc.response['Error']['Code'])
        if error_code == 403:
            err_msg = 'Access to bucket: "%s" is denied; check credentials'\
                      % bucket_name
            raise Exception(err_msg)
        elif error_code == 404:
            err_msg = 'Bucket: "%s" does not exist; check spelling and try '\
                      'again' % bucket_name
            raise Exception(err_msg)
        else:
            err_msg = 'Unable to connect to bucket: "%s". Error message:\n%s'\
                      % (bucket_name, exc)
    except Exception as exc:
        err_msg = 'Unable to connect to bucket: "%s". Error message:\n%s'\
                  % (bucket_name, exc)
        raise Exception(err_msg)

    # Return the bucket
    return bucket


# Function to return a RDS cursor
def return_cursor(creds_path):
    '''
    Method to return an Oracle DB cursor which is connected to database
    instance using credentials found in a local file.

    Parameters
    ----------
    creds_path : string (filepath)
        path to the csv file with 'DB_USER' as the header and the
        corresponding ASCII text for the user name underneath; same
        'DB_PASSWD', 'DB_HOST', 'DB_PORT', and 'DB_SID'headers and text

    Returns
    -------
    cursor : OracleCursor
        a cx_Oracle cursor object which is used to query and modify an
        Oracle database
    '''

    # Import packages
    import cx_Oracle

    # Get DB connection variables
    user, passwd, host, port, sid = return_rds_vars(creds_path)
    # Setup connection and cursor
    dsn = cx_Oracle.makedsn(host, port, sid)
    conn = cx_Oracle.connect(user, passwd, dsn)
    cursor = conn.cursor()

    # Return cursor
    return cursor


# Function to return RDS secure environment variables
def return_rds_vars(creds_path):
    '''
    Method to return user, password, host, port, and SID for an AWS RDS
    Oracle DB instance using credentials found in a local file.

    Parameters
    ----------
    creds_path : string (filepath)
        path to the csv file with 'db_user=' followed by oracle db
        username in the first row, 'db_passwd'= followed by oracle db
        password in the second row, 'db_host=' followed by oracle db
        hostname in the third row, 'db_port=' followed by oracle db
        portname in the fourth row, and 'db_sid' followed by oracle_db
        serviceID in the fifth row

    Returns
    -------
    db_user : string
        string of the username for the database connection
    db_passwd : string
        string of the password for the database connection
    db_host : string
        string of the host for the database connection
    db_port : string
        string of the port for the database connection
    db_sid : string
        string of the SID for the database connection
    '''

    # Import packages
    import csv

    # Init variables
    csv_reader = csv.reader(open(creds_path, 'r'))

    # Grab csv rows
    row1 = csv_reader.next()[0]
    row2 = csv_reader.next()[0]
    row3 = csv_reader.next()[0]
    row4 = csv_reader.next()[0]
    row5 = csv_reader.next()[0]

    # Get database credentials
    db_user = row1.split('=')[1]
    db_passwd = row2.split('=')[1]
    db_host = row3.split('=')[1]
    db_port = row4.split('=')[1]
    db_sid = row5.split('=')[1]

    # Return the DB variables
    return db_user,\
           db_passwd,\
           db_host,\
           db_port,\
           db_sid
