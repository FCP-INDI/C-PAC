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
        path to the csv file with 'AWSAccessKeyId=' followed by access
        key in the first row and 'AWSSecretAccessKey=' followed by
        secret access key in the second row

    Returns
    -------
    aws_access_key_id : string
        string of the AWS access key ID
    aws_secret_access_key : string
        string of the AWS secret access key
    '''

    # Import packages
    import csv

    # Init variables
    csv_reader = csv.reader(open(creds_path, 'r'))
    
    # Grab csv rows
    row1 = csv_reader.next()[0]
    row2 = csv_reader.next()[0]
    
    # And split out for keys
    aws_access_key_id = row1.split('=')[1]
    aws_secret_access_key = row2.split('=')[1]

    # Return keys
    return aws_access_key_id,\
           aws_secret_access_key,\


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
    import boto
    import boto.s3.connection

    # Get AWS credentials
    aws_access_key_id, aws_secret_access_key = return_aws_keys(creds_path)

    # Init connection
    cf = boto.s3.connection.OrdinaryCallingFormat()
    s3_conn = boto.connect_s3(aws_access_key_id, aws_secret_access_key,
                              calling_format=cf)
    # And fetch the bucket with the name argument
    bucket = s3_conn.get_bucket(bucket_name)

    # Return bucket
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