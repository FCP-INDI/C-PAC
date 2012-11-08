
#list of subjects that are to be excluded 
#exclusionSubjectList = '/home/data/settings/ex_subjects.txt'
exclusionSubjectList = None

#list of subjects that are included 
#if None extract data runs on all the subjects
subjectList = None
#subjectList = '/home/data/settings/include_subjects.txt'

#Put %s where site and subjects are in the path
#anatomicalTemplate = '/home/data/Incoming/cambridge_fcon/%s/%s/*/mprage_anonymized.nii.gz'
anatomicalTemplate = '/Users/sikkas01/sites/%s/%s/session_1/*/mprage.nii.gz'

#Functional Path
#Put  %s where site and subjects are in the path
#functionalTemplate = '/home/data/Incoming/cambridge_fcon/%s/%s/*/rest.nii.gz'
functionalTemplate = '/Users/sikkas01/sites/%s/%s/session_1/*/rest.nii.gz'

#list of sites
#if None extract data runs on all sites
siteList = None

#scan parameters csv file path
#This file is mandatory for slice timing correction
scanParametersCSV = None
#scanParametersCSV = /home/data/settings/scan_parameters.csv
