"""
list of subjects that are to be excluded, can be a text file or a list
"""
#exclusionSubjectList = '/home/data/settings/ex_subjects.txt'
#exclusionSubjectList = ['sub001', 'sub003']
exclusionSubjectList = None


"""
list of subjects that are included, can be a text file or a list
if None, extract data runs on all the subjects
"""
#subjectList = '/home/data/settings/include_subjects.txt'
#subjectList = ['sub002', 'sub003']
subjectList = None


"""
Put %s where site and subjects are in the path
"""
anatomicalTemplate = '/path/to/data/%s/%s/session_*/anat_*/mprage.nii.gz'


"""
Functional Path
Put %s where site and subjects are in the path
"""
functionalTemplate = '/path/to/data/%s/%s/session_*/rest_*/rest.nii.gz'


"""
list of sites, can be a text file or a list
if None, extract data runs on all sites
"""
#siteList = ['ABIDE', 'ADHD-200']
siteList = None


"""
Scan Parameters csv file path. This file is mandatory for 
slice timing correction. please use the right format for the csv, 
refer to user guide. If None, CPAC does not slice timing correction
"""
#scanParametersCSV = /home/data/settings/scan_parameters.csv
scanParametersCSV = None
