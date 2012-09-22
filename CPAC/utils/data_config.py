
#list of subjects that are to be excluded 
#exclusionSubjectList = '/Users/ranjeet.khanuja/Desktop/ex_subjects.txt'
exclusionSubjectList = None

#list of subjects that are included 
subjectList = None
#if None extract data runs on all the subjects
#subjectList = '/home/data/Projects/c-pac-subject-list-fix/rockland_subject_new.txt'



#Put %s where site and subjects are in the path
#anatomicalTemplate = '/home/data/Incoming/cambridge_fcon/%s/%s/*/mprage_anonymized.nii.gz'
anatomicalTemplate = '/home/data/Originals/DiscSci/NIFTI/T1/%s/anat.nii.gz'

#Functional Path
#Put  %s where site and subjects are in the path
#functionalTemplate = '/home/data/Incoming/cambridge_fcon/%s/%s/*/rest.nii.gz'
functionalTemplate = '/home/data/Originals/DiscSci/NIFTI/BOLD/%s/REST*.nii.gz'


#list of sites
#if None extract data runs on all sites
siteList = None

#slice timing parameters csv file path
sliceTimingParametersCSV = None