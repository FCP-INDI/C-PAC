
#list of subjects that are to be excluded 
#exclusionSubjectList = '/Users/ranjeet.khanuja/Desktop/ex_subjects.txt'
exclusionSubjectList = None

#list of subjects that are included 
#subjectList = None
#if None extract data runs on all the subjects
subjectList = '/home/data/Projects/c-pac-subject-list-fix/rockland_subject_new.txt'


#Anatomical file Path 
#Put %s where subjects are in the path
#anatomicalTemplate = '/home/data/Incoming/cambridge_fcon/*/%s/*/mprage_anonymized.nii.gz'
anatomicalTemplate = '/home/data/Originals/NYU_TRT/NYU_TRT_session1/%s/anat/anat.nii.gz'

#Functional Path
#Put  %s where subjects are in the path
#functionalTemplate = '/home/data/Incoming/cambridge_fcon/*/%s/*/rest.nii.gz'
functionalTemplate = '/home/data/Originals/NYU_TRT/NYU_TRT_session1/%s/func/lfo.nii.gz'
