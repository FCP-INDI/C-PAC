import os
import commands

"""

Welcome to the C-PAC configuation file.



Each setting is explained by comments immediately above it. Capitalized words (eg. THIS, ALL) indicate parts

of a setting string which should be replaced with user specific settings.



Possible setting options will be given (in parentheses).





1. Processing Settings



    a) runOnGrid : Run on local machine (False), or on compute cluster/grid (True)



    b) qsubArgs : For use only when runOnGrid = True



       -q QUEUE.q : Set the queue of nodes to use on the cluster/grid

				    The default setting is to use all nodes (-q all.q)

				    

				    If you do not know what value to enter, ask the system administrator for your cluster.

    

       -pe A-B : Specify the range of slots to use on each node. 

       			 A is the lower bound, B is the upper bound. A and B should both be numbers.

       				  

       			 This setting is optional, and should be added to the end of the qsubArgs string.

       			 For example, qsubArgs = -q all.q -pe 1-4 

       				

    c) numSubjectsAtOnce : Set how many subjects to run at once. This number will depend on computing resources.



    d) numCoresPerSubject : For use only when runOnGrid = False

			       			

			       			If your local machine has more than one processor, the number of cores to use per subject.

"""

runOnGrid = False

#only for multicore
numSubjectsAtOnce = 5

#for multicore and SGE/PBS
#for SGE/PBS signifies number of slots on a single node of a cluster
numCoresPerSubject = 2

#options are 'SGE' , 'PBS'
resourceManager = 'SGE'

queue = 'all.q'

#options for SGE only here,
#SGE users must set this enviroment,
#easy way to know your parallel environment is to execute the following on command on cluster
# $ qconf -spl

#A pipeline for each subject that needs preprocessing is spawned on single nodes of the cluster.
# To avoid I/O overhead the pipeline should only use the resources(cores) from that node.
# The users can enable this feature on Sun Grid Engine by modifying their parallel environment
# or adding a new parallel ennviroment using the exisiting environment parameters(some parameters tweaked)

#To create new environment using old environment follow ths  procedure

# 1. find the parallel environments u have on cluster
# $ qconf -spl

# 2. look through your parallel environments. Mine looks like this

#$ qconf -sp mpi
# pe_name            mpi
# slots              999
# user_lists         NONE
# xuser_lists        NONE
# start_proc_args    NONE
# stop_proc_args     NONE
# ---># allocation_rule    $fill_up
# control_slaves     TRUE
# job_is_first_task  FALSE
# urgency_slots      min
# accounting_summary TRUE

# 3. use old to create new environment
# $ qconf -sp mpi > cpac_mpi 

# 4. change the allocation_rule highlighted with the arrow to $pe_slots,
#    use your favourite text editor to accomplish this

# 5. Add your new envionment file to SGE 

# qconf -Ap cpac_mpi

# 6. Specify this new environment below
parallelEnvironment = 'mpi'





"""

2. Directory Setup



	NOTE: C-PAC currently requires that you create all directories manually before running.



	a) workingDirectory : Specify where C-PAC should store temporary and intermediate files.

						  Temporary cache files and intermediate data will both be put in this directory.

						  We recommended you delete this directory once C-PAC is finished running.



	b) crashLogDirectory : Specify where C-PAC should write crash logs if an error occurs.



						   NOTE: Due to a bug, C-PAC currently writes crash files to the parent folder of config.py.

	

	c) sinkDirectory : Specify where C-PAC should put processed data.

    d) dataPath : Specify where all the data is present for each site

"""


workingDirectory = '/home/bcheung/p_integration_test'

crashLogDirectory = '/home/bcheung/p_integration_test'

sinkDirectory = '/home/bcheung/p_integration_sink'

dataPath = '/home2/ssikka/nki_nyu_pipeline/testing/process'



"""

	d) Subject Directory Settings

		

		1) subjectDirectory : Specify the directory containing subject data sub-directories.

    	2) subjectList : This setting is optional. If no subject list is set, all subjects will be run.

						 

						 Set the full path to a .txt file containing a list of subjects to be processed.

						 This file should contain one subject identifier on each line.

						 Each subject identifier should reflect a folder in the subjectDirectory.



						 

						 If you are not using a subject list, set subjectList = None





		3) exclusionSubjectList : Similar to subjectList, but defines which subjects to exclude.

								  This file should be formatted in the same way as the subject list.

								  C-PAC will run all subjects not listed in this file.

								  This setting is useful if there are subjects with known bad data.



								  This setting is also optional.

								  If you are not using an exclude list, set exclusionSubjectList = None



								  If you have also specified a subject list, any subjects identified

								  in the exclusion list will not be run.



"""



## subjectDirectory = '/home2/data/Incoming/fcon_test/
## subjectList = None   
## the above two is not usable now...
exclusionSubjectList = None



"""

	e) Set Anatomical File Names and Locations



		1) anatTemplate : Tell C-PAC about your directory structure.

						  

						  In order to correctly process all your data, C-PAC needs to understand how

						  you have organized your files. This is done by defining a file template.



						  Using wildcards (%s) define the folder structure containing anatomical files.

						  Replace the subject folder name and file name with wildcards.



						  For example, if the path to the anatomical scan for subject 001 was:



						  /data/study_1/001/T1/session_1/brain.nii.gz



						  My template would be:



						  %s/T1/session_1/%s.nii.gz



						  Specify only the path within subjectDirectory.



						  C-PAC will automatically insert the subject identifier (parent folder name) in place of

						  the first wildcard, and insert the file name in place of the second wildcard.



		2) anatTemplateList : Leave the first value as 'subject'.

							  

							  Set the second value to the filename of your anatomical scans.

							  In the scenario described above, we would enter 'brain' as the second value.



"""



## DEPRECIATED anatTemplate = '*/%s/*/%s/%s.nii.gz'

## DEPRECIATED anatTemplateList = ['subject', 'scan', 'mprage']



"""

	    anatLogFile: In case of multiple anatomical scans per subject.

					 Specify the name of the log file which indicates which single anatomical scan to process.



	    			anatLogFilePath: Specify the location of the log file relative to the subject directory location.

					The pipeline substitutes 1st %s with the subject number/name and the 2nd %s is the name of the log file. 

					The User needs to specify only the structure as per his analysis. 





	    Ex: /sam/wave1/sub8000/

                        anat_1/mprage.nii.gz

                        anat_2/mprage.nii.gz

                        logs/log.txt

	    anatLogFile = 'log.txt'

	    anatLogFilePath = '%s/*/%s'

    

"""



anatLogFile = 'log.txt' #Don't need this, for IPN use only. -Yang

anatLogFilePath = '%s/*/*/%s'



"""

		Functional File Name and Location within Subject Directory

		Note: functionalFileName substituted into functionalDirectorySetup template

"""



#funcTemplate = '%s/%s/%s.nii.gz'

#funcTemplateList = ['subject', 'session', 'rest']

#funcSessionFile = '/home/ssikka/nki_nyu_pipeline/sessions.txt

# DEPRECIATED funcTemplate = '*/%s/*/%s/%s.nii.gz'

# DEPRECIATED funcTemplateList = ['subject','scan','rest']


"""

3. Pre-Processing Setup 



	a) FSLDIR : Set the directory where you have installed FSL



    b) priorDirectory : Set the directory where tissue priors are located

"""


#FSLDIR = '/usr/local/cmi/fsl/4.1' # for gelert

FSLDIR = commands.getoutput('echo $FSLDIR')

priorDirectory = '/home2/data/Projects/C-PAC/tissuepriors'



"""

	c) Optional Header and Timeseries Overrides



		1) startIdx : Starting time point (defaults to 0).



		2) stopIdx : Stopping time point (defaults to header specification)



		3) TR : Time of repetition (defaults to header specifications)



		

		NOTE: Set these values = None if you are not using overrides.

"""



startIdx = 0
stopIdx = None


TR = None



"""

4. Image Pre-Processing Settings

	

	NOTE: Users can specify multiple parameters for each setting by including them in an array (eg. [1, 2, 3]).

		  Setting multiple parameters will cause C-PAC to follow divergent analysis pipelines after that step.



		  For example, setting whiteMatterThreshold = [0,4] will result in two sets of processed data, one processed

		  with a white matter threshold of 0, and another with a threshold of 4.



	a) standardResolution : Enter the standard voxel size (in mm) for your images in standard space.



	b) MNI : Define your standard space. Options are MNI152 and MNI305



	c) fwhm : Specifies the width (in mm) for a low-pass 3D Gaussian filter.

			  To skip filtering, set this value to 0



			  NOTE: Recommended values are 1.5 or 2x voxel size.

			  		For example, if you are using 3x3x3mm voxels, set fwhm = [6]

"""
###Options are [1], [0], [1, 0]

runAnatomicalDataGathering = [1]
runAnatomicalPreprocessing = [1]
runRegistrationPreprocessing = [1]
runSegmentationPreprocessing = [1]

runFunctionalDataGathering = [1]
runFunctionalPreprocessing = [1]
runAnatomicalToFunctionalRegistration = [1]

runNuisance = [1,0]
runMedianAngleCorrection = [1,0]
runFrequencyFiltering = [1,0]
runRegisterFuncToMNI = [1]

standardResolution = '3mm'
MNI = 'MNI152'

fwhm = [4]

prior_path = os.path.join(priorDirectory, standardResolution)
PRIOR_CSF = os.path.join(prior_path, 'avg152T1_csf_bin.nii.gz')
PRIOR_GRAY = os.path.join(prior_path, 'avg152T1_gray_bin.nii.gz')
PRIOR_WHITE = os.path.join(prior_path, 'avg152T1_white_bin.nii.gz')
standardResolutionBrain = os.path.join(FSLDIR,
            'data/standard/MNI152_T1_%s_brain.nii.gz' % (standardResolution))
standard = os.path.join(FSLDIR,
            'data/standard/MNI152_T1_%s.nii.gz' % (standardResolution))
standardBrainMaskDiluted = os.path.join(FSLDIR,
            'data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' % (standardResolution))
configFile = os.path.join(FSLDIR,
            'etc/flirtsch/T1_2_MNI152_%s.cnf' % (standardResolution))
brainSymmetric = os.path.join(FSLDIR,
            'data/standard/MNI152_T1_2mm_brain_symmetric.nii.gz')
symmStandard = os.path.join(FSLDIR,
            'data/standard/MNI152_T1_2mm_symmetric.nii.gz')
twommBrainMaskDiluted = os.path.join(FSLDIR,
            'data/standard/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz')
configFileTwomm = os.path.join(FSLDIR,
            'etc/flirtsch/T1_2_MNI152_2mm.cnf')
identityMatrix = os.path.join(FSLDIR,
            'etc/flirtsch/ident.mat')


"""

	Set thresholds for use during tissue segmentation.



	C-PAC uses FSL for automatic tissue segmentation based on probabilistic atlases. As such, values entered here 

	correspond to probability thresholds for a given tissue type. For example, setting a value of .8 will result in 

	areas with a 80 percent probability of being a particular tissue type to be selected

"""

cerebralSpinalFluidThreshold = [0.4]

whiteMatterThreshold = [0.66]

grayMatterThreshold = [0.2]



"""

	Set whether to remove volumes with excessive movement ("Scrubbing" as described by Power et al. 2012)

"""

scrubData = [False]

scrubbingThreshold = [0.2]



"""

	d) Nuisance Signal Correction



		1) Corrections : Select which nuisance signals you would like to remove. 



				0 - Global mean signal regression

				1 - Compcor: A component based noise correction method (CompCor) for BOLD and perfusion based fMRI.

							 Extracts the principal components found in white matter and cerebral spinal fluid

						     Algorithm based on Behzadi et al. 2007 (NeuroImage)

				2 - White Matter

				3 - CSF regression (CSF)

				4 - GRAY Matter

				5 - First principal component

				6 - Linear trend

				7 - Motion Regression



				Using arrays, you can combine multiple techniques.

				For example, a setting of [2,3,4] would remove signal from white matter, gray matter, and CSF.

				Nested arrays allow for multiple pipelines. Example: [[1,2,3],[4,5,6]]



		2) nComponents : For use only when Compcor is selected

						 

						 Specify the number of components to regress (usually 5 or 6).



"""



Corrections = [{'compcor' : 1,
                'wm' : 1,
                'csf' : 1,
                'gm' : 0,
                'global' : 0,
                'pc1' : 0,
                'motion' : 1,
                'linear' : 1,
                'quadratic' : 0}]

nComponents = [5]



"""

    e) Median Angle Correction



		targetAngleDeg : If you would like to run median angle correction, enter the target angle in degrees.



						 If you do not wish to run this analysis, leave the array empty,[]



        Algorithm based on: H. He and T. T. Liu, A geometric view of global signal confounds in resting-state

        										 functional MRI, NeuroImage, Sep. 2011.

"""



targetAngleDeg = [90]



"""

	f) Temporal Filtering



		nuisanceBandpass : Select whether to apply temporal filters to your data. Options are True and False.



		nuisanceBandpassFreq : For use only when nuisanceBandpass = True



							   Specify which type of filtering to use, and at what frequencies.



							   These values are entered as a pair within parentheses.



							   The first value indicates the lower frequency bound, the second the higher bound.



							   To run a high-pass filter, set the first value to NONE.



							   To run a low-pass filter, set the second value to NONE.



							   You can specify multiple filtering options by nesting sets within the array.

							   Example: [(1,2),(2,3),(NONE,4)]



"""

nuisanceBandpass = True

nuisanceBandpassFreq =[(0.01, 0.1)]

"""





5. Analysis Setup



	a) derivatives : Select which analyses to run. Enter a value of True or False in the array position corresponding

					 to the analyses you wish to run (True to run, False to skip).

					



					 f/ALFF: (Fractional) Amplitude of Low Frequency Fluctuations. Runs both ALFF and fALFF.

			 

					 SCA: Seed-based Correlation Analysis, Biswal et al., (1995). doi:10.1002/mrm.1910340409



					 VMHC: Voxel-matched Homotopic Connectivity, Zuo, et al., (2010). doi:10.1523/JNEUROSCI.2612-10.2010

									   

					 ReHo : Regional homogeneity, Zang 2004 NeuroImage.

					 

					 tsE: Time Series Extraction



					 VerticesE: Vertices Time Series Extraction 



					 GA: FSL Group Analysis http://www.fmrib.ox.ac.uk/fsl/feat5/detail.html#higher

"""



#derivatives = [f/ALFF, SCA, VMHC, ReHo, tsE, VerticesE, GA] 

derivatives = [True, False, False, False, False, False, False]



"""

	b) ALFF/fALFF Options



    	For use only when f/ALFF = True



    	You must specify values for both filters.



    	1) highPassFreqALFF : Set frequency (in Hz) cutoff for the high-pass filter

		

		2) lowPassFreqALFF : Set frequency (in Hz) cutoff for the low-pass filter





    	NOTE: Regardless of user settings, data run through this analysis will never receive scrubbing.

"""



highPassFreqALFF = [0.01]

lowPassFreqALFF = [0.1]



"""

	c) Seed-Based Correlation Analysis Options



    	For use only when SCA = True



    	1) seedFile : Specify the full path to the file containing a list of seed regions.

    				  Each line of this file should contain the full path to one seed.

"""



seedFile = '/home2/data/Projects/ABIDE_MP/settings/seeds_list.txt' # yang

correlationSpace = 'mni' # what are other options?



"""

	d) Timeseries Extraction Options
	   For use only when tsE = True
	   1) unitTSOutputs : Specify output types for Unit Time Series data.
						  If you do not wish to run this analysis, enter False for both values.
	   2) unitDefinitionsDirectory : Specify the directory containing unit definitions.

		

									 NOTE: Definitions directory should contain one subdirectory for 

									 each set of units to be generated.



									 (e.g., Harvard-Oxford Atlas, AAL, Craddock, Dosenbach-160);



	   3) voxelTSOutputs : Specify output types for Voxel Time Series data.



						   If you do not wish to run this analysis, enter False for both values.



	   4) voxelMasksDirectory : Specify the directory containing mask definitions.



								NOTE: Definitions directory should contain one subdirectory for each

								mask or mask set to be used during voxel selection.



								This analysis will output one file per mask.

"""



# Output type: .csv, numPy

unitTSOutputs = [True, True]

unitDefinitionsDirectory = '/home2/data/Projects/NEO2012/mask_for_unitTS_extraction'

# Output type: .csv, numPy

voxelTSOutputs = [False, False]

voxelMasksDirectory = '/home2/data/Projects/NEO2012/mask_for_TS_extraction'



"""

    e) Vertices Time Series Extraction Options



	   For use only when VerticesE = True



	   1) verticesTSOutputs : Specify output types for Vertex Time Series data.



						   	  If you do not wish to run this analysis, enter False for both values.



	   2) runSurfaceRegistration : Registers vertex activation data to a surface model built by FreeSurfer.



	   3) reconSubjectsDirectory :

"""

# Output type: .csv, numPy

verticesTSOutputs = [False, False]

runSurfaceRegistraion = False

reconSubjectsDirectory = '/home/data/Projects/NEO2012/FS_outputs'





"""

	e) FSL Group Analysis

		

		NOTE: Separate group analysis is conducted for each derivative.



			- Separate group analysis conducted for each derivative

			- Not applicable to time series extraction derivatives *HUH?*



		1) derivativeList : Speicify one or more derivative outputs (e.g. SCA, VMHC, ReHo) on which

						    to run group analysis.



						    By default, some derivatives generated during the C-PAC run will automatically be included.



		2) modelsDirectory : Specify a model list file that contains one or more models to be executed for each 

							 derivative.

"""

derivativeList = ['ALFF_Z_FWHM_2standard', 'fALFF_Z_FWHM_2standard', 'sca_Z_FWHM']

modelsDirectory = '/home/data/Projects/abidehbm/group_models/'



"""

     	Templates for Group Analysis

     	

     	Design Files - These should be generated by the FSL Model Generation Utility

     	.con -> list of contrasts requested.

     	.fts -> list of F-tests requested.

     	.mat -> the actual values in the design matrix



     	model_name in the *TemplateList files will be fetched from the modelsDirectory specified above.

"""

matTemplateList = ['model_name']

conTemplateList = ['model_name']

ftsTemplateList = ['model_name']

grpTemplateList = ['model_name']



mat = '/Users/ranjeet.khanuja/Desktop/data2/models/%s.mat'

con = '/Users/ranjeet.khanuja/Desktop/data2/models/%s.con'

fts = '/Users/ranjeet.khanuja/Desktop/data2/models/%s.fts'

grp = '/Users/ranjeet.khanuja/Desktop/data2/models/%s.grp'



"""

    	Derivative Template

    	

    	The first argument is label which is actually the strategy name. The pipeline

    	automatically creates the lable-linkage file in the sym_links folder.

    	The second argument is derivative name. This will be fetched from 

    	the derivative list defined above

"""



dervTemplate = sinkDirectory+ '/sym_links/%s/%s/*/%s.nii.gz'

labelFile = sinkDirectory + '/sym_links/label_linkage.txt'

subList = '/home/data/Projects/abidehbm/settings/subject_list_group_analysis.txt'



dervTemplateList = ['label', 'derivative']



"""

	Statistical Options

"""

zThreshold = 2.3

pThreshold = 0.05

fTest = True

