import os
import commands

"""
=================
Computer Settings
=================
"""
# True = Run on compute cluster
# False = Run on local machine
runOnGrid = False

# Number of subjects to run simultaneously
# This number depends on computing resources
# Only applies when running on a local machine with multiple cores
numSubjectsAtOnce = 5

# Number of cores (local) or slots on a node (cluster) per subject
# Slots are cores on a cluster node
# This number depends on computing resources
# Only applies when local machine has multiple cores or runOnGrid = True
numCoresPerSubject = 2

# Options are 'SGE' (Sun Grid Engine) or 'PBS' (Portable Batch System)
# Only applies when runOnGrid = True
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
===============
Directory Setup
===============
"""
# The folder containing 
dataPath = '/home2/ssikka/nki_nyu_pipeline/testing/process'

workingDirectory = '/home/bcheung/p_integration_test'

crashLogDirectory = '/home/bcheung/p_integration_test'

sinkDirectory = '/home/bcheung/p_integration_sink'

"""
subjects
  """



## subjectDirectory = '/home2/data/Incoming/fcon_test/
## subjectList = None   
## the above two is not usable now...
exclusionSubjectList = None

"""
anat log
"""

???

anatLogFile = 'log.txt' #Don't need this, for IPN use only. -Yang

anatLogFilePath = '%s/*/*/%s'


"""
preproc setup
"""


#FSLDIR = '/usr/local/cmi/fsl/4.1' # for gelert

FSLDIR = commands.getoutput('echo $FSLDIR')

priorDirectory = '/home2/data/Projects/C-PAC/tissuepriors'



"""
==============================================
Optional Timeseries and Image Header Overrides
==============================================
"""
startIdx = 0

stopIdx = None

TR = None


"""
preproc
"""
###Options are [1], [0], [1, 0]

runAnatomicalDataGathering = [1]
runAnatomicalPreprocessing = [1]
runRegistrationPreprocessing = [1]


runFunctionalDataGathering = [1]
runFunctionalPreprocessing = [1]
runAnatomicalToFunctionalRegistration = [1]




runRegisterFuncToMNI = [1]



runVMHC = [1]

runGenerateMotionStatistics = [1]




runSymbolicLinks = [1]



standardResolution = '3mm'

fwhm = [4]




standardResolutionBrain = os.path.join('/usr/share/fsl/4.1/data/standard/MNI152_T1_%s_brain.nii.gz' % (standardResolution))
standard = os.path.join('/usr/share/fsl/4.1/data/standard/MNI152_T1_%s.nii.gz' % (standardResolution))
standardBrainMaskDiluted = os.path.join('/usr/share/fsl/4.1/data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' % (standardResolution))
configFile = os.path.join('/usr/share/fsl/4.1/etc/flirtsch/T1_2_MNI152_%s.cnf' % (standardResolution))
brainSymmetric = os.path.join('/usr/share/fsl/4.1/data/standard/MNI152_T1_2mm_brain_symmetric.nii.gz')
symmStandard = os.path.join('/usr/share/fsl/4.1/data/standard/MNI152_T1_2mm_symmetric.nii.gz')
twommBrainMaskDiluted = os.path.join('/usr/share/fsl/4.1/data/standard/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz')
configFileTwomm = os.path.join('/usr/share/fsl/4.1/etc/flirtsch/T1_2_MNI152_2mm.cnf')
identityMatrix = os.path.join('/usr/share/fsl/4.1/etc/flirtsch/ident.mat')


"""
=========================================
Probabilistic Tissue Segmentation Options
=========================================
"""
# Run automatic tissue segmentation
runSegmentationPreprocessing = [1]

# C-PAC uses FSL to automatically distinguish tissue types based on priors.

# Each prior represents the probability that a given voxel will be 
# of a particular tissue type (white matter, gray matter, or CSF).

# Please specify the location and name of your prior files.

prior_path = '/home/data/Projects/C-PAC/tissuepriors/%s' % standardResolution
PRIOR_CSF = os.path.join(prior_path, 'avg152T1_csf_bin.nii.gz')
PRIOR_GRAY = os.path.join(prior_path, 'avg152T1_gray_bin.nii.gz')
PRIOR_WHITE = os.path.join(prior_path, 'avg152T1_white_bin.nii.gz')

# Set thresholds for use during automatic tissue segmentation

# C-PAC uses FSL to automatically distinguish tissue types based on priors
# Make sure you have set the prior_path

cerebralSpinalFluidThreshold = [0.4]

whiteMatterThreshold = [0.66]

grayMatterThreshold = [0.2]

"""
==================================
Nusiance Signal Correction Options
==================================
"""
# Run nuisance signal corrections
runNuisance = [1]

# Select which nuisance signals to remove. 1 to run, 0 to skip
# compcor = CompCor
# wm = White Matter 
# csf = Cerebro Spinal Fluid
# gm = Gray Matter
# global = Global Mean Signal
# pc1 = First Principle Component
# motion = Motion
# linear = Linear Trend
# quadratic = Quadratic trend
Corrections = [{'compcor' : 1,
                'wm' : 1,
                'csf' : 1,
                'gm' : 0,
                'global' : 0,
                'pc1' : 0,
                'motion' : 1,
                'linear' : 1,
                'quadratic' : 0}]

# Number of Principle Components
# Only for use when 'compcor' : 1
nComponents = [5]

# Run median angle correction
runMedianAngleCorrection = [0]

# Target angle for median angle correction
# Only for use when runMedianAngleCorrection = [1] or [0, 1]
targetAngleDeg = [90]

# 1 to run, 0 to skip
runScrubbing = [1]

# Volumes with displacement greater than this value (in mm) will be removed.
# Only for use when runScrubbing = [1] or [0,1]
scrubbingThreshold = [0.2]

"""
==========================
Temporal Filtering Options
==========================
"""
# Apply Temporal Filtering
runFrequencyFiltering = [1]

nuisanceBandpass = True

nuisanceBandpassFreq =[(0.01, 0.1)]


"""
SCA
"""
runSCA = [1]

"""
==============================
Time Series Extraction Options
==============================
"""


runROITimeseries = [0]

# Output type: .csv, numPy
roiTSOutputs = [True, True]

roiDirectoryPath = '/home2/data/Projects/NEO2012/mask_for_unitTS_extraction'


runVoxelTimeseries = [1]

# Output type: .csv, numPy
voxelTSOutputs = [False, False]

maskDirectoryPath = seedDirPath


runSurfaceRegistraion = [0]
runVerticesTimeSeries = [0]

# Output type: .csv, numPy

verticesTSOutputs = [False, False]

reconSubjectsDirectory = '/home/data/Projects/NEO2012/FS_outputs'


"""
====================================
Regional Homogeneity (ReHo) Options) *
====================================
"""
# Run ReHo
runReHo = [1]

# Cluster size (number of voxels) to use.
# Options are 7, 19, and 27
clusterSize = 27

"""
==========================
Network Centrality Options
==========================
"""
# Calculate network centrality measures
runNetworkCentrality =[1]

# Select which centrality measures to calculate
# First value = Degree Centrality 
# Second value = Eigenvector Centrality
# Options are True/False
centralityMethodOptions = [True, True]

# Define how connections are defined during graph construction
# First value = Binarize (connection strenth is either 0 or 1)
# Second value = Weighted (connection strength is a correlation value)
# Options are True/False
centralityWeightOptions = [True, True]

# Select what type of threshold is applied to create an adjacency matrix
# 0 = Significance threshold (P-value)
# 1 = Sparsity threshold (Sparsity value)
# 2 = Correlation threshold (Pearson's r)
correlationThresholdOption = 1

# Based on the type of threshold selected above, enter the appropriate value
# Significance threshold = P-value
# Sparsity threshold = sparsity value
# Correlation threshold = Pearsons' r value
# examples: 0.05, 0.0744, 0.6
correlationThreshold = 0.0744

#path to mask/roi directory for netwrok centrality 
templateDirectoryPath = seedDirPath 

"""
====================================================
Bootstrap Analysis of Stable Clusters (BASC) Options
====================================================
"""

bascROIFile = '/home/data/Projects/nuisance_reliability_paper/seed_files/basil_ganglia/LEFT_BG_3_numbered+tlrc..nii.gz'

bascClusters = 6

# Number of bootstraps to apply individual stability matrices.
bascDatasetBootstraps = 100

# Number of bootstraps to apply original timeseries data.
bascTimeseriesBootstraps = 100

bascAffinityThresholdFile = '/home/bcheung/Dropbox/server_shares/CPAC_git/CPAC_main/C-PAC/CPAC/pipeline/subjects_affine.txt'

"""
================================================
Connectome-wide Association Study (CWAS) Options
================================================
"""

cwasROIFile = '/home/data/Projects/nuisance_reliability_paper/seed_files/basil_ganglia/LEFT_BG_3_numbered+tlrc..nii.gz'

cwasFSamples = 5000

cwasParallelNodes = 10

cwasRegressorFile = '/home/bcheung/Dropbox/server_shares/CPAC_git/CPAC_main/C-PAC/CPAC/pipeline/subject_regressors.txt'

"""
==========================================================================
Amplitude of Low Frequency Oscillations (ALFF) and fractional ALFF Options * *
==========================================================================
"""
# Calculate ALFF and fALFF
runALFF = [1]

# NOTE: Frequency filtering is not applied when calculating fALFF

# Frequency cutoff (in Hz) for a high-pass filter
# All frequencies lower than this value will be excluded from analysis
# To skip high-pass filtering, leave this array empty []
highPassFreqALFF = [0.01]

# Frequency cutoff (in Hz) for a low-pass filter
# All frequencies higher than this value will be excluded from analysis
# To skip low-pass filtering, leave this array empty []
lowPassFreqALFF = [0.1]

"""
======================
Group Analysis Options
======================
"""

derivativeList = ['ALFF_Z_FWHM_2standard', 'fALFF_Z_FWHM_2standard', 'sca_Z_FWHM']

modelsDirectory = '/home/data/Projects/abidehbm/group_models/'

"""
============================
Group Analysis File Template
============================
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
=========================
Derivatives File Template
=========================
"""
dervTemplate = sinkDirectory + '/sym_links/%s/%s/*/%s.nii.gz'

labelFile = sinkDirectory + '/sym_links/label_linkage.txt'

subList = '/home/data/Projects/abidehbm/settings/subject_list_group_analysis.txt'

dervTemplateList = ['label', 'derivative']

"""
===================
Statistical Options
===================
"""
zThreshold = 2.3

pThreshold = 0.05

fTest = True

