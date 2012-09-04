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
Directory Setup * *
===============
"""
# NOTE: Users must manually create these directories before running C-PAC
# Please specify the full path to each directory

# Directory where C-PAC should store temporary and intermediate files
workingDirectory = '/home/bcheung/p_integration_test'

# Directory where C-PAC should place crash logs
crashLogDirectory = '/home/bcheung/p_integration_test'

# Directory where C-PAC should put processed data
sinkDirectory = '/home/bcheung/p_integration_sink'


"""
==============================================
Optional Timeseries and Image Header Overrides
==============================================
"""
startIdx = 0

#if you specify none, CPAC calculates it on its own
stopIdx = None

#same here: based on headers
TR = None


"""
preproc
"""

standardResolution = '3mm'

fwhm = [4]


###Options are [1], [0], [1, 0]


runAnatomicalDataGathering = [1]
# anatomicalpreproc, segementation, nuisance won't run
# Anything that relies on those
runFunctionalDataGathering = [1]
# functional preproc wont run
# Anything that relies on that
runSymbolicLinks = [1]
#doesn't effect anything, but you wont get a pretty directory structure


runAnatomicalPreprocessing = [1]
# segmentation
# nuisance

runFunctionalPreprocessing = [1]
# nothing will run except segmentation and anatomical preproc

runRegistrationPreprocessing = [1]
# no derivatives will run


runAnatomicalToFunctionalRegistration = [1]
# nuisance and median angle won't run

runRegisterFuncToMNI = [1]
# no derivatives will run

runVMHC = [1]


# Generate motion statistics as 
# required for scrubbing
# can also be used as regressor in GA
runGenerateMotionStatistics = [1]



FSLDIR = commands.getoutput('echo $FSLDIR')
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
# required for nuisance and median angle

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
# Run nuisance signal correction
runNuisance = [1]

# Select which nuisance signal corrections to apply:

## compcor = CompCor
## wm = White Matter 
## csf = Cerebro Spinal Fluid
## gm = Gray Matter
## global = Global Mean Signal
## pc1 = First Principle Component
## motion = Motion
## linear = Linear Trend
## quadratic = Quadratic trend

# Options are 1 (correct) or 0 (ignore)
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
Temporal Filtering Options * *
==========================
"""
# Apply Temporal Filtering
runFrequencyFiltering = [1]

# First value = Lower bound for a band-pass filter
# Second value = Upper bound for a band-pass filter
# To use a high-pass filter, set the second value to NONE
# To use a low-pass filter, set the first value to NONE
nuisanceBandpassFreq =[(0.01, 0.1)]

"""
=======================================
Seed-based Correlation Analysis Options * *
=======================================
"""
# Run Seed-based Correlation Analysis
runSCA = [1]

# SCA will run on all ROI and voxel timeseries extracted below.

"""
==============================
Time Series Extraction Options *
==============================
"""
# Extract an average timeseries for each ROI
# Required if you wish to run ROI-based SCA
runROITimeseries = [0]

# Export ROI timeseries data
# First value = Output .csv
# Second value = Output numPy array
# Options are True/False
roiTSOutputs = [True, True]

# Directory containing ROI definitions
roiDirectoryPath = '/home2/data/Projects/NEO2012/mask_for_unitTS_extraction'

# Extract timeseries data for all individual voxels within a mask
# Required if you wish to run voxel-based SCA
runVoxelTimeseries = [1]

# Export voxel timeseries data
# First value = Output .csv
# Second value = Output numPy array
# Options are True/False
voxelTSOutputs = [False, False]

# Directory contaning masks
maskDirectoryPath = '/home2/data/Projects/NEO2012/mask_for_unitTS_extraction'

# Register timeseries data to a surface model built by FreeSurfer
# Required to run vertex timeseries extraction
runSurfaceRegistraion = [0]

# Extract timeseries data for surface vertices
runVerticesTimeSeries = [0]

# Export vertex timeseries data
# First value = Output .csv
# Second value = Output numPy array
# Options are True/False
verticesTSOutputs = [False, False]

# SOMETHING SOMETHING ABOUT FREESURFER. TALK TO SHARAD
reconSubjectsDirectory = '/home/data/Projects/NEO2012/FS_outputs'


"""
===================================
Regional Homogeneity (ReHo) Options * *
===================================
"""
# Calculate Regional Homogeneity
runReHo = [1]

# Cluster size (number of neighboring voxels)
# Options are 7, 19, and 27
clusterSize = 27

"""
==========================
Network Centrality Options * 
==========================
"""
# Calculate network centrality measures
runNetworkCentrality =[1]

# Select which centrality measures to calculate
# First value = Degree Centrality 
# Second value = Eigenvector Centrality
# Options are True/False
centralityMethodOptions = [True, True]

# Specify how connections are defined during graph construction
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

# Directory containing ROI definitions or masks
# Using ROIs will result in node-based centrality measures
# Using a mask will result in voxel-based centrality measures
templateDirectoryPath = '/home2/data/Projects/NEO2012/mask_for_unitTS' 

"""
====================================================
Bootstrap Analysis of Stable Clusters (BASC) Options
====================================================
"""

bascROIFile = '/home/data/Projects/nuisance_reliability_paper/seed_files/basil_ganglia/LEFT_BG_3_numbered+tlrc..nii.gz'

bascClusters = 6

# Number of bootstraps to apply to original timeseries data.
bascTimeseriesBootstraps = 100

# Number of bootstraps to apply to individual stability matrices.
bascDatasetBootstraps = 100

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
# Auto generated by extract_data
# Order of subjects in this list should match the order of the data in your model
groupAnalysisSubjectList = '/home/data/Projects/abidehbm/settings/subject_list_group_analysis.txt'
# Options come from list of resources
# Get list from Ranjeet, put here and in UG
derivativeList = ['alff_Z_standard', 'falff_Z_standard']
# SPecify path to FLS model(s)
# One path per model
# Generated from FSL, or through SHarads script (TALK TO HIM)
modelFile = '/home/data/Projects/abidehbm/setting/subject_list_model_list.txt'
# if you want to consider multiple scans at once
mixedScanAnalysis = False

zThreshold = 2.3
pThreshold = 0.05
fTest = True

