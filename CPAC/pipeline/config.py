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
====================
Data Directory Setup ***
====================
"""
# NOTE: Users must manually create these directories before running C-PAC

# Directory where C-PAC should store temporary and intermediate files
workingDirectory = '/path/to/working_directory'

# Directory where C-PAC should place crash logs
crashLogDirectory = '/path/to/crash_directory'

# Directory where C-PAC should put processed data
sinkDirectory = '/path/to/output_directory'

"""
========================
Resolution and Smoothing ***
========================
"""
# Set the resolution (in mm) to which images are transformed
# Transformation occurs during registration and is requried for many measures
standardResolution = '3mm'

# Width (FWHM, in mm) of the Gaussian kernel used for spatial smoothing
# To skip smoothing, set to 0
fwhm = [4]

"""
========================
Resource Directory Setup NOT FINISHED- NEED MORE INFO FOR FSL FILES
========================
"""
# Directory where FSL is located
# If you have added FSL to your .bashrc file, this will be set automatically
FSLDIR = commands.getoutput('echo $FSLDIR')

# The following options specify the path of various resources provided by FSL
# By default, C-PAC will automatically locate these files based on FSLDIR
# Most users will not need to modify these values

# For users wishing to use non-standard versions of these resources:
## 1) Delete the string in parentheses beginning with FSLDIR
## 2) Replace this value with the full path to the appropriate file

standardResolutionBrain = os.path.join(FSLDIR,'data/standard/MNI152_T1_%s_brain.nii.gz' % (standardResolution))

standard = os.path.join(FSLDIR,'data/standard/MNI152_T1_%s.nii.gz' % (standardResolution))

standardBrainMaskDiluted = os.path.join(FSLDIR,'data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' % (standardResolution))

configFile = os.path.join(FSLDIR,'etc/flirtsch/T1_2_MNI152_%s.cnf' % (standardResolution))

brainSymmetric = os.path.join(FSLDIR,'data/standard/MNI152_T1_2mm_brain_symmetric.nii.gz')

symmStandard = os.path.join(FSLDIR,'data/standard/MNI152_T1_2mm_symmetric.nii.gz')

twommBrainMaskDiluted = os.path.join(FSLDIR,'data/standard/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz')

configFileTwomm = os.path.join(FSLDIR,'etc/flirtsch/T1_2_MNI152_2mm.cnf')

identityMatrix = os.path.join(FSLDIR,'etc/flirtsch/ident.mat')

harvardOxfordMask = os.path.join(FSLDIR,'data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz')


"""
==================
Timeseries Options ***
==================
"""
# Ignore volumes before this timepoint
# Options are an integer or None (defaults to beginning of timeseries)
startIdx = 0

# Ignore volumes after this timepoint
# Options are an integer or None (defaults to end of timeseries)
stopIdx = None

# Specify a TR other than what is listen in image headers
# Options are an integer or None (defaults to header information)
TR = None

# Specify slice timing correction in functional preprocessing
# Make sure the input CPAC subject list have the scan parameters
# information, if slice timing correction option is set true
sliceTimingCorrection = True
"""
================================
Preprocessing Workflow Selection ***
================================
"""
# Set which preprocessing workflows to run.

# WARNING:
# Many measures and outputs require that these workflows be run.
# Please refer to the developer documentation before changing these settings.
# Options (here and for most other settings) are: 1 = run, 0 = skip

runAnatomicalDataGathering = [1]

runFunctionalDataGathering = [1]

runAnatomicalPreprocessing = [1]

runFunctionalPreprocessing = [1]

runRegistrationPreprocessing = [1]

runRegisterFuncToMNI = [1]

runAnatomicalToFunctionalRegistration = [1]

runSymbolicLinks = [1]

"""
=========================================
Probabilistic Tissue Segmentation Options *** NEED TO FIX PRIOR PATH
=========================================
"""
# Run automatic tissue segmentation
runSegmentationPreprocessing = [1]

# C-PAC uses FSL to automatically distinguish tissue types based on priors.
# Each prior represents the probability that a given voxel will be 
# of a particular tissue type (white matter, gray matter, or CSF).

# Please specify the location and name of your prior files.
# Priors distributed with FSL must be binarized to be used by C-PAC
# For information about how to do this, please see the User Guide
prior_path = '/path/to/tissuepriors/%s' % standardResolution

# These values will be set automatically based on prior_path
PRIOR_CSF = os.path.join(prior_path, 'avg152T1_csf_bin.nii.gz')
PRIOR_GRAY = os.path.join(prior_path, 'avg152T1_gray_bin.nii.gz')
PRIOR_WHITE = os.path.join(prior_path, 'avg152T1_white_bin.nii.gz')

# Set thresholds for use during automatic tissue segmentation.
# Values correspond to probability thresholds for a given tissue type.
# For example, setting a value of 0.8 will result in areas with a 80 percent 
# probability of being a particular tissue type to be classified as such

cerebralSpinalFluidThreshold = [0.4]

whiteMatterThreshold = [0.66]

grayMatterThreshold = [0.2]

"""
==================================
Nusiance Signal Correction Options *** 
==================================
"""
# Run nuisance signal correction
runNuisance = [1]

# Select which nuisance signal corrections to apply:
## compcor = CompCor
## wm = White Matter 
## csf = CSF
## gm = Gray Matter
## global = Global Mean Signal
## pc1 = First Principle Component
## motion = Motion
## linear = Linear Trend
## quadratic = Quadratic Trend

# Options are 1 (apply) or 0 (ignore)
Corrections = [{'compcor' : 1,
                'wm' : 1,
                'csf' : 1,
                'gm' : 0,
                'global' : 0,
                'pc1' : 0,
                'motion' : 1,
                'linear' : 1,
                'quadratic' : 0}]

# Number of Principle Components to calculate for CompCor (usually 5 or 6)
# Only for use when 'compcor' is set to 1
nComponents = [5]

# Run median angle correction
runMedianAngleCorrection = [0]

# Target angle for median angle correction
targetAngleDeg = [90]

# Run Scrubbing
runScrubbing = [1]

# Generate FD and DVARS motion statistics
# Required to run scrubbing, but can also be used as regressors in a GLM
runGenerateMotionStatistics = [1]

# Specify maximum acceptable Framewise Displacement (in mm)
# Any volume with displacement greater than this value will be removed.
# One volume before and two volumes after each over-threshold volume
# will also be removed
scrubbingThreshold = [0.2]

#number of preceding frames to the offending time 
#frames to be removed (i.e.,those exceeding FD threshold)
numRemovePrecedingFrames = 1

#number of following frames to the offending time 
#frames to be removed (i.e.,those exceeding FD threshold)
numRemoveSubsequentFrames = 2

"""
==========================
Temporal Filtering Options ***
==========================
"""
# Apply Temporal Filtering
runFrequencyFiltering = [1]

# First value = Lower bound for a band-pass filter
# Second value = Upper bound for a band-pass filter
# To use a high-pass filter, set the second value to None
# To use a low-pass filter, set the first value to None
nuisanceBandpassFreq =[(0.01, 0.1)]

"""
==============================
Timeseries Extraction Options ***
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
roiDirectoryPath = '/path/to/roi_definitions_directory'

# Extract timeseries data for all individual voxels within a mask
# Required if you wish to run voxel-based SCA
runVoxelTimeseries = [1]

# Export voxel timeseries data
# First value = Output .csv
# Second value = Output numPy array
# Options are True/False
voxelTSOutputs = [False, False]

# Directory contaning masks
maskDirectoryPath = '/path/to/mask_definitions_directory'

# Register timeseries data to a surface model built by FreeSurfer
# Required to run vertex timeseries extraction
runSurfaceRegistraion = [0]

# Directory where FreeSurfer outputs surface data
# This should be the same as SUBJECTS_DIR in .bashrc
reconSubjectsDirectory = '/path/to/fs_output_directory'

# Extract timeseries data for surface vertices
runVerticesTimeSeries = [0]

# Export vertex timeseries data
# First value = Output .csv
# Second value = Output numPy array
# Options are True/False
verticesTSOutputs = [False, False]

"""
=======================================
Seed-based Correlation Analysis Options ***
=======================================
"""
# Run Seed-based Correlation Analysis
runSCA = [1]

# IN ORDER TO RUN SCA, YOU MUST ALSO RUN TIMESERIES EXTRACTION.
# SCA will be run on all ROI and voxel timeseries extracted above.
# Seeds for SCA must be specified in roiDirectoryPath or maskDirectoryPath.

"""
===================================
Regional Homogeneity (ReHo) Options ***
===================================
"""
# Calculate Regional Homogeneity
runReHo = [1]

# Cluster size (number of neighboring voxels)
# Options are 7, 19, and 27
clusterSize = 27

"""
============================================
Voxel-mirrored Homotopic Connectivity (VMHC) ***
============================================
"""
# Calculate VMHC for all gray matter voxels
runVMHC = [1]

# There are no options for VMHC

"""
==========================================================================
Amplitude of Low Frequency Oscillations (ALFF) and fractional ALFF Options ***
==========================================================================
"""
# Calculate ALFF and fALFF
runALFF = [1]

# NOTE: Frequency filtering is not applied when calculating fALFF

# Frequency cutoff (in Hz) for a high-pass filter
highPassFreqALFF = [0.01]

# Frequency cutoff (in Hz) for a low-pass filter
lowPassFreqALFF = [0.1]

"""
==========================
Network Centrality Options ***
==========================
"""
# Calculate network centrality measures
runNetworkCentrality = [1]

# Select which centrality measures to calculate
# First value = Degree Centrality 
# Second value = Eigenvector Centrality
# Options are True/False
centralityMethodOptions = [True, True]

# Specify how connections are defined during graph construction
# First value = Binarized (connection strenth is either 0 or 1)
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
templateDirectoryPath = '/path/to/centrality_mask_roi_directory' 

#Option to generate adjacency matrix png image
# and adjacency matrix mat file
#Takes lot of memory. Do not turn it on for voxel based graph.
generateAdjacencyGraph = False
"""
====================================================
Bootstrap Analysis of Stable Clusters (BASC) Options **
====================================================
"""
# Run BASC
runBASC = [0]

# Path to a mask file. Voxels outside this mask will be excluded from BASC.
bascROIFile = '/path/to/basc_mask_file'

# Number of clusters at both the individual and group level.
bascClusters = 6

# Number of bootstraps to apply to original timeseries data.
bascTimeseriesBootstraps = 100

# Number of bootstraps to apply to individual stability matrices.
bascDatasetBootstraps = 100

# Path to a text file containing Affinity Thresholds for each subject.
# These are correlation thresholds applied prior to spectral clustering.
# Can be subject specific when subjects have differing numbers of timepoints.
# Subjects should be in the same order as in the main subject list.
bascAffinityThresholdFile = '/path/to/basc_affinity_threshold_file'

"""
================================================
Connectome-wide Association Study (CWAS) Options **
================================================
"""
# Run CWAS
runCWAS = [0]

# Path to a mask file. Voxels outside this mask will be excluded from CWAS.
cwasROIFile = '/path/to/cwas_mask_file'

# Number of permutation tests to run on the Psuedo-F statistic
cwasFSamples = 5000

# Number of NiPype nodes to be created while computing CWAS.
# This number depends on computing resources
cwasParallelNodes = 10

# Path to a text file containing phenotypic regressor.
cwasRegressorFile = '/path/to/cwas_regressor_file'
 
"""
============================
Group Statistics Options ***
============================
"""
# Calculate group statistics
runGroupAnalysis = [1]

# Path to list of subjects on which to run group statistics
# This file should be created automatically when you run extract_data.py
# The order of subjects in this list must match the order in your model
groupAnalysisSubjectList = '/path/to/subject_list_group_analysis.txt'

# Select which measures should be included in group analysis:
## sca_seed_Z_to_standard_smooth = voxel-based SCA
## sca_roi_Z_to_standard_smooth = ROI based SCA
## alff_Z_to_standard_smooth = ALFF
## falff_Z_to_standard_smooth  = fALFF
## vmhc_z_score_stat_map = VMHC
## reho_Z_to_standard_smooth = ReHo
derivativeList = ['sca_seed_Z_to_standard_smooth', \
                  'sca_roi_Z_to_standard_smooth', \
                  'alff_Z_to_standard_smooth', \
                  'falff_Z_to_standard_smooth', \
                  'vmhc_z_score_stat_map', \
                  'reho_Z_to_standard_smooth']

# Location of a text file contaning a list of FSL models
# Each line in this file should be the path to a model directory
# Each model directory should contain a .mat, .con, and .grp file
# If fTest = True (see below), model directories must also contain a .fts file
# These models can be generated through FSL, or using create_fsl_model.py
# For instructions on using create_fsl_model.py, see the user guide
modelFile = '/path/to/subject_list_model_list.txt'

# If a subjecs has multiple scans:
# False = Consdier only the first scan session during group analysis
# True = Consider all scan sessions
mixedScanAnalysis = False

zThreshold = 2.3

pThreshold = 0.05

# Run an F-test
# Options are True/False
fTest = True

