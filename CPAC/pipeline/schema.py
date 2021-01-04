from voluptuous import Schema, Required, All, Any, Length, Range, Match, In


schema = Schema({
    Required('pipelineName'): All(str, Length(min=1)),
    Required('workingDirectory'): str,
    Required('crashLogDirectory'): str,
    Required('logDirectory'): str,
    Required('outputDirectory'): str,

    'FSLDIR': str,
    'runOnGrid': bool,
    'resourceManager': str,
    'parallelEnvironment': str,
    'queue': str,

    'awsOutputBucketCredentials': str,
    's3Encryption': bool,  # check/normalize

    'maximumMemoryPerParticipant': float,
    'maxCoresPerParticipant': All(int, Range(min=1)),
    'numParticipantsAtOnce': All(int, Range(min=1)),
    'num_ants_threads': All(int, Range(min=1)),

    'write_func_outputs': bool,
    'write_debugging_outputs': bool,  # check/normalize
    'generateQualityControlImages': bool, # check/normalize
    'removeWorkingDir': bool,
    'run_logging': bool,
    'reGenerateOutputs': bool, # check/normalize
    'runSymbolicLinks': bool, # check/normalize

    'resolution_for_anat': All(str, Match(r'^[0-9]+mm$')),
    'template_brain_only_for_anat': str,
    'template_skull_for_anat': str,
    'template_brain_mask_for_anat': str,
    'template_symmetric_brain_only': str,
    'template_symmetric_skull': str,
    'dilated_symmetric_brain_mask': str,

    'already_skullstripped': bool,
    'skullstrip_option': In(['AFNI', 'BET']),
    'skullstrip_shrink_factor': float,
    'skullstrip_var_shrink_fac': bool,
    'skullstrip_shrink_factor_bot_lim':float,
    'skullstrip_avoid_vent': bool,
    'skullstrip_n_iterations': int,
    'skullstrip_pushout': bool,
    'skullstrip_touchup': bool,
    'skullstrip_fill_hole': int,
    'skullstrip_NN_smooth': int,
    'skullstrip_smooth_final': int,
    'skullstrip_avoid_eyes': bool,
    'skullstrip_use_edge': bool,
    'skullstrip_exp_frac': float,
    'skullstrip_push_to_edge': bool,  # check/normalize
    'skullstrip_use_skull': bool,  # check/normalize
    'skullstrip_perc_int': int,
    'skullstrip_max_inter_iter': int,
    'skullstrip_fac': int,
    'skullstrip_blur_fwhm': int,
    'bet_frac': float,
    'bet_mask_boolean': bool,  # check/normalize
    'bet_mesh_boolean': bool,  # check/normalize
    'bet_outline': bool,  # check/normalize
    'bet_padding': bool,  # check/normalize
    'bet_radius': int,
    'bet_reduce_bias': bool,  # check/normalize
    'bet_remove_eyes': bool,  # check/normalize
    'bet_robust': bool,  # check/normalize
    'bet_skull': bool,  # check/normalize
    'bet_surfaces': bool,  # check/normalize
    'bet_threshold': bool,  # check/normalize
    'bet_vertical_gradient': float,

    'regOption': [In(['FSL', 'ANTS'])],  # checck how to transform into a set
    'fnirtConfig': str,
    'ref_mask': str,
    'regWithSkull': [bool], # check/normalize

    'runSegmentationPreprocessing': [bool], # check/normalize
    'priors_path': str,
    'PRIORS_WHITE': str,
    'PRIORS_GRAY': str,
    'PRIORS_CSF': str,

    'slice_timing_correction': [bool], # check/normalize
    'TR': Any(None, float),
    'slice_timing_pattern': Any(str, int), # check for nifti header option
    'startIdx': All(int, Range(min=0)),
    'stopIdx': Any(None, All(int, Range(min=1))),

    'runEPI_DistCorr': [bool], # check/normalize
    'fmap_distcorr_skullstrip': [In(['BET', 'AFNI'])],
    'fmap_distcorr_frac': float, # check if it needs to be a list
    'fmap_distcorr_threshold': float,
    'fmap_distcorr_deltaTE': float, # check if it needs to be a list
    'fmap_distcorr_dwell_time': float, # check if it needs to be a list
    'fmap_distcorr_dwell_asym_ratio': float, # check if it needs to be a list
    'fmap_distcorr_pedir': In(["x", "y", "z", "-x", "-y", "-z"]),

    'runRegisterFuncToAnat': [bool], # check/normalize
    'runBBReg': [bool], # check/normalize
    'boundaryBasedRegistrationSchedule': str,

    'func_reg_input': In(['Mean Functional', 'Selected Functional Volume']),
    'func_reg_input_volume': All(int, Range(min=0)),
    'functionalMasking': [In(['3dAutoMask', 'BET'])],

    'runRegisterFuncToMNI': [bool], # check/normalize
    'resolution_for_func_preproc': All(str, Match(r'^[0-9]+mm$')),
    'resolution_for_func_derivative': All(str, Match(r'^[0-9]+mm$')),
    'template_brain_only_for_func': str,
    'template_skull_for_func': str,
    'identityMatrix': str,
    'configFileTwomm': str,

    'runICA': [bool], # check/normalize
    'aroma_denoise_type': In(['aggr', 'nonaggr']),
    
    'runNuisance': [bool], # check/normalize
    'lateral_ventricles_mask': str,
    'Regressors': [
        {
            'compcor': bool,
            'wm': bool,
            'csf': bool,
            'global': bool,
            'pc1': bool,
            'motion': bool,
            'linear': bool,
            'quadratic': bool,
            'gm': bool,
        }
    ],
    'nComponents': int, # check if list
    'runFristonModel': [bool], # check/normalize
    'runMotionSpike': [Any(None, In(['despiking', 'scrubbing']))], # check/normalize, check None

    'fdCalc': In(['power', 'jenkinson']), # check if it needs to be a list
    'spikeThreshold': float, # check if it needs to be a list
    'numRemovePrecedingFrames': int,
    'numRemoveSubsequentFrames': int,
    'runMedianAngleCorrection': [bool],
    'targetAngleDeg': float,
    'runFrequencyFiltering': [bool],
    'nuisanceBandpassFreq': [[float, float]], # how to check if [0] is > than [1]?

    'runROITimeseries': bool,
    'tsa_roi_paths': Any(None, {
        str: [In(['average', 'voxel', 'spatial_regression', 'pearson_correlation', 'partial_correlation'])],
        # normalize before running thrugh schema
    }),
    'roiTSOutputs': {
        'csv': bool,
        'numpy': bool
    }, # normalize before running thrugh schema

    'runSCA': bool,
    'sca_roi_paths': Any(None, {
        str: [In(['average', 'dual_regression', 'multiple_regression'])],
        # normalize before running thrugh schema
    }),
    'mrsNorm': bool,

    'runVMHC': bool,

    'runALFF': bool,
    'highPassFreqALFF': [float],
    'lowPassFreqALFF': [float],

    'runReHo': bool,
    'clusterSize': Any([7, 19, 27]),

    'runNetworkCentrality': bool,
    'templateSpecificationFile': str,
    'degWeightOptions': {
        'binarized': bool,
        'weighted': bool,
    },
    'degCorrelationThresholdOption': Any(In(["significance", "sparsity", "correlation"])),
    'degCorrelationThreshold': float,

    'eigWeightOptions': {
        'binarized': bool,
        'weighted': bool,
    },
    'eigCorrelationThresholdOption': Any(In(["significance", "sparsity", "correlation"])),
    'eigCorrelationThreshold': float,

    'lfcdWeightOptions': {
        'binarized': bool,
        'weighted': bool,
    },
    'lfcdCorrelationThresholdOption': Any(In(["significance", "sparsity", "correlation"])),
    'lfcdCorrelationThreshold': float,

    'memoryAllocatedForDegreeCentrality': float,

    'run_smoothing': [bool], # check/normalize
    'fwhm': float,
    'smoothing_order': In(['after', 'before']),

    'runZScoring': [bool], # check/normalize

    'run_fsl_feat': bool, # check/normalize
    'numGPAModelsAtOnce': int,
    'modelConfigs': [str],

    'run_basc': bool,
    'basc_resolution': All(str, Match(r'^[0-9]+mm$')),
    'basc_proc': int,
    'basc_memory': float,
    'basc_roi_mask_file': str,
    'basc_cross_cluster_mask_file': str,
    'basc_similarity_metric_list': [In(['correlation', 'euclidean', 'cityblock', 'cosine'])],
    'basc_timeseries_bootstrap_list': int,
    'basc_dataset_bootstrap_list': int,
    'basc_n_clusters_list': int,
    'basc_affinity_thresh': float,
    'basc_output_sizes': int,
    'basc_cross_cluster': bool,
    'basc_blocklength_list': int,
    'basc_group_dim_reduce': bool,
    'basc_inclusion': str,
    'basc_pipeline': str,
    'basc_scan_inclusion': str,

    'runMDMR': bool,
    'mdmr_inclusion': str,
    'mdmr_roi_file': str,
    'mdmr_regressor_file': str,
    'mdmr_regressor_participant_column': str,
    'mdmr_regressor_columns': str,
    'mdmr_permutations': int,
    'mdmr_parallel_nodes': int,

    'runISC': bool,
    'runISFC': bool,
    'isc_voxelwise': bool,
    'isc_roiwise': bool,
    'isc_permutations': int,
})