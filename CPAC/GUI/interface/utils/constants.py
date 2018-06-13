def enum(**enums):
    return type('Enum', (), enums)


control = enum(CHOICE_BOX = 0, 
               TEXT_BOX = 1, 
               COMBO_BOX = 2,
               INT_CTRL = 3,
               FLOAT_CTRL = 4,
               DIR_COMBO_BOX = 5,
               CHECKLIST_BOX = 6,
               LISTBOX_COMBO = 7,
               TEXTBOX_COMBO = 8,
               CHECKBOX_GRID = 9,
               GPA_CHECKBOX_GRID = 10,
               SPIN_BOX_FLOAT = 11)

dtype = enum(BOOL = 0,
             STR= 1,
             NUM= 2,
             LBOOL = 3,
             LSTR = 4,
             LNUM = 5,
             LOFL = 6,
             COMBO = 7,
             LDICT= 8 ) 

substitution_map = {'On': 1,
                    'Off': 0,
                    'On/Off': 10,
                    'ANTS & FSL': 11,
                    '3dAutoMask & BET': 12,
                    'AFNI & BET' : 12,
                    'ALFF':'alff_to_standard_zstd',
                    'ALFF (smoothed)':'alff_to_standard_zstd_smooth',
                    'f/ALFF':'falff_to_standard_zstd',
                    'f/ALFF (smoothed)':'falff_to_standard_zstd_smooth',
                    'ReHo':'reho_to_standard_zstd',
                    'ReHo (smoothed)':'reho_to_standard_smooth_zstd',
                    'ROI Average SCA':'sca_roi_files_to_standard_fisher_zstd',
                    'ROI Average SCA (smoothed)':'sca_roi_files_to_standard_smooth_fisher_zstd',
                    'Multiple Regression SCA':'sca_tempreg_maps_zstat_files',
                    'Multiple Regression SCA (smoothed)':'sca_tempreg_maps_zstat_files_smooth',
                    'VMHC':'vmhc_fisher_zstd_zstat_map',
                    'Network Centrality':'centrality_outputs_zstd',
                    'Network Centrality (smoothed)': 'centrality_outputs_smoothed_zstd',
                    'Dual Regression':'dr_tempreg_maps_zstat_files_to_standard',
                    'Dual Regression (smoothed)':'dr_tempreg_maps_zstat_files_to_standard_smooth',
                    'ROI Average Time Series Extraction': 'roi_average',
                    'ROI Voxelwise Time Series Extraction': 'roi_voxelwise',
                   }

multiple_value_wfs = ['runAnatomicalPreprocessing',
                      'runFunctionalPreprocessing',
                      'runRegistrationPreprocessing',
                      'runRegisterFuncToMNI',
                      'runAnatomicalToFunctionalRegistration',
                      'runSegmentationPreprocessing',
                      'runNuisance',
                      'runFrequencyFiltering',
                      'runMedianAngleCorrection',
                      'runScrubbing',
                      'runFristonModel'
                      'runEPI_DistCorr']
