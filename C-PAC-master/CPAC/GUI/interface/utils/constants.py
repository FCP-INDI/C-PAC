def enum(**enums):
    return type('Enum', (), enums)


control = enum(CHOICE_BOX=0, 
               TEXT_BOX=1, 
               COMBO_BOX=2,
               INT_CTRL=3,
               FLOAT_CTRL = 4,
               DIR_COMBO_BOX = 5,
               CHECKLIST_BOX =6,
               LISTBOX_COMBO = 7,
               TEXTBOX_COMBO = 8,
               CHECKBOX_GRID = 9)

dtype = enum(BOOL=0,
             STR=1,
             NUM=2,
             LBOOL=3,
             LSTR=4,
             LNUM=5,
             LOFL=6,
             COMBO=7,
             LDICT=8) 


substitution_map = {'On': 1,
                    'Off': 0,
                    'On/Off': 10,
                    'ANTS & FSL': 11,
                    '3dAutoMask & BET': 12,
                    'Voxelwise SCA': 'sca_seed_to_standard_smooth',
                    'Voxelwise SCA (z-scored)': 'sca_seed_Z_to_standard_smooth',
                    'ROI Average SCA':'sca_roi_to_standard_smooth',
                    'ROI Average SCA (z-scored)':'sca_roi_Z_to_standard_smooth',
                    'Multiple Regression SCA':'sca_tempreg_maps_z_files_smooth',
                    'ALFF':'alff_to_standard_smooth',
                    'ALFF (z-scored)':'alff_Z_to_standard_smooth',
                    'f/ALFF':'falff_to_standard_smooth',
                    'f/ALFF (z-scored)':'falff_Z_to_standard_smooth',
                    'VMHC':'vmhc_score_stat_map',
                    'VMHC (z-scored)':'vmhc_z_score_stat_map',
                    'ReHo':'reho_to_standard_smooth',
                    'ReHo (z-scored)':'reho_Z_to_standard_smooth',
                    'Network Centrality':'centrality_outputs_smoothed',
                    'Network Centrality (z-scored)': 'centrality_outputs_zscore',
                    'Dual Regression':'dr_tempreg_maps_z_files_smooth',
                    'End': 'None',
                    'ROI Average Time Series Extraction': 'roi_average',
                    'ROI Voxelwise Time Series Extraction': 'roi_voxelwise',
                    'Network Centrality': 'centrality_outputs_smoothed'
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
                      'runGenerateMotionStatistics',
                      'runScrubbing',
                      'runFristonModel']
