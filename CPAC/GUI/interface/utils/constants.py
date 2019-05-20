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
             STR = 1,
             NUM = 2,
             LBOOL = 3,
             LSTR = 4,
             LNUM = 5,
             LOFL = 6,
             COMBO = 7,
             LDICT = 8 ) 

substitution_map = {'On': 1,
                    'Off': 0,
                    'On/Off': 10,
                    'ANTS & FSL': 11,
                    '3dAutoMask & BET': 12,
                    'AFNI & BET' : 12,
                    'ALFF':'alff',
                    'f/ALFF':'falff',
                    'ReHo':'reho',
                    'ROI Average SCA':'sca_roi',
                    'Multiple Regression SCA':'sca_tempreg',
                    'VMHC':'vmhc',
                    'Network Centrality':'centrality',
                    'Dual Regression':'dr_tempreg',
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
                      'runMedianAngleCorrection',
                      'runScrubbing',
                      'runFristonModel'
                      'runEPI_DistCorr']
