"""
First Run
"""
# Specify which phenotypic variables and regressors to include in model
# These names should match column names in template_phenotypic.csv
columnsInModel = ['AgeAtScan', 'MeanFD', 'site', 'DxGroup']

# Set a type for each of the variables/regressors entered above.
# 1 = categorical, 0 = directional
categoricalVsDirectional = [0, 0, 1, 1]

# Specify whether to de-mean each column.
deMean = [1, 1, 0, 0]

# Specify the full path to template_phenotypic.csv
phenotypicFile = '/path/to/template_phenotypic.csv'

# Specify the full path to subject_list_group_analysis.txt
subjectListFile = '/path/to/subject_list_group_analysis.txt'

# IGNORE - not yet implemented
groupingVariable = ['DxGroup']

# IGNORE - not yet implemented
modelGroupVariancesSeparately = 0

outputModelFile = '/home/ssikka/nki_nyu_pipeline/model_fsl.csv'

# Specify a name for the model.
modelName = 'my_model'

# Full path to a contrast file based on the columns listed above.
contrastFile = '/home/ssikka/nki_nyu_pipeline/contrasts.csv'
