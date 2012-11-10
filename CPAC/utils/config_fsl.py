"""
First Run
"""
# Specify which phenotypic variables and regressors to include in model
# These names should match column names in template_phenotypic.csv
columnsInModel = ['AgeAtScan', 'MeanFD', 'site', 'DxGroup']

# Specify the full path to template_phenotypic.csv
phenotypicFile = '/path/to/template_phenotypic.csv'

# Specify the full path to subject_list_group_analysis.txt
subjectListFile = '/path/to/subject_list_group_analysis.txt'

"""
	Categorical :  1
	Directional :  0
"""

# Set a type for each of the variables/regressors entered above.
categoricalVsDirectional = [0, 0, 1, 1]

deMean = [1, 1, 0, 0]


groupingVariable = ['DxGroup']

# IGNORE - not yet implemented
modelGroupVariancesSeparately = 0

outputModelFile = '/home/ssikka/nki_nyu_pipeline/model_fsl.csv'

modelName = 'my_model'

contrastFile = '/home/ssikka/nki_nyu_pipeline/contrasts.csv'
