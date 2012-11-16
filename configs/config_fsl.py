# emacs =  -*- mode =  python; py-indent-offset =  4; indent-tabs-mode =  nil -*-
# vi =  set ft=python sts=4 ts=4 sw=4 et =
columnsInModel = ['AgeAtScan', 'MeanFD', 'site', 'DxGroup']

phenotypicFile = '/home/ssikka/nki_nyu_pipeline/phenotypic.csv'
subjectListFile = '/home/ssikka/nki_nyu_pipeline/subjects_gp.txt'

"""
	Categorical :  1
	Directional :  0
"""
categoricalVsDirectional = [0, 0, 1, 1]
deMean = [1, 1, 0, 0]

groupingVariable = ['DxGroup']
modelGroupVariancesSeparately = 0

outputModelFile = '/home/ssikka/nki_nyu_pipeline/model_fsl.csv'
modelName = 'my_model'
contrastFile = '/home/ssikka/nki_nyu_pipeline/contrasts.csv'
