"""
Configurable Pipeline for the Analysis of Connectomes
=====================================================

CPAC is a configurable, open-source, Nipype-based, automated processing pipeline
for resting state functional MRI (R-fMRI) data, for use by both novice and
expert users.
"""

import anat_preproc, \
       aroma, \
       EPI_DistCorr, \
       connectome, \
       func_preproc, \
       reho, \
       seg_preproc, \
       registration, \
       sca, \
       nuisance, \
       generate_motion_statistics, \
       alff, \
       qc, \
       seg_preproc, \
       vmhc, \
       median_angle, \
       timeseries, \
       network_centrality, \
       scrubbing, \
       group_analysis, \
       randomise, \
       easy_thresh,\
       utils, \
       pipeline, \
       cwas, \
       qpp, \
       GUI

__all__ = ['GUI', 'pipeline', 'anat_preproc', 'func_preproc', 'epi_distcorr',
           'registration', 'seg_preproc', 'reho', 'sca', 'nuisance',
           'alff', 'vmhc', 'median_angle', 'generate_motion_statistics',
           'timeseries', 'network_centrality', 'scrubbing', 'utils',
           'group_analysis','randomise', 'easy_thresh', 'aroma', 'qc', 'qpp']

from .info import __version__
version = __version__
