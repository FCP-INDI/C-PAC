"""The C-PAC pipeline and its underlying infrastructure"""
import os
import pkg_resources as p

AVAILABLE_PIPELINE_CONFIGS = os.listdir(
    p.resource_filename("CPAC", os.path.join("resources", "configs")))
AVAILABLE_PIPELINE_CONFIGS = [x.split('_')[2].replace('.yml', '') for
                              x in AVAILABLE_PIPELINE_CONFIGS if
                              'pipeline_config' in x]

__all__ = ['AVAILABLE_PIPELINE_CONFIGS']
