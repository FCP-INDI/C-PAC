"""The C-PAC pipeline and its underlying infrastructure"""
import os
import pkg_resources as p

ALL_PIPELINE_CONFIGS = os.listdir(
    p.resource_filename("CPAC", os.path.join("resources", "configs")))
ALL_PIPELINE_CONFIGS = [x.split('_')[2].replace('.yml', '') for
                        x in ALL_PIPELINE_CONFIGS if 'pipeline_config' in x]
ALL_PIPELINE_CONFIGS.sort()
AVAILABLE_PIPELINE_CONFIGS = [preconfig for preconfig in ALL_PIPELINE_CONFIGS
                              if not preconfig.startswith('benchmark-') and
                              not preconfig.startswith('regtest-')]

__all__ = ['ALL_PIPELINE_CONFIGS', 'AVAILABLE_PIPELINE_CONFIGS']
