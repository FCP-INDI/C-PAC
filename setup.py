
DISTNAME = 'CPAC'
DESCRIPTION = 'Configurable Pipeline for the Analysis of Connectomes'
LONG_DESCRIPTION = open('README.md').read()
MAINTAINER = ''
MAINTAINER_EMAIL = ''
URL = ''
LICENSE = ''
DOWNLOAD_URL = ''

import CPAC
VERSION = CPAC.__version__

from numpy.distutils.core import setup

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackage=True,
                       quiet=True)
    
    config.add_subpackage('CPAC')

    return config

if __name__ == "__main__":
     setup(configuration=configuration,
           name=DISTNAME,
           maintainer=MAINTAINER,
           maintainer_email=MAINTAINER_EMAIL,
           description=DESCRIPTION,
           license=LICENSE,
           url=URL,
           version=VERSION,
           download_url=DOWNLOAD_URL,
           long_description=LONG_DESCRIPTION,
           classifiers=['Intended Audience :: Science/Research',
                        'Programming Language :: Python'
                        ]
           )