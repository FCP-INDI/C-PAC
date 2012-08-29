import os
import sys


DISTNAME = 'CPAC'
DESCRIPTION = 'Configurable Pipeline for the Analysis of Connectomes'
LONG_DESCRIPTION = ''
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
                       delegate_options_to_subpackages=True,
                       quiet=True)
    
    config.add_subpackage('CPAC')

    return config

if __name__ == "__main__":
    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    
    os.chdir(local_path)
    sys.path.insert(0, local_path)
    
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
                      ])
    