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
 
print "Version: ", VERSION
 
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
    
    import site
    import shutil
    from shutil import rmtree
    from shutil import copytree

    # Temporary code to remove pre-existing CPAC directories
    # before re-installing - needs more elegant solution

    for sitePath in site.getsitepackages():

        for root,dirs,files in os.walk(sitePath):
            
            if 'CPAC-backup' in root:
                shutil.rmtree(root)


    for sitePath in site.getsitepackages():

        for root,dirs,files in os.walk(sitePath):
            
            if 'CPAC' in root:
                backupPath = sitePath + '/CPAC-backup'

                shutil.copytree(root,backupPath)
                shutil.rmtree(root)
                print "Backing up pre-existing CPAC directory into ", backupPath
                print "Removing directory ", root

 
    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
 
    
    os.chdir(local_path)
    sys.path.insert(0, local_path)
     
    setup(configuration=configuration,
          scripts=['scripts/log_py2js.py', 'scripts/cpac_gui'],
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

