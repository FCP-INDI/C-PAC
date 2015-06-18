#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

'''
Configurable Pipeline for the Analysis of Connectomes (C-PAC)
'''

# Package Information
DISTNAME = 'CPAC'
DESCRIPTION = 'Configurable Pipeline for the Analysis of Connectomes'
LONG_DESCRIPTION = ''
MAINTAINER = ''
MAINTAINER_EMAIL = ''
URL = 'http://fcp-indi.github.io/'
LICENSE = 'BSD License'
DOWNLOAD_URL = 'https://github.com/FCP-INDI/C-PAC/tarball/master'

# Import packages
import os, sys

# Import build helpers
try:
    from nisext.sexts import package_check, get_comrec_build
except ImportError:
    raise RuntimeError('Need nisext package from nibabel installation'
                       ' - please install nibabel first')

from build_helpers import INFO_VARS

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
 
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
        
    config.get_version('CPAC/__init__.py')
    config.add_subpackage('CPAC')
     
    # cython
    config.add_extension('CPAC.network_centrality.thresh_and_sum', 
                         sources=['CPAC/network_centrality/thresh_and_sum.pyx'], 
                         include_dirs=[get_numpy_include_dirs()])
    
    return config
 

################################################################################
# For some commands, use setuptools

if len(set(('develop', 'bdist_egg', 'bdist_rpm', 'bdist', 'bdist_dumb',
            'bdist_wininst', 'install_egg_info', 'egg_info', 'easy_install',
            )).intersection(sys.argv)) > 0:
    try:
        from setup_egg import extra_setuptools_args
    except ImportError:
        pass

# extra_setuptools_args can be defined from the line above, but it can
# also be defined here because setup.py has been exec'ed from
# setup_egg.py.
if not 'extra_setuptools_args' in globals():
    extra_setuptools_args = dict()

# Hard and soft dependency checking
package_check('matplotlib', INFO_VARS['MATPLOTLIB_MIN_VERSION'])
package_check('jinja2', INFO_VARS['JINJA_MIN_VERSION'])
#package_check('lockfile', INFO_VARS['PYLOCKFILE_MIN_VERSION']) # checking doesn't really work
package_check('yaml', INFO_VARS['PYYAML_MIN_VERSION'])

################################################################################


def main(**extra_args):
    from numpy.distutils.core import setup    
    from glob import glob
    
    # monkey-patch numpy distutils to use Cython instead of Pyrex
    from numpy.distutils.command.build_ext import build_ext
    from numpy.distutils.command.build_src import build_src
    from build_helpers import generate_a_pyrex_source
    build_src.generate_a_pyrex_source = generate_a_pyrex_source
    cmdclass = {
        'build_src': build_src, 
        'build_ext': build_ext
    }
    
    setup(name=INFO_VARS['NAME'],
          maintainer=INFO_VARS['MAINTAINER'],
          maintainer_email=INFO_VARS['MAINTAINER_EMAIL'],
          description=INFO_VARS['DESCRIPTION'],
          long_description=INFO_VARS['LONG_DESCRIPTION'],
          url=INFO_VARS['URL'],
          download_url=INFO_VARS['DOWNLOAD_URL'],
          license=INFO_VARS['LICENSE'],
          classifiers=INFO_VARS['CLASSIFIERS'],
          author=INFO_VARS['AUTHOR'],
          author_email=INFO_VARS['AUTHOR_EMAIL'],
          platforms=INFO_VARS['PLATFORMS'],
          version=INFO_VARS['VERSION'],
          requires=INFO_VARS['REQUIRES'],
          configuration = configuration,
          cmdclass = cmdclass,
          scripts = glob('scripts/*'),
          entry_points={
            'console_scripts': [
                'cpac_extract_parameters=CPAC.utils.extract_parameters:main'
                ]
            },
          #script_args = ['build_ext', '--inplace'], 
          **extra_args)

# Run main by default
if __name__ == "__main__":

    # Import packages
    import site, shutil
    from shutil import rmtree, copytree

    # Get in the right directory
    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0, local_path)

    # Run main function
    main(**extra_setuptools_args)
