import os, sys

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.get_version('CPAC/__init__.py')

    for root, _, _ in list(os.walk('CPAC')):
        if os.path.isfile(os.path.join(root, '__init__.py')):
            config.add_subpackage(root)

    config.add_data_dir('CPAC/GUI/resources')
    config.add_data_dir('CPAC/resources')

    config.add_define_macros([
        "NPY_NO_DEPRECATED_API",
        "NPY_1_7_API_VERSION"
    ])

    return config

def main(**extra_args):
    from numpy.distutils.core import setup
    from glob import glob
    from CPAC.info import (
        NAME,
        MAINTAINER,
        MAINTAINER_EMAIL,
        DESCRIPTION,
        LONG_DESCRIPTION,
        URL,
        DOWNLOAD_URL,
        LICENSE,
        CLASSIFIERS,
        AUTHOR,
        AUTHOR_EMAIL,
        PLATFORMS,
        VERSION,
        REQUIREMENTS,
    )

    setup(
        name=NAME,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        url=URL,
        download_url=DOWNLOAD_URL,
        license=LICENSE,
        classifiers=CLASSIFIERS,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        platforms=PLATFORMS,
        version=VERSION,
        install_requires=REQUIREMENTS,
        configuration=configuration,
        scripts=glob('scripts/*'),
        entry_points={
            'console_scripts': [
                'cpac = CPAC.__main__:main'
            ]
        },
        package_data={
            'CPAC': [
                'test_data/*',
                'test/templates/*',
                'qc/colors/*.txt'
            ]
        },
        **extra_args
    )

if __name__ == "__main__":
    main()