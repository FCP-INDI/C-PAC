# Copyright (C) 2022-2025  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
import os


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True,
    )

    config.get_version("CPAC/__init__.py")

    for root, _, _ in list(os.walk("CPAC")):
        if os.path.isfile(os.path.join(root, "__init__.py")):
            config.add_subpackage(root)

    config.add_data_dir("CPAC/resources")

    config.add_define_macros(["NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"])

    return config


def main(**extra_args):
    from glob import glob
    import sys

    from numpy.distutils.core import setup

    # prepend this file's parent directory to sys.path to find CPAC when building
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
    from CPAC.info import (
        AUTHOR,
        AUTHOR_EMAIL,
        CLASSIFIERS,
        DESCRIPTION,
        DOWNLOAD_URL,
        LICENSE,
        LONG_DESCRIPTION,
        MAINTAINER,
        MAINTAINER_EMAIL,
        NAME,
        PLATFORMS,
        REQUIREMENTS,
        URL,
        VERSION,
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
        extras_require={"graphviz": ["pygraphviz"]},
        configuration=configuration,
        scripts=glob("scripts/*"),
        entry_points={
            "console_scripts": [
                "cpac = CPAC.__main__:main",
                "C-PAC_nb_io = CPAC.pipeline.nb_io:main",
            ]
        },
        package_data={
            "CPAC": [
                "test_data/*",
                "test/templates/*",
                "qc/colors/*.txt",
                "qc/data/index.html",
            ]
        },
        **extra_args,
    )


if __name__ == "__main__":
    main()
