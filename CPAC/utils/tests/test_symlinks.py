import os
import tempfile
import pkg_resources as p

from CPAC.utils.symlinks import create_symlinks


mocked_outputs = \
    p.resource_filename(
        "CPAC",
        os.path.join(
            'utils',
            'tests',
            'test_symlinks-outputs.txt'
        )
    )


def test_symlinks():

    temp_dir = tempfile.mkdtemp(suffix='test_symlinks')

    paths = []
    with open(mocked_outputs, 'r') as f:
        for path in f.readlines():
            path = path.strip()
            if path:
                paths += [path]

    create_symlinks(
        temp_dir,
        'sym_links',
        'pipeline_benchmark-FNIRT', '1019436_1', paths
    )

    print("Links created at", temp_dir)

    # TODO test the generated links

    # Normal resource case
    # Several resources within same key case
    # QC case