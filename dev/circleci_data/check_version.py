import os
import subprocess

from sys import path


def main():
    """
    Function to get the previous version for dev versioning.

    Parameters
    ----------
    None

    Returns
    -------
    str
    """
    CPAC_PATH = os.path.abspath(os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    ))
    path.append(os.path.join(CPAC_PATH, 'CPAC'))

    from info import _version_major, _version_minor, _version_micro

    # get previous git tag for comparison
    o = subprocess.Popen('git describe HEAD~1 --always', shell=True,
                         cwd=os.path.join(CPAC_PATH, '.git'),
                         stdout=subprocess.PIPE).communicate(
                         )[0].decode().strip().split('-')[-1]

    return '.'.join([str(i) for i in [_version_major, _version_minor,
                                      _version_micro, f'{o}-dev']])


if __name__ == '__main__':
    print(main())
