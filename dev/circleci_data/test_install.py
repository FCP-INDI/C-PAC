import os

from pathlib import Path
from spython.main import Client


def test_AFNI_libraries():
    SINGULARITY_IMAGE_PATH = '/home/circleci/project/C-PAC-CI.simg'
    if not os.path.exists(SINGULARITY_IMAGE_PATH):
        try:
            SINGULARITY_IMAGE_PATH = [d for d in os.listdir(
                str(Path(__file__).parent.parent.parent)
            ) if (d.endswith('.simg') or d.endswith('.sif'))][0]
        except:
            raise Exception("Singularity image not in expected location.")
    if os.path.exists(SINGULARITY_IMAGE_PATH):
        afni_libraries = Client.execute(
            Client.instance(SINGULARITY_IMAGE_PATH),
            ['./dev/circleci_data/AFNI_libraries.sh']
        )
        assert "not found" not in afni_libraries, '\n'.join([
            line.strip() for line in afni_libraries.split(
                '\n'
            ) if "not found" in line
        ])