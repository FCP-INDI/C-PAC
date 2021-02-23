import os

from pathlib import Path
from spython.main import Client


def test_AFNI_libraries(singularity_image_path):
    if not os.path.exists(singularity_image_path):
        try:
            singularity_image_path = [d for d in os.listdir(
                str(Path(__file__).parent.parent.parent)
            ) if (d.endswith('.simg') or d.endswith('.sif'))][0]
        except Exception:
            raise FileNotFoundError(
                "Singularity image not in expected location.")
    if os.path.exists(singularity_image_path):
        afni_libraries = Client.execute(
            Client.instance(singularity_image_path),
            ['./dev/circleci_data/AFNI_libraries.sh']
        )
        assert "not found" not in afni_libraries, '\n'.join([
            line.strip() for line in afni_libraries.split(
                '\n'
            ) if "not found" in line
        ])
