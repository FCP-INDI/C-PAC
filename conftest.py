import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--singularity_image", action="store", default="/home/circleci/project/C-PAC-CI.simg", help="path to Singularity image"
    )


@pytest.fixture
def singularity_image_path(request):
    return request.config.getoption("--singularity_image")
