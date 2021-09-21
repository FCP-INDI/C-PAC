# install testing requirements
pip install -r /code/dev/circleci_data/requirements.txt

# run test with coverage as module
python -m coverage run --include */CPAC/*,*/run.py,*/dev/docker_data/* -m pytest --ignore-glob=*test_install.py --junitxml=test-results/junit.xml --doctest-modules dev/circleci_data /code/CPAC

echo "coverage saved to ${COVERAGE_FILE}"