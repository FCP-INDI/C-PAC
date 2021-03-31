# install testing requirements
pip install -r /code/dev/circleci_data/requirements.txt

# run test with coverage as module
python -m coverage run -m pytest --ignore-glob=*test_install.py --junitxml=test-results/junit.xml --doctest-ignore-import-errors --continue-on-collection-errors --doctest-modules dev/circleci_data /code/CPAC
