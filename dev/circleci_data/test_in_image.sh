# install testing requirements
pip install -r /code/dev/circleci_data/requirements.txt

# create conftest.py in root
printf "collect_ignore = ['run.py']\ncollect_ignore_glob = ['**/run.py', '**/test_mdmr_cython.py', '**/test_install.py']" > /code/conftest.py

# run test with coverage
coverage run -m pytest --junitxml=test-results/junit.xml --doctest-ignore-import-errors --continue-on-collection-errors --doctest-modules dev/circleci_data/ /code
