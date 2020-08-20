pip install -r /code/dev/circleci_data/requirements.txt
coverage run -m pytest --junitxml=test-results/junit.xml --doctest-ignore-import-errors --continue-on-collection-errors --doctest-modules --ignore=CPAC/cwas/tests/test_mdmr_cython.py --ignore=dev/circleci_data/test_install.py dev/circleci_data/ /code
