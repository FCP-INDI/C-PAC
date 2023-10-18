export PATH=$PATH:/home/$(whoami)/.local/bin

# install testing requirements
pip install -r /code/dev/circleci_data/requirements.txt

# run test with coverage as module
python -m coverage run --include $(circleci tests glob "*/CPAC/*,*/run.py,*/dev/docker_data/*" | circleci tests split --split-by=timings) -m pytest --ignore-glob="*test_install.py" --junitxml=test-results/junit.xml --doctest-modules dev/circleci_data /code/CPAC

echo "$?" > test-results/exitcode

echo "coverage saved to ${COVERAGE_FILE}"

exit $(cat test-results/exitcode)