#!/bin/bash

# Copyright (C) 2024  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.

# Update version comment strings
function wait_for_git_lock() {
    while [ -f "./.git/index.lock" ]; do
        echo "Waiting for the git lock file to be removed..."
        sleep 1
    done
}

cd CPAC || exit 1
VERSION=$(python -c "from info import __version__; print(('.'.join(('.'.join(__version__[::-1].split('-')[1].split('.')[1:])[::-1], __version__.split('-')[1])) if '-' in __version__ else __version__).split('+', 1)[0])")
cd ..
echo "v${VERSION}" > version
export _SED_COMMAND="s/^(# [Vv]ersion ).*$/# Version ${VERSION}/g"
if [[ "$OSTYPE" == "darwin"* ]]; then
    # Mac OSX
    find ./CPAC/resources/configs -name "*.yml" -exec sed -i '' -E "${_SED_COMMAND}" {} \;
else
    # Linux and others
    find ./CPAC/resources/configs -name "*.yml" -exec sed -i'' -r "${_SED_COMMAND}" {} \;
fi
wait_for_git_lock && git add version
VERSIONS=( `git show $(git log --pretty=format:'%h' -n 1 version | tail -n 1):version` `cat version` )
export PATTERN="(declare|typeset) -a"
if [[ "$(declare -p VERSIONS)" =~ $PATTERN ]]
then
  for DOCKERFILE in $(find ./.github/Dockerfiles -name "*.Dockerfile")
  do
    export IFS=""
    for LINE in $(grep "FROM ghcr\.io/fcp\-indi/c\-pac/.*\-${VERSIONS[0]}" ${DOCKERFILE})
    do
      echo "Updating stage tags in ${DOCKERFILE}"
      if [[ "$OSTYPE" == "darwin"* ]]; then
          # Mac OSX
          sed -i "" "s/\-${VERSIONS[0]}/\-${VERSIONS[1]}/g" ${DOCKERFILE}
      else
          # Linux and others
          sed -i "s/\-${VERSIONS[0]}/\-${VERSIONS[1]}/g" ${DOCKERFILE}
      fi
    done
  done
  unset IFS
fi
wait_for_git_lock && git add CPAC/resources/configs .github/Dockerfiles

# Overwrite top-level Dockerfiles with the CI Dockerfiles
wait_for_git_lock && cp .github/Dockerfiles/C-PAC.develop-jammy.Dockerfile Dockerfile
wait_for_git_lock && cp .github/Dockerfiles/C-PAC.develop-ABCD-HCP-bionic.Dockerfile variant-ABCD-HCP.Dockerfile
wait_for_git_lock && cp .github/Dockerfiles/C-PAC.develop-fMRIPrep-LTS-xenial.Dockerfile variant-fMRIPrep-LTS.Dockerfile
wait_for_git_lock && cp .github/Dockerfiles/C-PAC.develop-lite-jammy.Dockerfile variant-lite.Dockerfile
for DOCKERFILE in $(ls *Dockerfile)
do
  wait_for_git_lock && git add $DOCKERFILE
done
