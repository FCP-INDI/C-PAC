#!/bin/bash

antsBuildInstructions="https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS"

echo "
This script will download ANTs, build and install under the current directory. 

Developer tools including compilers, git and cmake must be installed

If you encounter errors, please see the installation instructions at

  $antsBuildInstructions
"

if [ $# -ne 1 ]
then 
    echo "
Takes one commandline argument, the version of ANTs to install.
"
exit 0
fi

echo "
Build will proceed in 5 seconds
"

sleep 5

mkdir -p /antsinstall

workingDir=/antsinstall
cd ${workingDir}
# Clone the repo
git clone https://github.com/ANTsX/ANTs.git

# If you want to build a particular release, do so here
cd ANTs
git checkout $1
cd -

# Number of threads used by make
buildThreads=4

# Where to build, should be an empty directory
buildDir=${workingDir}/build
installDir=/usr/lib/ants

mkdir -p $buildDir $installDir

cd $buildDir

cmake ${workingDir}/ANTs -DCMAKE_INSTALL_PREFIX=${installDir}
make 2>&1 | tee build.log
cd ANTS-build
make install 2>&1 | tee install.log
