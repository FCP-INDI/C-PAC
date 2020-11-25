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

workingDir=${PWD}

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
installDir=${workingDir}/install

mkdir $buildDir $installDir

cd $buildDir

# USE_VTK must be turned on to build antsSurf
cmake \
    -DCMAKE_INSTALL_PREFIX=$installDir \
    -DBUILD_SHARED_LIBS=OFF \
    -DUSE_VTK=OFF \
    -DSuperBuild_ANTS_USE_GIT_PROTOCOL=OFF \
    -DBUILD_TESTING=OFF \
    -DRUN_LONG_TESTS=OFF \
    -DRUN_SHORT_TESTS=OFF \
    ${workingDir}/ANTs 2>&1 | tee cmake.log

if [[ $? -ne 0 ]]; then
  echo "ANTs SuperBuild configuration failed. Please review documentation at

    $antsBuildInstructions

  If opening an issue, please attach
  
  ${buildDir}/cmake.log
  ${buildDir}/CMakeCache.txt
  ${buildDir}/CMakeFiles/CMakeError.log
  ${buildDir}/CMakeFiles/CMakeOutput.log
  
"
  exit 1 
fi

make -j $buildThreads 2>&1 | tee build.log

if [[ ! -f "CMakeFiles/ANTS-complete" ]]; then
  echo "ANTs compilation failed. Please review documentation at

    $antsBuildInstructions

  If opening an issue, please attach

  ${buildDir}/build.log
  ${buildDir}/cmake.log
  ${buildDir}/CMakeCache.txt
  ${buildDir}/CMakeFiles/CMakeError.log
  ${buildDir}/CMakeFiles/CMakeOutput.log
  
"
  exit 1
fi

cd ANTS-build
make install 2>&1 | tee install.log

antsRegExe="${installDir}/bin/antsRegistration"

if [[ ! -f ${antsRegExe} ]]; then
  echo "Installation failed. Please review documentation at

    $antsBuildInstructions

  If opening an issue, please attach

  ${buildDir}/build.log
  ${buildDir}/cmake.log
  ${buildDir}/CMakeCache.txt
  ${buildDir}/ANTS-build/install.log
  ${buildDir}/CMakeFiles/CMakeError.log
  ${buildDir}/CMakeFiles/CMakeOutput.log

"
  exit 1
fi

echo "Installation complete, running ${antsRegExe}"

${antsRegExe} --version

echo "
Binaries and scripts are located in 

  $installDir

Please see post installation instructions at 

  $antsBuildInstructions

"