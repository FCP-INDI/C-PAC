#! /bin/bash
#
# Usage:
# chmod +x /path/to/cpac_install.sh
# sudo /path/to/cpac_install.sh
#
##### Get system info and set variables #####
# --- Create env vars instantiation script ---
# Source env vars
CPAC_ENV=/etc/profile.d/cpac_env.sh
touch $CPAC_ENV
. ~/.bashrc
. $CPAC_ENV
#
echo '########## CHECKING SYSTEM INFO AND SOFTWARE... ##########'
# --- Check for pre-existing software ---
FSLFLAG=`which flirt | grep flirt -c`
AFNIFLAG=`which afni | grep afni -c`
ANTSFLAG=`which ANTS | grep ANTS -c`
C3DFLAG=`which c3d | grep c3d -c`
NIPYPEFLAG=`python -c 'import nipype; print nipype.__version__' | grep 0.9.2 -c`
# --- Get Ubuntu version and architecture ---
. /etc/lsb-release
DIST_CODENAME=$DISTRIB_CODENAME
ARCH64=`uname -a | grep x86_64 -c`
ARCH386=`uname -a | grep i386 -c`
ARCH686=`uname -a |grep i686 -c`
AFNI_DOWNLOAD=linux_openmp
if [ $ARCH64 -ne 0 ]
then
    AFNI_DOWNLOAD=linux_openmp_64
    C3D_DOWNLOAD=c3d-0.8.2-Linux-x86_64
elif [ $ARCH386 -ne 0 ]
then
    C3D_DOWNLOAD=c3d-0.8.2-Linux-i386
elif [ $ARCH686 -ne 0 ]
then
    C3D_DOWNLOAD=c3d-0.8.2-Linux-i686
fi
# FSL repo only has 10.04 LTS distro, anything between 10.04 and 11.10 gets 10.04 distro
if [ "$DISTRIB_CODENAME" == "lucid" ] || [ "$DISTRIB_CODENAME" == "maverick" ] ||  [ "$DISTRIB_CODENAME" == "natty" ] || [ "$DISTRIB_CODENAME" == "oneiric" ]
then
    DIST_CODENAME="lucid"
fi
##### Acquire aptitude-supported packages and dependencies #####
echo '########## UPDATING SOFTWARE VIA APTITUDE... ##########'
# --- Insure software is up-to-date ---
apt-get update
apt-get upgrade -y
# --- Install python package ---
echo '---------- ACQUIRING PYTHON DEPENDENCIES... ----------'
apt-get install -y python-numpy python-scipy python-matplotlib python-networkx python-traits python-wxgtk2.8 python-yaml python-jinja2 python-lockfile python-pygraphviz python-nibabel python-nose cython ipython
# --- Install other utilities/libraries ---
# Git is needed for version control of source code (ANTs)
# Make is needed to compile source code (cmake, ANTs)
# Unzip is needed to unarchive zip files
# CPAC GUI needs libcanberra-gtk module
# AFNI (afni itself) needs: libxp6, netpbm. AFNI tools (e.g. 3dSkullStrip) need: libglu1, gsl-bin
# CMAKE and ANTs require zlib1g-dev to build latest ANTs from source
echo '---------- ACQUIRING NEEDED UTILITIES/LIBRARIES... ----------'
apt-get install -y cmake git make unzip libcanberra-gtk-module libxp6 netpbm libglu1 gsl-bin zlib1g-dev
# --- Sun Grid Engine compatibility (uncomment if you want SGE compatibility) ---
#apt-get install -y libmotif4 nfs-common nfs-kernel-server
#
##### Install needed packages #####
echo '########## INSTALLING SOFTWARE TOOLS... ##########'
# --- Create CPAC-Downloads directory for files storage ---
if [ ! -d CPAC-Downloads ]; then mkdir -p CPAC-Downloads; fi
cd CPAC-Downloads
DOWNLOADS_DIR=`pwd`
#--- Install Nipype 0.9.2 ---
if [ $NIPYPEFLAG -eq 0 ]
then
    echo '---------- INSTALLING NIPYPE... ----------'
    git clone -b 0.10.0 https://github.com/nipy/nipype.git
    cd nipype
    python setup.py install
    cd ..
fi
# --- Install FSL ---
if [ $FSLFLAG -eq 0 ]
then
    echo '---------- INSTALLING FSL... ----------'
    wget -O- http://neuro.debian.net/lists/${DIST_CODENAME}.us-nh.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
    apt-key adv --recv-keys --keyserver pgp.mit.edu 2649A5A9
    apt-get update
    apt-get install -y fsl-5.0-complete
    FSLDIR=/usr/share/fsl/5.0
    echo '# Path to FSL' >> $CPAC_ENV
    echo 'FSLDIR=/usr/share/fsl/5.0' >> $CPAC_ENV
    echo '. ${FSLDIR}/etc/fslconf/fsl.sh' >> $CPAC_ENV
    echo 'PATH=${FSLDIR}/bin:${PATH}' >> $CPAC_ENV
    echo 'export FSLDIR PATH' >> $CPAC_ENV
fi
# Source fsl.sh via symbolic link
#ln -s ${FSLDIR}/etc/fslconf/fsl.sh /etc/profile.d/fsl.sh
# --- Install AFNI ---
if [ $AFNIFLAG -eq 0 ]
then
    echo '---------- INSTALLING AFNI... ----------'
    wget http://afni.nimh.nih.gov/pub/dist/tgz/${AFNI_DOWNLOAD}.tgz
    tar xfz ${AFNI_DOWNLOAD}.tgz
    mv $AFNI_DOWNLOAD /opt/afni
    echo '# Path to AFNI' >> $CPAC_ENV
    echo 'export PATH=/opt/afni:$PATH' >> $CPAC_ENV
    echo 'export DYLD_FALLBACK_LIBRARY_PATH=/opt/afni' >> $CPAC_ENV
fi
# --- Install ANTs ---
if [ $ANTSFLAG -eq 0 ]
then
    echo '---------- INSTALLING ANTs... ----------'
    #wget http://downloads.sourceforge.net/project/advants/ANTS/ANTS_1_9_x/ANTs-1.9.x-Linux.tar.gz
    #tar xfz ANTs-1.9.x-Linux.tar.gz
    #mv ANTs-1.9.x-Linux /opt/ants
    cd $DOWNLOADS_DIR
    git clone -b v2.1.0rc2 https://github.com/stnava/ANTs.git
    mkdir /opt/ants
    cd /opt/ants
    cmake -c -g ${DOWNLOADS_DIR}/ANTs
    make
    ANTSPATH=/opt/ants/bin
    cp ${DOWNLOADS_DIR}/ANTs/Scripts/antsIntroduction.sh ${ANTSPATH}
    cp ${DOWNLOADS_DIR}/ANTs/Scripts/antsAtroposN4.sh ${ANTSPATH}
    cp ${DOWNLOADS_DIR}/ANTs/Scripts/antsBrainExtraction.sh ${ANTSPATH}
    cp ${DOWNLOADS_DIR}/ANTs/Scripts/antsCorticalThickness.sh ${ANTSPATH}
    echo '# Path to ANTS' >> $CPAC_ENV
    echo 'export ANTSPATH=/opt/ants/bin/' >> $CPAC_ENV
    echo 'export PATH=/opt/ants/bin:$PATH' >> $CPAC_ENV
    cd ${DOWNLOADS_DIR}
fi
# Replace antsIntroduction.sh
#wget https://raw.github.com/stnava/ANTs/master/Scripts/antsIntroduction.sh
#mv antsIntroduction.sh $ANTSPATH
#chmod 755 ${ANTSPATH}/antsIntroduction.sh
# --- Install C3D ---
if [ $C3DFLAG -eq 0 ]
then
    echo '---------- INSTALLING C3D... ----------'
    wget http://sourceforge.net/projects/c3d/files/c3d/c3d-0.8.2/${C3D_DOWNLOAD}.tar.gz
    tar xfz ${C3D_DOWNLOAD}.tar.gz
    mv $C3D_DOWNLOAD /opt/c3d
    echo '# Path to C3D' >> $CPAC_ENV
    echo 'export PATH=/opt/c3d/bin:$PATH' >> $CPAC_ENV
fi
#
##### Download and CPAC image resources #####
echo '########## ACQUIRING CPAC IMAGE RESOURCES... ##########'
wget http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tgz
tar -xzvf cpac_resources.tgz
cd cpac_image_resources
./install_resources.sh $FSLDIR
#
##### Install CPAC #####
echo '########## INSTALLING CPAC... ##########'
cd ../..
wget https://github.com/FCP-INDI/C-PAC/archive/master.zip
unzip master.zip
cd C-PAC-master
python setup.py install
cd ../..
#
##### Cleanup and resource #####
echo '########## CLEANING UP... ##########'
# --- Remove unnexessary files ---
apt-get autoremove -y
# --- Re-source env vars and exit ---
if [ -d ~/.matplotlib ]; then
    chmod -R 777 ~/.matplotlib
fi
if [ -d ~/.config/matplotlib ]; then
    chmod -R 777 ~/.config/matplotlib
fi
exit
. /etc/profile.d/cpac_env.sh
#
echo '########## DONE! ##########'
echo 'TO BEGIN USING CPAC, OPEN A NEW TERMINAL WINDOW AND EXECUTE \"cpac_gui\"'
