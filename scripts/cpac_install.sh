#! /bin/bash

# cpac_provision.sh
# =================================================================================================
# Version: 0.3.9
# Author(s): John Pellman, Daniel Clark
# Based off of cpac_install.sh by Daniel Clark.
# Description: Will perform specific operations to install C-PAC dependencies and C-PAC.
# Checks for user privileges and performs installation either locally or system-wide.
# Can be customized using flags.
# =================================================================================================
# Flags:
# -s : System-level dependencies only.
# -p : Python dependencies only
# -n : Install specific neuroimaging packages.  Accepts any number of the following as arguments:
#	afni, fsl, c3d, ants, cpac
#	will issue warnings if dependencies for these neuroimaging packages are not fulfilled.
#	If multiple packages are to be specified, they must be surrounded by quotation marks.
# -a : Install all neuroimaging suites not already installed.  Will also tell you if all neuroimaging suites are already installed and on the path.
# -l : Local install. Equivalent to -pa ; will not run FSL installer, but will issue a warning if running on Ubuntu. 
# -r : Root install.  Equivalent to -spa
# -h : Bring up the help dialog.
# =================================================================================================
# Example usage:
#	cpac_provision.sh -n "fsl afni"
#	Will install FSL and AFNI.  The list of neuroimaging suites to install is iterated through sequentially.
#	In this case, FSL would first be installed before AFNI.
# TODO: Use Juju for local installations. Prompt user to ask if they would like to do an entirely local install.

function install_system_dependencies {
	echo "Installing C-PAC system dependencies..."
	system_dependencies_installed ; if [ $? -eq 0 ]; then
		echo System dependencies are already installed!
		echo Moving on...
		echo '[ '$(date)' ] : C-PAC system dependencies are already installed, do not need to be re-installed.' >> ~/cpac.log
		return
	fi
	if [ $LOCAL -eq 0 ]; then
		if [ $DISTRO == 'CENTOS' ]; then
			yum update -y
			cd /tmp && wget http://dl.fedoraproject.org/pub/epel/7/x86_64/e/epel-release-7-5.noarch.rpm && rpm -Uvh epel-release-7-5.noarch.rpm 
			yum install -y cmake git make unzip netpbm gcc python-devel gcc-gfortran gcc-c++ libgfortran lapack lapack-devel blas libcanberra-gtk2 libXp.x86_64 mesa-libGLU-9.0.0-4.el7.x86_64 gsl-1.15-13.el7.x86_64 wxBase wxGTK wxGTK-gl wxPython graphviz graphviz-devel.x86_64 zlib-devel
			yum autoremove -y
		elif [ $DISTRO == 'UBUNTU' ]; then
			apt-get update
			apt-get upgrade -y
			apt-get install -y cmake git make unzip libcanberra-gtk-module libxp6 netpbm libglu1-mesa gsl-bin zlib1g-dev graphviz graphviz-dev pkg-config build-essential
			apt-get autoremove -y
		else
			echo Linux distribution not recognized.  System-level dependencies cannot be installed.
			echo '[ '$(date)' ] : C-PAC system dependencies could not be installed (Linux distribution not recognized).' >> ~/cpac.log
			cd $INIT_DIR
			exit 1
		fi	
	elif [ $LOCAL -eq 1 ]; then
		echo System-level dependencies cannot be installed since you do not have root privileges.
		echo Re-run this script as root or have your system administrator run it.
		cd $INIT_DIR
		echo '[ '$(date)' ] : C-PAC system dependencies could not be installed (not root).' >> ~/cpac.log
		exit 1
	else
		echo Invalid value for variable 'LOCAL'.
		echo This script is unable to determine whether or not you are running it as root.
		echo '[ '$(date)' ] : C-PAC system dependencies could not be installed (unable to determine if root).' >> ~/cpac.log
		cd $INIT_DIR
		exit 1
	fi
	echo '[ '$(date)' ] : C-PAC system dependencies succesfully installed.' >> ~/cpac.log
}

function system_dependencies_installed {
	if [ $DISTRO == 'CENTOS' ]; then
		for package in git make unzip netpbm gcc python-devel gcc-gfortran gcc-c++ libgfortran lapack lapack-devel blas libcanberra-gtk2 libXp.x86_64 mesa-libGLU-9.0.0-4.el7.x86_64 gsl-1.15-13.el7.x86_64 wxBase wxGTK wxGTK-gl wxPython graphviz graphviz-devel.x86_64; do
			yum list installed ${package} > /dev/null 2>&1
		done
	elif [ $DISTRO == 'UBUNTU' ]; then
		dpkg -s cmake git make unzip libcanberra-gtk-module libxp6 netpbm libglu1-mesa gsl-bin zlib1g-dev graphviz graphviz-dev pkg-config > /dev/null 2>&1
	fi
	return $?
}

function install_python_dependencies {
	echo "Installing C-PAC Python dependencies..."
	python_dependencies_installed ; if [ $? -eq 0 ]; then
		echo Python dependencies are already installed!
		echo Moving on...
		echo '[ '$(date)' ] : Python dependencies are already installed - do not need to be re-installed.' >> ~/cpac.log
		return
	fi
	system_dependencies_installed ; if [ $? -ne 0 ]; then
		echo Python dependencies cannot be installed unless system-level dependencies are installed first.
		echo Have your system administrator install system-level dependencies as root.
		echo Exiting now...
		echo '[ '$(date)' ] : Python dependencies could not be installed (system-level dependencies not installed.' >> ~/cpac.log
		cd $INIT_DIR
		exit 1
	fi
	cd /tmp 
	wget http://repo.continuum.io/miniconda/Miniconda-3.8.3-Linux-x86_64.sh 
	chmod +x Miniconda-3.8.3-Linux-x86_64.sh
	if [ $LOCAL -eq 0 ]; then
		./Miniconda-3.8.3-Linux-x86_64.sh -b -p /usr/local/bin/miniconda
		chmod -R 775 /usr/local/bin/miniconda
		chmod g+s /usr/local/bin/miniconda
		export PATH=/usr/local/bin/miniconda/bin:${PATH}
		echo 'export PATH=/usr/local/bin/miniconda/bin:${PATH}' >> ~/cpac_env.sh
	elif [ $LOCAL -eq 1 ] && [ ! -d ~/miniconda ]; then
		./Miniconda-3.8.3-Linux-x86_64.sh -b
		export PATH=~/miniconda/bin:${PATH}
		echo 'export PATH=~/miniconda/bin:${PATH}' >> ~/cpac_env.sh
	fi
	if [ ! -d ~/miniconda/envs/cpac ] || [ ! -d /usr/local/bin/miniconda/envs/cpac ]; then
		conda create -y -n cpac python
		source activate cpac
		conda install -y cython numpy scipy matplotlib networkx traits pyyaml jinja2 nose ipython pip wxpython
 		pip install lockfile pygraphviz nibabel nipype patsy memory_profiler psutil
		source deactivate
		echo 'source activate cpac' >> ~/cpac_env.sh
	fi
}

function python_dependencies_installed {
	if [ ! -d ~/miniconda/envs/cpac ] || [ ! -d /usr/local/bin/miniconda/envs/cpac ]; then
		return 1
	fi
	source activate cpac &> /dev/null
	python -c "import cython, numpy, scipy, matplotlib, networkx, traits, yaml, jinja2, nose, pip, lockfile, pygraphviz, nibabel, nipype, wx" 2> /dev/null && which ipython &> /dev/null
	status=$?
	source deactivate &> /dev/null
	return $status
}

function install_fsl {
	echo "Installing FSL."
	which fsl &> /dev/null ; if [ $? -eq 0 ]; then
		echo FSL is already installed!
		echo Moving on...
		echo '[ '$(date)' ] : FSL is already installed - does not need to be re-installed.' >> ~/cpac.log
		return
	fi
	system_dependencies_installed ; if [ $? -ne 0 ]; then
		echo FSL cannot be installed unless system-level dependencies are installed first.
		echo Have your system administrator install system-level dependencies as root.
		echo Exiting now...
		echo '[ '$(date)' ] : FSL installation failed - system-level dependencies are not installed.' >> ~/cpac.log
		cd $INIT_DIR
		exit 1
	fi
	if [ $DISTRO == 'CENTOS' ]; then
    		cd /tmp
		wget fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py
	fi
	if [ $LOCAL -eq 0 ]; then
		if [ $DISTRO == 'CENTOS' ]; then
			python fslinstaller.py -d /usr/share
			FSLDIR=/usr/share/fsl/
			mkdir $FSLDIR/5.0
			mv $FSLDIR/bin $FSLDIR/5.0/bin
			ln -s $FSLDIR/data $FSLDIR/5.0/data
			mv $FSLDIR/doc $FSLDIR/5.0/doc
			mv $FSLDIR/etc $FSLDIR/5.0/etc
			mv $FSLDIR/tcl $FSLDIR/5.0/tcl
		# Debian-based distros must use NeuroDebian instead of the installer.
		elif [ $DISTRO == 'UBUNTU' ]; then
			wget -O- http://neuro.debian.net/lists/$(lsb_release -cs).us-nh.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
			apt-key adv --recv-keys --keyserver pgp.mit.edu 2649A5A9
			apt-get update
			apt-get install -y fsl-5.0-complete
		fi
		FSLDIR=/usr/share/fsl/5.0
		. ${FSLDIR}/etc/fslconf/fsl.sh
		PATH=${FSLDIR}/bin:${PATH}
		export FSLDIR PATH
		echo '# Path to FSL' >> ~/cpac_env.sh
		echo 'FSLDIR=/usr/share/fsl/5.0' >> ~/cpac_env.sh
		echo '. ${FSLDIR}/etc/fslconf/fsl.sh' >> ~/cpac_env.sh
		echo 'PATH=${FSLDIR}/bin:${PATH}' >> ~/cpac_env.sh
		echo 'export FSLDIR PATH' >> ~/cpac_env.sh
	elif [ $LOCAL -eq 1 ]; then
		if [ $DISTRO == 'CENTOS' ]; then
                        python fslinstaller.py -d ~ 
                        FSLDIR=~/fsl/
                        mkdir $FSLDIR/5.0
                        mv $FSLDIR/bin $FSLDIR/5.0/bin
                        ln -s $FSLDIR/data $FSLDIR/5.0/data
                        mv $FSLDIR/doc $FSLDIR/5.0/doc
                        mv $FSLDIR/etc $FSLDIR/5.0/etc
                        mv $FSLDIR/tcl $FSLDIR/5.0/tcl
			FSLDIR=~/fsl/5.0
			. ${FSLDIR}/etc/fslconf/fsl.sh
			PATH=${FSLDIR}/bin:${PATH}
			export FSLDIR PATH
			echo '# Path to FSL' >> ~/cpac_env.sh
                	echo 'FSLDIR=~/fsl/5.0' >> ~/cpac_env.sh
                	echo '. ${FSLDIR}/etc/fslconf/fsl.sh' >> ~/cpac_env.sh
                	echo 'PATH=${FSLDIR}/bin:${PATH}' >> ~/cpac_env.sh
                	echo 'export FSLDIR PATH' >> ~/cpac_env.sh
		elif [ $DISTRO == 'UBUNTU' ]; then
			echo FSL cannot be installed without root privileges on Ubuntu Linux.
			echo '[ '$(date)' ] : FSL installation failed - need root privileges on Ubuntu.' >> ~/cpac.log
			cd $INIT_DIR
			install_cpac_env
			exit 1
		fi
	else
		echo Invalid value for variable 'LOCAL'.
		echo This script is unable to determine whether or not you are running it as root.
		echo '[ '$(date)' ] : FSL could not be installed (unable to determine if root).' >> ~/cpac.log
		cd $INIT_DIR
		exit 1
	fi
}

function install_afni {
	echo "Installing AFNI."
	which afni &> /dev/null ; if [ $? -eq 0 ]; then
		echo AFNI is already installed!
		echo Moving on...
		echo '[ '$(date)' ] : AFNI is already installed - does not need to be re-installed.' >> ~/cpac.log
		return
	fi
	system_dependencies_installed ; if [ $? -ne 0 ]; then
		echo AFNI cannot be installed unless system-level dependencies are installed first.
		echo Have your system administrator install system-level dependencies as root.
		echo Exiting now...
		echo '[ '$(date)' ] : AFNI installation failed - system-level dependencies are not installed.' >> ~/cpac.log
		cd $INIT_DIR
		exit 1
	fi
	cd /tmp
	if [ $(uname -p) == 'x86_64' ]; then
		AFNI_DOWNLOAD=linux_openmp_64
	else
		AFNI_DOWNLOAD=linux_openmp
	fi
    	wget http://afni.nimh.nih.gov/pub/dist/tgz/${AFNI_DOWNLOAD}.tgz
    	tar xfz ${AFNI_DOWNLOAD}.tgz
	if [ $LOCAL -eq 0 ]; then
    		mv ${AFNI_DOWNLOAD} /opt/afni
		export PATH=/opt/afni:$PATH
		export DYLD_FALLBACK_LIBRARY_PATH=/opt/afni
    		echo '# Path to AFNI' >> ~/cpac_env.sh
   		echo 'export PATH=/opt/afni:$PATH' >> ~/cpac_env.sh
    		echo 'export DYLD_FALLBACK_LIBRARY_PATH=/opt/afni' >> ~/cpac_env.sh

	elif [ $LOCAL -eq 1 ]; then
    		mv ${AFNI_DOWNLOAD} ~/afni 
		export PATH=~/afni:$PATH
		export DYLD_FALLBACK_LIBRARY_PATH=~/afni
    		echo '# Path to AFNI' >> ~/cpac_env.sh
    		echo 'export PATH=~/afni:$PATH' >> ~/cpac_env.sh
    		echo 'export DYLD_FALLBACK_LIBRARY_PATH=~/afni' >> ~/cpac_env.sh
	else
		echo Invalid value for variable 'LOCAL'.
		echo This script is unable to determine whether or not you are running it as root.
		echo '[ '$(date)' ] : AFNI could not be installed (unable to determine if root).' >> ~/cpac.log
		cd $INIT_DIR
		exit 1
	fi
}

function install_ants {
	echo "Installing ANTS."
	which ANTS &> /dev/null ; if [ $? -eq 0 ]; then
		echo ANTS is already installed!
		echo Moving on...
		echo '[ '$(date)' ] : ANTS is already installed - does not need to be re-installed.' >> ~/cpac.log
		return
	fi
	system_dependencies_installed ; if [ $? -ne 0 ]; then
		echo ANTS cannot be installed unless system-level dependencies are installed first.
		echo Have your system administrator install system-level dependencies as root.
		echo Exiting now...
		echo '[ '$(date)' ] : ANTS installation failed - system-level dependencies are not installed.' >> ~/cpac.log
		cd $INIT_DIR
		exit 1
	fi
	which c3d &> /dev/null ; if [ $? -ne 0 ]; then
		echo ANTS cannot be installed unless c3d is installed first.
		echo Install c3d and then try again.
		echo Exiting now...
		echo '[ '$(date)' ] : ANTS installation failed - C3D is not installed.' >> ~/cpac.log
		cd $INIT_DIR
		install_cpac_env
		exit 1
	fi
    	cd /tmp
    	git clone https://github.com/stnava/ANTs.git
	if [ $LOCAL -eq 0 ]; then
		mkdir /opt/ants
		cd /opt/ants
		cmake -c -g /tmp/ANTs
		make -j 4
		ANTSPATH=/opt/ants/bin
		cp /tmp/ANTs/Scripts/antsIntroduction.sh ${ANTSPATH}
		cp /tmp/ANTs/Scripts/antsAtroposN4.sh ${ANTSPATH}
		cp /tmp/ANTs/Scripts/antsBrainExtraction.sh ${ANTSPATH}
		cp /tmp/ANTs/Scripts/antsCorticalThickness.sh ${ANTSPATH}
		export ANTSPATH
		export PATH=/opt/ants/bin:$PATH
		echo '# Path to ANTS' >> ~/cpac_env.sh
		echo 'export ANTSPATH=/opt/ants/bin/' >> ~/cpac_env.sh
		echo 'export PATH=/opt/ants/bin:$PATH' >> ~/cpac_env.sh
	elif [ $LOCAL -eq 1 ]; then
		mkdir ~/ants
		cd ~/ants
		cmake -c -g /tmp/ANTs
		make -j 4
		ANTSPATH=~/ants/bin
		cp /tmp/ANTs/Scripts/antsIntroduction.sh ${ANTSPATH}
		cp /tmp/ANTs/Scripts/antsAtroposN4.sh ${ANTSPATH}
		cp /tmp/ANTs/Scripts/antsBrainExtraction.sh ${ANTSPATH}
		cp /tmp/ANTs/Scripts/antsCorticalThickness.sh ${ANTSPATH}
		export ANTSPATH
                export PATH=/opt/ants/bin:$PATH
		echo '# Path to ANTS' >> ~/cpac_env.sh
		echo 'export ANTSPATH=~/ants/bin/' >> ~/cpac_env.sh
		echo 'export PATH=~/ants/bin:$PATH' >> ~/cpac_env.sh
	else
		echo Invalid value for variable 'LOCAL'.
		echo This script is unable to determine whether or not you are running it as root.
		echo '[ '$(date)' ] : ANTS could not be installed (unable to determine if root).' >> ~/cpac.log
		cd $INIT_DIR
		exit 1
	fi
}

function install_c3d {
	echo "Installing C3D."
	which c3d &> /dev/null ; if [ $? -eq 0 ]; then
		echo c3d is already installed!
		echo Moving on...
		echo '[ '$(date)' ] : C3D is already installed - does not need to be re-installed.' >> ~/cpac.log
		return
	fi
	ARCHITECTURE=$(uname -p)
	case $ARCHITECTURE in
    		x86_64 )
        		C3D_DOWNLOAD=c3d-0.8.2-Linux-x86_64
        		;;
    		i386 )
        		C3D_DOWNLOAD=c3d-0.8.2-Linux-i386
        		;;
   		i686 )
        		C3D_DOWNLOAD=c3d-0.8.2-Linux-i686
     			;;
	esac
	cd /tmp
	wget http://sourceforge.net/projects/c3d/files/c3d/c3d-0.8.2/${C3D_DOWNLOAD}.tar.gz
	tar xfz ${C3D_DOWNLOAD}.tar.gz
	if [ $LOCAL -eq 0 ]; then
		mv $C3D_DOWNLOAD /opt/c3d
		export PATH=/opt/c3d/bin:$PATH
		echo '# Path to C3D' >> ~/cpac_env.sh
		echo 'export PATH=/opt/c3d/bin:$PATH' >> ~/cpac_env.sh
	elif [ $LOCAL -eq 1 ]; then
		mv $C3D_DOWNLOAD ~/c3d
		export PATH=~/c3d/bin:$PATH
		echo '# Path to C3D' >> ~/cpac_env.sh
		echo 'export PATH=~/c3d/bin:$PATH' >> ~/cpac_env.sh
	else
		echo Invalid value for variable 'LOCAL'.
		echo This script is unable to determine whether or not you are running it as root.
		echo '[ '$(date)' ] : C3D could not be installed (unable to determine if root).' >> ~/cpac.log
		cd $INIT_DIR
		exit 1
	fi
}

function install_cpac_resources {
	echo "Installing C-PAC Image Resources."
	# Determines if C-PAC image resources are all already installed.
	RES_PRES=1
	for res in MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz MNI152_T1_2mm_brain_symmetric.nii.gz MNI152_T1_2mm_symmetric.nii.gz MNI152_T1_3mm_brain_mask_dil.nii.gz MNI152_T1_3mm_brain_mask.nii.gz MNI152_T1_3mm_brain_mask_symmetric_dil.nii.gz MNI152_T1_3mm_brain.nii.gz MNI152_T1_3mm_brain_symmetric.nii.gz MNI152_T1_3mm.nii.gz MNI152_T1_3mm_symmetric.nii.gz; do
		[ ! -f $FSLDIR/data/standard/$res ] && RES_PRES=0
	done
	[ ! -d $FSLDIR/data/standard/tissuepriors/2mm ] || [ ! -d $FSLDIR/data/standard/tissuepriors/3mm ] || [ ! -f $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz ] && RES_PRES=0 
	if [ $RES_PRES -eq 1 ]; then
		echo CPAC Resources are already present!
		echo Moving on...
		echo '[ '$(date)' ] : C-PAC resources are already installed - do not need to be re-installed.' >> ~/cpac.log
		return
	fi
	which fsl &> /dev/null ; if [ $? -ne 0 ]; then
		echo CPAC templates cannot be copied unless FSL is installed first.
		echo Install FSL and then try again.
		echo Exiting now...
		echo '[ '$(date)' ] : C-PAC resources installation failed - FSL is not installed.' >> ~/cpac.log
		cd $INIT_DIR
		install_cpac_env
		exit 1
	fi
	cd /tmp
	wget http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tgz
	tar xfz cpac_resources.tgz 2> /dev/null
	cd cpac_image_resources
	cp -n MNI_3mm/* $FSLDIR/data/standard
	cp -n symmetric/* $FSLDIR/data/standard
	cp -nr tissuepriors/2mm $FSLDIR/data/standard/tissuepriors
	cp -nr tissuepriors/3mm $FSLDIR/data/standard/tissuepriors
	cp -n HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz $FSLDIR/data/atlases/HarvardOxford
}

function install_cpac {
	echo "Installing C-PAC."
	python -c "import CPAC" 2> /dev/null ; if [ $? -eq 0 ]; then
		echo CPAC is already installed!
		echo Moving on...
		echo '[ '$(date)' ] : C-PAC is already installed - does not need to be re-installed.' >> ~/cpac.log
		return
	fi
	which fsl &> /dev/null ; if [ $? -ne 0 ]; then
		echo CPAC cannot be installed unless FSL is installed first.
		echo Install FSL and then try again.
		echo Exiting now...
		echo '[ '$(date)' ] : C-PAC installation failed - FSL is not installed.' >> ~/cpac.log
		cd $INIT_DIR
		install_cpac_env
		exit 1
	fi
	which afni &> /dev/null ; if [ $? -ne 0 ]; then
		echo CPAC cannot be installed unless AFNI is installed first.
		echo Install AFNI and then try again.
		echo Exiting now...
		echo '[ '$(date)' ] : C-PAC installation failed - AFNI is not installed.' >> ~/cpac.log
		cd $INIT_DIR
		install_cpac_env
		exit 1
	fi
	python_dependencies_installed ; if [ $? -ne 0 ]; then
		echo CPAC cannot be installed unless Python dependencies are installed first.
		echo Install Python dependencies and then try again.
		echo Exiting now...
		echo '[ '$(date)' ] : C-PAC installation failed - Python dependencies are not installed.' >> ~/cpac.log
		cd $INIT_DIR
		install_cpac_env
		exit 1
	fi
	source activate cpac
	cd /tmp
	git clone https://github.com/FCP-INDI/C-PAC.git
	cd C-PAC
	python setup.py install
	source deactivate
}

function install_cpac_env {
	echo "Installing C-PAC environmental variables"
	if [ -f ~/cpac_env.sh ]; then
		# Append cpac_env.sh to end of bashrc and remove if this is not root.  Otherwise move cpac_env.sh to /etc/profile.d
		if [ $LOCAL -eq 1 ]; then
			cat ~/cpac_env.sh >> ~/.bashrc
			rm ~/cpac_env.sh
		elif [ $LOCAL -eq 0 ]; then
			if [ -f /etc/profile.d/cpac_env.sh ]; then
				# Since functions will not re-install already installed software, this should only append
				# packages that weren't already in cpac_env.sh.
				cat ~/cpac_env.sh >> /etc/profile.d/cpac_env.sh
				rm ~/cpac_env.sh
			else
				mv ~/cpac_env.sh /etc/profile.d/
			fi
		fi
	fi
}

# Check to see if user has root privileges.  If not, perform local install.
[ $EUID -eq 0 ] && LOCAL=0 || LOCAL=1

# Check to see whether the distribution is CentOS or Ubuntu.
[ -f /etc/redhat-release ] && DISTRO=CENTOS
which lsb_release &> /dev/null && [ $(lsb_release -si) == 'Ubuntu' ] && DISTRO=UBUNTU

INIT_DIR=$(pwd)
: ${LOCAL:? "LOCAL needs to be set and non-empty."}
: ${DISTRO:? "DISTRO needs to be set and non-empty."}
while getopts ":spn:alrh" opt; do
	case $opt in
		s)
			install_system_dependencies
			install_cpac_env
      			;;
    		p) 
      			install_python_dependencies
			install_cpac_env
      			;;
    		n) 
      			suites=($OPTARG)
      			for suite in ${suites[@]}; do
				case $suite in 
					afni)
						install_afni
						install_cpac_env
						;;
					fsl)
						install_fsl
						install_cpac_env
						;;
					c3d)
						install_c3d
						install_cpac_env
						;;
					ants)
						install_ants
						install_cpac_env
						;;
					cpac)
						install_cpac_resources
						install_cpac
						install_cpac_env
						;;
					*)
						echo Invalid neuroimaging suite: $suite
						echo CPAC provisioning script will continue.
						echo '[ '$(date)' ] : Unexpected neuroimaging suite: ' $suite  >> ~/cpac.log
						;;
				esac
      			done
      			;;
    		a) 
			install_afni
			if [ $LOCAL -eq 1 ] && [ $DISTRO == 'UBUNTU' ]; then
				echo FSL cannot be installed locally on Ubuntu.
				echo Contact your system administrator to install FSL.
				echo Continuing the installation...
				echo '[ '$(date)' ] : FSL installation failed - need root privileges on Ubuntu.' >> ~/cpac.log
			else
				install_fsl
			fi
			install_c3d
			install_ants
			install_cpac_env
			;;
		l) 
			install_python_dependencies
			install_afni
			if [ $LOCAL -eq 1 ] && [ $DISTRO == 'UBUNTU' ]; then
				echo FSL cannot be installed locally on Ubuntu.
				echo Contact your system administrator to install FSL.
				echo Continuing the installation...
				echo '[ '$(date)' ] : FSL installation failed - need root privileges on Ubuntu.' >> ~/cpac.log
			else
				install_fsl
			fi
			install_c3d
			install_ants
			install_cpac_resources
			install_cpac
			install_cpac_env
      			;;
   		 r)
			install_system_dependencies 
     			install_python_dependencies
			install_afni
			install_fsl
			install_c3d
			install_ants
			install_cpac_resources
			install_cpac
			install_cpac_env
			;;
           	 h)
			echo "
cpac_provision.sh 
 =================================================================================================
 Version: 0.3.9
 Author(s): John Pellman, Daniel Clark
 Based off of cpac_install.sh by Daniel Clark.
 Description: Will perform specific operations to install C-PAC dependencies and C-PAC.
 Checks for user privileges and performs installation either locally or system-wide.
 Can be customized using flags.
 =================================================================================================
 Flags:
 -s : System-level dependencies only.
 -p : Python dependencies only
 -n : Install specific neuroimaging packages.  Accepts any number of the following as arguments:
	afni, fsl, c3d, ants, cpac
	will issue warnings if dependencies for these neuroimaging packages are not fulfilled.
	If multiple packages are to be specified, they must be surrounded by quotation marks.
 -a : Install all neuroimaging suites not already installed.  Will also tell you if all neuroimaging suites are already installed and on the path.
 -l : Local install. Equivalent to -pa ; will not run FSL installer, but will issue a warning if running on Ubuntu. 
 -r : Root install.  Equivalent to -spa
 -h : Bring up the help dialog.
=================================================================================================
 Example usage:
	cpac_provision.sh -n \"fsl afni\"
	Will install FSL and AFNI.  The list of neuroimaging suites to install is iterated through sequentially.
	In this case, FSL would first be installed before AFNI.
					"
			;;
   		\?)
     			echo "Invalid option: -$OPTARG" >&2
     			cd $INIT_DIR
			exit 1
    			;;
		:)
      			echo "Option -$OPTARG requires an argument." >&2
			cd $INIT_DIR
     			exit 1
      			;;
	esac
done
cd $INIT_DIR

