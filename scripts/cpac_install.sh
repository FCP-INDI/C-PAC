#! /bin/bash

# these are packages that are common to centos 5, 6, and 7
centos_packages=("git" "make" "unzip" "netpbm" "gcc" "python-devel"\
    "gcc-gfortran" "gcc-c++" "libgfortran" "lapack" "lapack-devel" "blas"\
    "libcanberra-gtk2" "libXp.x86_64" "mesa-libGLU-9.0.0-4.el7.x86_64"\
    "gsl-1.15-13.el7.x86_64" "wxBase" "wxGTK" "wxGTK-gl" "wxPython" "graphviz"\
    "graphviz-devel.x86_64")

# configuration options that are specific to centos 5
centos5_epel_url="http://dl.fedoraproject.org/pub/epel/5/x86_64/epel-release-5-4.noarch.rpm"
centos5_epel_rpm="epel-release-5-4.noarch.rpm"
centos5_packages=("mesa-libGLU-6.5.1-7.11.el5_9.i386" "gsl-1.13-3.el5.x86_64"\
    "libxml2-devel libpng-1.2.10-17.el5_8.i386")

# configuration options that are specific to centos 6
centos6_epel_url="http://dl.fedoraproject.org/pub/epel/6/x86_64/e/epel-release-6-8.noarch.rpm"
centos6_epel_rpm="epel-release-6-8.noarch.rpm"
centos6_packages=("mesa-libGLU-11.0.7-4.el6.x86_64" "gsl-1.13-1.el6.x86_64"\
    "libcanberra-gtk2" "libxml2-devel" "libpng-1.2.49-2.el6_7.i686")

# configuration options that are specific to centos 7
centos7_epel_url="http://dl.fedoraproject.org/pub/epel/7/x86_64/e/epel-release-7-5.noarch.rpm"
centos7_epel_rpm="epel-release-7-5.noarch.rpm"
centos7_packages=("mesa-libGLU-9.0.0-4.el7.x86_64" "gsl-1.15-13.el7.x86_64"\
    "libcanberra-gtk2" "libxml-devel" "libpng12.x86_64")

# are all of the ubuntu packages the same regardless of the version?
ubuntu_packages=("cmake" "git" "graphviz" "graphviz-dev" "gsl-bin" "libcanberra-gtk-module" \
    "libexpat1-dev" "libgiftiio-dev" "libglib2.0-dev" "libglu1-mesa" "libglu1-mesa-dev" \
    "libgsl0-dev" "libjpeg-progs" "libmotif-dev" "libxml2" "libxml2-dev" "libxext-dev" \
    "libxft2" "libxft-dev" "libxi-dev" "libxmu-headers" "libxmu-dev" "libxpm-dev" "libxslt1-dev" \
    "libxp6" "libxp-dev" "make" "mesa-common-dev" "mesa-utils" "netpbm" "pkg-config" \
    "build-essential" "xvfb" "xauth" "libgl1-mesa-dri" "tcsh" "unzip" "zlib1g-dev" "m4")

conda_packages=("pandas" "cython" "numpy" "scipy" "matplotlib" "networkx" "traits" "pyyaml" "jinja2" "nose" "ipython" "pip" "wxpython")

pip_packages=("future" "prov" "simplejson" "lockfile" "pygraphviz" "nibabel" "nipype" "patsy" "memory_profiler" "psutil" "configparser" "indi_tools")

function set_system_deps {
    system_pkgs=''
    epel_rpm=''
    epel_url=''

    if [ $DISTRO == 'CENTOS' ]
    then
        # add in the packages that are common to all
        system_pkgs=${centos_packages[@]}

        yum update -y && yum install -y wget

        # add in the packages that are specific to the redhat-release
        version=$(rpm -q --queryformat '%{VERSION}' centos-release)
        case ${version} in
            5)
                epel_url=centos5_epel_url
                epel_rpm=centos5_epel_rpm
                system_pkgs+=(${centos5_packages})
                ;;
            6)
                epel_url=centos6_epel_url
                epel_rpm=centos6_epel_rpm
                system_pkgs+=(${centos6_packages})
                ;;
            7)
                epel_url=centos7_epel_url
                epel_rpm=centos7_epel_rpm
                system_pkgs+=(${centos7_packages})
                ;;
            *)
                echo "Unknown version ${version}"
        esac
    elif [ $DISTRO == 'UBUNTU' ]
    then
        # take care of initing apt-get and installing wget
        echo "!!!!!! CC"
        apt-get update && apt-get upgrade -y && apt-get install -y wget
        system_pkgs=${ubuntu_packages[@]}
    else
        echo "Unknown distribution ${DISTRO}"
        exit 1
    fi
}

# CC - reformatted this to have better control of the output
function print_usage {
    echo ""
    echo "Usage: cpac_install.sh -[spnalrh]"
    echo "========================================================================="
    echo "Version: 0.4.0"
    echo "Author(s): John Pellman, Daniel Clark"
    echo "Based off of cpac_install.sh by Daniel Clark."
    echo "Description: Will perform specific operations to install C-PAC"
    echo "  dependencies and C-PAC. Checks for user privileges and performs"
    echo "  installation either locally or system-wide."
    echo "========================================================================="
    echo "One or more command line options are required:"
    echo "  -s : System-level dependencies only."
    echo "  -p : Python dependencies only"
    echo "  -n : Install specific neuroimaging packages.  Accepts any number of the"
    echo "       following as arguments: afni, fsl, c3d, ants, cpac"
    echo "       will issue warnings if dependencies for these neuroimaging packages"
    echo "       are not fulfilled. If multiple packages are to be specified, they"
    echo "       must be surrounded by quotation marks."
    echo "  -a : Install all neuroimaging suites not already installed.  Will also"
    echo "       tell you if all neuroimaging suites are already installed and on"
    echo "       the path."
    echo "  -l : Local install. Equivalent to -pa ; will not run FSL installer, but"
    echo "       will issue a warning if running on Ubuntu."
    echo "  -r : Root install.  Equivalent to -spa"
    echo "  -h : Print this help message."
    echo "========================================================================="
    echo "Example usage: cpac_install.sh -n \"fsl afni\""
    echo "  Will install FSL and AFNI. The list of neuroimaging suites to install"
    echo "  is iterated through sequentially. In this case, FSL would first be"
    echo "  installed before AFNI."
    echo ""
}

function install_system_dependencies {

    echo "Installing C-PAC system dependencies... [${missing_system_dependencies[@]}][${#missing_system_dependencies[@]}]"

    if [ ${#missing_system_dependencies[@]} -eq 0 ]
    then
        echo "sys packages to be installed ${system_packages[@]}"
        echo "System dependencies are already installed!"
        echo "Moving on..."
        echo "[ $(date) ] : C-PAC system dependencies are already installed, do" \
            "not need to be re-installed." >> ~/cpac.log
        return
    fi
    if [ $LOCAL -eq 0 ]
    then
        system_dependencies_installed=1
        if [ $DISTRO == 'CENTOS' ]
        then
            version=$(rpm -q --queryformat '%{VERSION}' centos-release)

            # update the repositories
            #yum update -y
            cd /tmp && wget ${epel_url} && rpm -Uvh ${epel_rpm}

            yum install -y ${missing_system_dependencies[@]} 
            if [ $? -ne 0 ]
            then
                system_dependencies_installed=0
                echo "[ $(date) ] yum failed to install packages: ${missing_system_dependencies[@]}"
                echo "[ $(date) ] yum failed to install packages: ${missing_system_dependencies[@]}" \
                    >> ~/cpac.log
            else	
                echo "[ $(date) ] : yum Installed C-PAC system dependency"\
                    "${missing_system_dependencies[@]}"
                echo "[ $(date) ] : yum Installed C-PAC system dependency"
                    "${missing_system_dependencies[@]}" \
                    >> ~/cpac.log
            fi
            #for p in ${missing_system_dependencies[@]}
            #do
                #echo "[ $(date) ] : Installing C-PAC system dependency $p"
                #yum install -y ${p}
                #if [ $? -ne 0 ]
                #then
                    #system_dependencies_installed=0
                    #echo "failed to install package: ${p}"
                    #echo "[ $(date) ] : Failed to install C-PAC system dependency $p" \
                        #>> ~/cpac.log
		#else
                    #echo "[ $(date) ] : Installed C-PAC system dependency $p"
                    #echo "[ $(date) ] : Installed C-PAC system dependency $p" \
                        #>> ~/cpac.log
                #fi
            #done
        elif [ $DISTRO == 'UBUNTU' ]
        then
            #apt-get update
            #apt-get upgrade -y

            apt-get install -y ${missing_system_dependencies[@]} 
            if [ $? -ne 0 ]
            then
                system_dependencies_installed=0
                echo "[ $(date) ] apt-get failed to install packages: ${missing_system_dependencies[@]}"
                echo "[ $(date) ] apt-get failed to install packages: ${missing_system_dependencies[@]}" \
                    >> ~/cpac.log
            else	
                echo "[ $(date) ] : apt-get Installed C-PAC system dependency"\
                    "${missing_system_dependencies[@]}"
                echo "[ $(date) ] : apt-get Installed C-PAC system dependency"
                    "${missing_system_dependencies[@]}" \
                    >> ~/cpac.log
            fi
           # for p in ${missing_system_dependencies[@]}
           # do
           #     echo "[ $(date) ] : Installing C-PAC system dependency $p"
           #     apt-get install -y ${p}
           #     if [ $? -ne 0 ]
           #     then
           #         system_dependencies_installed=0
           #         echo "failed to install package: ${p}"
           #         echo "[ $(date) ] : Failed to install C-PAC system dependency $p" \
           #             >> ~/cpac.log
	   #     else
           #         echo "[ $(date) ] : Installed C-PAC system dependency $p" 
           #         echo "[ $(date) ] : Installed C-PAC system dependency $p" \
           #             >> ~/cpac.log
           #     fi
           # done
            # finish up
            apt-get autoremove -y
        else
            echo "Linux distribution not recognized.  System-level dependencies cannot" \
                 "be installed."
            echo "[ $(date) ] : C-PAC system dependencies could not be installed (Linux" \
                "distribution not recognized)." >> ~/cpac.log
            cd $INIT_DIR
            exit 1
        fi
    elif [ $LOCAL -eq 1 ]
    then
        echo "System-level dependencies cannot be installed since you do not have"\
             "root privileges."
        echo "Re-run this script as root or have your system administrator run it."
        cd $INIT_DIR
        echo "[ $(date) ] : C-PAC system dependencies could not be installed (not root)."\
            >> ~/cpac.log
        exit 1
    else
        echo "Invalid value for variable 'LOCAL'."
        echo "This script is unable to determine whether or not you are running it as root."
        echo "[ $(date) ] : C-PAC system dependencies could not be installed (unable to"\
            "determine if root)." >> ~/cpac.log
        cd $INIT_DIR
        exit 1
    fi
    if [ ${system_dependencies_installed} -eq 0 ]
    then
        echo "[ $(date) ] : C-PAC system dependencies not fully installed." >> ~/cpac.log
    else
        echo "[ $(date) ] : C-PAC system dependencies succesfully installed." >> ~/cpac.log
    fi
}

function get_missing_system_dependencies()
{
    missing_system_dependencies=()
    system_dependencies_installed=1

    if [ $DISTRO == 'CENTOS' ]
    then
        for package in ${system_pkgs[@]}
        do
            yum list installed ${package} > /dev/null 2>&1
            if [ $? -ne 0 ]
            then
                system_dependencies_installed=0
                ARRAY+=(${package})
                echo "[ $(date) ] : Missing system dependency ${package}" >> ~/cpac.log
            fi
        done
    elif [ $DISTRO == 'UBUNTU' ]
    then
        for package in ${system_pkgs[@]}
        do
            dpkg -s ${package} > /dev/null 2>&1
            if [ $? -ne 0 ]
            then
                system_dependencies_installed=0
                missing_system_dependencies+=(${package})
                echo "[ $(date) ] : Missing system dependency ${package}" >> ~/cpac.log
            fi
        done
    else
        echo "[ $(date) ] : Do not know how to check for packages installed on ${DISTRO}" >> ~/cpac.log
    fi
    echo "missing ${missing_system_dependencies[@]}"
}

function install_python_dependencies {

    if [ ${python_dependencies_installed} -eq 1 ]
    then
        echo "[ $(date) ] C-PAC Python dependencies installed!"
        echo "[ $(date) ] C-PAC Python dependencies installed!" >> ~/cpac.log
        return
    fi

    if [ ${system_dependencies_installed} -ne 1 ]
    then
        echo "Python dependencies cannot be installed unless system-level dependencies are installed first."
        echo "Have your system administrator install system-level dependencies as root."
        echo "Exiting now..."
        echo "[ $(date) ] : Python dependencies could not be installed (system-level" \
            "dependencies not installed." >> ~/cpac.log
        cd $INIT_DIR
        exit 1
    fi

    # for now always install miniconda, in the future should only install 
    # if not there
    echo "[ $(date) ] Installing miniconda!"
    echo "[ $(date) ] Installing miniconda!" >> ~/cpac.log

    cd /tmp
    if [ ! -f Miniconda-3.8.3-Linux-x86_64.sh ]
    then
        wget http://repo.continuum.io/miniconda/Miniconda-3.8.3-Linux-x86_64.sh
        if [ $? -ne 0 ]
        then
            echo "[ $(date) ] Could not download miniconda installation script!"
            echo "[ $(date) ] Could not download miniconda installation script!" >> ~/cpac.log
            return 
        fi
    fi
    chmod +x Miniconda-3.8.3-Linux-x86_64.sh
    if [ $LOCAL -eq 0 ]
    then
        ./Miniconda-3.8.3-Linux-x86_64.sh -b -p /usr/local/bin/miniconda
        if [ $? -ne 0 ]
        then
            echo "[ $(date) ] Miniconda installation failed!"
            echo "[ $(date) ] Miniconda installation failed!" >> ~/cpac.log
            #return 
        fi
        chmod -R 775 /usr/local/bin/miniconda
        chmod g+s /usr/local/bin/miniconda
        export PATH=/usr/local/bin/miniconda/bin:${PATH}
        echo 'export PATH=/usr/local/bin/miniconda/bin:${PATH}' >> ~/cpac_env.sh
    elif [ $LOCAL -eq 1 ] && [ ! -d ~/miniconda ]
    then
        ./Miniconda-3.8.3-Linux-x86_64.sh -b
        if [ $? -ne 0 ]
        then
            echo "[ $(date) ] Miniconda installation failed!"
            echo "[ $(date) ] Miniconda installation failed!" >> ~/cpac.log
            return 
        fi
        export PATH=~/miniconda/bin:${PATH}
        echo 'export PATH=~/miniconda/bin:${PATH}' >> ~/cpac_env.sh
    fi

    # for docker dont install virtualenv
    #conda create -y -n cpac python
    #source activate cpac
    conda install -y ${missing_conda_dependencies[@]}
    if [ $? -ne 0 ]
    then
        echo "[ $(date) ] Conda install ${p} failed!"
        echo "[ $(date) ] Conda install ${p} failed!" >> ~/cpac.log
        exit 1 
    fi
    #for p in ${missing_conda_dependencies[@]}
    #do
        #echo "[ $(date) ] Conda install ${p}!"
        #conda install -y ${p}
        #if [ $? -ne 0 ]
        #then
            #echo "[ $(date) ] Conda install ${p} failed!"
            #echo "[ $(date) ] Conda install ${p} failed!" >> ~/cpac.log
            #exit 1 
        #fi
    #done

    pip install ${missing_pip_dependencies[@]}
    if [ $? -ne 0 ]
    then
        echo "[ $(date) ] Pip install ${missing_pip_dependencies[@]} failed!"
        echo "[ $(date) ] Pip install ${p} failed!" >> ~/cpac.log
        exit 1 
    fi
    #for p in ${missing_pip_dependencies[@]}
    #do
        #echo "[ $(date) ] Pip install ${p}!"
        #pip install ${p}
        #if [ $? -ne 0 ]
        #then
            #echo "[ $(date) ] Pip install ${p} failed!"
            #echo "[ $(date) ] Pip install ${p} failed!" >> ~/cpac.log
            #exit 1 
        #fi
    #done

    #echo 'source activate cpac' >> ~/cpac_env.sh
    cd /tmp
    git clone https://github.com/FCP-INDI/INDI-Tools.git
    cd INDI-Tools/
    python setup.py install
    #source deactivate
    cd $INIT_DIR
}

function get_missing_python_dependencies {

    python_dependencies_installed=0
    missing_pip_dependencies=()
    missing_conda_dependencies=()

    # first we check to make sure that we have python
    if [ ! -f /usr/local/bin/miniconda/bin/python ]
    then
        python_installed=0
    else
        python_installed=1
    fi

    if [ ${python_installed} -eq 0 ]
    then
        echo "[ $(date) ] : Python is not installed, need to install all"\
             "Python dependencies." >> ~/cpac.log
        missing_pip_dependencies=${pip_packages[@]}
        missing_conda_dependencies=${conda_packages[@]}
    else
        # if we find an enviroment, then enable it
        if [ -d ~/miniconda/envs/cpac ] || [ -d /usr/local/bin/miniconda/envs/cpac ]
        then
            echo "[ $(date) ] : Found C-PAC virtual environment, activating" >> ~/cpac.log
            source activate cpac &> /dev/null
        fi

        python_dependencies_installed=1
        for p in ${pip_packages[@]}
        do
            if [ ${p} == "indi_tools" ]
            then
                /usr/local/bin/miniconda/bin/python -c "import indi_aws" 2> /dev/null
                if [ $? -ne 0 ]
                then
                    echo "[ $(date) ] : Python package $p not installed" >> ~/cpac.log
                    missing_pip_dependencies+=($p)
                    python_dependencies_installed=0
                else
                    echo "[ $(date) ] : Python package $p installed" >> ~/cpac.log
                fi
            else
                /usr/local/bin/miniconda/bin/python -c "import ${p}" 2> /dev/null
                if [ $? -ne 0 ]
                then
                    echo "[ $(date) ] : Python package $p not installed" >> ~/cpac.log
                    missing_pip_dependencies+=($p)
                    python_dependencies_installed=0
                else
                    echo "[ $(date) ] : Python package $p installed" >> ~/cpac.log
                fi
            fi
        done

        for p in ${conda_packages[@]}
        do
            if [ ${p} == "wxpython" ]
            then
                /usr/local/bin/miniconda/bin/python -c "import wx" 2> /dev/null
                retval=$?
            elif [ ${p} == "pyyaml" ]
            then
                /usr/local/bin/miniconda/bin/python -c "import yaml" 2> /dev/null
                retval=$?
            elif [ ${p} == "ipython" ]
            then
                if [ -f /usr/local/bin/miniconda/bin/ipython ]
                then
                    retval=0
                else
                    retval=1
                fi
            else
                /usr/local/bin/miniconda/bin/python -c "import ${p}" 2> /dev/null
                retval=$?
            fi
            if [ $retval -ne 0 ]
            then
                echo "[ $(date) ] : Python package $p not installed" >> ~/cpac.log
                missing_conda_dependencies+=($p)
                python_dependencies_installed=0
            else
                echo "[ $(date) ] : Python package $p installed" >> ~/cpac.log
            fi
        done

        # if we find an enviroment, then disable it
        if [ -d ~/miniconda/envs/cpac ] || [ -d /usr/local/bin/miniconda/envs/cpac ]
        then
            echo "[ $(date) ] : Found C-PAC virtual environment, de-activating" >> ~/cpac.log
            source deactivate &> /dev/null
        fi
    fi
}

function install_fsl {
    echo "Installing FSL."
    which fsl &> /dev/null ; if [ $? -eq 0 ]; then
        echo FSL is already installed!
        echo Moving on...
        echo '[ '$(date)' ] : FSL is already installed - does not need to be re-installed.' >> ~/cpac.log
        return
    fi
    if [ $system_dependencies_installed -ne 1 ]
    then
        echo "FSL cannot be installed unless system-level dependencies are installed first."
        echo "Have your system administrator install system-level dependencies as root."
        echo "Exiting now..."
        echo "[ $(date) ] : FSL installation failed - system-level dependencies are not installed." >> ~/cpac.log
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
            if [ $? -ne 0 ]
            then
                echo "FSL Install failed!"
                exit 1
            fi

            FSLDIR=/usr/share/fsl/
            mkdir $FSLDIR/5.0
            mv $FSLDIR/bin $FSLDIR/5.0/bin
            ln -s $FSLDIR/data $FSLDIR/5.0/data
            mv $FSLDIR/doc $FSLDIR/5.0/doc
            mv $FSLDIR/etc $FSLDIR/5.0/etc
            mv $FSLDIR/tcl $FSLDIR/5.0/tcl
        # Debian-based distros must use NeuroDebian instead of the installer.
        elif [ $DISTRO == 'UBUNTU' ]; then
            wget -O- http://neuro.debian.net/lists/trusty.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
            apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9
            apt-get update
            apt-get install -y fsl-5.0-complete
            if [ $? -ne 0 ]
            then
                echo "FSL Install failed!"
                exit 1
            fi

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
    elif [ $LOCAL -eq 1 ]
    then
        if [ $DISTRO == 'CENTOS' ]
        then
            python fslinstaller.py -d ~
            if [ $? -ne 0 ]
            then
                echo "FSL Install failed!"
                exit 1
            fi

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
        elif [ $DISTRO == 'UBUNTU' ]
        then
            echo FSL cannot be installed without root privileges on Ubuntu Linux.
            echo "[ $(date) ] : FSL installation failed - need root privileges" \
                "on Ubuntu." >> ~/cpac.log 
            cd $INIT_DIR
            install_cpac_env
            exit 1
        fi
    else
        echo "Invalid value for variable 'LOCAL'."
        echo "This script is unable to determine whether or not you are running it as root."
        echo "[ $(date) ] : FSL could not be installed (unable to determine " \
            "if root)." >> ~/cpac.log
        cd $INIT_DIR
        exit 1
    fi
}

function install_afni {
    echo "Installing AFNI."
    which afni &> /dev/null ; if [ $? -eq 0 ]; then
        echo AFNI is already installed!
        echo Moving on...
        echo "[ $(date) ] : AFNI is already installed - does not need to be" \
            " re-installed." >> ~/cpac.log
        return
    fi
    if [ $system_dependencies_installed -ne 1 ]
    then
        echo "AFNI cannot be installed unless system-level dependencies are installed first."
        echo "Have your system administrator install system-level dependencies as root."
        echo "Exiting now..."
        echo "[ $(date) ] : AFNI installation failed - system-level dependencies are" \
            "not installed." >> ~/cpac.log
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

    #wget http://fcp-indi.s3.amazonaws.com/resources/cc_afni_trusty_openmp_64.tar.gz
    #tar xfz cc_afni_trusty_openmp_64.tar.gz
    #rm cc_afni_trusty_openmp_64.tar.gz
    #AFNI_DOWNLOAD=afni

    #git clone https://github.com/ccraddock/afni.git
    #cd /tmp/afni/src
    #cp Makefile.linux_openmp_64_trusty Makefile
    #make vastness

    if [ $? -ne 0 ]
    then
        echo "AFNI Install failed!"
        exit 1
    fi

    if [ $LOCAL -eq 0 ]
    then
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
    if [ ${system_dependencies_installed} -ne 1 ]
    then
        echo ANTS cannot be installed unless system-level dependencies are installed first.
        echo Have your system administrator install system-level dependencies as root.
        echo Exiting now...
        echo '[ '$(date)' ] : ANTS installation failed - system-level dependencies are not installed.' >> ~/cpac.log
        cd $INIT_DIR
        exit 1
    fi
    which c3d &> /dev/null ; if [ $? -ne 0 ]; then
        echo "ANTS cannot be installed unless c3d is installed first."
        echo "Install c3d and then try again."
        echo "Exiting now..."
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
	# go slow, -j 4 causes seg fault w/ building containers
        make
        if [ $? -ne 0 ]
        then
            echo "ANTS compile failed."
            echo "Exiting now..."
            echo "[ $(date) ] : ANTS installation failed - compile failed." >> ~/cpac.log
            cd $INIT_DIR
            exit 1
        fi
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
	# go slow, -j 4 causes seg fault w/ building containers
        make
        if [ $? -ne 0 ]
        then
            echo "ANTS compile failed."
            echo "Exiting now..."
            echo "[ $(date) ] : ANTS installation failed - compile failed." >> ~/cpac.log
            cd $INIT_DIR
            exit 1
        fi
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

cpac_resources=("$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz" \
    "$FSLDIR/data/standard/MNI152_T1_2mm_brain_symmetric.nii.gz" \
    "$FSLDIR/data/standard/MNI152_T1_2mm_symmetric.nii.gz" \
    "$FSLDIR/data/standard/MNI152_T1_3mm_brain_mask_dil.nii.gz" \
    "$FSLDIR/data/standard/MNI152_T1_3mm_brain_mask.nii.gz" \
    "$FSLDIR/data/standard/MNI152_T1_3mm_brain_mask_symmetric_dil.nii.gz" \
    "$FSLDIR/data/standard/MNI152_T1_3mm_brain.nii.gz" \
    "$FSLDIR/data/standard/MNI152_T1_3mm_brain_symmetric.nii.gz" \
    "$FSLDIR/data/standard/MNI152_T1_3mm.nii.gz" \
    "$FSLDIR/data/standard/MNI152_T1_3mm_symmetric.nii.gz" \
    "$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz")

cpac_resdirs=("$FSLDIR/data/standard/tissuepriors/2mm" \
              "$FSLDIR/data/standard/tissuepriors/3mm")

function install_cpac_resources {
    echo "Installing C-PAC Image Resources."
    # Determines if C-PAC image resources are all already installed.
    RES_PRES=1
    for res in ${cpac_resources[@]}
    do
        if [ ! -f $FSLDIR/data/standard/$res ]
        then
            RES_PRES=0
        fi
    done
    for resdir in ${cpac_resdirs[@]}
    do
        if [ ! -d ${resdir} ]
        then
            RES_PRES=0
        fi
    done

    if [ ${RES_PRES} -eq 1 ]
    then
        echo "CPAC Resources are already present!"
        echo "Moving on..."
        echo "[ $(date) ] : C-PAC resources are already installed - do not need to be re-installed." >> ~/cpac.log
        return
    fi

    #which fsl &> /dev/null ; if [ $? -ne 0 ]; then
    if [ ! -d "$FSLDIR/data" ]
    then
        echo "CPAC templates cannot be copied unless FSL is installed first."
        echo "Install FSL and then try again."
        echo "Exiting now..."
        echo "[ $(date) ] : C-PAC resources installation failed - FSL is not installed." >> ~/cpac.log
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
    if [ ${python_dependencies_installed} -ne 1 ]
    then
        echo CPAC cannot be installed unless Python dependencies are installed first.
        echo Install Python dependencies and then try again.
        echo Exiting now...
        echo "missing python dependencies"
        echo ${missing_conda_dependencies[@]}
        echo ${missing_pip_dependencies[@]}
        echo '[ '$(date)' ] : C-PAC installation failed - Python dependencies are not installed.' >> ~/cpac.log
        cd $INIT_DIR
        install_cpac_env
        exit 1
    fi
    #source activate cpac
    cd /tmp
    #wget https://github.com/FCP-INDI/C-PAC/archive/v1.0.0.tar.gz
    #tar xzvf v1.0.0.tar.gz
    git clone https://github.com/FCP-INDI/C-PAC.git C-PAC-1.0.0
    cd C-PAC-1.0.0
    python setup.py install
    rm -rf /tmp/C-PAC-1.0.0
    #source deactivate
}

function install_cpac_env {
    echo "Installing C-PAC environmental variables"
    if [ -f ~/cpac_env.sh ]
    then
        # Append cpac_env.sh to end of bashrc and remove if this is not root.
        # Otherwise move cpac_env.sh to /etc/profile.d
        if [ $LOCAL -eq 1 ]
        then
            cat ~/cpac_env.sh >> ~/.bashrc
            rm ~/cpac_env.sh
        elif [ $LOCAL -eq 0 ]
        then
            if [ -f /etc/profile.d/cpac_env.sh ]
            then
                # Since functions will not re-install already installed 
                # software, this should only append
                # packages that weren't already in cpac_env.sh.
                cat ~/cpac_env.sh >> /etc/profile.d/cpac_env.sh
                rm ~/cpac_env.sh
            else
                mv ~/cpac_env.sh /etc/profile.d/
            fi
        fi
    fi
}


##### MAIN ENTRY POINT



# Check to see if user has root privileges.  If not, perform local install.
# CC undid the obfuscation
if [ $EUID -eq 0 ]
then
    # user is superuser
    LOCAL=0
else
    # user is not superuser
    LOCAL=1
fi

# Check to see whether the distribution is CentOS or Ubuntu.
#  CC: broke this out to make it understandable, and to make it work with
#      bare-bone installations that do not have lsb_release installed
if [ -f /etc/redhat-release ]
then
    DISTRO=CENTOS
elif [ -f /etc/lsb-release ]
then
    source /etc/lsb-release
    DISTRO=${DISTRIB_ID^^}
fi

INIT_DIR=$(pwd)

# again easier to read
if [ -z ${LOCAL} ]
then
    echo "LOCAL needs to be set and non-empty."
    exit 1
fi

if [ -z ${DISTRO} ]
then
    echo "DISTRO needs to be set and non-empty. Check that /etc/redhat-release\n"
    echo "or /etc/lsb-release exist."
    exit 1
fi

# tell the user what we are doing
if [ ${LOCAL} -eq 1 ]
then
    echo "Installing the C-CPAC ecosystem locally on ${DISTRO} with $@"
else
    echo "Installing the C-CPAC ecosystem system-wide on ${DISTRO} with $@"
fi

set_system_deps

# get an accounting of the missing dependencies
get_missing_system_dependencies
get_missing_python_dependencies

echo "missing python dependencies"
echo ${missing_conda_dependencies[@]}
echo ${missing_pip_dependencies[@]}


if [ ${system_dependencies_installed} -eq 1 ]
then
    echo "All required system dependencies are installed."
elif [ ${LOCAL} -eq 1 ]
then
    echo "The following system dependences need to be installed as super"\
        "user before the C-PAC installation can continue:"
    for p in ${missing_system_dependencies}
    do
        echo "  $p"
    done
    exit 1
fi

# CC if user doesn't provide any command line arguments, install everything
if [ $# -eq 0 ]
then
    if [ ${LOCAL} -eq 0 ]
    then
        install_system_dependencies
    else
        echo "Installing system dependencies requires you to be superuser, skipping ..."
    fi
    install_python_dependencies
    install_afni
    if [ ${LOCAL} -eq 0 ] || [ ${DISTRO} == "CENTOS" ]
    then
        install_fsl
    else
        echo "Installing FSL on Ubuntu requires you to be superuser, skipping ..."
    fi
    install_c3d
    install_ants
    install_cpac_resources
    install_cpac
    install_cpac_env
fi

while getopts ":spn:alrh" opt
do
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
            for suite in ${suites[@]}
            do
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
                        echo "Invalid neuroimaging suite: $suite"
                        echo "CPAC provisioning script will continue."
                        echo "[ $(date) ] : Unexpected neuroimaging suite: $suite" \
                            >> ~/cpac.log
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
                echo "[ $(date) ] : FSL installation failed - need root privileges" \
                    "on Ubuntu." >> ~/cpac.log
            else
                install_fsl
            fi
            install_c3d
            install_ants
            install_cpac_resources
            install_cpac
            install_cpac_env
            ;;
        l)
            install_python_dependencies
            install_afni
            if [ $LOCAL -eq 1 ] && [ $DISTRO == 'UBUNTU' ]
            then
                echo "FSL cannot be installed locally on Ubuntu."
                echo "Contact your system administrator to install FSL."
                echo "Continuing the installation..."
                echo "[$(date)] : FSL installation failed - need root" \
                    "privileges on Ubuntu." >> ~/cpac.log
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
            print_usage
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
