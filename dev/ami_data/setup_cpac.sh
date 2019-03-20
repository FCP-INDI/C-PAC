#!/bin/bash

set -e

apt-get update
apt-get install -y x2goserver lubuntu-desktop lxde-icon-theme xvfb
apt-get remove -y lxlock xscreensaver xscreensaver-data gnome-screensaver

# Disable asking for user password
groupadd -r autologin
gpasswd -a ubuntu autologin
cat <<EOT > /etc/lightdm/lightdm.conf
[Seat:*]
pam-service=lightdm
pam-autologin-service=lightdm-autologin
autologin-user=ubuntu  
autologin-user-timeout=0
session-wrapper=/etc/X11/Xsession
greeter-session=lightdm-greeter
EOT

rm -f /etc/xdg/autostart/gnome-screensaver.desktop
rm -f /etc/xdg/autostart/org.gnome.SettingsDaemon.ScreensaverProxy.desktop
rm -f /etc/xdg/autostart/light-locker.desktop

Xvfb :99 & export DISPLAY=:99
su -c 'lxsession' ubuntu &  # to create configs

sleep 10

# sed -z 's/\s*Button\s*{\s*id=lxde-screenlock.desktop\s*}//g' /home/ubuntu/.config/lxpanel/LXDE/panels/panel

# Disable screen lock
cat <<EOT > /home/ubuntu/.config/xfce4/xfconf/xfce-perchannel-xml/xfce4-power-manager.xml
<?xml version="1.0" encoding="UTF-8"?>

<channel name="xfce4-power-manager" version="1.0">
  <property name="xfce4-power-manager" type="empty">
    <property name="power-button-action" type="empty"/>
    <property name="dpms-enabled" type="bool" value="false"/>
  </property>
</channel>
EOT

chown -R ubuntu: /home/ubuntu/.config

apt-get install -y \
    build-essential \
    cmake \
    graphviz \
    graphviz-dev \
    gsl-bin \
    libcanberra-gtk-module \
    libexpat1-dev \
    libgiftiio-dev \
    libglib2.0-dev \
    libglu1-mesa \
    libglu1-mesa-dev \
    libjpeg-progs \
    libgl1-mesa-dri \
    libglw1-mesa \
    libxml2 \
    libxml2-dev \
    libxext-dev \
    libxft2 \
    libxft-dev \
    libxi-dev \
    libxmu-headers \
    libxmu-dev \
    libxpm-dev \
    libxslt1-dev \
    m4 \
    make \
    mesa-common-dev \
    mesa-utils \
    netpbm \
    pkg-config \
    tcsh \
    unzip \
    xvfb \
    xauth \
    zlib1g-dev \
    dh-autoreconf \
    libgsl-dev \
    libmotif-dev \
    libtool \
    libx11-dev \
    libxext-dev \
    x11proto-xext-dev \
    x11proto-print-dev \
    xutils-dev

git clone git://anongit.freedesktop.org/xorg/lib/libXp /tmp/libXp && \
    cd /tmp/libXp && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    cd / && \
    rm -rf /tmp/libXp

mkdir -p /opt/c3d && \
    curl -sSL "http://downloads.sourceforge.net/project/c3d/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz" \
    | tar -xzC /opt/c3d --strip-components 1

export C3DPATH=/opt/c3d/
export PATH=$C3DPATH/bin:$PATH

libs_path=/usr/lib/x86_64-linux-gnu && \
    if [ -f $libs_path/libgsl.so.19 ]; then \
        ln $libs_path/libgsl.so.19 $libs_path/libgsl.so.0; \
    fi && \
    mkdir -p /opt/afni && \
    curl -sO https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz && \
    tar zxv -C /opt/afni --strip-components=1 -f linux_openmp_64.tgz && \
    rm -rf linux_openmp_64.tgz

export PATH=/opt/afni:$PATH


apt-get install -y --no-install-recommends \
    fsl-core \
    fsl-atlases \
    fsl-mni152-templates

export FSLDIR=/usr/share/fsl/5.0
export FSLOUTPUTTYPE=NIFTI_GZ
export FSLMULTIFILEQUIT=TRUE
export POSSUMDIR=/usr/share/fsl/5.0
export LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH
export FSLTCLSH=/usr/bin/tclsh
export FSLWISH=/usr/bin/wish
export PATH=/usr/lib/fsl/5.0:$PATH


cd /tmp && \
curl -sO http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tar.gz && \
tar xfz cpac_resources.tar.gz && \
cd cpac_image_resources && \
cp -n MNI_3mm/* $FSLDIR/data/standard && \
cp -n MNI_4mm/* $FSLDIR/data/standard && \
cp -n symmetric/* $FSLDIR/data/standard && \
cp -nr tissuepriors/2mm $FSLDIR/data/standard/tissuepriors && \
cp -nr tissuepriors/3mm $FSLDIR/data/standard/tissuepriors && \
cp -n HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz $FSLDIR/data/atlases/HarvardOxford


apt-get install -y ants


mkdir -p /opt/ICA-AROMA
curl -sL https://github.com/rhr-pruim/ICA-AROMA/archive/v0.4.3-beta.tar.gz | tar -xzC /opt/ICA-AROMA --strip-components 1
chmod +x /opt/ICA-AROMA/ICA_AROMA.py
export PATH=/opt/ICA-AROMA:$PATH


curl -sO https://repo.continuum.io/miniconda/Miniconda-3.8.3-Linux-x86_64.sh && \
bash Miniconda-3.8.3-Linux-x86_64.sh -b -p /usr/local/miniconda && \
rm Miniconda-3.8.3-Linux-x86_64.sh


export PATH=/usr/local/miniconda/bin:$PATH


conda install -y \
        blas

conda install -y  \
        cython==0.26 \
        jinja2==2.7.2 \
        matplotlib=2.0.2 \
        networkx==1.11 \
        nose==1.3.7 \
        numpy==1.13.0 \
        pandas==0.23.4 \
        scipy==0.19.1 \
        traits==4.6.0 \
        wxpython==3.0.0.0 \
        pip

pip install --upgrade pip==9.0.1
pip install -r /opt/C-PAC/requirements.txt
pip install xvfbwrapper

mkdir /ndmg_atlases && \
curl https://s3.amazonaws.com/mrneurodata/data/resources/ndmg_atlases.zip -o /ndmg_atlases/ndmg_atlases.zip && \
cd /ndmg_atlases && unzip /ndmg_atlases/ndmg_atlases.zip && \
rm /ndmg_atlases/ndmg_atlases.zip

apt-get clean && \
apt-get autoremove -y && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

pip install -e /opt/C-PAC