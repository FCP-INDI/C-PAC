#!/bin/bash

set -e

wget -O /etc/apt/sources.list.d/neurodebian.sources.list http://neuro.debian.net/lists/bionic.us-ca.full
for i in {1..5}; do 
    apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9 && break || sleep 5;
done

apt-get update
apt-get install -y x2goserver lubuntu-desktop lxde-icon-theme xvfb
apt-get remove -y lxlock xscreensaver xscreensaver-data gnome-screensaver
apt-get upgrade

# Cleaning up NetworkManager, use netplan
systemctl stop NetworkManager.service
systemctl disable NetworkManager.service
rm -Rf /etc/NetworkManager

# Cleaning up NetworkManager, use netplan
systemctl stop NetworkManager.service
systemctl disable NetworkManager.service
rm -Rf /etc/NetworkManager

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

# Xvfb :99 & export DISPLAY=:99
# su -c 'lxsession' ubuntu &  # to create configs
# sleep 5

# sed -z 's/\s*Button\s*{\s*id=lxde-screenlock.desktop\s*}//g' /home/ubuntu/.config/lxpanel/LXDE/panels/panel

# Disable screen lock
mkdir -p /home/ubuntu/.config/xfce4/xfconf/xfce-perchannel-xml/
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

mkdir -p /opt/c3d
curl -sSL "http://downloads.sourceforge.net/project/c3d/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz" | tar -xzC /opt/c3d --strip-components 1

echo 'C3DPATH=/opt/c3d/' >> /etc/bash.bashrc
echo 'PATH=$C3DPATH/bin:$PATH' >> /etc/bash.bashrc

# libpng12 for afni
wget -q -O /tmp/libpng12.deb http://mirrors.kernel.org/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1_amd64.deb \
  && dpkg -i /tmp/libpng12.deb \
  && rm /tmp/libpng12.deb

libs_path=/usr/lib/x86_64-linux-gnu 
if [ -f $libs_path/libgsl.so.19 ]; then \
    ln $libs_path/libgsl.so.19 $libs_path/libgsl.so.0; \
fi
mkdir -p /opt/afni
curl -sO https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz
tar zxv -C /opt/afni --strip-components=1 -f linux_openmp_64.tgz
rm -rf linux_openmp_64.tgz

echo 'PATH=/opt/afni:$PATH' >> /etc/bash.bashrc

apt-get install -y --no-install-recommends \
    fsl-core \
    fsl-atlases \
    fsl-mni152-templates

source /etc/fsl/5.0/fsl.sh
echo 'source /etc/fsl/5.0/fsl.sh' >> /etc/bash.bashrc

curl -sL http://fcon_1000.projects.nitrc.org/indi/cpac_resources.tar.gz -o /tmp/cpac_resources.tar.gz && \
tar xfz /tmp/cpac_resources.tar.gz -C /tmp && \
cp -n /tmp/cpac_image_resources/MNI_3mm/* $FSLDIR/data/standard && \
cp -n /tmp/cpac_image_resources/MNI_4mm/* $FSLDIR/data/standard && \
cp -n /tmp/cpac_image_resources/symmetric/* $FSLDIR/data/standard && \
cp -nr /tmp/cpac_image_resources/tissuepriors/2mm $FSLDIR/data/standard/tissuepriors && \
cp -nr /tmp/cpac_image_resources/tissuepriors/3mm $FSLDIR/data/standard/tissuepriors && \
cp -n /tmp/cpac_image_resources/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz $FSLDIR/data/atlases/HarvardOxford

apt-get install -y ants

mkdir -p /opt/ICA-AROMA
curl -sL https://github.com/rhr-pruim/ICA-AROMA/archive/v0.4.3-beta.tar.gz | tar -xzC /opt/ICA-AROMA --strip-components 1
chmod +x /opt/ICA-AROMA/ICA_AROMA.py
echo 'PATH=/opt/ICA-AROMA:$PATH' >> /etc/bash.bashrc

curl -sO https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
bash Miniconda2-latest-Linux-x86_64.sh -b -p /usr/local/miniconda && \
rm Miniconda2-latest-Linux-x86_64.sh

echo 'PATH=/usr/local/miniconda/bin:$PATH' >> /etc/bash.bashrc
export PATH=/usr/local/miniconda/bin:$PATH

conda install -y  \
    pip \
    wxpython==3.0.0.0

pip install --upgrade pip==9.0.1
pip install -r /opt/C-PAC/requirements.txt
pip install xvfbwrapper

curl -L https://github.com/neurodata/neuroparc/archive/master.zip -o /tmp/neuroparc.zip && \
unzip /tmp/neuroparc.zip -d /tmp/neuroparc 'neuroparc-master/atlases/*' && \
cp -r /tmp/neuroparc/neuroparc-master/atlases /ndmg_atlases && rm /tmp/neuroparc.zip

pip install /opt/C-PAC

conda install -y -c anaconda glib
conda install -y -c conda-forge harfbuzz

cpac
