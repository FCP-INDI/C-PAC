Bootstrap: docker
From: fcpindi/c-pac:latest
IncludeCmd: yes

%environment
FREESURFER_HOME=

%post

ln -s /usr/lib/x86_64-linux-gnu/libgsl.so.23.0.0 \
      /usr/lib/x86_64-linux-gnu/libgsl.so.23 || echo "Link exists"

ln -s /usr/lib/x86_64-linux-gnu/libgsl.so.23.0.0 \
      /usr/lib/x86_64-linux-gnu/libgsl.so.0 || echo "Link exists"

if [ -e /usr/lib/x86_64-linux-gnu/libGL.so.1 ] ; then
    rm -f /usr/lib/x86_64-linux-gnu/libGL.so.1 || echo "Link exists"
fi

ln -s /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1.2.0 \
   /usr/lib/x86_64-linux-gnu/libGL.so.1 || echo "Link exists"
