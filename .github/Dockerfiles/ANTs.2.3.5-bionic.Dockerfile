FROM ubuntu:bionic-20200112 as builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    software-properties-common \
                    build-essential \
                    apt-transport-https \
                    ca-certificates \
                    curl \
                    gnupg \
                    software-properties-common \
                    wget \
                    ninja-build \
                    git \
                    unzip \
                    zlib1g-dev

RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null \
    | apt-key add - \
  && apt-add-repository -y 'deb https://apt.kitware.com/ubuntu/ bionic main' \
  && apt-get update \
  && apt-get -y install cmake=3.18.3-0kitware1 cmake-data=3.18.3-0kitware1 

RUN mkdir -p /tmp/ants/build \
    && cd /tmp/ants/build \
    && mkdir -p /opt/ants \
    && git config --global url."https://".insteadOf git:// \
    && git clone -b v2.3.5 --depth 1 https://github.com/ANTsX/ANTs.git /tmp/ants/ANTs \
    && cmake \
        -DCMAKE_INSTALL_PREFIX=/opt/ants \
        -DBUILD_SHARED_LIBS=OFF \
        -DUSE_VTK=OFF \
        -DSuperBuild_ANTS_USE_GIT_PROTOCOL=OFF \
        -DBUILD_TESTING=OFF \
        -DRUN_LONG_TESTS=OFF \
        -DRUN_SHORT_TESTS=OFF \
      ../ANTs 2>&1 | tee cmake.log \
    && make -j 4 2>&1 | tee build.log \
    && cd ANTS-build \
    && make install 2>&1 | tee install.log

# download OASIS templates for niworkflows-ants skullstripping
RUN mkdir /ants_template && \
    curl -sL https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/3133832/Oasis.zip -o /tmp/Oasis.zip && \
    unzip /tmp/Oasis.zip -d /tmp &&\
    mv /tmp/MICCAI2012-Multi-Atlas-Challenge-Data /ants_template/oasis && \
    rm -rf /tmp/Oasis.zip /tmp/MICCAI2012-Multi-Atlas-Challenge-Data

FROM ubuntu:bionic-20200112
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD: ANTs 2.3.5 stage"
COPY --from=builder /opt/ants /opt/ants
COPY --from=builder /ants_template /ants_template

ENV ANTSPATH="/opt/ants/bin" \
    PATH="/opt/ants/bin:$PATH" \
    LD_LIBRARY_PATH="/opt/ants/lib:$LD_LIBRARY_PATH"
RUN apt-get update \
    && apt install -y --no-install-recommends zlib1g-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /data

CMD ["/bin/bash"]
