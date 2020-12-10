FROM neurodebian:bionic-non-free

USER root

ARG DEBIAN_FRONTEND=noninteractive

ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ANTSPATH="/usr/lib/ants/bin" \
    PATH="/usr/lib/ants/bin:$PATH" \
    LD_LIBRARY_PATH="/usr/lib/ants/lib:$LD_LIBRARY_PATH"

RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           apt-utils \
           bzip2 \
           ca-certificates \
           curl \
           locales \
           unzip \
           cmake \
           g++ \
           gcc \
           git \
           make \
           zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /tmp/ants/build \
    && git clone https://github.com/ANTsX/ANTs.git /tmp/ants/source \
    && cd /tmp/ants/source \
    && git fetch --tags \
    && git checkout v2.3.4 \
    && cd /tmp/ants/build \
    && cmake -DBUILD_SHARED_LIBS=ON /tmp/ants/source \
    && make -j1 \
    && mkdir -p /usr/lib/ants \
    && mv bin lib /usr/lib/ants/ \
    && mv /tmp/ants/source/Scripts/* /usr/lib/ants/bin \
    && rm -rf /tmp/ants
