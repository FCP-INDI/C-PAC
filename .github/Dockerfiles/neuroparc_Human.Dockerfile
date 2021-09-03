# neuroparc:v1.0-human

# using neurodebian runtime as parent image
FROM neurodebian:bionic-non-free

LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y git

RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install git-lfs
RUN git lfs install

# Get atlases
RUN mkdir -p /ndmg_atlases/label && \
    GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 --branch v1.0 https://github.com/neurodata/neuroparc.git /tmp/neuroparc && \
    cd /tmp/neuroparc && \
    git lfs install --skip-smudge && \
    git lfs pull -I "atlases/label/Human/*" && \
    cp -r /tmp/neuroparc/atlases/label/Human /ndmg_atlases/label && \
    cd -
