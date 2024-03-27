#!/bin/bash
# Copyright (C) 2021-2024  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
FROM ghcr.io/fcp-indi/c-pac/stage-base:fmriprep-lts-v1.8.8.dev1
LABEL org.opencontainers.image.description "Full C-PAC image with software dependencies version-matched to [fMRIPrep LTS](https://reproducibility.stanford.edu/fmriprep-lts#long-term-support-lts)"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
USER root

# install C-PAC & set up runscript
COPY dev/circleci_data/pipe-test_ci.yml /cpac_resources/pipe-test_ci.yml
COPY . /code
RUN pip cache purge && pip install -e /code
# set up runscript
COPY dev/docker_data /code/docker_data
RUN rm -Rf /code/docker_data/checksum && \
    mv /code/docker_data/* /code && \
    rm -Rf /code/docker_data && \
    chmod +x /code/CPAC/_entrypoints/run.py && \
    chmod +x /code/CPAC/_entrypoints/run-with-freesurfer.sh
ENTRYPOINT ["/code/CPAC/_entrypoints/run-with-freesurfer.sh"]

# link libraries & clean up
RUN sed -i 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && locale-gen \
    && apt-get clean \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/pip/* \
    && ldconfig \
    && chmod 777 / \
    && chmod -R 777 /home/c-pac_user \
    && chmod 777 $(ls / | grep -v sys | grep -v proc)
ENV PYTHONUSERBASE=/home/c-pac_user/.local
ENV PATH=$PATH:/home/c-pac_user/.local/bin \
    PYTHONPATH=$PYTHONPATH:$PYTHONUSERBASE/lib/python3.10/site-packages

# set user
WORKDIR /home/c-pac_user
# USER c-pac_user
