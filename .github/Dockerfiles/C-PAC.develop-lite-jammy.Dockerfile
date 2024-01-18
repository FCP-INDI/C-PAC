# Copyright (C) 2022-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
FROM ghcr.io/fcp-indi/c-pac/stage-base:lite-v1.8.7dev1
LABEL org.opencontainers.image.description "Full C-PAC image without FreeSurfer"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
USER root

# install C-PAC
COPY dev/circleci_data/pipe-test_ci.yml /cpac_resources/pipe-test_ci.yml
COPY . /code
COPY --from=ghcr.io/fcp-indi/c-pac_templates:latest /cpac_templates /cpac_templates
RUN pip cache purge && pip install -e /code
# set up runscript
COPY dev/docker_data /code/docker_data
RUN rm -Rf /code/docker_data/checksum && \
    mv /code/docker_data/* /code && \
    rm -Rf /code/docker_data && \
    chmod +x /code/run.py && \
    rm -Rf /code/run-with-freesurfer.sh
ENTRYPOINT ["/code/run.py"]

# link libraries & clean up
# link libraries & clean up
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/* \
    && find / -type f -print0 | sort -t/ -k2 | xargs -0 rdfind -makehardlinks true \
    && rm -rf results.txt \
    && apt-get remove rdfind -y \
    && apt-get clean \
    && apt-get autoremove -y \
    && ldconfig \
    && chmod 777 / \
    && chmod 777 $(ls / | grep -v sys | grep -v proc)
ENV PYTHONUSERBASE=/home/c-pac_user/.local
ENV PATH=$PATH:/home/c-pac_user/.local/bin \
    PYTHONPATH=$PYTHONPATH:$PYTHONUSERBASE/lib/python3.10/site-packages

# set user
WORKDIR /home/c-pac_user
# USER c-pac_user
