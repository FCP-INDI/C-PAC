# Copyright (C) 2023  C-PAC Developers

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
FROM neurodebian:bionic-non-free as builder
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
COPY dev/docker_data/checksum/Python.3.10.12.sha384 /tmp/checksum.sha384
RUN apt-get update \
    && apt-get install -y \
      build-essential \
      libbz2-dev \
      libffi-dev \
      libgdbm-dev \
      liblzma-dev \
      libncurses5-dev \
      libreadline-dev \
      libsqlite3-dev \
      libssl-dev \
      tk-dev \
      uuid-dev \
      wget \
      zlib1g-dev \
    && wget https://www.python.org/ftp/python/3.10.12/Python-3.10.12.tgz \
    && sha384sum --check /tmp/checksum.sha384 \
    && tar xvf Python-3.10.12.tgz \
    && cd Python-3.10.12 \
    && ./configure --prefix=/usr/local/python3.10 --enable-optimizations \
    && make && make install

# only keep what we need
FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
Ubuntu Bionic base image"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=builder /usr/local/python3.10 /usr/local/python3.10
