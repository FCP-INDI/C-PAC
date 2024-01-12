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
FROM ghcr.io/fcp-indi/c-pac/ubuntu:jammy-non-free as base

USER root

RUN apt-get update \
    && apt-get install -y connectome-workbench

USER c-pac_user

# FROM scratch
# LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
# connectome-workbench 1.5.0 stage"
# LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
# COPY --from=base /lib/x86_64-linux-gnu/ld-linux-x86-64.so.2 /lib/x86_64-linux-gnu/ld-linux-x86-64.so.2
# COPY --from=base /lib/x86_64-linux-gnu/libGL.so.1 /lib/x86_64-linux-gnu/libGL.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libGLU.so.1 /lib/x86_64-linux-gnu/libGLU.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libGLX.so.0 /lib/x86_64-linux-gnu/libGLX.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libGLdispatch.so.0 /lib/x86_64-linux-gnu/libGLdispatch.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libLLVM-15.so.1 /lib/x86_64-linux-gnu/libLLVM-15.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libOSMesa.so.8 /lib/x86_64-linux-gnu/libOSMesa.so.8
# COPY --from=base /lib/x86_64-linux-gnu/libOpenGL.so.0 /lib/x86_64-linux-gnu/libOpenGL.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libQt5Core.so.5 /lib/x86_64-linux-gnu/libQt5Core.so.5
# COPY --from=base /lib/x86_64-linux-gnu/libQt5Gui.so.5 /lib/x86_64-linux-gnu/libQt5Gui.so.5
# COPY --from=base /lib/x86_64-linux-gnu/libQt5Network.so.5 /lib/x86_64-linux-gnu/libQt5Network.so.5
# COPY --from=base /lib/x86_64-linux-gnu/libQt5Xml.so.5 /lib/x86_64-linux-gnu/libQt5Xml.so.5
# COPY --from=base /lib/x86_64-linux-gnu/libX11.so.6 /lib/x86_64-linux-gnu/libX11.so.6
# COPY --from=base /lib/x86_64-linux-gnu/libXau.so.6 /lib/x86_64-linux-gnu/libXau.so.6
# COPY --from=base /lib/x86_64-linux-gnu/libXdmcp.so.6 /lib/x86_64-linux-gnu/libXdmcp.so.6
# COPY --from=base /lib/x86_64-linux-gnu/libbrotlicommon.so.1 /lib/x86_64-linux-gnu/libbrotlicommon.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libbrotlidec.so.1 /lib/x86_64-linux-gnu/libbrotlidec.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libbsd.so.0 /lib/x86_64-linux-gnu/libbsd.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libc.so.6 /lib/x86_64-linux-gnu/libc.so.6
# COPY --from=base /lib/x86_64-linux-gnu/libcom_err.so.2 /lib/x86_64-linux-gnu/libcom_err.so.2
# COPY --from=base /lib/x86_64-linux-gnu/libdouble-conversion.so.3 /lib/x86_64-linux-gnu/libdouble-conversion.so.3
# COPY --from=base /lib/x86_64-linux-gnu/libedit.so.2 /lib/x86_64-linux-gnu/libedit.so.2
# COPY --from=base /lib/x86_64-linux-gnu/libffi.so.8 /lib/x86_64-linux-gnu/libffi.so.8
# COPY --from=base /lib/x86_64-linux-gnu/libfreetype.so.6 /lib/x86_64-linux-gnu/libfreetype.so.6
# COPY --from=base /lib/x86_64-linux-gnu/libftgl.so.2 /lib/x86_64-linux-gnu/libftgl.so.2
# COPY --from=base /lib/x86_64-linux-gnu/libgcc_s.so.1 /lib/x86_64-linux-gnu/libgcc_s.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libglapi.so.0 /lib/x86_64-linux-gnu/libglapi.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libglib-2.0.so.0 /lib/x86_64-linux-gnu/libglib-2.0.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libgomp.so.1 /lib/x86_64-linux-gnu/libgomp.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libgraphite2.so.3 /lib/x86_64-linux-gnu/libgraphite2.so.3
# COPY --from=base /lib/x86_64-linux-gnu/libgssapi_krb5.so.2 /lib/x86_64-linux-gnu/libgssapi_krb5.so.2
# COPY --from=base /lib/x86_64-linux-gnu/libharfbuzz.so.0 /lib/x86_64-linux-gnu/libharfbuzz.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libicudata.so.70 /lib/x86_64-linux-gnu/libicudata.so.70
# COPY --from=base /lib/x86_64-linux-gnu/libicui18n.so.70 /lib/x86_64-linux-gnu/libicui18n.so.70
# COPY --from=base /lib/x86_64-linux-gnu/libicuuc.so.70 /lib/x86_64-linux-gnu/libicuuc.so.70
# COPY --from=base /lib/x86_64-linux-gnu/libk5crypto.so.3 /lib/x86_64-linux-gnu/libk5crypto.so.3
# COPY --from=base /lib/x86_64-linux-gnu/libkeyutils.so.1 /lib/x86_64-linux-gnu/libkeyutils.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libkrb5.so.3 /lib/x86_64-linux-gnu/libkrb5.so.3
# COPY --from=base /lib/x86_64-linux-gnu/libkrb5support.so.0 /lib/x86_64-linux-gnu/libkrb5support.so.0
# COPY --from=base /lib/x86_64-linux-gnu/liblzma.so.5 /lib/x86_64-linux-gnu/liblzma.so.5
# COPY --from=base /lib/x86_64-linux-gnu/libm.so.6 /lib/x86_64-linux-gnu/libm.so.6
# COPY --from=base /lib/x86_64-linux-gnu/libmd.so.0 /lib/x86_64-linux-gnu/libmd.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libmd4c.so.0 /lib/x86_64-linux-gnu/libmd4c.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libpcre.so.3 /lib/x86_64-linux-gnu/libpcre.so.3
# COPY --from=base /lib/x86_64-linux-gnu/libpcre2-16.so.0 /lib/x86_64-linux-gnu/libpcre2-16.so.0
# COPY --from=base /lib/x86_64-linux-gnu/libpng16.so.16 /lib/x86_64-linux-gnu/libpng16.so.16
# COPY --from=base /lib/x86_64-linux-gnu/libresolv.so.2 /lib/x86_64-linux-gnu/libresolv.so.2
# COPY --from=base /lib/x86_64-linux-gnu/libstdc++.so.6 /lib/x86_64-linux-gnu/libstdc++.so.6
# COPY --from=base /lib/x86_64-linux-gnu/libtinfo.so.6 /lib/x86_64-linux-gnu/libtinfo.so.6
# COPY --from=base /lib/x86_64-linux-gnu/libxcb.so.1 /lib/x86_64-linux-gnu/libxcb.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libxml2.so.2 /lib/x86_64-linux-gnu/libxml2.so.2
# COPY --from=base /lib/x86_64-linux-gnu/libz.so.1 /lib/x86_64-linux-gnu/libz.so.1
# COPY --from=base /lib/x86_64-linux-gnu/libzstd.so.1 /lib/x86_64-linux-gnu/libzstd.so.1
# COPY --from=base /lib64/ld-linux-x86-64.so.2 /lib64/ld-linux-x86-64.so.2
# COPY --from=base /usr/bin/wb_* /usr/bin
# COPY --from=base /usr/share/bash-completion/completions/wb_command /usr/share/bash-completion/completions/wb_command
# COPY --from=base /usr/share/doc/connectome-workbench /usr/share/doc/connectome-workbench
# COPY --from=base /usr/share/bash-completion/completions/wb_shortcuts /usr/share/bash-completion/completions/wb_shortcuts
