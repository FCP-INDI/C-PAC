FROM ghcr.io/fcp-indi/c-pac/ubuntu:xenial-20200114
USER root

# Installing and setting up c3d
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    convert3d=0.0.20190204-1~nd16.04+1
ENV C3DPATH /usr/bin
ENV PATH $C3DPATH/bin:$PATH

ENTRYPOINT ["/bin/bash"]

# Link libraries for Singularity images
RUN ldconfig

RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

FROM scratch
LABEL org.opencontainers.image.description "NOT INTENDED FOR USE OTHER THAN AS A STAGE IMAGE IN A MULTI-STAGE BUILD \
c3d 1.0.0 (Xenial) stage"
LABEL org.opencontainers.image.source https://github.com/FCP-INDI/C-PAC
COPY --from=c3d /usr/lib/libITK* /usr/lib/
COPY --from=c3d /usr/lib/libitk* /usr/lib/
COPY --from=c3d /usr/lib/x86_64-linux-gnu/libgdcm* /usr/lib/x86_64-linux-gnu/
COPY --from=c3d /usr/lib/x86_64-linux-gnu/libfftw3* /usr/lib/x86_64-linux-gnu/
COPY --from=c3d /lib/x86_64-linux-gnu/lib*.so* /lib/x86_64-linux-gnu/
COPY --from=c3d /usr/lib/x86_64-linux-gnu/lib*.so* /usr/lib/x86_64-linux-gnu/
COPY --from=c3d /lib64/ld-linux-x86-64.so.2 /lib64/ld-linux-x86-64.so.2
COPY --from=c3d /usr/bin/c*d /usr/bin/
COPY --from=c3d /usr/bin/c3d_* /usr/bin/
COPY --from=c3d /usr/share/doc/convert3d /usr/share/doc/convert3d
COPY --from=c3d /usr/lib/c3d_gui-1.1.0/Convert3DGUI /usr/lib/c3d_gui-1.1.0/Convert3DGUI
