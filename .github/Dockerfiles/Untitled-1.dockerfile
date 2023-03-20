FROM ghcr.io/fcp-indi/c-pac/fsl:6.0.4-python3.10-bionic as FSL
FROM ghcr.io/fcp-indi/c-pac/ubuntu:python3.10-bionic-non-free

ENV FSLDIR=/usr/share/fsl/6.0 \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=/usr/share/fsl/6.0 \
    LD_LIBRARY_PATH=/usr/lib/fsl/6.0:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    PATH=/usr/share/fsl/6.0/bin:$PATH \
    TZ=America/New_York
COPY --from=FSL /usr/bin/tclsh /usr/bin/tclsh
COPY --from=FSL /usr/bin/wish /usr/bin/wish
COPY --from=FSL /usr/share/fsl /usr/share/fsl
COPY --from=FSL /usr/lib /usr/lib
COPY --from=FSL /lib/x86_64-linux-gnu/lib*so* /lib/x86_64-linux-gnu/
COPY --from=FSL /usr/share/fsl/5.0/data/standard/tissuepriors/*mm /usr/share/fsl/6.0/data/standard/tissuepriors/

ENTRYPOINT [ "/bin/bash" ]
USER c-pac_user
