FROM docker.io/fcpindi/c-pac:release-v1.3.0
LABEL org.opencontainers.image.authors="The C-PAC Team <CNL@childmind.org>"

# install cpac
COPY . /code
RUN pip install -e /code
