FROM docker.io/fcpindi/c-pac:release-v1.3.0
LABEL org.opencontainers.image.authors="The C-PAC Team <CNL@childmind.org>" \
      org.opencontainers.image.version="1.3.0.post3.dev1"
ARG BIDS_VALIDATOR_VERSION="1.2.3"

# Install the validator
RUN apt-get update && \
    apt-get install -y curl && \
    curl -sL https://deb.nodesource.com/setup_8.x | bash - && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    npm install -g "bids-validator@${BIDS_VALIDATOR_VERSION}"

# install cpac
COPY . /code
RUN pip install -e /code
