#!/bin/bash

VERSION=$1

# if ! [[ "$VERSION" =~ "[0-9]+\.[0-9]+\.[0-9]+" ]]; then
#     echo "Invalid version: $VERSION"
#     echo "Use format NN.NN.NN"
#     exit
# fi

docker build -t fcpindi/c-pac:latest .
docker tag fcpindi/c-pac:latest fcpindi/c-pac:v$VERSION
# docker push fcpindi/c-pac:latest
# docker push fcpindi/c-pac:v$VERSION