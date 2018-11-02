#!/bin/bash

VERSION=$1

# if ! [[ "$VERSION" =~ "[0-9]+\.[0-9]+\.[0-9]+" ]]; then
#     echo "Invalid version: $VERSION"
#     echo "Use format NN.NN.NN"
#     exit
# fi

docker build -t childmind/c-pac:latest .
docker tag childmind/c-pac:latest childmind/c-pac:v$VERSION
docker push childmind/c-pac:latest
docker push childmind/c-pac:v$VERSION