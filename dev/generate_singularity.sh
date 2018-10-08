#!/bin/bash

OUT=`mktemp -d`

docker run \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${OUT}:/output \
    --privileged -t --rm \
    singularityware/docker2singularity \
    ${1}

find ${OUT} -name "*.img"