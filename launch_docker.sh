#!/usr/bin/env bash

IMAGE_NAME=tcell-activation

docker build -t ${IMAGE_NAME} .
docker run --rm --name tcell \
    -v /Users/natalie/Documents/projects/guido/RP:/app/data/rp:ro \
    -v /Users/natalie/Documents/projects/guido/RNA:/app/data/rna:ro \
    -it ${IMAGE_NAME} /bin/bash
