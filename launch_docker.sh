#!/usr/bin/env bash

IMAGE_NAME=tcell-activation

docker build -t ${IMAGE_NAME} . && docker run \
    -v tcell_data:/app/data \
    -it ${IMAGE_NAME} /bin/bash
