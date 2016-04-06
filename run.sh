#!/bin/bash

docker run -it --rm -v ~/Desktop/:/mnt -v ~/git/zymo/mirror-seq/:/scripts/mirror-seq -v ~/.aws/:/root/.aws/ hunterchung/zymo-docker bash
