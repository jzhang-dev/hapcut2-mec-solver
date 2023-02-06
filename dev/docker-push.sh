#!/usr/bin/env bash

docker build -t hapcut2_mec_solver . &&
echo $DOCKER_HUB_TOKEN | docker login -u jzhang0246 --password-stdin &&
docker tag hapcut2_mec_solver:latest jzhang0246/hapcut2-mec-solver:$1
docker tag hapcut2_mec_solver:latest jzhang0246/hapcut2-mec-solver:latest
docker push jzhang0246/hapcut2-mec-solver:$1
docker push jzhang0246/hapcut2-mec-solver:latest
