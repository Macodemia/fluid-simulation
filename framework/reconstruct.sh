#!/usr/bin/env bash
./build_release.sh

PWD=$(pwd)

FILES=$@

SORTED_FILES=$(echo $FILES | tr " " "\n"| sort -V | sed "s|^|$PWD/|")

cd build/release/reconstruction

./learnSPH_reconstruction $SORTED_FILES
