#!/usr/bin/env bash
cd "$(dirname "$0")"

./build.sh
cd build/debug
cd tests
./learnSPH_tests
