#!/usr/bin/env bash
cd "$(dirname "$0")"

./build_release.sh
cd build/release
cd tests
./learnSPH_tests
