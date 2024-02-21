#!/usr/bin/env bash
cd "$(dirname "$0")"

mkdir -p build/debug
cd build/debug
cmake ../../ -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_BUILD_TYPE=Debug

cp compile_commands.json ../../ 

core_count=$(getconf _NPROCESSORS_ONLN)
make -j$core_count
