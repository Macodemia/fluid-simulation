#!/usr/bin/env bash
cd "$(dirname "$0")"

mkdir -p build/release
cd build/release
cmake ../../ -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_BUILD_TYPE=Release

cp compile_commands.json ../../ 

core_count=$(getconf _NPROCESSORS_ONLN)
make -j$core_count
