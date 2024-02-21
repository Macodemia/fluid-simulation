#!/usr/bin/env bash
cd "$(dirname "$0")"

./build_debug.sh

echo ""
echo "-------------------------------------------------------------------"
echo ""

if [ -z "$2" ]
then
    SCENARIO_FILE=""
else
    SCENARIO_FILE="$(pwd)/$2"
fi

if [ -z "$3" ]
then
    THREAD_COUNT=""
else
    THREAD_COUNT="$3"
fi

cd build/debug
cd app
./learnSPH_app $1 $SCENARIO_FILE $THREAD_COUNT
