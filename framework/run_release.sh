#!/usr/bin/env bash
cd "$(dirname "$0")"

./build_release.sh

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

cd build/release
cd app

OUTPUT_FILE="$(echo "../res/$2/$2_$1_output" | sed 's/scenarios\///g' | sed 's/.txt//g')"
OUTPUT_FILE+=".txt"
mkdir -p "$(dirname $OUTPUT_FILE)"
echo "Output is also written to $OUTPUT_FILE"

./learnSPH_app $1 $SCENARIO_FILE $THREAD_COUNT | tee "$OUTPUT_FILE"
