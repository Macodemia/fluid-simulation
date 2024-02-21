#!/usr/bin/env bash

# Check for correct amount of parameters
if [ "$#" -ne 3 ]; then
    echo "USAGE: ./create_side_by_side_video.sh video1.mp4 video2.mp4 output.mp4"
    exit
fi

input1=$1
input2=$2
output=$3

# Validate parameters (check if file exists)
if [[ ! -f "$1" ]]; then
    echo "Passed file ($input1) does not exist!"
    exit
fi
if [[ ! -f "$2" ]]; then
    echo "Passed file ($input2) does not exist!"
    exit
fi

ffmpeg -i $input1 -i $input2 -filter_complex hstack $output
