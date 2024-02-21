#!/usr/bin/env bash

# Check for correct amount of parameters
if [ "$#" -ne 2 ]; then
    echo "USAGE: ./create_video.sh input_directory output.mp4"
    exit
fi

input1=$1
output=$2

echo "Reading pic.%04d.png files from: $input1"
echo "Writing output video to: $output"
echo "Running command: ffmpeg -r 30 -i $1/pic.%04d.png $output -y"

ffmpeg -r 30 -i $1/pic.%04d.png $output -y
