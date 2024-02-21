#!/usr/bin/env bash
#./create_video.sh /run/media//Users/Macodemia/Desktop/FluidSim/presentation/water_droplet_without_gravity_with_surface_tension/raytraced_implicit/ videos/water_droplet/water_droplet_i

for raytraced_explicit in /run/media/Macodemia/Users/Macodemia/Desktop/FluidSim/presentation/**/raytraced_explicit/; do
    echo "$raytraced_explicit"
    path=$(dirname $raytraced_explicit | sed 's|/run/media/Macodemia/Users/Macodemia/Desktop/FluidSim/presentation/||')
    mkdir -p "videos/$path"
    file="$path""_explicit.mp4"
    ./create_video.sh "$raytraced_explicit" "videos/$path/$file"
done

for raytraced_implicit in /run/media/Macodemia/Users/Macodemia/Desktop/FluidSim/presentation/**/raytraced_implicit/; do
    echo "$raytraced_implicit"
    path=$(dirname $raytraced_implicit | sed 's|/run/media/Macodemia/Users/Macodemia/Desktop/FluidSim/presentation/||')
    mkdir -p "videos/$path"
    file="$path""_implicit.mp4"
    ./create_video.sh "$raytraced_implicit" "videos/$path/$file"
done
