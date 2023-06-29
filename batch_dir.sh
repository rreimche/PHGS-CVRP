#!/bin/bash

num_subdirectories=$1

# Create subdirectories if they don't exist
for ((i = 1; i <= num_subdirectories; i++)); do
  dir="$i"
  if [ ! -d "$dir" ]; then
    mkdir "$dir"
  fi
done

# Execute the desired command in each subdirectory
for dir in */; do
  if [ -d "$dir" ]; then
    cd "$dir" || exit
    ../../batch.sh ../../build/hgs ../../Instances/CVRP 1
    cd ..
  fi
done