#!/bin/bash

EXECUTABLE_FILE="$1"  # Use the first command-line argument as the executable file path
DIRECTORY="$2"  # Replace with the actual directory path
LIMIT="$3"  # Use the third command-line argument as the limit for the loop

count=0  # Initialize a counter variable


for FILENAME in "$DIRECTORY"/*; do
  if [ -f "$FILENAME" ]; then
    N_NODES=$(grep 'DIMENSION' "$FILENAME" | awk -F' ' '{print $3}')
    TIME_LIMIT=$(awk "BEGIN {printf \"%.0f\n\", $N_NODES * 240/100}")
    FILE_PATH=$(dirname "$FILENAME")
    FILE_NAME=$(basename "$FILENAME")
    #echo "FILE: "$FILE_PATH"/"$FILE_NAME", CAPACITY: $N_NODES, TIME_LIMIT: $TIME_LIMIT"
    echo "" "EXECUTING:$EXECUTABLE_FILE" "$FILE_PATH"/"$FILE_NAME" "$FILE_NAME.sol" -seed 1 -t "$TIME_LIMIT" ""
    "$EXECUTABLE_FILE" "$FILE_PATH"/"$FILE_NAME" "$FILE_NAME.sol" -seed 1 -t "$TIME_LIMIT"

    count=$((count + 1))  # Increment the counter

    # Check if the counter exceeds the limit
    if [ ! -z "$LIMIT" ] && [ "$count" -ge "$LIMIT" ]; then
      break  # Exit the loop if the limit is reached
    fi

  fi
done