#!/bin/bash

# Function to convert symlinked subpaths to real paths
convert_symlinks() {
  # Iterate over all items (files and directories) in the current directory
  for ITEM in "$1"/*; do
    echo $ITEM
    if [[ -L "$ITEM" ]]; then
      # If the item is a symlink, convert it to a real path
      TARGET=$(readlink -f "$ITEM")
      rm "$ITEM"
      cp -r "$TARGET" "$ITEM"
      echo "Converted symlink to real path: $ITEM"
    elif [[ -d "$ITEM" ]]; then
      # If the item is a directory, recursively call the function on it
      convert_symlinks "$ITEM"
    fi
  done
}

convert_symlinks $1
