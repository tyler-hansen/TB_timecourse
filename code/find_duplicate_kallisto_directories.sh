#!/bin/bash

# Function to list all directories in a given folder
list_directories() {
  find "$1" -type d -exec basename {} \; 2>/dev/null
}

# Function to find duplicate directory names
find_duplicate_directories() {
  declare -A directory_count
  declare -A directory_locations

  for folder in "$1" "$2" "$3" "$4"; do
    dir_list=($(list_directories "$folder"))

    for dir in "${dir_list[@]}"; do
      ((directory_count["$dir"]++))
      directory_locations["$dir"]+=" $folder"
    done
  done

  # Find directories with count greater than 1
  duplicates=()
  for dir in "${!directory_count[@]}"; do
    if [ "${directory_count[$dir]}" -gt 1 ]; then
      duplicates+=("$dir")
    fi
  done

  # Print the duplicate directory names and their locations
  if [ ${#duplicates[@]} -eq 0 ]; then
    echo "No duplicate directories found."
  else
    echo "Duplicate directories found:"
    for duplicate in "${duplicates[@]}"; do
      locations="${directory_locations[$duplicate]}"
      echo "$duplicate found in:$locations"
    done
  fi
}

# Check if four folder paths are provided
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 folder1 folder2 folder3 folder4"
  exit 1
fi

# Call the function to find duplicate directory names
find_duplicate_directories "$1" "$2" "$3" "$4"
