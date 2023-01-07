#!/usr/bin/env bash

# Set the root directory to search
root_dir="$1"

# Set the list of file extensions to exclude
exclude_extensions="$2"

# Build the find command
find_cmd="find $root_dir -type f"

# Loop through the list of file extensions and add a -not -name clause for each one
for ext in $(echo $exclude_extensions | sed "s/,/ /g")
do
  find_cmd="$find_cmd -not -name \*.$ext"
done

# Execute the find command and compute the MD5 hash of each file
output=$(eval $find_cmd -exec md5sum {} +)

# Loop through the output and generate the desired text for each file
while read -r line; do
  path=$(echo $line | awk '{print $2}')
  hash=$(echo $line | awk '{print $1}')
  printf "assert new File(\"%s\").exists()\n" "$path"
  printf "assert path(\"%s\").md5 == \"%s\"\n" "$path" "$hash"
  echo
done <<< "$output"
