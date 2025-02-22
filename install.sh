#!/bin/bash

# Ensure a prefix is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <install-prefix>"
    exit 1
fi

# Set installation prefix
prefix=$(realpath "$1")

# Run fpm install with the provided prefix
fpm install --prefix="$prefix"

# Ensure installation was successful
if [ $? -ne 0 ]; then
    echo "fpm install failed"
    exit 1
fi

# Define paths
project_root=$(pwd)
include_src="$project_root/include"
include_dest="$prefix/include"

# Create destination include directory if it doesn't exist
mkdir -p "$include_dest"

# Copy C and C++ headers
cp "$include_src"/*.h* "$include_dest" 2>/dev/null

# Verify copy operation
if [ $? -eq 0 ]; then
    echo "Headers copied successfully to $include_dest"
else
    echo "Warning: No header files found in $include_src"
fi

