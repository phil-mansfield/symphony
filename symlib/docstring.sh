#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 2 ]; then
	echo "Usage: ./docstring.sh <target_file> [target_class.]<target_function>"
	exit 1
fi

target_file="$1"
target_func="$2"

python3 ../docs/scripts/docstring.py $target_file $target_func
