#!/bin/bash

# Get the directory of the script
TPCLINES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set an environment variable with the script's directory
export TPCLINES_DIR="$TPCLINES_DIR"

alias runtpclines='sh $TPCLINES_DIR/runtpclines.sh'
