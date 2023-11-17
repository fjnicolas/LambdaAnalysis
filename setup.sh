#!/bin/bash

# Get the directory of the script
LAMBDAANA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

LAMBDAANA_DIR="$LAMBDAANA_DIR"
LAMBDAANA_SRC="$LAMBDAANA_DIR""/src"
LAMBDAANA_BUILD="$LAMBDAANA_DIR""/build"

# Set an environment variable with the script's directory
export LAMBDAANA_DIR="$LAMBDAANA_DIR"
export LAMBDAANA_SRC="$LAMBDAANA_SRC"
export LAMBDAANA_BUILD="$LAMBDAANA_BUILD"

echo "Top directory: $LAMBDAANA_DIR"
echo "Source directory: $LAMBDAANA_SRC"
echo "Build directory: $LAMBDAANA_BUILD"

alias runtpclines='$LAMBDAANA_BUILD/src/RunAlgoTPCLines'
alias runfrans='$LAMBDAANA_BUILD/src/RunAlgoFRANS'
alias runHandScan='sh $LAMBDAANA_DIR/scripts/runHandScanEvents.sh'
alias runHandScan2='sh $LAMBDAANA_DIR/scripts/runHandScanEvents2.sh'
#alias runLAMBDAANA='sh $LAMBDAANA_DIR/runLAMBDAANA.sh'
#alias runFRANS='sh $LAMBDAANA_DIR/runFRANS.sh'

#Search Paths
export LAMBDAANA_SEARCHPATH=$LAMBDAANA_DIR/scripts:.:$LAMBDAANA_DIR/src/config/
