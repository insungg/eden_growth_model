#!/bin/bash

PARAM_FILE="parameters_rarb13.txt"
OUTPUT_DIR="../data"
EXECUTABLE="./eden_simulator"

for i in {2..20}; do
    # Modify fileName in parameters_normal.txt
    sed -i.bak "s/^fileName=.*/fileName=rarb13hs${i}/" $PARAM_FILE

    # Run the simulation
    $EXECUTABLE $PARAM_FILE $OUTPUT_DIR

    # Restore the original parameters_normal.txt
    mv ${PARAM_FILE}.bak $PARAM_FILE
done
