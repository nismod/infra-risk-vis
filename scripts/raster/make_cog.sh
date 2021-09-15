#!/bin/env bash

# TODO validate arguments

INPUT_FILE=$1
OUTPUT_FILE=$2


# TODO tune COG conversion parameters

gdal_translate "$INPUT_FILE" "$OUTPUT_FILE" -of COG