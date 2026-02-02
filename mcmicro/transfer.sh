#!/bin/bash

# Destination folder where the files will be copied
DESTINATION_FOLDER="~/projects/jhu/rawcsv/"

# Find and copy the files
find . -type f -path "./LSP*/quantification/*cellRing.csv" -exec cp {} "~/projects/jhu/rawcsv/" \;

