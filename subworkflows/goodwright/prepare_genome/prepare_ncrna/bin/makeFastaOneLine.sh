#!/bin/bash

# Make a fasta into one-line version

file=$1
new_name=$2

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${file} | tail -n +2 > ${new_name}

