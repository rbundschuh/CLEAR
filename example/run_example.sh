#!/bin/bash

# Example job script for CLEAR
# Python scripts must be present in $PATH

# Assumes that the bed files have already been generated and are in current
# working directory

for i in *.bed;
do
  echo "python2 ../make_dat.py ncbiRefSeq.txt $i; python2 ../fitter.py $i.dat > $i.dat.txt"
done
