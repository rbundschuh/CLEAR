#!/bin/bash

# Example job script for CLEAR
# Python scripts must be present in $PATH

# Assumes that the bed files have already been generated and are in current
# working directory

for i in *.bed
do
  python make_dat.py ncbiRefSeq.txt $i # generates `dat` files
  python fitter.py $i.dat > $i.dat.txt # finds passing transcripts
done

python grouper.py *.dat.txt > CLEAR_passed.txt # generates passed gene list
python make_violin_plots.py # generates CLEAR_violins.pdf