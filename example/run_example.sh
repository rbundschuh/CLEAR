#!/bin/bash

# Example job script for CLEAR
# Python scripts must be present in $PATH

# Assumes that the bed files have already been generated and are in current
# working directory

#for i in *.bed;
#do
#  python2 ../make_dat.py ncbiRefSeq.txt $i # generates `dat` files
#  python2 ../fitter.py $i.dat > $i.dat.txt # finds passing transcripts
#done

python2 ../grouper.py --require-samples 3 *.dat.txt > CLEAR_passed.txt # generates passed gene list
python2 ../make_violin_plots.py # generates CLEAR_violins.pdf
