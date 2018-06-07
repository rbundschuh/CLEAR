# CLEAR: Coverage-based Limiting-cell Experiment Analysis for RNA-seq

## System Requirements
CLEAR is built to run on any system with Anaconda Python (see: https://www.anaconda.com/download/) properly installed.

The following packages should additionally be installed:
- `matplotlib`
- `numpy`

## Installation
1. Ensure that python is properly installed and available on the system path.
2. Clone the CLEAR repository into a working folder for installation.
3. Retreive needed reference files in NCBI table format
  - Example download locations:
    - NCBI RefSeq GRCh38: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/ncbiRefSeq.txt.gz
    - NCBI RefSeq GRCm38: http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/ncbiRefSeq.txt.gz
4. Ensure that reference files are uncompressed.

## Creating BED Coverage Maps
1. Use an external tool such as `BedTools` to create a `.bed` file of coverages from each aligned sample to be interrogated. 

## Running CLEAR
For each `.bed` file, do the following:
  1. A
  2. B


## Generating Visualization Violin Plots
1. Open the folder with all previously-generated `.dat` files.
2. Run `python make_violin_plots.py` to create a file `violins.py` containing violin plots of all samples in the folder.
