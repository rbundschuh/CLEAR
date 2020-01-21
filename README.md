# CLEAR: Coverage-based Limiting-cell Experiment Analysis for RNA-seq

## System Requirements
CLEAR is built to run on any system with Anaconda Python (see: https://www.anaconda.com/download/) properly installed.

The following packages should additionally be installed:
- `matplotlib`
- `numpy`

Processing samples from Human genome requires approximately 1GB of RAM, but this can vary based on the size of the reference genome and complexity of the coverage map.

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
  1. Run `python make_dat.py [Name of reference file] [Name of the .bed file from above]`, which will generate a file ending in `.dat` with the mu callings.
  2. Run `python fitter.py [Name of dat file generated in (1)]` which will print the passing transcripts to the terminal
    - This can be printed to a file using `python fitter.py [Name of dat file generated in (1)] > file_name`

## Grouping Transcript Lists
For each transcript name file produced above, run the `grouper.py` command as follows:
	`python grouper.py [name of file 1] [name of file 2] [name of...]`

You can also add the `--require-samples [#]` parameter, where `[#]` is the number of samples a transcript
must apper in to be included in the output. This can be used to relax the "passing in all samples" requirement
used in the manuscript.

## Wrapper script
See `wrapper.sh` for a complete wrapper for running the steps outlined above (`bash wrapper.sh`).

## Generating Visualization Violin Plots
1. Open the folder with all previously-generated `.dat` files.
2. Run `python make_violin_plots.py` to create a file `CLEAR_violins.pdf` containing violin plots of all samples in the folder.

## Example Case
An example case containing 6 cells' data from Zeisel et al. [1], allowing you to test your installation.
To run, simply run `cd example & bash run_example.sh`.
The expected results are in the `result` folder


[1] Zeisel A, Muñoz-Manchado AB, Codeluppi S, Lönnerberg P et al. Brain structure. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science 2015 Mar 6;347(6226):1138-42. PMID: 25700174x

