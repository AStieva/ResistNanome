# ResistNanome
A pipeline-tool to easily extract the resistome information of a metagenomic sample sequenced using the Nanopore technique.

This pipeline should be reasonably fast in giving the antibiotic resistance genes present in the metagenomic Nanopore dataset, combined with the bacteria. 
Getting this data is the main objective, but a multitude of tools are build in to improve on the output.

## Usage

When calling the program, two arguments are mandatory, so should always be given: `--inp` or `-i` and `--outdir` or `-o`. After putting down these arguments, the (path to the) input data and the path to the output folder should be put in, respectively. for example when only calling the resistome  analasys: 

`./Resistnanomewrapper.py --inp file.fastq.gz --outdir /path/to/output/folder/ --resistome`

The output from this will be a pdf file which gives the name of the bacteria, which read and the antibiotic resistance genes. However, the output doesn't give any information about the quality of the reads and doesn't filter out anything.

## Options

```
usage: ./wrapper.py --inp file.fastq.gz --outdir /path/to/output/folder/ [options]

-h, --help        Show this help message and exit
-i, --inp         Place here the input files, in fastq or fastq.gz format
-o, --outdir      Set the output directory, please make sure to input the whole path
-p, --prefix      Set the prefix for all output files, default is the name of the input files. This name can't contain 'temp_' as this is used to define temporary data
-t, --threads     Set the number of threads used, default is maximum available up to 16 threads
--keep            If called, all the intermediate data is also kept. Otherwise only outputs with new information wil be kept
--demux           Execute demultiplexing with Porechop, the verbose will be saved to a pdf file
-fl, --filtlong   Execute filtlong, add --minlen [bp] to add a custom minimum length of reads, default for this is 1000
--minlen          The option to put in the minimum length of the reads saved in bp, default is 1000. Needs filtlong to be used
-qc, --QC         Execute the quality control
--host            Execute the host contamination screening, the database has a lot of vertebrate, but if yours isn't in it, please feel free to add it
-ar, --resistome  Execute the resistome analysis/antibiotic resistance screening
-cs, --taxonomy   Execute the bacterial community screening
--gz              gzip-ing the fastq file(s) left at the end

## Implemented tools



##### This page is still under construction
