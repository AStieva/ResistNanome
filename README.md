# ResistNanome
A pipeline-tool to easily extract the resistome information of a metagenomic sample sequenced using the Nanopore technique.

This pipeline should be reasonably fast in giving the antibiotic resistance genes present in the metagenomic Nanopore dataset, combined with the bacteria. 
Getting this data is the main objective, but a multitude of tools are build in to improve on the output.

-------------------------------------------

## Installation

ResitNanome contains ost of what is needed. Only bokeh is to be pre-installed for the use of nanoQC.

To install Resistnanome through Unix:

```
git clone https://github.com/AStieva/ResistNanome
```

To include the databases:
```
cd ResistNanome/database
./getDB.sh
```
To delete what's not needed after downloading the databases
```
rm *.tar.gz
```
If you want to use nanoQC, please install bokeh like this: `pip3 install bokeh`

## Usage

When calling the program, two arguments are mandatory, so should always be given: `--inp` or `-i` and `--outdir` or `-o`. After putting down these arguments, the (path to the) input data and the path to the output folder should be put in, respectively. For example when only calling the resistome  analysis: 

`./Resistnanomewrapper.py --inp file.fastq.gz --outdir /path/to/output/folder/ --resistome`

The output from this will be a pdf file which gives the top 10 reads and a tab seperated values file with all information. However, the output doesn't give any information about the quality of the reads and doesn't filter out anything.

## Options

```
usage: ./wrapper.py --inp file.fastq.gz --outdir /path/to/output/folder/ [options]

-h, --help        Show this help message and exit
-i, --inp         Set the input files, in fastq or fastq.gz format
-o, --outdir      Set the output directory, please make sure to input the whole path
-p, --prefix      Set the prefix for all output files, default is the name of the input files. This name can't contain 'temp_' as this is used to define temporary data
-t, --threads     Set the number of threads used, default is maximum available up to 16 threads
--keep            If called, all the intermediate data is also kept. Otherwise only outputs with new information wil be kept
--demux           Execute demultiplexing with Porechop, the verbose will be saved to a pdf file
-fl, --filtlong   Execute filtlong, add --minlen [bp] to add a custom minimum length of reads, default for this is 1000
--minlen          The option to put in the minimum length of the reads saved in bp, default is 1000. Needs filtlong to be used
-qc, --QC         Execute the quality control
--host            Execute the host contamination screening, the database has a lot of vertebrate, but if yours isn't in it, please feel free to add it
--lvl             The classification level for determining the abundance levels with Bracken. Default is all. Options are: K (kingdom level), P (phylum), C (class), O (order), F (family), G (genus), and S (species)
-ar, --resistome  Execute the resistome analysis/antibiotic resistance screening
-N, --repN        Replace the resistancy gene with a series of N. This can take a long time (up to multiple days for a few GB data)
--phred           The option of giving a minimum phred score for resistome filtering
-cs, --taxonomy   Execute the bacterial community screening
-gz, --gz         gzip-ing the fastq file(s) left at the end
```

## Implemented tools

### NanoQC
Wouter De Coster, Svenn D’Hert, Darrin T Schultz, Marc Cruts, Christine Van Broeckhoven, "NanoPack: visualizing and processing long-read sequencing data", Bioinformatics, Volume 34, Issue 15, 01 August 2018, Pages 2666–2669, [paper](https://doi.org/10.1093/bioinformatics/bty149)

Tool for creating Quality Controll graphs.

[GitHub](https://github.com/wdecoster/nanoQC)

>#### bokeh
 >[Information page](https://docs.bokeh.org/en/latest/index.html)

 >Nanoqc uses bokeh for the visualisation

 >[GitHub](https://github.com/bokeh/bokeh)

### Filtlong
Ryan Wick filtlong [GitHub](https://github.com/rrwick/Filtlong)

Tool for filtering Nanopore data for length and quality.

### Porechop
Ryan R. Wick, Louise M. Judd, Claire L. Gorrie, Kathryn E. Holt, "Completing bacterial genome assemblies with multiplex MinION sequencing", MICROBIAL GENOMICS Volume 3, Issue 10, First Published: 14 September 2017 [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000132)

Tool for demultiplexing and cutting off barcodes.

[GitHub](https://github.com/rrwick/Porechop)

### Minimap2
Heng Li, "Minimap2: pairwise alignment for nucleotide sequences", Bioinformatics, Volume 34, Issue 18, 15 September 2018, Pages 3094–3100, [paper](https://doi.org/10.1093/bioinformatics/bty191)

Tool for alignment.

[GitHub](https://github.com/lh3/minimap2)

### Kraken2
Derrick E. Wood, Jennifer Lu, Ben Langmead, "Improved metagenomic analysis with Kraken 2", bioRxiv 762302 [paper](https://doi.org/10.1101/762302)
 
Tool for community screening.
 
[GitHub](https://github.com/DerrickWood/kraken2)
 
### Bracken
Jennifer Lu, Florian P. Breitwieser, Peter Thielen, Steven L. Salzberg, "Bracken: estimating species abundance in metagenomics data", Article in Computer Science, published January 2017, [paper](https://peerj.com/articles/cs-104/)

Tool for giving an abundance estimate based of a kraken(2) report.

[GitHub](https://github.com/jenniferlu717/Bracken)

### PyFPDF
Olivier Plathey [FPDF](http://www.fpdf.org/) Original FPDF for PHP
Max Pat, Mariano Reingart, Roman Kharin [readthedocs](https://pyfpdf.readthedocs.io/en/latest/) Python version of FPDF. Original python version by Max Pat, forked version used.

Tool for making PDF-files with code. Used to write output directly to PDF.

[GitHub FPDF](https://github.com/Setasign/FPDF)

[GitHub Python](https://github.com/reingart/pyfpdf)

### pysam
Andreas Heger, Kevin Jacobs et al. 2009 [readthedocs](https://pysam.readthedocs.io/en/latest/#)

Tool for manipulating SAM files.

[GitHub](https://github.com/pysam-developers/pysam)

### KMA
Philip T.L.C. Clausen, Frank M. Aarestrup & Ole Lund, "Rapid and precise alignment of raw reads against redundant databases with KMA", BMC Bioinformatics, 2018;19:307, [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2336-6)

Tool for k-mer alignment. Often used for determenation of the resistome.

[Bitbucket](https://bitbucket.org/genomicepidemiology/kma/src/master/)

-------------------------------------------
##### This page is still under construction
