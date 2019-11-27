#!/usr/bin/python3

# The big one, used for calling all the different parts and making sure everything runs smoothly

import argparse, gzip, sys, os, multiprocessing, re, shutil, pysam, csv, datetime, threading
from numpy import median
from Bio import SeqIO
from PyPDF2 import PdfFileMerger, PdfFileReader

# Setting a default maximum of threads, even if maximum is higher
def cpu_threads(max_threads):
  if multiprocessing.cpu_count() > max_threads:
    return max_threads
  else:
    return multiprocessing.cpu_count()

# Arguments to be put in by user. They may not be in the correct order, some may need to be split in different options
parser = argparse.ArgumentParser(description = "A pipeline to make it easy to extract the resistome (and other information) from a set of metagenomic Nanopore data")
parser.add_argument("--inp", "-i",
                    action = "store",
                    required = True,
                    help = "Place here the input files, in fastq or fastq.gz format")
parser.add_argument("--outdir", "-o",
                    action = "store",
                    required = True,
                    help = "Set the output directory, please make sure to input the whole path")
parser.add_argument("--prefix", "-p",
                    action = "store",
                    help = "Set the prefix for all output files, default is the name of the input files. This name can't start with 'temp_' as this is used to define temporary data")
parser.add_argument("--threads", "-t",
                    default = cpu_threads(16),
                    required = False,
                    type = int,
                    help = "Set the number of threads used, default is maximum available up to 16 threads")
parser.add_argument("--keep",
                    action = "store_true",
                    default = False,
                    help = "If called, all the intermediate data is also kept. Otherwise only outputs with new information wil be kept")
parser.add_argument("--demux",
                    action = "store_true",
                    default = False,
                    help = "Execute demultiplexing with Porechop, the verbose will be saved to a pdf file")
parser.add_argument("--filtlong", "-fl",
                    action = "store_true",
                    default = False,
                    required = False,
                    help = "Execute filtlong, add --minlen [bp] to add a custom minimum length of reads, default for this is 1000")
parser.add_argument("--minlen",
                    action = "store",
                    default = 1000,
                    required = False,
                    help = "The option to put in the minimum length of the reads saved in bp, default is 1000. Needs filtlong to be used")
parser.add_argument("--QC", "-qc",
                    action = "store_true",
                    default = False,
                    help = "Execute the quality control")
parser.add_argument("--host",
                    action = "store_true",
                    default = False,
                    help = "Execute the host contamination screening, the database has a lot of vertebrate, but if yours isn't in it, please feel free to add it")
parser.add_argument("--lvl",
                    action = "store",
                    default = ["S", "G", "F", "O", "C", "P", "K"],
                    required = False,
                    help = "The classification level for determining the abundance levels with Bracken. Default is all."
                           "Options are: K (kingdom level), P (phylum), C (class), O (order), F (family), G (genus), and S (species)")
parser.add_argument("--resistome", "-ar",
                    action = "store_true",
                    default = False,
                    help = "Execute the resistome analysis")
parser.add_argument("--taxonomy", "-cs",
                    action = "store_true",
                    default = False,
                    help = "Execute the bacterial community screening, output on genome level")
parser.add_argument("--gz", "-gz",
                    action = "store_true",
                    default = False,
                    help = "gzip-ing the last input file")
args = parser.parse_args()

# If filtlong is not called, but the minlen is, it won't just run without filtlong
if args.minlen is not 1000 and args.filtlong == False:
    print("A minimum length has been passed, but the filtlong step is not called. Would you like to continue without filtering?")
    filt = input("(Y/N) ")
    if filt == "N" or filt == "n" or filt == "No" or filt == "no":
        print("To implement filtlong, please add --filtlong or -fl to your arguments")
        sys.exit()
    else:
        print("Continuing without filtlong")

# Demux disclamer
if args.demux and (args.filtlong or args.host or args.resistome or args.taxonomy):
    print("You've selected demultilexing, only this and possibly QC will be ran, the output will be zipped. Please run the demultiplexed output files seperately")
    print("Would you like to continue?")
    dem = input("(Y/N) ")
    if dem == "N" or dem == "n" or dem == "No" or dem == "no":
        print("Please take the --demux argument out to run without demultiplexing")
        sys.exit()
    else:
        print("Continue with demultiplexing")

# Make sure the directory path exists, if not, create it
if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)

# Setting the input as indata
indata = args.inp

# The prefix will be specified here, so it doesn't has to happen later
pre = os.path.basename(args.inp)
prefix = re.sub("\.gz|\.fastq|\.fq", "", pre)
if args.prefix is not None:
    prefix = args.prefix

# Creating a log of the pipeline
RNlog = os.path.join(args.outdir, "{}_ResistNanome-info.log".format(prefix))
f = open(RNlog, "a")
dt = datetime.datetime.now()
f.write(str(dt) + "\tPipeline started\n")
f.close()

# Creating a 'shortcut' to directories
directory = os.path.abspath(os.path.dirname(sys.argv[0]))
tool_dir = os.path.join(directory, "tools/")
lib_dir = os.path.join(directory, "database/")
script_dir = os.path.join(directory, "scripts/")

# Making sure the scripts are found
sys.path.append(script_dir)
sys.path.append(os.path.join(tool_dir, "pyfpdf"))
from fpdf import FPDF

# Function for running kraken2 and bracken
def profiling(db, inp, part, id):
    print("Running Kraken2 for {}".format(part))
    os.system("{} --db {} --report {}_{}kreport.txt --threads {} --use-names --output {}_{}kraken.txt {}".format(
        os.path.join(tool_dir, "Kraken/kraken2"), db, os.path.join(args.outdir, prefix), id, args.threads,
        os.path.join(args.outdir, prefix), id, inp))
def abundance(db, id, len, level):
    if isinstance(level, list):
        for lvl in level:
            os.system(
                "{} -d {} -i {}_{}kreport.txt -o {}_{}bracken.txt -r {} -l {}".format(os.path.join(tool_dir, "bracken"),
                                                                                      db,
                                                                                      os.path.join(args.outdir, prefix),
                                                                                      lvl + "_" + id,
                                                                                      os.path.join(args.outdir, prefix),
                                                                                      id, len, lvl))
    else:
        os.system(
            "{} -d {} -i {}_{}kreport.txt -o {}_{}bracken.txt -r {} -l {}".format(os.path.join(tool_dir, "bracken"), db,
                                                                                  os.path.join(args.outdir, prefix), id,
                                                                                  os.path.join(args.outdir, prefix), id,
                                                                                  len, level))
def med_round(data):
    len_read = []
    for seq in SeqIO.parse(data, "fastq"):
        len_read.append(int(len(seq)))
    med = int(median(len_read))
    smed = str(med)
    if len(smed) >= 4 and med > 2549:
        x = -3
    elif len(smed) == 3 or med <= 2549:
        x = -2
    return int(round(med, x))

# Running first QC, if asked for QC
if args.QC:
    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tStart NanoQC1\n")
    f.close()
    from QC import nanoqc
    nanoqc("QC1", "./tools/nanoQC/nanoQC.py", indata)
    print("Done, the html report can be found in the map 'QC' in the outdir")

# Running demultiplexing, if asked for
deminput = os.path.normpath(os.path.join(args.outdir, "{}_demulti/".format(prefix)))
if args.demux:
    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tStart Demultiplexing\n")
    f.close()
    from Demultiplexing import demux
    demux("./tools/porechop-runner.py", indata)
    indata = deminput

#Running Filtlong if asked for. Possibly with different min-length
if args.filtlong and not args.demux:
    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tStart filtlong\n")
    f.close()
    from Filtlong import filtlong
    filtlong(os.path.join(tool_dir, "filtlong"), indata, args.minlen)
    indata = os.path.join(args.outdir, "temp_filtlong.fastq")

# Running second QC, if asked for QC and asked for demultiplexing or filtlong, for comparison
if args.QC and (args.demux or args.filtlong):
    from QC import nanoqc
    if os.path.isdir(indata):
        lst = sorted(os.listdir(indata))
        qclst = []
        for fl in lst:
            fl = os.path.join(indata, fl)
            qclst.append(fl)
        for qc in qclst:
            f = open(RNlog, "a")
            dt = datetime.datetime.now()
            f.write(str(dt) + "\tStart " + re.sub("\.fastq|\.gz", "", os.path.basename(qc)) + "\n")
            f.close()
            nanoqc("QC_" + re.sub("\.fastq|\.gz", "", os.path.basename(qc)), "./tools/nanoQC/nanoQC.py", qc)
        print("Done, the html reports can be found in the map 'QC' in the outdir")
    else:
        f = open(RNlog, "a")
        dt = datetime.datetime.now()
        f.write(str(dt) + "\tStart NanoQC2\n")
        f.close()
        nanoqc("QC2", "./tools/nanoQC/nanoQC.py", indata)
        print("Done, the html report can be found in the map 'QC' in the outdir")

# Running screening for/extracting of host DNA, if asked for
if args.host and not args.demux:
    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tStart host filtering\n")
    f.close()
    from Minifilter import unmapped
    unmapped(os.path.abspath("/mnt/db/db/mash_db/"), indata)
    indata = os.path.join(args.outdir, "temp_novert.fastq")

# Running resistome analysis and/or community profiling, in multithreading
def resistome():
    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tStart determination of resistome\n")
    f.close()
    from Minifilter import resist
    resist(os.path.abspath("/mnt/db/db/resistome/resfinder.fasta"), indata)
    resist_indata = os.path.join(args.outdir, "temp_resistome.fastq")
    if not os.path.isfile(resist_indata):
        print("There was'n any antibiotic resistance found!")
def taxonomy():
    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tStart taxonomy\n")
    f.close()
    profiling("/mnt/docker/ResistNanome/db/Kraken2_db/", indata, "taxonomy", "t")
    med = int(med_round(indata))
    abundance("/mnt/docker/ResistNanome/db/Kraken2_db/", "t", med, "G")

if args.resistome and args.taxonomy and not args.demux:
    # Creating threads
    t1 = threading.Thread(target = resistome())
    t2 = threading.Thread(target = taxonomy())
    # Starting threads
    t1.start()
    t2.start()
    # Nothing else will start 'till both are done
    t1.join()
    t2.join()
elif args.resistome and not args.demux and not args.taxonomy:
    resistome()
elif args.taxonomy and not args.demux and not args.resistome:
    taxonomy()

f = open(RNlog, "a")
dt = datetime.datetime.now()
f.write(str(dt) + "\tStart cleaning output\n")
f.close()

# Merging all pdf-files or changing the name if there's only one
if not args.demux:
    outlist = sorted(os.listdir(args.outdir))
    Pdf = []
    for i in outlist:
        if re.search(".+\.pdf", i) and not re.search(".+_ResistNanome\.pdf", i):
            ii = os.path.join(args.outdir, i)
            Pdf.append(ii)
    if Pdf:
        out = os.path.join(args.outdir, "{}_ResistNanome.pdf".format(prefix))
        try:
            Pdf[1]
        except IndexError:
            os.rename(Pdf[0], out)
        else:
            if not os.path.exists(out):
                open(out, "x").close()
            merge = PdfFileMerger()
            for path in Pdf:
                merge.append(PdfFileReader(path))
            merge.append(out)
            merge.close()
#            for path in Pdf:
#                os.remove(path)

# The temporary output directory can be deleted, if keep argument is called, this argument will be false, so all data will be kept
if not args.keep:
    outlist = sorted(os.listdir(args.outdir))
    rmlist = []
    for i in outlist:
        if re.search("temp_.+", i):
            ii = os.path.join(args.outdir, i)
            rmlist.append(ii)
    for r in rmlist:
        if os.path.isdir(r):
            rm = sorted(os.listdir(r))
            for f in rm:
                os.remove(os.path.join(r, f))
            os.rmdir(r)
        else:
            os.remove(r)

# Zipping fastq files to .gz
fq_files = []
if args.gz:
    print("Zipping files")
    outlist = sorted(os.listdir(args.outdir))
    for thing in outlist:
        if re.search(".+\.fastq$", str(thing)) or re.search(".+\.fq$", str(thing)):
            thing_path = os.path.join(args.outdir, thing)
            fq_files.append(thing_path)
# Zipping demux files to .gz, this will be done when demux is called
elif args.demux:
    outdir = sorted(os.listdir(deminput))
    for bc in outdir:
        if re.search(".+\.fastq", str(bc)) or re.search(".+\.fq", str(bc)):
            bc_path = os.path.join(deminput, bc)
            fq_files.append(bc_path)
# The actual zipping part
if args.gz or args.demux:
    for fq in fq_files:
        outp =  fq + ".gz"
        with open(fq, "rb") as full:
            with gzip.open(outp, "wb")  as output:
                shutil.copyfileobj(full, output)
        os.remove(fq)
        
f = open(RNlog, "a")
dt = datetime.datetime.now()
f.write(str(dt) + "\tFinished!\n")
f.close()
                
