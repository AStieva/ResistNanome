#!/usr/bin/python3

# The big one, used for calling all the different parts and making sure everything runs smoothly

import argparse, gzip, sys, os, multiprocessing, re, shutil, pysam, csv, datetime, threading
from numpy import median
from Bio import SeqIO

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
parser.add_argument("--repN", "-N",
                    action="store_true",
                    default=False,
                    required=False,
                    help="Replace the resistancy gene with a series of N. This can take a long time")
parser.add_argument("--phred",
                    action="store",
                    default=0,
                    required=False,
                    help="The option of giving a minimum phred score for resistome filtering")
parser.add_argument("--taxonomy", "-cs",
                    action = "store_true",
                    default = False,
                    help = "Execute the bacterial community screening, output set with --lvl")
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
        
# If KMA is not called, but the repN is, it won't just run without KMA
if args.repN is True and args.resistome is False:
    print(
        "Called for replacing resistome genes, but resistome analysis isn't called. Continue without resistome analysis?")
    rep = input("(Y/N) ")
    if rep == "N" or rep == "n" or rep == "No" or rep == "no":
        print("To implement resistome analysis, please add --resistome or -ar to your arguments")
        sys.exit()
    else:
        print("Continuing without resistome analysis")

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
sys.path.append(os.path.join(tool_dir, "PyFPDF", "fpdf", "fpdf.py"))
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
    med = round(median(len_read))
    smed = str(med)
    if len(smed) >= 4 and med > 2549:
        step = -3
    elif len(smed) == 3 or med <= 2549:
        step = -2
    rounding = round(med, step)
    return int(rounding)

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
    unmapped(os.path.abspath(os.path.join(lib_dir, "/mash_db/")), indata)
    indata = os.path.join(args.outdir, "temp_novert.fastq")

# Running resistome analysis and/or community profiling, in multithreading
def resistome():
    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tStart determination of resistome (res)\n")
    f.close()
    from KMA import resistome
    resistome(indata, os.path.abspath(os.path.join(lib_dir, "KMA_ResFinder")), args.phred)
    resist_indata = os.path.join(args.outdir, "temp_resistome.fasta")
    if not os.path.isfile(resist_indata):
        print("There wasn't any antibiotic resistance found!")

def taxonomy():
    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tStart taxonomy (tax)\n")
    f.close()
    profiling(os.path.abspath(os.path.join(lib_dir, "/Kraken2_Nanodb/")), indata, "taxonomy", "t")
    med = int(med_round(indata))
    abundance(os.path.abspath(os.path.join(lib_dir, "/Kraken2_Nanodb/")), "t", med, "G")

    com = []
    kra = "{}_tkraken.txt".format(os.path.join(args.outdir, prefix))
    with open(kra, "rt") as csvf:
        reader = csv.reader(csvf, delimiter="\t")
        for read in reader:
            if read[0] == "C":
                com.append(read[1] + ":" + read[2])

    tax = []
    taxduo = []
    if type(args.lvl) is list:
        bra = "{}_S_tbracken.txt".format(os.path.join(args.outdir, prefix))
    else:
        bra = "{}_tbracken.txt".format(os.path.join(args.outdir, prefix))
    for b in com:
        tt = b.split(":")
        ID = tt[1].split(" ")
        lID = list(ID[-1])
        lID.pop()
        sID = "".join(lID)
        with open(bra, "rt") as csvfi:
            reader = csv.reader(csvfi, delimiter="\t")
            for col in reader:
                if col[1] == sID:
                    taxduo.append("{}:{}".format(b, col[-1]))
    for du in taxduo:
        if du not in tax:
            tax.append(du)

    tax.sort(reverse = True)
    
    print("Writing taxonomy output")
    toutput = os.path.join(args.outdir, "{}_taxonomy_output.txt".format(prefix))
    tt = open(toutput, "a")
    tt.write("Read\tName (tax ID)\tPercentage of bacteria\n")
    tt.close()
    taxinfo = []
    for e, inf in enumerate(tax):
        o = inf.split(":")
        perc = float(o[-1]) / 100
        peround = "{0:.6f}".format(perc)
        if e <= 9:
            taxinfo.append("{}\n\t\t{} - {}%\n".format(o[1], o[2], peround))
        tt = open(toutput, "a")
        tt.write("{}\t{}\t{}\n".format(o[1], o[2], peround))
        tt.close()            

    tlog = os.path.join(args.outdir, "temp_taxonomy_topout.txt")
    t = open(tlog, "a")
    t.write("Taxonomy")
    t.write("\nRead\n\t\tName (tax ID) - Percentage out of bacteria\n\n")
    for i in taxinfo:
        t.write(i)
    t.close()

    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tTaxonomy finished\n")
    f.close()

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

# Merging all output-files to one PDF
if not args.demux:
    outlist = sorted(os.listdir(args.outdir))
    Pdf = []
    for i in outlist:
        if re.search("output", i):
            ii = os.path.join(args.outdir, i)
            Pdf.append(ii)

    pdf = FPDF()
    out = os.path.join(args.outdir, "{}_ResistNanome.pdf".format(prefix))
    if not os.path.exists(out):
        open(out, "x").close()
    textlist = []
    try:
        Pdf[1]
    except IndexError:
        with open(Pdf[0], "r") as f:
            pdf.add_page(orientation="P")
            for row in f.readlines():
                if row == f.readlines[0]:
                    pdf.set_font("Arial", "BU", size=10)
                    pdf.write(4, row)
                elif re.search("Read|Name", row):
                    pdf.set_font("Arial", "B", size=10)
                    pdf.write(4, row)
                else:
                    pdf.set_font("Arial", size=10)
                    pdf.write(4, row)
        pdf.output(out, "F")
    else:
        for path in Pdf:
            p = open(path, "r")
            pdf.add_page(orientation="P")
            lines = p.readlines()
            for row in lines:
                if row == lines[0]:
                    pdf.set_font("Arial", "BU", size=10)
                    pdf.write(4, row)
                elif re.search("Read|Name", row):
                    pdf.set_font("Arial", "B", size=10)
                    pdf.write(4, row)
                else:
                    pdf.set_font("Arial", size=10)
                    pdf.write(4, row)
            p.close()
        pdf.output(out, "F")

    if not args.keep:
        for path in Pdf:
            os.remove(path)

f = open(RNlog, "a")
dt = datetime.datetime.now()
f.write(str(dt) + "\tStart cleaning output\n")
f.close()

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
                
