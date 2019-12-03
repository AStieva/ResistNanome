#!/usr/bin/python3

# Filtering the input so only the reads with the resistome are left

import argparse, gzip, sys, os, multiprocessing, re, shutil, pysam, csv, datetime, threading
from numpy import median
from fpdf import FPDF
from Bio import SeqIO
from PyPDF2 import PdfFileMerger, PdfFileReader
from __main__ import *

def resistome(inp, db, phred):
    print("Filtering for resistome")

    # Filtering with KMA
    if phred == 0:
        os.system("{} -i {} -o {} -t_db {} -bcNano".format(os.path.join(tool_dir, "KMA"), inp,
                                                           os.path.join(args.outdir, prefix + "_kma"), db))
    else:
        os.system("{} -i {} -o {} -t_db {} -bcNano -mp {}".format(os.path.join(tool_dir, "KMA"), inp,
                                                                  os.path.join(args.outdir, prefix + "_kma"), db,
                                                                  phred))

    # Creating lists with needed information, reads are in each, to link them
    resreads = [] # reads
    stome = [] # reads, gene, score
    info = [] # reads, gene, start, stop
    res = os.path.join(args.outdir, prefix + "_kma.frag.gz")
    with gzip.open(res, "rt") as csvfile:
        resread = csv.reader(csvfile, delimiter = "\t")
        for row in resread:
            resreads.append(row[6])
            stome.append("{}:{}:{}".format(row[6], row[5].rpartition("_")[0], row[2]))
            info.append("{}:{}:{}:{}".format(row[6], row[5], row[3], row[4]))
    reads = []
    for x in resreads:
        if x not in reads:
            reads.append(x)

    # Combining resistome and read lists to one; the best matching resistome for each read
    tome = []
    resmatch = []
    for r in reads:
        for t in stome:
            v = t.split(":")
            if re.search(r, t):
                tome.append(v[2])
            mval = max(tome)
            if mval == v[2]:
                resmatch.append("{}:{}".format(r, v[1]))
    Nlist = []
    for R in resmatch:
        T = R.split(":")
        for N in info:
            if re.search(T[1], N):
                Nlist.append(N)

    # Writing a fastq file with only resistome-files, taking out the resistome genes. output in fasta
    filtered = os.path.join(args.outdir, "temp_resistome.fastq")
    for record in SeqIO.parse(inp, "fastq"):
        filt = {record.id : record.seq}
    for n in Nlist:
        i = n.split(":")
        start = int(i[2]) - 1
        end = int(i[3]) - 1
        for key, val in filt:
            if i[0] == key:
                for x, y in enumerate(val):
                    if x == range(start, end):
                        val[x] = "N"
    for key, val in filt:
        with open(filtered, "a") as handle:
            SeqIO.write(val, key, handle, "fastq")

    # # Writing a fastq file with only resistome-files
    # filtered = os.path.join(args.outdir, "temp_resistome.fastq")
    # filt = []
    # for n in Nlist: # herschrijven
    #     N = n.split(":")
    #     for i in SeqIO.parse(inp, "fastq"):
    #         if N[0] == i.id:
    #             filt.append(i)
    # for f in filt:
    #     with open(filtered, "a") as handle:
    #         SeqIO.write(f, handle, "fastq")


    # Running Kraken2
    profiling("/mnt/docker/ResistNanome/db/Kraken2_db/", filtered, "resistome identification", "r")

    # Bracken addition
    med = int(med_round(os.path.join(args.outdir, "temp_resistome.fastq")))

    abnudance("/mnt/docker/ResistNanome/db/Kraken2_db/", "r", med, args.lvl)

    resis = []
    kra = "{}_rkreport.txt".format(os.path.join(args.outdir, prefix))
    with open(kra, "rt") as csvf:
        reader = csv.reader(csvf, delimiter = "\t")
        for row in reader:
            if row[0] == "C":
                for l in resmatch:
                    if re.search(row[1], l):
                        resis.append("{}:{}".format(l, row[2]))

    resinfo = []
    resdub = []
    if type(args.lvl) is list:
        bra = "{}_S_rkreport_bracken.txt".format(os.path.join(args.outdir, prefix))
    else:
        bra = "{}_rkreport_bracken.txt".format(os.path.join(args.outdir, prefix))
    with open(bra, "rt") as csvfi:
        reader = csv.reader(csvfi, delimiter = "\t")
        for row in reader:
            for t in resis:
                ts = t.split(":")
                idn = re.sub("\D", "", ts[2])
                if idn == row[4]:
                    resdub.append("{}:{}".format(t, row[0]))
    for dub in resdub:
        if dub not in resinfo:
            resinfo.append(dub)

    print("Writing resistome output")
    Rlog = os.path.join(args.outdir, "{}_resistome_output.txt".format(prefix))
    f = open(Rlog, "a")
    f.write("Read\t|\tName (taxonomy ID) Percentage bacteria of this type\nResistancy gene\n\n")
    for inf in resinfo:
        o = inf.split(":")
        resi = "{}\t|\t{} {}%\n{}\n".format(o[0], o[2], o[3], o[1])
        f.write(resi)
    f.close()