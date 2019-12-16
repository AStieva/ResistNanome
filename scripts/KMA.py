#!/usr/bin/python3

# Filtering the input so only the reads with the resistome are left

import argparse, gzip, sys, os, multiprocessing, re, shutil, pysam, csv, datetime, threading
from numpy import median
from fpdf import FPDF
from Bio import SeqIO
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
    resreads = []
    info = []  # reads, score, gene, start, stop
    res = os.path.join(args.outdir, prefix + "_kma.frag.gz")
    with gzip.open(res, "rt") as csvfile:
        resread = csv.reader(csvfile, delimiter="\t")
        for row in resread:
            resreads.append(row[6])
            info.append("{}:{}:{}:{}:{}".format(row[6], row[2], row[5].rpartition("_")[0], row[3], row[4]))
    reads = []
    for Rid in resreads:
        if Rid not in reads:
            reads.append(Rid)

    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\t\t(res) Creating KMA-output\n")
    f.close()

    # Combining resistome and read lists to one; the best matching resistome for each read
    tome = []
    stome = {}
    resmatch = {}
    match = []
    for rid in reads:
        for t in info:
            u = t.split(":")
            if re.search(rid, t):
                tome.append(u[1])
        stome[rid] = max(tome)
        tome.clear()
    for k, v in stome.items():
        for t in info:
            u = t.split(":")
            if re.search(v, t) and k == u[0]:
                resmatch[u[0]] = u[2:]
    for k, v in resmatch.items():
        vv = ":".join(v)
        match.append(k + ":" + vv)

    # Writing a fastq file with only resistome-files, possibly taking out the resistant genes
    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\t\t(res) Writing fasta output file\n")
    f.close()
    matchseq = []
    for t in match:
        for record in SeqIO.parse(inp, "fastq"):
            if re.search(record.id, t):
                matchseq.append(t + ":" + str(record.seq))
    if args.repN:
        # Replacing the bases of the resistome gene with 'N'
        filtered = os.path.join(args.outdir, "temp_resistome.fasta")
        filt = []
        for n in matchseq:
            i = n.split(":")
            start = int(i[2])
            stop = int(i[3])
            seql = list(i[4])
            for x, y in enumerate(seql):
                if start <= x < stop:
                    seql[x] = "N"
                    valtoo = "".join(seql)
                    filt.append(valtoo + ":" + i[0])
        with open(filtered, "a") as handle:
            for key in filt:
                val = key.split(":")
                handle.write(">{}\n{}\n".format(val[1], val[0]))
    else:
        # Writing a fastq file with only resistome-files
        filtered = os.path.join(args.outdir, "temp_resistome.fastq")
        filt = []
        for n in matchseq:
            N = n.split(":")
            for i in SeqIO.parse(inp, "fastq"):
                if N[0] == i.id:
                    filt.append(i)
        for f in filt:
            with open(filtered, "a") as handle:
                SeqIO.write(f, handle, "fastq")


    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\t\t(res) Taxonomisation on KMA-output\n")
    f.close()

    # Running Kraken2
    profiling("/mnt/docker/ResistNanome/db/Kraken2_db/", filtered, "resistome identification", "r")

    # Bracken addition
    med = med_round(inp)

    abundance("/mnt/docker/ResistNanome/db/Kraken2_db/", "r", med, args.lvl)

    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\t\t(res) Creating information output\n")
    f.close()

    clas = []
    kra = "{}_rkraken.txt".format(os.path.join(args.outdir, prefix))
    with open(kra, "rt") as csvf:
        reader = csv.reader(csvf, delimiter="\t")
        for row in reader:
            if row[0] == "C":
                clas.append(row[1] + ":" + row[2])
    resis = []
    for l in match:
        lm = l.split(":")
        for c in clas:
            row = c.split(":")
            if re.search(row[0], l):
                resis.append("{}:{}:{}".format(lm[0], lm[1], row[1]))

    resinfo = []
    resduo = []
    if type(args.lvl) is list:
        bra = "{}_S_rbracken.txt".format(os.path.join(args.outdir, prefix))
    else:
        bra = "{}_rbracken.txt".format(os.path.join(args.outdir, prefix))
    for t in resis:
        tt = t.split(":")
        ID = tt[2].split(" ")
        lID = list(ID[-1])
        lID.pop()
        sID = "".join(lID)
        with open(bra, "rt") as csvfi:
            reader = csv.reader(csvfi, delimiter="\t")
            for row in reader:
                if row[1] == sID:
                    resduo.append("{}:{}".format(t, row[-1]))
    for du in resduo:
        if du not in resinfo:
            resinfo.append(du)

    print("Writing resistome output")
    resi = []
    for inf in resinfo:
        o = inf.split(":")
        perc = float(o[-1]) / 100
        peround = "{0:.6f}".format(perc)
        resi.append("{}\t|\t{}\n\t\t{} - {}%\n".format(o[0], o[1], o[2], peround))

    Rlog = os.path.join(args.outdir, "{}_resistome_output.txt".format(prefix))
    t = open(Rlog, "a")
    t.write("Resistome")
    t.write("\nRead\t|\tResistancy gene\n\t\tName (tax ID) - Percentage out of bacteria\n\n")
    for r in resi:
        t.write(r)
    t.close()

    f = open(RNlog, "a")
    dt = datetime.datetime.now()
    f.write(str(dt) + "\tResistome finished\n")
    f.close()
