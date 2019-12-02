#!/usr/bin/python3

# The file for filtering reads against a database using minimap2

import argparse, gzip, sys, os, multiprocessing, re, shutil, pysam, csv
from numpy import median
from fpdf import FPDF
from Bio import SeqIO
from PyPDF2 import PdfFileMerger, PdfFileReader
from __main__ import *

def unmapped(ref, inp):
    print("Taking out reads with vertebrate DNA")

    sam = os.path.join(args.outdir, "temp_novert.sam")    
    filtered = os.path.join(args.outdir, "temp_novert.fastq")

    # Using minimap2 to align reads to vertebrate db
    # Using pysam to filter and translate to fastq
    os.system("{} -ax map-ont {} {} | grep -v '^@' > {}".format(os.path.join(tool_dir, "Minimap2"), ref, inp, sam))

    with pysam.AlignmentFile(sam, 'r', check_sq = False) as samfile:
        for read in samfile.fetch():
            if read.flag == 4 and read.qname is not None and read.seq is not None and read.qual is not None:
                with open(filtered, "a") as handle:
                    handle.write("@" + read.qname + "\n" + read.seq +"\n+\n" + read.qual + "\n")

    print("Host DNA is filtered out")


def resist(db, inp):
    print("Filtering for resistome")

    # Using minimap2 to filter anything resistome out
    sam = os.path.join(args.outdir, "temp_resistome.sam")
    filtered = os.path.join(args.outdir, "temp_resistome.fastq")
    os.system("{} -ax map-ont {} --sam-hit-only {} > {}".format(os.path.join(tool_dir, "Minimap2"), db, inp, sam))
       
    samfile = pysam.AlignmentFile(sam, 'rb', check_sq=False)

    # Making a list for writing output and translate sam to fastq
    unid = []
    res = []
    for read in samfile.fetch():
        if read.seq is not None:
            with open(filtered, "a") as handle:
                handle.write("@" + read.qname + "\n" + read.seq + "\n+\n" + read.qual + "\n")
            unid.append(read.qname)
        res.append("{}:{}".format(read.qname, samfile.get_reference_name(read.rname)))
    
    print("Only resistome data left")

    # Using the kraken2 formula from wrapper for taxonomisation
    profiling("/mnt/docker/ResistNanome/db/Kraken2_db/", filtered, "resistome identification", "r")

    # Addition of bracken
    med = int(med_round(os.path.join(args.outdir, "temp_resistome.fastq")))
        
    abundance("/mnt/docker/ResistNanome/db/Kraken2_db/", "r", med, args.lvl)

    # Update the output-list
    tres = []
    name = {}
    for i in unid:
        for r in res:
            v = r.split(":")
            if re.search(i, r):
                tres.append(re.sub("resfinder~~~|~~~.+", "", v[1]) + "\n")
        RES = "".join(tres)
        name[i] = RES
        tres.clear()
 
    # Make the output look pretty 
    names = {}
    kra = "{}_rkraken.txt".format(os.path.join(args.outdir, prefix))
    with open(kra) as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row[0] == "C":
                for key, value in name.items():
                    if row[1] == key:
                        names[row[2] + "\t\t|\t\t" + key] = value

    # Writing the output list to a pdf file
    pdf = FPDF()
    pdf_out = os.path.join(args.outdir, "{}_resistome_info.pdf".format(prefix))
    pdf.add_page(orientation="P")
    pdf.set_font("Arial", style="B", size=10)
    pdf.cell(200, 10, txt="Resistome", ln=1, align="C")
    pdf.set_font("Arial", size = 10)
    pdf.write(4, "Bacteria (taxID)\t\t|\t\tRead ID\nResistancy\n\n")
    for key, value in names.items():
        pdf.write(4, key + "\n" + value)
    pdf.output(pdf_out, "F")

    samfile.close()
