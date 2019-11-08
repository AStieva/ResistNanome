#!/usr/bin/python

# The file for filtering reads against a database using minimap2

import argparse, gzip, sys, os, multiprocessing, re, shutil, pysam
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
            if read.flag == 4:
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

    names = []
    for read in samfile.fetch():
        names.append(read.qname + "\t\t|\t\t" + samfile.get_reference_name(read.rname))
        if read.seq is not None:
            with open(filtered, "a") as handle:
                handle.write("@" + read.qname + "\n" + read.seq + "\n+\n" + read.qual + "\n")

    print("Only resistome data left")

    pdf = FPDF()
    pdf_out = os.path.join(args.outdir, "{}_resistome_info.pdf".format(prefix))
    pdf.add_page(orientation="P")
    pdf.set_font("Arial", style="B", size=10)
    pdf.cell(200, 10, txt="Resistome", ln=1, align="C")
    pdf.set_font("Arial", style = "U", size = 10)
    pdf.write(4, "Read ID\t\t|\t\tResistancy\n\n")
    for l in names:
        pdf.write(4, re.sub("resfinder~~~|~~~.+", "", l) + "\n")
    pdf.output(pdf_out, "F")

    samfile.close()