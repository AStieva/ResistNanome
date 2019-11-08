#!/usr/bin/python

# This is the file for Demultiplexing data

import argparse, gzip, sys, os, multiprocessing, re, shutil, pysam
from fpdf import FPDF
from Bio import SeqIO
from PyPDF2 import PdfFileMerger, PdfFileReader
from __main__ import *

def demux(tool, inp):
    # Make a directory for the demultiplex output
    deminput = os.path.normpath(os.path.join(args.outdir, "{}_demulti/".format(prefix)))
    if not os.path.exists(deminput):
        os.mkdir(deminput)

    # To run Porechop, the correct tool needs to be found, the input has to be put in, format makes sure the output is the correct format and output leads to the directory the output should go
    print("Running Porechop")

    # Running porechop and writing the verbose output to a text file, to be replaced by a pdf file without loading information
    log = os.path.join(args.outdir, "temp_demultiplex_info.txt")
    if not os.path.exists(log):
        open(log, "x").close()

    os.system("{} -i {} --format fastq -b {} > {}".format(tool, inp, deminput, log))

    file = []
    with open(log, 'r') as l:
        for line in l:
            if not re.search(".+%\)", line):
                file.append(line)
    pdf = FPDF()
    pdf_out = os.path.join(args.outdir, "{}_demultiplex_info.pdf".format(prefix))
    pdf.add_page(orientation="P")
    pdf.set_font("Arial", size=10)
    pdf.cell(200, 10, txt="demultiplex", ln=1, align="C")

    for i in file:
        pdf.set_font("Arial", size=10)
        if "[4m" in i:
            pdf.set_font("Arial", style="U", size=10)
            pdf.write(3, i + "\n")
        elif "[32m" in i:
            pdf.set_text_color(255, 0, 0)
            pdf.write(3, i + "\n")
        elif "[31m" in i:
            pdf.set_text_color(0, 204, 0)
            pdf.write(4, i + "\n")
        else:
            pdf.set_text_color(1)
            pdf.write(3, i + "\n")
    pdf.output(pdf_out, "F")

    print("Done with running Porechop")
