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

