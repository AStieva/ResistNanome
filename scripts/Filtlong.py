#!/usr/bin/python

# This is the file for filtering and cutting long data

import argparse, gzip, sys, os, multiprocessing, re, shutil, pysam
from fpdf import FPDF
from Bio import SeqIO
from PyPDF2 import PdfFileMerger, PdfFileReader
from __main__ import *

def filtlong(tool, inp, len):
    print("Running Filtlong")

    # The min length could be changed, or is 1000. The min mean q stands for the minimum quality, the 1 percentile worst quality reads are cut out.
    os.system("{} --min_length {} --min_mean_q 1 {} > {}".format(tool, len, inp, os.path.join(args.outdir, "temp_{}_filtlong.fastq".format(prefix))))

    print("Filtlong finished")
