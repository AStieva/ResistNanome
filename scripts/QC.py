#!/usr/bin/python

# This is the file for the quality control

import argparse, gzip, sys, os, multiprocessing, re, shutil, pysam
from fpdf import FPDF
from Bio import SeqIO
from PyPDF2 import PdfFileMerger, PdfFileReader
from __main__ import *

def nanoqc(name, tool, inp):
    QCinput = os.path.normpath(os.path.join(args.outdir, "temp_{}_{}/".format(prefix, name)))
    if not os.path.exists(QCinput):
        os.mkdir(QCinput)

    QCollect = os.path.normpath(os.path.join(args.outdir, "QC/"))
    if not os.path.exists(QCollect):
        os.mkdir(QCollect)

    print("start {}".format(name))

    # NanoQC tool with a length of 400, so X = 200 per graph
    os.system("{} --outdir {} -l 400 {}".format(tool, QCinput, inp))

    os.rename(os.path.join(QCinput, "nanoQC.html"), os.path.join(QCollect, "{}_{}.html".format(prefix, name)))



def nanoqcimage(name, tool, inp):

    QCinput = os.path.normpath(os.path.join(args.outdir, "temp_{}_{}/".format(prefix, name)))
    if not os.path.exists(QCinput):
        os.mkdir(QCinput)

    print("start QC")

    # Just the nanoQC tool, it's output will be in pdf format.
    os.system("{} --outdir {} -f png {}".format(tool, QCinput, inp))

    # Putting the output in a PDF-file
    pdf_path = os.path.join(args.outdir, "{}_{}.pdf".format(prefix, name))
    images = sorted(os.listdir(QCinput), reverse = True)
    img = []
    for i in images:
        if ".png" in i:
            im = QCinput + i
            img.append(im)
    pdf = FPDF()

    for i in img:
        pdf.add_page(orientation = "P")
        pdf.set_font("Arial", size = 12)
        pdf.cell(200, 10, txt = "{}{}".format(os.path.basename(i).replace(".png", ""), name), ln = 1, align = "C")
        pdf.ln(2)
        pdf.image(i, type = "PNG", w = 210)
    pdf.output(pdf_path, "F")

    print("This QC is done and could already be looked at in the outputfolder. This will be part of the final output")
