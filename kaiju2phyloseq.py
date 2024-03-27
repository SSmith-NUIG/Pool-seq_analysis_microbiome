#!/usr/bin/env python
# Made using Python 3.7
# V1.00 Written by Gisle Vestergaard (gislevestergaard@gmail.com)
# This takes a directory with Kaiju (http://kaiju.binf.ku.dk/) 
# output files and creates Phyloseq otu and tax tables to be imported 
# in R for statistical analysis. This uses the program kaiju-addTaxonNames
# which is part of the kaiju installation to add taxa names to each OTU.
# The provided names and nodes dmp has to be those used by Kaiju. and 
# kaiju-addTaxonNames has to be in your PATH.

import argparse
import Bio
import glob
import pandas
import subprocess
import sys
import os

parser = argparse.ArgumentParser(description='Prints fasta entries lenght')
parser.add_argument("-i", "--indir", required = True, 
					help="Input directory containing the wanted *.kaiju files")
parser.add_argument("-n", "--names", required = True, 
					help="Name of Kaiju names.dmp file")
parser.add_argument("-m", "--nodes", required = True, 
					help="Name of kaiju nodes.dmp file")
parser.add_argument("-t", "--threads", default=1, required = False, 
					help="Input amount of available cores")
parser.add_argument("-o", "--output", required = True, 
					help="Basename for output files")
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

# First, we will add taxonomical names to the kaiju output using kaiju-addTaxonNames
kaijus = []
kaijus = glob.glob(args.indir + "/*.out")
for file in kaijus:
	output = file + ".names"
	subprocess.call(['kaiju-addTaxonNames', '-n', args.names, '-t', args.nodes, '-r', 'superkingdom,phylum,class,order,family,genus,species'
, '-i', file, '-o', output])

# Secondly, we can now read all samples with taxonomic annotation and merge to one otumatrix
kaijunames = []
kaijunames = glob.glob(args.indir + "/*.out.names")
otumatrix  = {}

# Extract and count reads from each sample. Merge Unknowns as well.
for file in kaijunames:
	sample = file.split("/")[-1].split(".out.names")[0]
	otumatrix[sample] = {}
	with open(file) as file_in:
		lines = []
		for line in file_in:
			type = line.split("\t")[0]
			taxa = line.split("\t")[-1].rstrip()
			if type == "U":
				otumatrix[sample]["Unknown"] = otumatrix[sample].get("Unknown", 0) + 1			
			elif type == "C":
				if taxa == "NA; NA; NA; NA; NA; NA; NA;":
					otumatrix[sample]["Unknown"] = otumatrix[sample].get("Unknown", 0) + 1
				else:
					otumatrix[sample][taxa] = otumatrix[sample].get(taxa, 0) + 1

# Create list of all otus
otus = []
for sample in otumatrix:
	for otu in otumatrix[sample] :
		if otu not in otus:
			otus.append(otu)
# Create Phyloseq otutable & taxtable
otutable   = {}
order_otu  = []
otu_file   = open(args.output + ".otu.tab", "w")
tax_file   = open(args.output + ".tax.tab", "w")
counter	   = 1
# Adds sample names to first row
for sample in otumatrix:
	otu_file.write('\t' + sample)
otu_file.write('\n')
# Add header to tax_file
tax_file.write('\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n')
# Adds OTUs rows
for otu in otus:
	otu_file.write("OTU" + str(counter).zfill(5))
	tax_file.write("OTU" + str(counter).zfill(5) + "\t" + "\t".join(otu.split(";")).rstrip() + "\n")
	counter = counter+1
	for sample in otumatrix:	
		if otu in otumatrix[sample]:
			otu_file.write("\t" + str(otumatrix[sample][otu]))
		else:
			otu_file.write("\t" + "0")
	otu_file.write("\n")
otu_file.close
tax_file.close
