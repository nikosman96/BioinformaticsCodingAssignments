#BINF6410 Assignment 3: FASTA Parser 

import sys


import numpy as np

sequence = open("seq.fasta", "r")
fprimer = list("GTGCCAGCMGCCGCGGTAA")
#Hardcoded the filename and primers into variables that can be used.

#Reverse complement the reverse primer using Biopython
from Bio.Seq import Seq

reverseprimer = Seq("ACAGCCATGCANCACCT")
str(reverseprimer)
rprimer = list(reverseprimer.reverse_complement())

#Read the file to determine the raw number of sequences and the total sequence length

numLines = 0
numChars = 0
with open("seq.fasta") as file:
  for line in file:
    li=line.strip()
    if not li.startswith(">"):
        numLines += 1
        numChars += len(line)
print('Number of Sequences: ', numLines)
print("Total Sequence Length: ", numChars)

sequences = {}
good_reads = []
bad_reads = []

with open("paired.fasta", "r") as file:
        for line in file:
                line=line.rstrip()
                if (line[0] == ">"):
                        header = line
                        sequences[header] = ""
                else:
                        data = line
                        sequences[header] += data

min_length = 150

# figure out which reads are good/bad
for header in sequences.keys():
        if (len(sequences[header]) > min_length):
                good_reads.append(header)
        else:
                bad_reads.append(header)

# write good reads
with open("paired.fasta", "w+") as good_out:
        for header in good_reads:
                good_out.write("{}\n{}\n".format(header, sequences[header]))

#Write bad reads to a seperate log file.
with open("process.log", "w+") as bad_out:
        for header in bad_reads:
                bad_out.write("{}\nExcluded because too short\n".format(header))

