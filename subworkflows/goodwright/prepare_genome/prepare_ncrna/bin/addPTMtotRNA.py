from Bio import SeqIO
import sys

# Author: Charlotte Capitanchik
# adds ptms to tRNA file

fasta = sys.argv[1]
fasta_sequences = SeqIO.parse(open(fasta), 'fasta')
new_file = sys.argv[2]
f_file = open(new_file,"w")

# Read the fasta file as a dictionary

for f in fasta_sequences:
    name, sequence = f.id, f.seq
    sequence = str(sequence).upper()
    if sequence[-3:] != "CCA":
        sequence = sequence + "CCA"
    if name.find("His") != -1:
        sequence = "G" + sequence
    f_file.write(">"+name+"\n"+sequence+"\n")

f_file.close()