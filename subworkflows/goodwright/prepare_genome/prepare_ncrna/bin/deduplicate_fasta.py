# Charlotte Capitanchik #
# 21.03.2018 #
# Get the unique sequences in a fasta file and merge the duplicates into the new fasta header so you know what has been collapsed #

from Bio import SeqIO
import sys

fasta =  sys.argv[1]


fasta_sequences = SeqIO.parse(open(fasta), 'fasta')
seq_dict = {}
total_count = 0

# Read the fasta file as a dictionary

for f in fasta_sequences:
    name, sequence = f.id, f.seq
    sequence = str(sequence).upper()
    # Make the sequence DNA if it is RNA (e.g from gtRNAdb)
    sequence = sequence.replace("U","T")
    seq_dict[name] = sequence
    total_count +=1

# Flip the dictionary in reverse so the key is the sequence
# and the value is all the RNAs that have the same identical sequence

rev_seq_dict = {}
for key, value in seq_dict.items():
    rev_seq_dict.setdefault(value, set()).add(key)

total_non_redundant = len(rev_seq_dict)
for key in rev_seq_dict.keys():
    rev_seq_dict[key] = ";".join(rev_seq_dict[key])

report_name = fasta.replace(".fa","-deduplicated-report.txt")
report_file = open(report_name, "w")
report_file.write("There were " + str(total_count) + " sequences in the original file: " + str(fasta) + "\n")
report_file.write("If we remove identical sequences there are " + str(total_non_redundant) + " sequences remaining.\n")
report_file.close()

new_ff_name = fasta.replace(".fa", "-deduplicated.fa")
new_fasta_file = open(new_ff_name, "w")
for key,value in rev_seq_dict.items():
    new_fasta_file.write(">" + value + "\n" + key + "\n")