import sys
import os
from Bio import SeqIO

sample = sys.argv[1] #current sample number being run 

#initilizing all my paths
input_file = f"results/{sample}_assembly/contigs.fasta"
longest_out = f"results/{sample}_longest_contig.fasta"
blast_out = f"results/{sample}_blast_results.txt"
db_path = "results/hcmv_db"

#finds the longest contig
records = list(SeqIO.parse(input_file, "fasta"))
longest_record = max(records, key=lambda x: len(x.seq))
SeqIO.write(longest_record, longest_out, "fasta")

#constructs the BLAST command

#tab-delimited format 6 with specific columns
outfmt = "6 sacc pident length qstart qend sstart send bitscore evalue stitle"
blast_command = f"blastn -query {longest_out} -db {db_path} -max_target_seqs 5 -max_hsps 1 -outfmt '{outfmt}' -out {blast_out}"

os.system(blast_command)