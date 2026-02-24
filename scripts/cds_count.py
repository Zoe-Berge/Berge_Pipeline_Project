#importing dependancies 
import sys
from Bio import SeqIO

def extract_cds(input_gbk, output_fasta): #function with the input of the genbank file and the fasta file 
    #read the GenBank file
    with open(input_gbk, "r") as handle:
        record = SeqIO.read(handle, "genbank") #loads genbank into and assigns it to record 
    
    #scans and filters for CDS features
    cds_features = [f for f in record.features if f.type == "CDS"]
    
    with open(output_fasta, "w") as out:
        for feature in cds_features:
            protein_id = feature.qualifiers.get("protein_id", [None])[0] #pulls protein ID and uses it as the header 
            
            if protein_id:
                sequence = feature.extract(record.seq) #cuts and pulls coding sequences, assinging them to a new variable  
                out.write(f">{protein_id}\n{sequence}\n") #writes out the info
    
if __name__ == "__main__": #checks if its being run by main program (snakemake) 
    if len(sys.argv) == 3:
        extract_cds(sys.argv[1], sys.argv[2]) #executes function using genbank file and output fasta file from terminal list 

