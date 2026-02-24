import sys
from Bio import SeqIO

def extract_cds(input_gbk, output_fasta):
    #read the GenBank file
    with open(input_gbk, "r") as handle:
        record = SeqIO.read(handle, "genbank")
    
    #filter for CDS features
    cds_features = [f for f in record.features if f.type == "CDS"]
    
    with open(output_fasta, "w") as out:
        for feature in cds_features:
            protein_id = feature.qualifiers.get("protein_id", [None])[0]
            
            if protein_id:
                sequence = feature.extract(record.seq)
                out.write(f">{protein_id}\n{sequence}\n")
    
if __name__ == "__main__":
    if len(sys.argv) == 3:
        extract_cds(sys.argv[1], sys.argv[2])

