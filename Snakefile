import sys
executable = sys.executable
#added the line above to make sure it executes through by virtual environment and not just python3

#list of sample SSR number to be looked at 
SAMPLES = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]
#SAMPLES = ["SRR5660030"]
#Above is a test data option 
DATA_DIR = "data"
#defined seperately so it can be redirected to test_data

rule all:
    input: 
        "Berge_PipelineReport.txt"

rule download_patient_data: #download seperately !!!! ????
    output:
        r1 = DATA_DIR + "/{sample}_1.fastq",
        r2 = DATA_DIR + "/{sample}_2.fastq"
    shell:
        # This only runs if the files aren't already in test_data/
        "fasterq-dump --split-files {wildcards.sample} -O " + DATA_DIR

rule download_hcmv_genbank: #downloads the refrence genome 
    output: "data/hcmv_ref.gbk"
    shell:
        "efetch -db nucleotide -id NC_006273.2 -format gbwithparts > {output}"

rule download_hcmv_genome: #downloads the refrence genome 
    output:
        "data/hcmv_genome.fasta"
    shell:
        "esearch -db nucleotide -query 'NC_006273.2' | efetch -format fasta > {output}"

#extract coding sequences using my python script
rule extract_cds:
    input: "data/hcmv_ref.gbk"
    output: 
        fasta = "data/hcmv_cds.fasta",
        statement = "results/cds_statement.txt"
    run:
        shell("python scripts/cds_count.py {input} {output.fasta}")
        with open(output.fasta, "r") as f:
            count = sum(1 for line in f if line.startswith(">")) #counts lines with > because each one is a sequence
        with open(output.statement, "w") as out:
            out.write(f"The HCMV genome (GCF_000845245.1) has {count} CDS.\n") #writes the number of sequences

#uses the coding sequences and indexs them using kallisto
rule kallisto_index:
    input: "data/hcmv_cds.fasta"
    output: "results/hcmv_index.idx"
    shell: "kallisto index -i {output} {input}"

#quantifys TPM for each sample
rule kallisto_quant:
    input:
        idx = "results/hcmv_index.idx",
        r1 = DATA_DIR + "/{sample}_1.fastq", 
        r2 = DATA_DIR + "/{sample}_2.fastq" #r1 and r2 are my paired end fastq files 
    output:
        "results/{sample}/abundance.tsv", #has estimated counts and TPM for every gene.
        "results/{sample}/abundance.h5"  #stores the bootstrab data to be used by sleuth 
    shell:
        "kallisto quant -i {input.idx} -o results/{wildcards.sample} -b 10 {input.r1} {input.r2}"
        #-b 1 re-runs the estimation for the paired end fastqs 10 times for better accuracy 

#runs the R script for Sleuth
rule run_sleuth:
    input:
        #this tells snakemake to wait until all of the samples are quantified
        #it was giving me errors withough this
        expand("results/{sample}/abundance.tsv", sample=SAMPLES)
    output:
        "results/significant_transcripts.temp"
    shell:
        "Rscript scripts/sleuth_analysis.R"

rule bowtie2_index: #builds the bowtie index 
    input:
        "data/hcmv_genome.fasta"
    output:
        # Bowtie2 creates multiple files starting with this prefix
        multiext("results/hcmv_bowtie_idx", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    shell:
        "bowtie2-build {input} results/hcmv_bowtie_idx"

#Bowtie2 Mapping
rule bowtie2_mapping_filter:
    input:
        r1 = DATA_DIR + "/{sample}_1.fastq",
        r2 = DATA_DIR + "/{sample}_2.fastq",
        idx = rules.bowtie2_index.output #calls the previous rules output 
    output:
        mapped_1 = "results/{sample}_mapped_1.fastq",
        mapped_2 = "results/{sample}_mapped_2.fastq",
        counts = "results/{sample}_counts.txt"
    shell:
        """
        #maps reads and pipes to samtools to keep only mapped pairs (-f 1 -F 12)

        bowtie2 -x results/hcmv_bowtie_idx -1 {input.r1} -2 {input.r2} | \
        samtools fastq -f 1 -F 12 -1 {output.mapped_1} -2 {output.mapped_2} -
        
        #count read pairs (total lines / 4)

        BEFORE=$(echo $(cat {input.r1} | wc -l) / 4 | bc)
        AFTER=$(echo $(cat {output.mapped_1} | wc -l) / 4 | bc)
        
        #write it to report 

        echo "Sample {wildcards.sample} had $BEFORE read pairs before and $AFTER read pairs after Bowtie2 filtering." > {output.counts}
        """

#SPAdes Assembly
rule assemble:
    input: 
        r1 = "results/{sample}_mapped_1.fastq",
        r2 = "results/{sample}_mapped_2.fastq"
    output: "results/{sample}_assembly/contigs.fasta"
    shell: "spades.py -k 127 -1 {input.r1} -2 {input.r2} -o results/{wildcards.sample}_assembly"

rule download_subfamily: #downloads the betaherpesvirinae data for a local database  
    output:
        "data/betaherpesvirinae.fasta"
    shell:
        #searches and gets all the results from the search 
        "esearch -db nucleotide -query 'Betaherpesvirinae[Organism]' | efetch -format fasta > {output}"

rule make_blast_db: #this makes the local database 
    input:
        "data/betaherpesvirinae.fasta"
    output:
        "results/hcmv_db.nsq" 
    shell:
        #makes the database with the nucleotide results 
        "makeblastdb -in {input} -dbtype nucl -out results/hcmv_db"

#run the blast script
rule run_blast_analysis:
    input:
        db = "results/hcmv_db.nsq",
        contigs = "results/{sample}_assembly/contigs.fasta"
    output:
        longest = "results/{sample}_longest_contig.fasta",
        results = "results/{sample}_blast_results.txt"
    shell:
        "{executable} scripts/find_and_blast.py {wildcards.sample}"

rule build_final_report:
    input:
        #assigns all of the gathered information
        cds = "results/cds_statement.txt",
        sleuth = "results/significant_transcripts.temp",
        #expand inputs each of the sample numbers
        counts = expand("results/{sample}_counts.txt", sample=SAMPLES),
        blast = expand("results/{sample}_blast_results.txt", sample=SAMPLES)
    output:
        "Berge_PipelineReport.txt"
    run:
        with open(output[0], "w") as out:
            #write CDS info
            with open(input.cds) as f: out.write(f.read() + "\n")
            
            #write Sleuth info
            out.write("target_id\ttest_stat\tpval\tqval\n")
            with open(input.sleuth) as f: out.write(f.read() + "\n")
            
            #write bowtie2 and blast info for all the samples 
            for sample in SAMPLES:
                #add bowtie2 counts
                count_file = f"results/{sample}_counts.txt"
                with open(count_file) as f: out.write(f.read())
                
                #add blast hits
                out.write(f"\n{sample}:\n")
                out.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
                blast_file = f"results/{sample}_blast_results.txt"
                with open(blast_file) as f: out.write(f.read() + "\n")