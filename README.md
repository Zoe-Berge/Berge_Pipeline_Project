# Berge_Pipeline_Project
COMP 483 COMP BIO PYTHON PIPELINE PROJECT

This pipeline takes raw sequence reads from Human cytomegalovirus(HCMV) and performs differential expression analysis and genome assembly.

# Dependencies:
For the pipeline to be successful, download the following tools in your environment:

- Python 3.x (with Biopython)
  --> https://biopython.org/wiki/Documentation 
- Snakemake
  --> https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-1-mapping-reads 
- R (sleuth and tidyverse)
  --> Sleuth walkthrough: https://pachterlab.github.io/sleuth_walkthroughs/pval_agg/analysis.html
  --> Dplyr(Tidyverse): https://dplyr.tidyverse.org/   
- Kallisto
  --> https://pachterlab.github.io/kallisto/manual 
- Bowtie2
  --> https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
- SPAdes
  --> https://ablab.github.io/spades/ 
- NCBI blastn
  --> https://blast.ncbi.nlm.nih.gov/Blast.cgi
- SRA Toolkit
  --> https://github.com/ncbi/sra-tools/wiki/HowTo:-binary-installation
- SAMtools
  --> http://www.htslib.org/
- Entrez Direct
  --> https://www.ncbi.nlm.nih.gov/books/NBK179288/

# Data:
The required sequencing data is automated to download inside the Snakemake from the NCBI Sequence Read Archive(SRA) using the fastq-dump tool
Includes downloading the following: SRR5660030, SRR5660033, SRR5660044, SRR5660045

If you choose to download manually, use the following commands, and Snakemake will skip over the downloading step:

fasterq-dump --split-files SRR5660030

fasterq-dump --split-files SRR5660033

fasterq-dump --split-files SRR5660044

fasterq-dump --split-files SRR5660045

# Overview:

1. Extracts coding sequence(CDS) features from the HCMV genome (NC_006273.2/GCF_000845245.1) and uses Kallisto to index the results
2. Uses Kallisto to estimate transcript abundance(TPM) for each sample
3. Uses Sleuth(R package) to identify significant genes between 2dpi(2 days post-infection) and 6dpi(6 days post-infection) conditions
4. Maps reads against the HCMV genome using Bowtie2 and extracts mapped pairs with SAMtools to remove host background(human data)
5. Performs de novo assembly on the filtered viral reads using SPAdes
6. Finds the longest contig for each sample and uses BLAST+ to query a local database of Betaherpesvirinae(retrieved via Entrez Direct) to identify the top matching strains

# How to run:

1. Clone the repository using the following commands: 

  git clone https://github.com/Zoe-Berge/Berge_Pipeline_Project.git

2. Ensure all tools listed in dependencies are installed
3. Ensure you are in the pipeline directory
   
   cd Berge_Pipeline_Project

5. Run the snakemake command:

   snakemake --cores 2
