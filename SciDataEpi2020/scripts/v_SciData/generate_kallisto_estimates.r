# data is from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85331
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5576565/

# Preprocess
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# FASTQC
# fastqc the files to check for complications..       
  system("data/raw/fastqc SRR4011*")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# FASTP
# Trimming before alignment
  system("mkdir trimmed")

  system("fastp 
    -i data/raw/SRR4011898_H9_0_1_1.fastq.gz 
    -I data/raw/SRR4011898_H9_0_1_2.fastq.gz 
    -o data/raw/trimmed/SRR4011898_H9_0_1_1_trimmed.fastq.gz 
    -O data/raw/trimmed/SRR4011898_H9_0_1_2_trimmed.fastq.gz 
    -w 8 
    -h data/raw/trimmed/fastp_report_H9_0_1
  ")
  system("fastp 
    -i data/raw/SRR4011899_H9_0_2_1.fastq.gz 
    -I data/raw/SRR4011899_H9_0_2_2.fastq.gz 
    -o data/raw/trimmed/SRR4011899_H9_0_2_1_trimmed.fastq.gz 
    -O data/raw/trimmed/SRR4011899_H9_0_2_2_trimmed.fastq.gz 
    -w 8 
    -h data/raw/trimmed/fastp_report_H9_0_2
  ")
  system("fastp 
    -i data/raw/SRR4011904_H9_30_1_1.fastq.gz 
    -I data/raw/SRR4011904_H9_30_1_2.fastq.gz 
    -o data/raw/trimmed/SRR4011904_H9_30_1_1_trimmed.fastq.gz 
    -O data/raw/trimmed/SRR4011904_H9_30_1_2_trimmed.fastq.gz 
    -w 8 
    -h data/raw/trimmed/fastp_report_H9_30_1
  ")
  system("fastp 
    -i data/raw/SRR4011905_H9_30_2_1.fastq.gz 
    -I data/raw/SRR4011905_H9_30_2_2.fastq.gz 
    -o data/raw/trimmed/SRR4011905_H9_30_2_1_trimmed.fastq.gz 
    -O data/raw/trimmed/SRR4011905_H9_30_2_2_trimmed.fastq.gz 
    -w 8 -h data/raw/trimmed/fastp_report_H9_30_2
  ")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# FASTQC       
# fastqc the trimmed files to check for further complications..       
  system("data/raw/trimmed/fastqc SRR4011*")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# KALLISTO
# The first  step is to create a kallisto index. This uses a series of transcript FASTAs. 
#   I download the fasta file from "ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/"
#   Giving the following file: Homo_sapiens.GRCh38.cdna.all.fa.gz
# We use this compressed fasta format to generate a kmer index of sorts (i think)


# Place the cDNA into the data directory
  system("kallisto index -i data/HS_GRCh38.cdna.all.test.idx 
  ~/Documents/genomes/human_hg38/Homo_sapiens.GRCh38.cdna.all.fa.gz")

# The second step is to use this index and perform pseudoalignment to is using our sequenced fasta files...
# I did this manually as we only have a few samples...
# Use this index and perform pseudoalignment of our trimmed sequenced fasta files...
  system("kallisto quant 
    -i data/HS_GRCh38.cdna.all.idx -o data/kallisto_out/d0_1/trimmed 
    -b 100 
    -t 12 
    data/raw/trimmed/SRR4011898_H9_0_1_1_trimmed.fastq.gz 
    data/raw/trimmed/SRR4011898_H9_0_1_2_trimmed.fastq.gz
  ")
  system("kallisto quant 
    -i data/HS_GRCh38.cdna.all.idx 
    -o data/kallisto_out/d0_2/trimmed 
    -b 100 
    -t 12 
    data/raw/trimmed/SRR4011899_H9_0_2_1_trimmed.fastq.gz 
    data/raw/trimmed/SRR4011899_H9_0_2_2_trimmed.fastq.gz
  ")
  system("kallisto quant 
    -i data/HS_GRCh38.cdna.all.idx 
    -o data/kallisto_out/d30_1/trimmed 
    -b 100 
    -t 12 
    data/raw/trimmed/SRR4011904_H9_30_1_1_trimmed.fastq.gz 
    data/raw/trimmed/SRR4011904_H9_30_1_2_trimmed.fastq.gz
  ")
  system("kallisto quant 
    -i data/HS_GRCh38.cdna.all.idx 
    -o data/kallisto_out/d30_2/trimmed 
    -b 100 
    -t 12 
    data/raw/trimmed/SRR4011905_H9_30_2_1_trimmed.fastq.gz 
    data/raw/trimmed/SRR4011905_H9_30_2_2_trimmed.fastq.gz
  ")
  
