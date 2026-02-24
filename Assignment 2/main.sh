############################################################
# RNA-seq Alignment and Quantification Pipeline
# Author: [Your Name]
# Environment: Remote HPC Cluster (SLURM-managed)
# Purpose: Align RNA-seq reads, assemble transcripts, and quantify gene expression
############################################################


############################################################
# Step 1: Set up working directory and download reference data
############################################################

# Request an interactive SLURM session with 4 CPU cores and 100 GB RAM
salloc -c 4 --mem=100g

# Create main project folder
mkdir A_ToxidoAlignment
cd A_ToxidoAlignment

# Create a subdirectory for input data
mkdir data
cd data

# Download example RNA-seq dataset (contains HBR and UHR samples)
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar

# Extract the compressed archive
tar -xvf HBR_UHR_ERCC_ds_5pc.tar

# Decompress all FASTQ files
gunzip -v *.gz

# Download the human genome (GRCh38 primary assembly)
wget ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip -v Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38.fa

# Give read and execute permissions for the reference genome
chmod 755 hg38.fa

# Download Ensembl GTF annotation file (release 110)
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip -v Homo_sapiens.GRCh38.110.gtf.gz
mv Homo_sapiens.GRCh38.110.gtf GRCh38.110.gtf
chmod 755 GRCh38.110.gtf

# Go back to the main directory
cd ~/A_ToxidoAlignment


############################################################
# Step 2: Quality control using FastQC
############################################################

# Check current directory
pwd
cd ~/A_ToxidoAlignment

# Create folder to store FastQC reports
mkdir fastqc_reports
cd fastqc_reports

# Load FastQC module from HPC environment
module load fastqc

# Verify FastQC installation and options
fastqc --help

# Run FastQC on all FASTQ files and output results to fastqc_reports/
cd ~/A_ToxidoAlignment
fastqc -o fastqc_reports data/*.fastq

# Return to main directory
cd ~/A_ToxidoAlignment


############################################################
# Step 3: Build HISAT2 genome index
############################################################

pwd
cd ~/A_ToxidoAlignment
mkdir hg38_index
cd hg38_index

# Allocate resources for index building
salloc -c 4 --mem=100g

# Load HISAT2 module
module load hisat2

# Build genome index for HISAT2 aligner
hisat2-build ../data/hg38.fa hg38_ref_ind

# Return to main project folder
cd ~/A_ToxidoAlignment


############################################################
# Step 4: Align paired-end reads using HISAT2
############################################################

pwd
cd ~/A_ToxidoAlignment
mkdir aligned_rna_seq
cd aligned_rna_seq

# Allocate resources for alignment
salloc -c 4 --mem=100g
module load hisat2

# Loop through all read1 files and align with paired read2
for R1 in ../data/*read1.fastq; do
  # Define matching read2 file
  R2=${R1/read1.fastq/read2.fastq}
  
  # Extract base name to use for output SAM files
  base=$(basename "$R1" _chr22.read1.fastq)
  
  # Run HISAT2 with 4 threads and output SAM file
  hisat2 -x ../hg38_index/hg38_ref_ind \
         -1 "$R1" -2 "$R2" \
         -p 4 -S "${base}_aligned.sam"
done

# (Output: 6 SAM alignment files, one per sample)

# Preview the first 500 lines of a SAM file
head -500 HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq_aligned.sam

cd ~/A_ToxidoAlignment


############################################################
# Step 5: Convert and sort alignments using SAMtools
############################################################

pwd
cd ~/A_ToxidoAlignment
mkdir sorted_rna-seq
cd sorted_rna-seq

# Request resources for SAMtools operations
salloc -c 4 --mem=100g
module load samtools

# Sort each SAM file into a BAM file
for file in ../aligned_rna_seq/*.sam; do
  base=$(basename "$file" _aligned.sam)
  samtools sort "$file" -o "${base}_sorted.bam"
done

# Alternative sorting loop with indexing
for sam in ../aligned_rna_seq/*.sam; do
  b=$(basename "$sam")
  b=${b%.sam}
  b=${b%-chr22.read1.fastq_aligned}
  b=${b%.read1.fastq_aligned}
  out="${b}_sorted.bam"
  samtools sort -@4 -o "$out" "$sam"
  samtools index "$out"
done

# Clean up filenames by removing redundant labels
for f in *_sorted.bam; do
  mv "$f" "${f/_Build37-ErccTranscripts/}"
done

# Inspect BAM file contents
cat HBR_Rep1_ERCC-Mix2_sorted.bam | head
samtools view HBR_Rep1_ERCC-Mix2_sorted.bam | head

# Index all BAM files for random access
for bam in *.bam; do
  samtools index "$bam"
done

cd ~/A_ToxidoAlignment


############################################################
# Step 6: Transcript assembly using StringTie
############################################################

pwd
cd ~/A_ToxidoAlignment
salloc -c 4 --mem=100g
mkdir assembly
cd assembly
module load stringtie

# Assemble transcripts for each sample BAM file
for bam in ../sorted_rna-seq/*.bam; do
  base=$(basename "$bam" _sorted.bam)
  stringtie "$bam" -o "${base}.gtf" \
    -G /global/teaching-home/sa3317041/A_ToxidoAlignment/data/GRCh38.110.gtf -e
done

# List output GTF files and preview one
cd ~/A_ToxidoAlignment/assembly
ls -lh
head -5 HBR_Rep1_ERCC-Mix2.gtf
cd ~/A_ToxidoAlignment


############################################################
# Step 7: Merge transcript assemblies across all samples
############################################################

cd ~/A_ToxidoAlignment/assembly
ls -lh *.gtf

# Create a text file listing all sample GTFs to merge
nano mergelist.txt
# (Add the following lines manually)
# assembly/HBR_Rep1_ERCC-Mix2.gtf
# assembly/HBR_Rep2_ERCC-Mix2.gtf
# assembly/HBR_Rep3_ERCC-Mix2.gtf
# assembly/UHR_Rep1_ERCC-Mix1.gtf
# assembly/UHR_Rep2_ERCC-Mix1.gtf
# assembly/UHR_Rep3_ERCC-Mix1.gtf

cd ~/A_ToxidoAlignment
module load stringtie

# Merge all sample assemblies into one reference transcriptome
stringtie --merge \
  -G /global/teaching-home/sa3317041/A_ToxidoAlignment/data/GRCh38.110.gtf \
  -o /global/teaching-home/sa3317041/A_ToxidoAlignment/assembly/stringtie_merged.gtf \
  /global/teaching-home/sa3317041/A_ToxidoAlignment/assembly/mergelist.txt

cd ~/A_ToxidoAlignment


############################################################
# Step 8: Quantify expression levels for each sample (Ballgown input)
############################################################

pwd
cd ~/A_ToxidoAlignment
salloc -c 4 --mem=100g
module load stringtie
mkdir -p ballgown

# Generate normalized transcript abundances (FPKM/TPM) and Ballgown tables
for bam in sorted_rna-seq/*.bam; do
  base=$(basename "$bam" _sorted.bam)
  stringtie "$bam" \
    -p 4 \
    -G assembly/stringtie_merged.gtf \
    -e -B \
    --rf \
    -o ballgown/${base}/transcripts.gtf \
    -A ballgown/${base}_norm_cnts.tsv
done

# Preview normalized count table
head -5 ballgown/HBR_Rep1_ERCC-Mix2_norm_cnts.tsv


############################################################
# Step 9: Generate gene and transcript count matrices using prepDE.py
############################################################

cd ~/A_ToxidoAlignment/ballgown

# Download official prepDE Python script from StringTie
wget https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
module load python
chmod +x prepDE.py3

# Prepare a sample sheet (CSV) linking sample names to transcript paths
printf "sample,path\n" > samples.csv
for d in */; do
  s=${d%/}
  if [ -f "$s/transcripts.gtf" ]; then
    echo "$s,$PWD/$s/transcripts.gtf" >> samples.csv
  fi
done

# Cleaning up any extra spaces in the CSV
sed -i 's/, /,/g' samples.csv

# Run prepDE to generate count matrices
python prepDE.py -i samples.csv -g gene_count_matrix.csv -t transcript_count_matrix.csv

# List and preview results
ls -lh gene_count_matrix.csv transcript_count_matrix.csv
head -5 gene_count_matrix.csv
head -5 transcript_count_matrix.csv


############################################################
# Step 10: Combine normalized counts using MATLAB
############################################################

pwd
cd ~/A_ToxidoAlignment
salloc -c 4 --mem=100g

# Check available MATLAB modules and load one
module avail matlab
module load matlab/R2023a

# Run MATLAB script in non-interactive mode
matlab -nodisplay -nodesktop < script_combine_norm_counts.m

# Preview combined normalized counts output
head -20 counts_combined_tpm.csv
