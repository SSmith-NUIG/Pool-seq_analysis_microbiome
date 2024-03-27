#!/bin/sh

#SBATCH --job-name="align_col"
#SBATCH -o /data/ssmith/logs/VCF_pipeline_%A_%a.out
#SBATCH -e /data/ssmith/logs/VCF_pipeline_%A_%a.err
#SBATCH --array=210-214
#SBATCH -N 1
#SBATCH -n 8

module load Anaconda3
source activate wgs_env


# get our sample name which corresponds with the SLURM_ARRAY_TASK_ID we are processing.
# This will pull the second column value (sample name) from the novo2name CSV file
sample_name=$(awk -v pattern="$SLURM_ARRAY_TASK_ID"  -F',' '$1 == pattern { print $2 }' /data2/ssmith/novo2name.csv)

# get the basenames of the sample we are processing (forward and reverse files)
file_name_1=$(basename /data2/ssmith/fastqs/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_1.fq.gz)
file_name_2=$(basename /data2/ssmith/fastqs/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_2.fq.gz)

echo BEGINNING FASTQC SCRIPT

# the directory path we need for fastQC
dirPath=/data2/ssmith/fastqs/

fastqc "$dirPath""$SLURM_ARRAY_TASK_ID"_*.fq.gz --noextract -t 6 -a /data2/ssmith/adapters.txt --outdir=/data2/ssmith/QC/initial_QC/

echo TRIMMOMATIC ADAPTER TRIMMING

# run trimmomatic, see the doc for sepcific arguments
trimmomatic \
PE \
/data2/ssmith/fastqs/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_1.fq.gz \
/data2/ssmith/fastqs/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_2.fq.gz \
/data2/ssmith/trimmed/"$file_name_1" /data2/ssmith/discarded/"$file_name_1" \
/data2/ssmith/trimmed/"$file_name_2" /data2/ssmith/discarded/"$file_name_2" \
ILLUMINACLIP:/data2/ssmith/adapters.fa:2:30:10 \
MINLEN:50


echo BEGINNING FASTQC SCRIPT
# run fastQC on our trimmed reads
dirPath=/data2/ssmith/trimmed/
fastqc "$dirPath""$SLURM_ARRAY_TASK_ID"_*.fq.gz --noextract -t 6 -a /data2/ssmith/adapters.txt --outdir=/data2/ssmith/QC/after_QC/

echo BWA-MEM2 ALIGNMENT

# align the trimmed reads using bwa-mem2 and sort using samtools
bwa-mem2 mem -t 4 /data/ssmith/c_l_genome/apis_c_l_genome.fa  \
/data2/ssmith/trimmed/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_1.fq.gz \
/data2/ssmith/trimmed/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_2.fq.gz \
| samtools sort -@ 4 -o /data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.bam

echo SAMTOOLS CREATE UNMAPPED READS BAM FILE

# output all unmapped reads for microbiome analysis
samtools view -b -f 4 /data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.bam \
> /data2/ssmith/unmapped_bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_unmapped.bam

unmapped_name=$(basename /data2/ssmith/unmapped_bams/"$SLURM_ARRAY_TASK_ID"_*_unmapped.bam _unmapped.bam)

# sort the unmapped bam
samtools sort \
-n /data2/ssmith/unmapped_bams/"$unmapped_name"_unmapped.bam \
-o /data2/ssmith/unmapped_bams/"$unmapped_name"_unmapped.sorted.bam

rm /data2/ssmith/unmapped_bams/"$unmapped_name"_unmapped.bam

# turn the unmapped bams into unmapped fastqs for humann and metaphlan analysis
bamToFastq \
-i /data2/ssmith/unmapped_bams/"$unmapped_name"_unmapped.sorted.bam \
-fq /data2/ssmith/unmapped_fastqs/"$unmapped_name"_1.fastq.gz \
-fq2 /data2/ssmith/unmapped_fastqs/"$unmapped_name"_2.fastq.gz \

rm /data2/ssmith/unmapped_bams/"$unmapped_name"_unmapped.sorted.bam

# combine forward and reverse unmapped fastqs into one file
cat /data2/ssmith/unmapped_fastqs/"$unmapped_name"_1.fastq.gz /data2/ssmith/unmapped_fastqs/"$unmapped_name"_2.fastq.gz \
> /data2/ssmith/unmapped_fastqs/"$unmapped_name".fastq
