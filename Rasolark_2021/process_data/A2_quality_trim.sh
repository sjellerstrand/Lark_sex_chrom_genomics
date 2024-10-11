#!/bin/bash -l

#SBATCH -A snic2021-5-469
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J quality_trim
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
METADATA=$WORKDIR/metadata;

### Load modules
module load bioinfo-tools trimmomatic/0.36 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/A2_quality_trim;
OUTDIR=$OUTDIR/A2_quality_trim;

### Set up file info
find $WORKDIR/raw -name "*.fastq.gz" > $OUTDIR/INDS1.txt;
cat $OUTDIR/INDS1.txt | awk -F"_" '{NF-=2; gsub(" ", "_"); print}' | \
sort -u > $OUTDIR/INDS2.txt;
cat $OUTDIR/INDS2.txt | rev | cut -d'/' -f1 | rev \
> $OUTDIR/INDS3.txt;
paste -d '\t' $OUTDIR/INDS2.txt \
$OUTDIR/INDS3.txt > $OUTDIR/INDS.txt;

### Remove adapters and quality trim reads

## Define function
trim() {

# Input parameters
FILE=$1;
ID=$2;

# Remove adapters and quality trim reads
trimmomatic PE $FILE\_R1_001.fastq.gz $FILE\_R2_001.fastq.gz \
$OUTDIR/$ID\_R1_paired.fq.gz $OUTDIR/$ID\_R1_unpaired.fq.gz \
$OUTDIR/$ID\_R2_paired.fq.gz $OUTDIR/$ID\_R2_unpaired.fq.gz \
ILLUMINACLIP:$METADATA/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:120;
};

## Excecute function in parallell
export OUTDIR METADATA;
export -f trim;
parallel --colsep '\t' 'trim {1} {2}' :::: $OUTDIR/INDS.txt;
