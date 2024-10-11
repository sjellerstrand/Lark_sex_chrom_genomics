#!/bin/bash -l

#SBATCH -A snic2021-5-469
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J fastqc_raw
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT/;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;

### Load modules
module load bioinfo-tools FastQC/0.11.8 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/A1_fastqc_raw;
OUTDIR=$OUTDIR/A1_fastqc_raw;

### Set up file info
find $WORKDIR/raw -name "*.fastq.gz" > $OUTDIR/INDS.txt;

### Create quality reports
export OUTDIR;
parallel 'fastqc {} -o $OUTDIR' :::: $OUTDIR/INDS.txt;

