#!/bin/bash -l

#SBATCH -A snic2021-5-469
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J fastqc_trim
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;

### Load modules
module load bioinfo-tools FastQC/0.11.8 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/A3_fastqc_quality_trim;
OUTDIR=$OUTDIR/A3_fastqc_quality_trim;

### Set up file info
find $WORKDIR/quality_trim -name "*.fq.gz" > $OUTDIR/INDS.txt;

### Create quality reports
export OUTDIR;
parallel 'fastqc {} -o $OUTDIR' :::: $OUTDIR/INDS.txt;
