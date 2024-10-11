#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 30:00:00
#SBATCH -J ancestral_generation_times
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
ITERATIONS=15000000;
BURNIN=5000000;

MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
METADATA=$WORKDIR/metadata;
RESOURCES=$METADATA/resources;

### Load modules
module load R_packages/4.0.0;

### Create folders
mkdir $OUTDIR/ancestral_generation_times;
OUTDIR=$OUTDIR/ancestral_generation_times;

### Plot statistics
Rscript $MAINDIR/scripts/$PROJECT/analyses/ancestral_generation_times.r \
--args ITERATIONS=$ITERATIONS BURNIN=$BURNIN OUTDIR=$OUTDIR METADATA=$METADATA RESOURCES=$RESOURCES;

