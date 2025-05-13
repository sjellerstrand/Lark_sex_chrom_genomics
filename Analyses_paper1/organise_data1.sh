#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J organise_data1
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/cfs/klemming/projects/supr/snic2020-2-25/user_data/simon/Sylvioidea;
PROJECT1=Skylark_2021;
PROJECT2=Rasolark_2021;
WORKDIR1=$MAINDIR/data/$PROJECT1;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/$PROJECT1;
GENES=$WORKDIR1/D4_align_shared_genes/metadata/genes_info.tsv;

### Load modules
module load PDC R;

### Create folders
mkdir $OUTDIR/organise_data1;
OUTDIR=$OUTDIR/organise_data1;

### Plot statistics
Rscript $MAINDIR/scripts/$PROJECT1/analyses/organise_data1.r \
--args PROJECT1=$PROJECT1 PROJECT2=$PROJECT2 WORKDIR1=$WORKDIR1 GENES=$GENES OUTDIR=$OUTDIR;

