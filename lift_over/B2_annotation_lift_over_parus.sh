#!/bin/bash -l

#SBATCH -A naiss2023-5-169
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 2-00:00:00
#SBATCH -J annotation_lift_over
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/poster;
REF_name=Pmaj;
REF0=$MAINDIR/data/reference/Parus_major\
/GCF_001522545.3_Parus_major1.1_genomic.fasta;
SPEC_name=Aarv;
SPEC0=$WORKDIR/B1_prepare_reference\
/GCA_902810485.1_skylark_genome_genomic_no_IUPAC.fasta;
kraken=/crex/proj/snic2020-2-25/bin/kraken-accessed-2021-12-21/kraken/bin;

### Load modules
module load bioinfo-tools samtools/1.14 satsuma2/2016-12-07 BEDTools/2.29.2;

### Activate conda environment
source /sw/apps/bioinfo/CGAT/0.3.3/rackham/conda-install/etc/profile.d/conda.sh;
conda activate cgat-s;

### Create folders
mkdir $OUTDIR/B2_annotation_lift_over_parus;
OUTDIR=$OUTDIR/B2_annotation_lift_over_parus;

### Perform annotation lift-over
NAMES=$REF_name\_$SPEC_name;
cat $REF0 | awk '{if($0 ~ /^.*>/) print $1; else print;}' \
> $OUTDIR/$REF_name.fa;
REF=$OUTDIR/$REF_name.fa;
samtools faidx $REF;
cat $SPEC0 | awk '{if($0 ~ /^.*>/) print $1; else print;}' \
> $OUTDIR/$SPEC_name.fa;
SPEC=$OUTDIR/$SPEC_name.fa;

# Create synteny with Satsuma2
SatsumaSynteny2 -t $REF -q $SPEC \
-o $OUTDIR/satsuma_$NAMES -threads 20;
