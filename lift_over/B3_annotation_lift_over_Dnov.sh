#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 12:00:00
#SBATCH -J annotation_lift_over
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
REF_name=Tgut;
REF0=$MAINDIR/data/reference/Taeniopygia_guttata\
/GCF_003957565.2_bTaeGut1.4.pri_genomic.no_W.fasta;
SPEC_name=Dnov;
SPEC0=$MAINDIR/data/reference/Dromaius_novaehollandiae\
/GCF_036370855.1_bDroNov1.hap1_genomic.no_W_chrnr.fasta;

### Load modules
module load bioinfo-tools BEDTools/2.29.2;

### Activate conda environment
conda activate liftover;

### Create folders
mkdir $OUTDIR/B3_annotation_lift_over_Dnov;
OUTDIR=$OUTDIR/B3_annotation_lift_over_Dnov;

### Perform annotation lift-over
NAMES=$REF_name\_$SPEC_name;
cat $REF0 | awk '{if($0 ~ /^.*>/) print $1; else print;}' \
> $OUTDIR/$REF_name.fa;
REF=$OUTDIR/$REF_name.fa;
cat $SPEC0 | awk '{if($0 ~ /^.*>/) print $1; else print;}' \
> $OUTDIR/$SPEC_name.fa;
SPEC=$OUTDIR/$SPEC_name.fa;

# Modify GTF to remove unwanted whitespaces in attributes
cat $(echo "${REF0%.*}".gtf) | awk -F'\t' -v OFS='\t' '
{gsub(" ","_",$9); gsub("_\""," \"",$9); gsub(";_","; ",$9); gsub("; $",";",$9); print}' \
> $OUTDIR/$REF_name\_modified.gtf;
GTF=$OUTDIR/$REF_name\_modified.gtf;

# Perform annotation lift-over
liftoff -g $GTF -o $OUTDIR/$SPEC_name\_$REF_name\_liftover.gtf -u $OUTDIR/unmapped_features.txt \
-dir $OUTDIR/intermediate_files -polish -flank 0.5 -p 20 -overlap 1.0 $SPEC $REF -cds;
