#!/bin/bash -l

#SBATCH -A naiss2023-5-169
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 20:00:00
#SBATCH -J annotation_lift_over
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
REF_name=Tgut; ## Pmaj ## Falb
REF0=$MAINDIR/data/reference/Taeniopygia_guttata\
/GCF_003957565.2_bTaeGut1.4.pri_genomic.no_W.fasta;
SPEC_name=Aarv;
SPEC0=$WORKDIR/B1_prepare_reference\
/GCA_902810485.1_skylark_genome_genomic_no_IUPAC.fasta;

### Load modules
module load bioinfo-tools BEDTools/2.29.2;

### Activate conda environment
conda activate liftover;

### Create folders
mkdir $OUTDIR/B3_annotation_lift_over;
OUTDIR=$OUTDIR/B3_annotation_lift_over;

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

# Remove duplicates (will remove these where different exons matches different scaffolds)
cgat gtf2gtf --method=remove-duplicates --duplicate-feature gene  \
-I $OUTDIR/$SPEC_name\_$REF_name\_liftover.gtf_polished \
-S $OUTDIR/$SPEC_name\_$REF_name\_liftover_polished_rmdup.gtf;

## Create region files

# Genomic regions
cat $SPEC0.fai | cut -f1,2 > $OUTDIR/$SPEC_name\_genome.bed;

# CDS regions
cat $OUTDIR/$SPEC_name\_$REF_name\_liftover.gtf_polished | \
awk -F'\t' '$3 == "CDS" {print}' | cgat gff2bed | grep -v "^#" | \
awk '{print $1"\t"$2"\t"$3}' | sort -k 1,1 -k 2n,2 | \
bedtools merge > $OUTDIR/$SPEC_name\_cds.bed;

# Non-CDS regions
bedtools complement -i $OUTDIR/$SPEC_name\_cds.bed -g $OUTDIR/$SPEC_name\_genome.bed \
> $OUTDIR/$SPEC_name\_noncds.bed;

# Exonic regions
cat $OUTDIR/$SPEC_name\_$REF_name\_liftover.gtf_polished | \
awk -F'\t' '$3 == "exon" {print}' | cgat gff2bed | grep -v "^#" | \
awk '{print $1"\t"$2"\t"$3}' | sort -k 1,1 -k 2n,2 | \
bedtools merge > $OUTDIR/$SPEC_name\_exonic.bed;

# Non-exonic regions
bedtools complement -i $OUTDIR/$SPEC_name\_exonic.bed -g $OUTDIR/$SPEC_name\_genome.bed \
> $OUTDIR/$SPEC_name\_nonexonic.bed;