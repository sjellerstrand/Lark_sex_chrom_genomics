#!/bin/bash -l

#SBATCH -A snic2021-5-469
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
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
REF_name=Tgut;
REF0=$MAINDIR/data/reference/Taeniopygia_guttata\
/GCF_003957565.2_bTaeGut1.4.pri_genomic.no_W.fasta;
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
mkdir $OUTDIR/B2_annotation_lift_over;
OUTDIR=$OUTDIR/B2_annotation_lift_over;

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

## Annotation lift-over with kraken

# Create config file
echo -e "[genomes]\n\
$REF_name\t$REF\n\
$SPEC_name\t$SPEC\n\
\n\
[pairwise-maps]\
\n$REF_name\t$SPEC_name\t\
$OUTDIR/satsuma_$NAMES/satsuma_summary.chained.out"\
> $OUTDIR/$NAMES.config;

# Modify GTF to remove unwanted whitespaces in attributes
cat $(echo "${REF0%.*}".gtf) | awk -F'\t' -v OFS='\t' '
{gsub(" ","_",$9); gsub("_\""," \"",$9); gsub(";_","; ",$9); gsub("; $",";",$9); print}' \
> $OUTDIR/$REF_name\_modified.gtf;

# Annotation lift-over
$kraken/RunKraken -c $OUTDIR/$NAMES.config -s $OUTDIR/$REF_name\_modified.gtf \
-S $REF_name -T $SPEC_name -o $OUTDIR/$SPEC_name\_$REF_name\_liftover.gtf;

# Sort GTF file
cgat gtf2gtf --method=sort --sort-order=contig+gene \
-I $OUTDIR/$SPEC_name\_$REF_name\_liftover.gtf \
-S $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted.gtf;

# Remove duplicates (will remove these where different exons matches different scaffolds)
cgat gtf2gtf --method=remove-duplicates --duplicate-feature gene  \
-I $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted.gtf \
-S $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted_rmdup.gtf;

# Get longest transcript
cgat gtf2gtf --method=filter --filter longest-transcript \
-I $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted_rmdup.gtf \
-S $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted_rmdup_longesttrans.gtf;

## Create region files

# Genomic regions
cat $SPEC0.fai | cut -f1,2 > $OUTDIR/$SPEC_name\_genome.bed;

# Get features with information
cat $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted_rmdup.gtf | awk -F'\t' '{ print $1"\t"$4-1"\t"$5"\t"$3" "$9}' | \
awk -F' ' '{for(i=4; i<=NF; i++) {if($i ~ /^gene_id/) print $1"\t"$2"\t"$3"\t"$4"\t"$(i+1)"\t"$0}}' | \
awk -F' ' '{for(i=5; i<=NF; i++) {if($i ~ /^exon_number/) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\texon_"$(i+1)"\t"$0}}'| \
awk -F' ' '{for(i=5; i<=NF; i++) {if($i ~ /^transcript_id/) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$(i+1)}}' | \
tr -d '"' | tr -d ';' > $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted_rmdup_features.bed;

# CDS regions
cat $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted_rmdup.gtf | \
awk -F'\t' '$3 == "CDS" {print}' | cgat gff2bed | grep -v "^#" | \
awk '{print $1 "\t" $2 "\t" $3}' | sort -k 1,1 -k 2n,2 | \
bedtools merge > $OUTDIR/$SPEC_name\_cds.bed;

# CDS regions of longest transcript
cat $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted_rmdup_longesttrans.gtf | \
awk -F'\t' '$3 == "CDS" {print}' | cgat gff2bed | grep -v "^#" | \
awk '{print $1 "\t" $2 "\t" $3}' | sort -k 1,1 -k 2n,2 | \
bedtools merge > $OUTDIR/$SPEC_name\_cds_longesttrans.bed;

# Non-CDS regions (takes removed duplicate regions into account)
cat $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted.gtf | \
awk -F'\t' '$3 == "CDS" {print}' | cgat gff2bed | grep -v "^#" | \
awk '{print $1 "\t" $2 "\t" $3}' | sort -k 1,1 -k 2n,2 | \
bedtools merge > $OUTDIR/$SPEC_name\_cdstemp.bed;
bedtools complement -i $OUTDIR/$SPEC_name\_cdstemp.bed -g $OUTDIR/$SPEC_name\_genome.bed \
> $OUTDIR/$SPEC_name\_noncds.bed;
rm $OUTDIR/$SPEC_name\_cdstemp.bed;

# Exonic regions
cat $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted_rmdup.gtf | \
awk -F'\t' '$3 == "exon" {print}' | cgat gff2bed | grep -v "^#" | \
awk '{print $1 "\t" $2 "\t" $3}' | sort -k 1,1 -k 2n,2 | \
bedtools merge > $OUTDIR/$SPEC_name\_exonic.bed;

# Non-exonic regions (takes removed duplicate regions into account)
cat $OUTDIR/$SPEC_name\_$REF_name\_liftover_sorted_rmdup.gtf | \
awk -F'\t' '$3 == "exon" {print}' | cgat gff2bed | grep -v "^#" | \
awk '{print $1 "\t" $2 "\t" $3}' | sort -k 1,1 -k 2n,2 | \
bedtools merge > $OUTDIR/$SPEC_name\_exonictemp.bed;
bedtools complement -i $OUTDIR/$SPEC_name\_exonictemp.bed -g $OUTDIR/$SPEC_name\_genome.bed \
> $OUTDIR/$SPEC_name\_nonexonic.bed;
rm $OUTDIR/$SPEC_name\_exonictemp.bed;

### Remove temporary files
rm $OUTDIR/$REF_name.fa $OUTDIR/$SPEC_name.fa;
