#!/bin/bash -l

#SBATCH -A snic2022-5-484
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 50:00:00
#SBATCH -J coordinate_lift_over
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
PROJECT2=Skylark_2021
WORKDIR=$MAINDIR/data/$PROJECT;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
REF_name=Pmaj;
REF0=$MAINDIR/data/reference/Parus_major\
/GCF_001522545.3_Parus_major1.1_genomic.fasta;
SPEC_name=Aarv;
SPEC0=$WORKDIR2/B1_prepare_reference\
/GCA_902810485.1_skylark_genome_genomic_no_IUPAC.fasta;
ALIGNMENT=$WORKDIR2/B2_annotation_lift_over_parus/satsuma_Pmaj_Aarv/satsuma_summary.chained.out;
kraken=/crex/proj/snic2020-2-25/bin/kraken-accessed-2021-12-21/kraken/bin;

### Load modules
module load bioinfo-tools samtools/1.14 BEDTools/2.29.2;

### Activate conda environment
#source $conda/etc/profile.d/conda.sh;
source /sw/apps/bioinfo/CGAT/0.3.3/rackham/conda-install/etc/profile.d/conda.sh;
conda activate cgat-s;

### Create folders
mkdir $OUTDIR/coordinate_lift_over;
OUTDIR=$OUTDIR/coordinate_lift_over;

### Perform annotation lift-over
NAMES=$REF_name\_$SPEC_name;
cat $REF0 | awk '{if($0 ~ /^.*>/) print $1; else print;}' \
> $OUTDIR/$REF_name.fa;
REF=$OUTDIR/$REF_name.fa;
samtools faidx $REF;
cat $SPEC0 | awk '{if($0 ~ /^.*>/) print $1; else print;}' \
> $OUTDIR/$SPEC_name.fa;
SPEC=$OUTDIR/$SPEC_name.fa;

## Annotation lift-over with kraken

# Create config file
echo -e "[genomes]\n\
$REF_name\t$REF\n\
$SPEC_name\t$SPEC\n\
\n\
[pairwise-maps]\
\n$REF_name\t$SPEC_name\t\
$ALIGNMENT"\
> $OUTDIR/$NAMES.config;

cat $WORKDIR/A8_PhaseWY/final_output/stats/genome_summary/focal/$PROJECT\_Genome_summary.bed | \
grep -v "Missing data" | \
awk -F'\t' '{if($4=="Autosomal") print;\
else if($4=="Sex phase & depth difference") print $1"\t"$2"\t"$3"\tboth"; \
else if($4=="Sex depth difference") print $1"\t"$2"\t"$3"\thetgamdrop"; \
else if($4=="Sex phase difference") print $1"\t"$2"\t"$3"\tphase"}' \
> $OUTDIR/$PROJECT\_Genome_summary.bed
 # Convert bed file to gtf
cgat bed2gff --as-gtf -I  $OUTDIR/$PROJECT\_Genome_summary.bed -S $OUTDIR/$PROJECT\_Genome_summary.gtf;

# lift over
$kraken/RunKraken -c $OUTDIR/$NAMES.config -s $OUTDIR/$PROJECT\_Genome_summary.gtf -S $SPEC_name -T $REF_name -o $OUTDIR/$PROJECT\_$REF_name\_Genome_summary_all.gtf -m 100000000 -M 100000000 -a;
cat $OUTDIR/$PROJECT\_$REF_name\_Genome_summary_all.gtf | grep -v "Kraken_mapped \"FALSE\"" > $OUTDIR/$PROJECT\_$REF_name\_Genome_summary_mapped.gtf;

# Convert gtf to bed
cat $OUTDIR/$PROJECT\_$REF_name\_Genome_summary_mapped.gtf  | cut -f1,4,5,9 | \
awk -F'[\t"]' '{print $1"\t"$2-1"\t"$3"\t"$5 }' | \
bedtools sort -faidx $OUTDIR/$REF_name.fa.fai \
> $OUTDIR/$PROJECT\_$REF_name\_Genome_summary_mapped.bed;
cat $OUTDIR/$REF_name.fa.fai | cut -f1,2 > $OUTDIR/$REF_name.fai

bedtools complement -i $OUTDIR/$PROJECT\_$REF_name\_Genome_summary_mapped.bed -g $OUTDIR/$REF_name.fai | \
awk -F'\t' '{print $1"\t"$2"\t"$3"\tMissing data" }' \
> $OUTDIR/$PROJECT\_$REF_name\_missing_region.bed;
cat $OUTDIR/$PROJECT\_$REF_name\_Genome_summary_mapped.bed $OUTDIR/$PROJECT\_$REF_name\_missing_region.bed | \
bedtools sort -faidx $OUTDIR/$REF_name.fa.fai | \
awk -F'\t' '{print $1"\t"$2+1"\t"$3"\t"$4}' \
> $OUTDIR/$PROJECT\_$REF_name\_Genome_summary_final.txt;
