#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J genome_scans_neutral
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

# Set parameters
SEX_SYSTEM=ZW;
NONDIPLOIDS=No;

## Fixed number of callable sites for windows
WINDOW=10000;
STEP=$WINDOW;

## Distance from exons
EXON_DIST=100000;

## Minimum callable sites per genome type (autosomal, homogametic, heterogametic) for scaffold to be analysed
MIN_SCAFFOLD_SIZE=1000000;

MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
PROJECT2=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
WORKDIR2=$MAINDIR/data/$PROJECT2;
VCFS=$WORKDIR/A8_PhaseWY/final_output/vcfs;
BEDS=$WORKDIR/A9_define_regions;
EXONS=$WORKDIR2/B3_annotation_lift_over/*_exonic.bed;
OUTDIR0=$MAINDIR/working/Simon/$PROJECT;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
METADATA=$WORKDIR/metadata;
FUNCTIONS=$MAINDIR/scripts/$PROJECT/analyses;
genomics_general=/crex/proj/snic2020-2-25/bin/genomics_general-accessed-2021-12-09;

# Load modules
module load bioinfo-tools bcftools/1.17 vcflib/1.0.1 BEDTools/2.29.2 python/3.9.5 R_packages/4.0.0;

# Define functions
genome_scans_windows=$FUNCTIONS/genome_scans_windows.r;
merge_scan_data=$FUNCTIONS/merge_scan_data.r;

# Create folders
mkdir $OUTDIR0/E2a_genome_scans_neutral;
OUTDIR0=$OUTDIR0/E2a_genome_scans_neutral;
mkdir $OUTDIR0/Reuseable_files \
$OUTDIR0/windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST \
$OUTDIR0/windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST/windows;
OUTDIR1=$OUTDIR0/windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST;

# Perform window based scans across the genome

## Determine datasets to run
if [ $SEX_SYSTEM != NA ]; then
  DATASETS=$(echo autosomal homogametic heterogametic)
else
  DATASETS=autosomal;
fi;

## Run through datasets
for DATA in $DATASETS; do

  ### Create subfolders
  mkdir $OUTDIR1/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST;
  OUTDIR2=$OUTDIR1/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST;

  ### Prepare BED files
  bedtools slop -i $EXONS -g $REF.fai -b $EXON_DIST | bedtools merge \
  > $OUTDIR0/Reuseable_files/exons_with_flank_$EXON_DIST.bed;
  bedtools subtract -a $BEDS/$PROJECT\_$DATA.bed \
  -b $OUTDIR0/Reuseable_files/exons_with_flank_$EXON_DIST.bed \
  > $OUTDIR0/Reuseable_files/$DATA\_callable_$EXON_DIST.bed;
  MASK=$OUTDIR0/Reuseable_files/$DATA\_callable_$EXON_DIST.bed;

  if [ $(find $OUTDIR1 | grep "$OUTDIR1/windows/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST\_infile.txt" | wc -l) == 0 ]; then

    ### Calculate windows with consideration to missing data
    Rscript $genome_scans_windows --args WINDOW=$WINDOW STEP=$STEP MIN_SCAFFOLD_SIZE=$MIN_SCAFFOLD_SIZE \
    EXON_DIST=$EXON_DIST PROJECT=$PROJECT MASK=$MASK REF=$REF OUTDIR1=$OUTDIR1 OUTDIR2=$OUTDIR2 DATA=$DATA;
    cat $OUTDIR1/windows/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST.txt | tail -n+2 | cut -f1,2,3 \
    > $OUTDIR1/windows/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST\_infile.txt;

  fi;

  ### Calculate populations statistics

  #### Get ploidy info
  if [ $DATA == "autosomal" ]; then

    if [ $NONDIPLOIDS == "Yes" ]; then
      # Get ploidy file from metadata folder
      ploidy=$METADATA/ploidy.txt;
    else
      # Set all individuals as diploid
      cat $METADATA/sample_info.txt | cut -f1 | tail -n +2 | awk -F'\t' '{print $1"\t2"}' \
      > $OUTDIR2/ploidy.txt;
      ploidy=$OUTDIR2/ploidy.txt;
    fi;
  else
    if [ $SEX_SYSTEM != NA ]; then

      if [ $SEX_SYSTEM == ZW ]; then
        HOMGAM=Male;
        HETGAM=Female;
      elif [ $SEX_SYSTEM == XY ]; then
        HOMGAM=Female;
        HETGAM=Male;
      fi;
      if [ $DATA == "homogametic" ]; then
        ### Set all homogametes as diploids and all heterogametes as haploid
        cat $METADATA/sample_info.txt | awk -F'\t' -v HETGAM=$HETGAM '$3==HETGAM {print $1"\t1"}' \
        > $OUTDIR2/ploidy1.txt;
        cat $METADATA/sample_info.txt | awk -F'\t' -v HOMGAM=$HOMGAM '$3==HOMGAM {print $1"\t2"}' \
        > $OUTDIR2/ploidy2.txt;
        cat $OUTDIR2/ploidy1.txt $OUTDIR2/ploidy2.txt > $OUTDIR2/ploidy.txt;
        rm $OUTDIR2/ploidy1.txt $OUTDIR2/ploidy2.txt;
        ploidy=$OUTDIR2/ploidy.txt;

      elif [ $DATA == "heterogametic" ]; then
        ### Set all heterogametes as haploid
        cat $METADATA/sample_info.txt | awk -F'\t' -v HETGAM=$HETGAM '$3==HETGAM {print $1"\t1"}' \
        > $OUTDIR2/ploidy.txt;
        ploidy=$OUTDIR2/ploidy.txt;
      fi;
    fi;
  fi;

  #### Convert input vcf to correct format
  if [ $(find $OUTDIR0 | grep "$OUTDIR0/Reuseable_files/$PROJECT\_$DATA.geno.gz" | wc -l) == 0 ]; then
    bcftools view $VCFS/$PROJECT\_$DATA.vcf.gz | \
    vcfclassify - | \
    vcffilter -s -f "!( INS | DEL | MNP )" -f "AC > 0 & AF < 1" | \
    bgzip -c > $OUTDIR0/Reuseable_files/$PROJECT\_$DATA\_snps.vcf.gz;
    python3 $genomics_general/VCF_processing/parseVCF.py -i $OUTDIR0/Reuseable_files/$PROJECT\_$DATA\_snps.vcf.gz \
    --ploidyFile $ploidy | bgzip -c > $OUTDIR0/Reuseable_files/$PROJECT\_$DATA.geno.gz;
    rm $OUTDIR0/Reuseable_files/$PROJECT\_$DATA\_snps.vcf.gz;
  fi;

  #### Set window parameters
  PARAMETERS=$(echo --windType predefined --windCoords $OUTDIR1/windows/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST\_infile.txt);

  #### Calculate Pi and Tajima's D
  python3 $genomics_general/popgenWindows.py -g $OUTDIR0/Reuseable_files/$PROJECT\_$DATA.geno.gz \
  -f phased --writeFailedWindows --ploidyFile $ploidy -o $OUTDIR2/$PROJECT\_$DATA.pi_tajimas_D  \
  --analysis popDist popFreq $PARAMETERS;

  #### Merge data
  Rscript $merge_scan_data --args WINDOW=$WINDOW STEP=$STEP \
  PROJECT=$PROJECT REF=$REF OUTDIR1=$OUTDIR1 OUTDIR2=$OUTDIR2 DATA=$DATA EXON_DIST=$EXON_DIST;

done;
