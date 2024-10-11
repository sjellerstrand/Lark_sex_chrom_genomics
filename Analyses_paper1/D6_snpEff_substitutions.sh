#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -J snpEff_substitutions
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT1=Skylark_2021;
PROJECT2=Rasolark_2021;
WORKDIR1=$MAINDIR/data/$PROJECT1;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/Simon/$PROJECT1;
GENES=$WORKDIR1/D4_align_shared_genes/metadata/genes_info.tsv;
FEATURES1=$WORKDIR1/D1_filter_genes/$PROJECT1\_features.tsv;
FEATURES2=$WORKDIR2/D1_filter_genes/$PROJECT2\_features.tsv;
DB_NAME=Tgut;
DB_REF=$MAINDIR/data/Skylark_2021/B3_annotation_lift_over/$DB_NAME.fa;
DB_GTF=$MAINDIR/data/reference/Taeniopygia_guttata/GCF_003957565.2_bTaeGut1.4.pri_genomic.no_W.gtf;
METADATA1=$WORKDIR1/metadata;
METADATA2=$WORKDIR2/metadata;
OUTGROUPS=$(echo -e \
"Pmaj \
${MAINDIR}/data/reference/Parus_major/GCF_001522545.3_Parus_major1.1_genomic.fasta \
${WORKDIR}/D1_filter_genes/Pmaj_uniqID.gtf
:\
Falb \
${MAINDIR}/data/reference/Ficedula_albicollis/GCA_000247815.2_FicAlb1.5_genomic.fasta \
${WORKDIR}/D1_filter_genes/Falb_uniqID.gtf");

### Load modules
module load bioinfo-tools SeqKit/0.15.0 bcftools/1.17 vcflib/1.0.1 BEDTools/2.29.2 snpEff/4.3t python/3.9.5 gnuparallel/20180822;
msa2vcf=/crex/proj/snic2020-2-25/bin/msa2vcf-accessed-2023-11-14/msa2vcf.py;
diploid2haploid=/crex/proj/snic2020-2-25/nobackup/simon/bin_tools/diploid2haploid.py;
snpeff_config=/sw/bioinfo/snpEff/4.3t/rackham/snpEff.config.orig;

### Activate conda environment
conda activate liftover;

### Create folders
mkdir $OUTDIR/D6_snpEff_substitutions;
OUTDIR=$OUTDIR/D6_snpEff_substitutions;
mkdir $OUTDIR/metadata \
$OUTDIR/shared_genes \
$OUTDIR/shared_genes/autosomal_autosomal \
$OUTDIR/shared_genes/sex_dropout_sex_dropout \
$OUTDIR/shared_genes/sex_phase_sex_phase \
$OUTDIR/shared_genes/sex_phase_sex_dropout \
$OUTDIR/shared_genes/sex_dropout_sex_phase \
$OUTDIR/shared_genes/sex_phase_autosomal \
$OUTDIR/merge_files \
$OUTDIR/merge_files/autosomal_autosomal \
$OUTDIR/merge_files/sex_dropout_sex_dropout \
$OUTDIR/merge_files/sex_phase_sex_phase \
$OUTDIR/merge_files/sex_phase_sex_dropout \
$OUTDIR/merge_files/sex_dropout_sex_phase \
$OUTDIR/merge_files/sex_phase_autosomal \
$OUTDIR/final_input \
$OUTDIR/data \
$OUTDIR/data/$DB_NAME \
$OUTDIR/data/genomes;

# Set up metadata files
OUTGROUPS=$(echo $OUTGROUPS | tr ':' '\n' | cut -d' ' -f1 | tr '\n' '|');
OUTGROUPS=${OUTGROUPS::-1};

cat $WORKDIR1/D4_align_shared_genes/metadata/female_A_A_sampleID.txt | grep -v $DB_NAME \
> $OUTDIR/metadata/female_A_A_sampleID.txt;
cat $WORKDIR1/D4_align_shared_genes/metadata/female_ZW_ZW_sampleID.txt | grep -v $DB_NAME \
> $OUTDIR/metadata/female_ZW_ZW_sampleID.txt;
cat $WORKDIR1/D4_align_shared_genes/metadata/female_ZW_Z_sampleID.txt | grep -v $DB_NAME \
> $OUTDIR/metadata/female_ZW_Z_sampleID.txt;
cat $WORKDIR1/D4_align_shared_genes/metadata/female_Z_ZW_sampleID.txt | grep -v $DB_NAME \
> $OUTDIR/metadata/female_Z_ZW_sampleID.txt;
cat $WORKDIR1/D4_align_shared_genes/metadata/female_Z_Z_sampleID.txt | grep -v $DB_NAME \
> $OUTDIR/metadata/female_Z_Z_sampleID.txt;
cat $WORKDIR1/D4_align_shared_genes/metadata/female_ZW_A_sampleID.txt | grep -v $DB_NAME \
> $OUTDIR/metadata/female_ZW_A_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_A_sampleID.txt | grep -v $DB_NAME | grep -vE $OUTGROUPS \
> $OUTDIR/metadata/$PROJECT1\_female_A_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_Z_sampleID.txt | grep -v $DB_NAME | grep -vE $OUTGROUPS \
> $OUTDIR/metadata/$PROJECT1\_female_Z_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_W_sampleID.txt | grep -v $DB_NAME | grep -vE $OUTGROUPS \
> $OUTDIR/metadata/$PROJECT1\_female_W_sampleID.txt;
cat $WORKDIR1/D4_align_shared_genes/metadata/female_A_A_sampleID.txt | \
grep -vFf $WORKDIR1/D3_codon_aware_alignment/metadata/female_A_sampleID.txt \
> $OUTDIR/metadata/$PROJECT2\_female_A_sampleID.txt;
cat $WORKDIR1/D4_align_shared_genes/metadata/female_Z_Z_sampleID.txt | \
grep -vFf $WORKDIR1/D3_codon_aware_alignment/metadata/female_Z_sampleID.txt \
> $OUTDIR/metadata/$PROJECT2\_female_Z_sampleID.txt;
cat $WORKDIR1/D4_align_shared_genes/metadata/female_Z_ZW_sampleID.txt | \
grep -vFf $WORKDIR1/D3_codon_aware_alignment/metadata/female_Z_sampleID.txt | \
grep -vFf $OUTDIR/metadata/$PROJECT2\_female_Z_sampleID.txt \
> $OUTDIR/metadata/$PROJECT2\_female_W_sampleID.txt;

### Set up gene and feature information
cat $GENES | awk -F'\t' '{if($10 > 0 && $12 > 0) print}' | grep $PROJECT1\_$PROJECT2 \
> $OUTDIR/metadata/genes_in.tsv;
GENES=$OUTDIR/metadata/genes_in.tsv;
cat $FEATURES1 | sort -k13 \
> $OUTDIR/metadata/features1.tsv;
cat $FEATURES2 | sort -k13 \
> $OUTDIR/metadata/features2.tsv;
join -j13 -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.9,1.10,2.10,1.11,2.11,1.12,2.12,1.13 \
$OUTDIR/metadata/features1.tsv $OUTDIR/metadata/features2.tsv | tr ' ' '\t' \
> $OUTDIR/metadata/features_shared.tsv;
FEATURES=$OUTDIR/metadata/features_shared.tsv;
rm $OUTDIR/metadata/features1.tsv \
$OUTDIR/metadata/features2.tsv;

export WORKDIR1 WORKDIR2 OUTDIR PROJECT1 PROJECT2 DB_NAME METADATA1 METADATA2 msa2vcf diploid2haploid \
GENES FEATURES OUTGROUPS;

### Convert sequences to database species reference
shared_to_DB_ref() {

  GENE=$1;
  TRANS=$2;
  REGION0_1=$3;
  REGION0_2=$4;

  if [ $REGION0_1 == "sex_phase" ] || [ $REGION0_1 == "partial_dropout" ]; then
    REGION1="sex_phase";
  else
    REGION1=$REGION0_1;
  fi;
  if [ $REGION0_2 == "sex_phase" ] || [ $REGION0_2 == "partial_dropout" ]; then
    REGION2="sex_phase";
  else
    REGION2=$REGION0_2;
  fi;

  # Choose regions
  if [ $REGION1 == "autosomal" ] && [ $REGION2 == "autosomal" ]; then
    SUFFIX="A_A";
  elif [ $REGION1 == "sex_dropout" ] && [ $REGION2 == "sex_dropout" ]; then
    SUFFIX="Z_Z";
  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "sex_phase" ]; then
    SUFFIX="ZW_ZW";
  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "sex_dropout" ]; then
    SUFFIX="ZW_Z";
  elif [ $REGION1 == "sex_dropout" ] && [ $REGION2 == "sex_phase" ]; then
    SUFFIX="Z_ZW";
  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "autosomal" ]; then
    SUFFIX="ZW_A";
  fi;

  echo $GENE;
  mkdir $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE;

  # Find gene CDS features
  cat $FEATURES | grep $GENE\_$TRANS \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list_temp.tsv;

  # Sort CDS sequences by exon order
  for EXON in $(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list_temp.tsv | cut -f17); do
    echo $(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list_temp.tsv | grep $EXON) $(echo $EXON | rev | cut -d'_' -f2 | rev);
  done | tr ' ' '\t' | sort -nk18,18 | cut -f1-17 > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv;
  rm $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list_temp.tsv;

  # Look for insertions compared to database species
  cat $WORKDIR1/D4_align_shared_genes/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females_$SUFFIX.fasta \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females.fasta;
  cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females.fasta | \
  head -n2 | tail -n1 | grep -ob "[-\!]" | cut -f1 | awk '{print $0+1}' \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_insertions.txt;

  # Seperate gene into CDS sequences
  POS_START=1;
  for EXON_NR in $(seq 1 1 $(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | wc -l)); do
    DB_SCAFFOLD=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f1);
    DB_START=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f2);
    DB_END=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f3);
    DB_STRANDEDNESS=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f4);
    MISSING1=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f9);
    MISSING2=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f10)
    EXON=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f17);

    # Find end of CDS sequence
    POS_END=$(echo $POS_START + $DB_END - $DB_START | bc);

    # Look for insertions relative to database species and extend if necessary to end position
    X=$POS_START;
    Y=$POS_END;
    Z=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_insertions.txt | \
    awk -v X=$X -v Y=$Y 'BEGIN{sum=0} {if($1>=X && $1<=Y) {sum+=1}} END {print sum}');
    while [ $Z -gt 0 ]; do
      X=$(echo $Y + 1 | bc);
      Y=$(echo $Y + $Z | bc);
      Z=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_insertions.txt | \
      awk -v X=$X -v Y=$Y 'BEGIN{sum=0} {if($1>=X && $1<=Y) {sum+=1}} END {print sum}');
    done;
    POS_END=$Y;

    if [ $MISSING1 == "callable" ] && [ $MISSING2 == "callable" ]; then

      mkdir $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON;
      cd $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON;

      # Extract CDS sequence
      if [ $DB_STRANDEDNESS == "+" ]; then

        # Extract CDS sequence and add dummy bases in start of alignment
        seqkit subseq -r $POS_START:$POS_END \
        $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females.fasta | \
        awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | \
        awk '{if($0 ~ /^>/) print; else print "XXXXXXXXXXXXXXXXXXXX"$0}' | tail -n+2 | tr '\!' '-' \
        > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON\_alignment_NT_snpeff_females_leadingdummy.fasta;

      else

        # Extract CDS sequence, add dummy bases in start of alignment
        seqkit subseq -r $POS_START:$POS_END \
        $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females.fasta | \
        awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | \
        while read alignment; do
          if [[ $alignment =~ ^'>' ]]; then
            echo $alignment;
          else
            echo $alignment | tr ACGTacgt TGCAtgca | rev;
          fi;
        done | awk '{if($0 ~ /^>/) print; else print "XXXXXXXXXXXXXXXXXXXX"$0}' | tail -n+2 | tr '\!' '-' \
        > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON\_alignment_NT_snpeff_females_leadingdummy.fasta;

      fi;

      # Convert alignment to vcf
      python $msa2vcf $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON\_alignment_NT_snpeff_females_leadingdummy.fasta | \
      awk -F'\t' -v DB_SCAFFOLD=$DB_SCAFFOLD -v DB_START=$DB_START \
      '{OFS="\t"; if($0 ~ /^##contig/) print "##contig=<ID="DB_SCAFFOLD">"; else if($0 ~ /^#/) print; \
      else if($2<=20) next; else {$1=DB_SCAFFOLD; $2=$2+DB_START-20-1; print}}' | \
      bcftools view -s ^$DB_NAME | \
      bcftools view -S $OUTDIR/metadata/female_$SUFFIX\_sampleID.txt | \
      bgzip -c > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON.vcf.gz;

      # Remove temporary files
      rm $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON\_alignment_NT_snpeff_females_leadingdummy.fasta \
      $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/snptable.txt;

      cd -;

   fi;

    # Update start position for next CDS sequence
    POS_START=$(echo $POS_END + 1 | bc);

  done;

  # Merge CDS back to gene
  find $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE -name "*.vcf.gz" \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_vcf_list.txt;
  bcftools concat -f $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_vcf_list.txt -Ou | \
  bcftools sort -T $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/temp_merge | \
  vcfallelicprimitives --keep-info --keep-geno | \
  awk -F $'\t' 'BEGIN {OFS = FS} /^[#]/ {print; next} {for (i=10; i<=NF; i++) { gsub("\\|","/",$i)} print}' | \
  bcftools norm -m -any --multi-overlaps "." | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  awk -F'\t' '{if($0 ~ /^#/) {print} else if($0 !~ /^#/ && $5 !~ /N/) {print}}' | \
  bgzip -c > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp1.vcf.gz;
  tabix $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp1.vcf.gz;

  # Filter alternate alleles present in outgroup lineages
  bcftools view $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp1.vcf.gz \
  -s$(echo $OUTGROUPS | tr '|' ',') | \
  vcffixup - | \
  vcffilter -f "AC > 0 " | \
  bgzip -c > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp_outgroup.vcf.gz;
  tabix $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp_outgroup.vcf.gz;
  bcftools isec -C $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp1.vcf.gz \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp_outgroup.vcf.gz -w 1 | \
  bcftools view -s^$(echo $OUTGROUPS | tr '|' ',') | \
  bcftools norm -m +any | \
  vcffixup - | \
  bgzip -c > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp2.vcf.gz;

  ## Convert genotypes to haploid
  python $diploid2haploid \
  -i $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp2.vcf.gz \
  -o $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE.vcf;
  bgzip $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE.vcf;
  tabix $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE.vcf.gz;

  # Remove temporary files
  rm -r $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_* \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females.fasta* \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp*;

};

## Excecute function in parallell
export -f shared_to_DB_ref;
parallel --colsep '\t' 'shared_to_DB_ref {}' :::: $GENES;

## Merge all genes across regions
for REGION1 in $(echo autosomal sex_phase sex_dropout); do

  if [ $REGION1 == "autosomal" ]; then

    REGION2=autosomal;
    find $OUTDIR/shared_genes/$REGION1\_$REGION2 -name "*.vcf.gz" \
    > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt;
    bcftools concat -a -f $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt | \
    bcftools sort --temp-dir $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp | \
    bcftools norm --rm-dup all | \
    vcffixup - | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;
    tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;

    bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -h \
    > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
    bedtools sort -i $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -faidx $DB_REF.fai \
    >> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
    bgzip $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
    tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz;

    # Create lineage files with lineage specific substitutions
    for PROJECT in $(echo $PROJECT1 $PROJECT2); do

        if [ $PROJECT == $PROJECT1 ]; then
          REGION=$REGION1;
        else
          REGION=$REGION2;
        fi;

      bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
      -S $OUTDIR/metadata/$PROJECT\_female_A_sampleID.txt | \
      vcfallelicprimitives --keep-info --keep-geno | \
      vcffixup - | \
      vcffilter -f "( AC / NS ) = 1" | \
      bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_A.vcf.gz;
      tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_A.vcf.gz;

    done;

  else

    for REGION2 in $(echo sex_phase sex_dropout); do

      find $OUTDIR/shared_genes/$REGION1\_$REGION2 -name "*.vcf.gz" \
      > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt;
      bcftools concat -a -f $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt | \
      bcftools sort --temp-dir $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp | \
      bcftools norm --rm-dup all | \
      vcffixup - | \
      vcffilter -f "AC > 0" | \
      bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;
      tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;

      bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -h \
      > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
      bedtools sort -i $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -faidx $DB_REF.fai \
      >> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
      bgzip $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
      tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz;

      # Create lineage files with lineage specific substitutions
      for PROJECT in $(echo $PROJECT1 $PROJECT2); do

        if [ $PROJECT == $PROJECT1 ]; then
            REGION=$REGION1;
        else
            REGION=$REGION2;
        fi;

        bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
        -S $OUTDIR/metadata/$PROJECT\_female_Z_sampleID.txt | \
        vcfallelicprimitives --keep-info --keep-geno | \
        vcffixup - | \
        vcffilter -f "( AC / NS ) = 1" | \
        bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_Z.vcf.gz;
        tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_Z.vcf.gz;

        if [ $REGION == "sex_phase" ]; then

          bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
          -S $OUTDIR/metadata/$PROJECT\_female_W_sampleID.txt | \
          vcfallelicprimitives --keep-info --keep-geno | \
          vcffixup - | \
          vcffilter -f "( AC / NS ) = 1" | \
          bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_W.vcf.gz;
          tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_W.vcf.gz;

        fi;

      done;

    done;

  fi;

done;

REGION1=sex_phase;
REGION2=autosomal;

find $OUTDIR/shared_genes/$REGION1\_$REGION2 -name "*.vcf.gz" \
> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt;
bcftools concat -a -f $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt | \
bcftools sort --temp-dir $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp | \
bcftools norm --rm-dup all | \
vcffixup - | \
vcffilter -f "AC > 0" | \
bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;

bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -h \
> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
bedtools sort -i $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -faidx $DB_REF.fai \
>> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
bgzip $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz;

bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
-S $OUTDIR/metadata/$PROJECT1\_female_Z_sampleID.txt | \
vcfallelicprimitives --keep-info --keep-geno | \
vcffixup - | \
vcffilter -f "( AC / NS ) = 1" | \
bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT1\_$REGION1\_Z.vcf.gz;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT1\_$REGION1\_Z.vcf.gz;

bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
-S $OUTDIR/metadata/$PROJECT1\_female_W_sampleID.txt | \
vcfallelicprimitives --keep-info --keep-geno | \
vcffixup - | \
vcffilter -f "( AC / NS ) = 1" | \
bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT1\_$REGION1\_W.vcf.gz;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT1\_$REGION1\_W.vcf.gz;

bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
-S $OUTDIR/metadata/$PROJECT2\_female_A_sampleID.txt | \
vcfallelicprimitives --keep-info --keep-geno | \
vcffixup - | \
vcffilter -f "( AC / NS ) = 1" | \
bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT2\_$REGION2\_A.vcf.gz;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT2\_$REGION2\_A.vcf.gz;

### Create database
cp $snpeff_config $OUTDIR/metadata/snpEff.config;
echo $DB_NAME genome, version $DB_NAME >> $OUTDIR/metadata/snpEff.config;
echo $DB_NAME.genome : $DB_NAME >> $OUTDIR/metadata/snpEff.config;

# Create transcripts only gtf
cat $GENES | cut -f2 > $OUTDIR/metadata/$PROJECT1\_$PROJECT2\_transcripts.txt;
cgat gtf2gtf --method=filter --filter-method=transcript --map-tsv-file $OUTDIR/metadata/$PROJECT1\_$PROJECT2\_transcripts.txt -I $DB_GTF \
-S $OUTDIR/data/$DB_NAME/genes.gtf;

# Copy files to the correct folders and names
cp $DB_REF $OUTDIR/data/genomes/$DB_NAME.fa;
cd $OUTDIR;

# Build database
snpEff build -gtf22 -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -treatAllAsProteinCoding -v $DB_NAME;

# Merge lineage specific files and annotate fixed substitutions
for REGION in $(echo A Z W); do

  for PROJECT in $(echo $PROJECT1 $PROJECT2); do

    find $OUTDIR/merge_files/ -name "$PROJECT\_*_$REGION.vcf.gz" \
    > $OUTDIR/final_input/$PROJECT\_$REGION\_vcf_list.txt;

    bcftools concat -a -f $OUTDIR/final_input/$PROJECT\_$REGION\_vcf_list.txt | \
    bcftools sort --temp-dir $OUTDIR/final_input/$PROJECT\_$REGION\_temp | \
    bcftools norm --rm-dup all | \
    vcfallelicprimitives --keep-info --keep-geno | \
    vcffixup - | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIR/final_input/$PROJECT\_$REGION\_temp1.vcf.gz;
    tabix $OUTDIR/final_input/$PROJECT\_$REGION\_temp1.vcf.gz;

    bcftools view $OUTDIR/final_input/$PROJECT\_$REGION\_temp1.vcf.gz -h \
    > $OUTDIR/final_input/$PROJECT\_$REGION.vcf;
    bedtools sort -i $OUTDIR/final_input/$PROJECT\_$REGION\_temp1.vcf.gz -faidx $DB_REF.fai \
    >> $OUTDIR/final_input/$PROJECT\_$REGION.vcf;
    bgzip $OUTDIR/final_input/$PROJECT\_$REGION.vcf;
    tabix $OUTDIR/final_input/$PROJECT\_$REGION.vcf.gz;

    mkdir $OUTDIR/$PROJECT\_$REGION;

    snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION.html \
    -csvStats $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
    $OUTDIR/final_input/$PROJECT\_$REGION.vcf.gz \
    > $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION\_impacts.vcf;
    bgzip $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION\_impacts.vcf;
    tabix $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION\_impacts.vcf.gz;

    # Extract information from impact annotation
    bcftools view $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION\_impacts.vcf.gz -H | \
    cut -f8 | cut -d'|' -f4,16 | cut -d',' -f1 | cut -d';' -f1 | grep WARNING | tr '|' '\t' | sort -u \
    > $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION\_warnings.tsv;
    cat $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION.genes.txt | \
    awk -F'\t' -v PROJECT=$PROJECT -v REGION=$REGION 'BEGIN {COL=0} {if(NR==2) {for(i=1; i<=NF; i++) {if($i=="variants_impact_HIGH") COL=i}} \
    else if(COL > 0) {print $1"\t"$3"\t"$COL"\t"PROJECT"_"REGION}}' \
    > $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION\_high_impacts.tsv;
    bcftools view $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION\_impacts.vcf.gz -H | \
    cut -f8 | awk -F'LOF=\\(' '{print $2}' | awk -F ')' '{print $1}' | \
    cut -d'|' -f1 | tr -s '\n' | sort -u | awk -F'\t' -v PROJECT=$PROJECT -v REGION=$REGION '{print $1"\t"PROJECT"_"REGION}' \
    > $OUTDIR/$PROJECT\_$REGION/$PROJECT\_$REGION\_LOF.tsv;

    # Remove temporary files
    rm $OUTDIR/final_input/$PROJECT\_$REGION\_vcf_list.txt \
    $OUTDIR/final_input/$PROJECT\_$REGION\_temp*;

  done;

done;

### Annotate sex-linked substitutions
mkdir $OUTDIR/$PROJECT1\_Z_sexlinked \
$OUTDIR/$PROJECT1\_W_sexlinked \
$OUTDIR/$PROJECT2\_Z_sexlinked \
$OUTDIR/$PROJECT2\_W_sexlinked;

### Project1 Z
bcftools view $OUTDIR/final_input/$PROJECT1\_Z.vcf.gz \
-T ^$OUTDIR/final_input/$PROJECT1\_W.vcf.gz | \
bcftools view -T ^$OUTDIR/final_input/$PROJECT2\_W.vcf.gz | \
bgzip -c > $OUTDIR/final_input/$PROJECT1\_Z_sexlinked.vcf.gz;
tabix $OUTDIR/final_input/$PROJECT1\_Z_sexlinked.vcf.gz;

snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked.html \
-csvStats $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
$OUTDIR/final_input/$PROJECT1\_Z_sexlinked.vcf.gz \
> $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked_impacts.vcf;
bgzip $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked_impacts.vcf;
tabix $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked_impacts.vcf.gz;

# Extract information from impact annotation
bcftools view $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked_impacts.vcf.gz -H | \
cut -f8 | cut -d'|' -f4,16 | cut -d',' -f1 | cut -d';' -f1 | grep WARNING | tr '|' '\t' | sort -u \
> $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked_warnings.tsv;
cat $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked.genes.txt | \
awk -F'\t' -v PROJECT=$PROJECT1 'BEGIN {COL=0} {if(NR==2) {for(i=1; i<=NF; i++) {if($i=="variants_impact_HIGH") COL=i}} \
else if(COL > 0) {print $1"\t"$3"\t"$COL"\t"PROJECT"_Z_sexlinked"}}' \
> $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked_high_impacts.tsv;
bcftools view $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked_impacts.vcf.gz -H | \
cut -f8 | awk -F'LOF=\\(' '{print $2}' | awk -F ')' '{print $1}' | \
cut -d'|' -f1 | tr -s '\n' | sort -u | awk -F'\t' -v PROJECT=$PROJECT1 '{print $1"\t"PROJECT"_Z_sexlinked"}' \
> $OUTDIR/$PROJECT1\_Z_sexlinked/$PROJECT1\_Z_sexlinked_LOF.tsv;

### Project1 W
bcftools view $OUTDIR/final_input/$PROJECT1\_W.vcf.gz \
-T ^$OUTDIR/final_input/$PROJECT1\_Z.vcf.gz | \
bcftools view -T ^$OUTDIR/final_input/$PROJECT2\_Z.vcf.gz | \
bgzip -c > $OUTDIR/final_input/$PROJECT1\_W_sexlinked.vcf.gz;
tabix $OUTDIR/final_input/$PROJECT1\_W_sexlinked.vcf.gz;

snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked.html \
-csvStats $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
$OUTDIR/final_input/$PROJECT1\_W_sexlinked.vcf.gz \
> $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked_impacts.vcf;
bgzip $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked_impacts.vcf;
tabix $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked_impacts.vcf.gz;

# Extract information from impact annotation
bcftools view $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked_impacts.vcf.gz -H | \
cut -f8 | cut -d'|' -f4,16 | cut -d',' -f1 | cut -d';' -f1 | grep WARNING | tr '|' '\t' | sort -u \
> $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked_warnings.tsv;
cat $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked.genes.txt | \
awk -F'\t' -v PROJECT=$PROJECT1 'BEGIN {COL=0} {if(NR==2) {for(i=1; i<=NF; i++) {if($i=="variants_impact_HIGH") COL=i}} \
else if(COL > 0) {print $1"\t"$3"\t"$COL"\t"PROJECT"_W_sexlinked"}}' \
> $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked_high_impacts.tsv;
bcftools view $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked_impacts.vcf.gz -H | \
cut -f8 | awk -F'LOF=\\(' '{print $2}' | awk -F ')' '{print $1}' | \
cut -d'|' -f1 | tr -s '\n' | sort -u | awk -F'\t' -v PROJECT=$PROJECT1 '{print $1"\t"PROJECT"_W_sexlinked"}' \
> $OUTDIR/$PROJECT1\_W_sexlinked/$PROJECT1\_W_sexlinked_LOF.tsv;

### Project2 Z
bcftools view $OUTDIR/final_input/$PROJECT2\_Z.vcf.gz \
-T ^$OUTDIR/final_input/$PROJECT2\_W.vcf.gz | \
bcftools view -T ^$OUTDIR/final_input/$PROJECT1\_W.vcf.gz | \
bgzip -c > $OUTDIR/final_input/$PROJECT2\_Z_sexlinked.vcf.gz;
tabix $OUTDIR/final_input/$PROJECT2\_Z_sexlinked.vcf.gz;

snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked.html \
-csvStats $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
$OUTDIR/final_input/$PROJECT2\_Z_sexlinked.vcf.gz \
> $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked_impacts.vcf;
bgzip $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked_impacts.vcf;
tabix $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked_impacts.vcf.gz;

# Extract information from impact annotation
bcftools view $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked_impacts.vcf.gz -H | \
cut -f8 | cut -d'|' -f4,16 | cut -d',' -f1 | cut -d';' -f1 | grep WARNING | tr '|' '\t' | sort -u \
> $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked_warnings.tsv;
cat $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked.genes.txt | \
awk -F'\t' -v PROJECT=$PROJECT2 'BEGIN {COL=0} {if(NR==2) {for(i=1; i<=NF; i++) {if($i=="variants_impact_HIGH") COL=i}} \
else if(COL > 0) {print $1"\t"$3"\t"$COL"\t"PROJECT"_Z_sexlinked"}}' \
> $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked_high_impacts.tsv;
bcftools view $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked_impacts.vcf.gz -H | \
cut -f8 | awk -F'LOF=\\(' '{print $2}' | awk -F ')' '{print $1}' | \
cut -d'|' -f1 | tr -s '\n' | sort -u | awk -F'\t' -v PROJECT=$PROJECT2 '{print $1"\t"PROJECT"_Z_sexlinked"}' \
> $OUTDIR/$PROJECT2\_Z_sexlinked/$PROJECT2\_Z_sexlinked_LOF.tsv;

### Project2 W
bcftools view $OUTDIR/final_input/$PROJECT2\_W.vcf.gz \
-T ^$OUTDIR/final_input/$PROJECT2\_Z.vcf.gz | \
bcftools view -T ^$OUTDIR/final_input/$PROJECT1\_Z.vcf.gz | \
bgzip -c > $OUTDIR/final_input/$PROJECT2\_W_sexlinked.vcf.gz;
tabix $OUTDIR/final_input/$PROJECT2\_W_sexlinked.vcf.gz;

snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked.html \
-csvStats $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
$OUTDIR/final_input/$PROJECT2\_W_sexlinked.vcf.gz \
> $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked_impacts.vcf;
bgzip $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked_impacts.vcf;
tabix $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked_impacts.vcf.gz;

# Extract information from impact annotation
bcftools view $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked_impacts.vcf.gz -H | \
cut -f8 | cut -d'|' -f4,16 | cut -d',' -f1 | cut -d';' -f1 | grep WARNING | tr '|' '\t' | sort -u \
> $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked_warnings.tsv;
cat $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked.genes.txt | \
awk -F'\t' -v PROJECT=$PROJECT2 'BEGIN {COL=0} {if(NR==2) {for(i=1; i<=NF; i++) {if($i=="variants_impact_HIGH") COL=i}} \
else if(COL > 0) {print $1"\t"$3"\t"$COL"\t"PROJECT"_W_sexlinked"}}' \
> $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked_high_impacts.tsv;
bcftools view $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked_impacts.vcf.gz -H | \
cut -f8 | awk -F'LOF=\\(' '{print $2}' | awk -F ')' '{print $1}' | \
cut -d'|' -f1 | tr -s '\n' | sort -u | awk -F'\t' -v PROJECT=$PROJECT2 '{print $1"\t"PROJECT"_W_sexlinked"}' \
> $OUTDIR/$PROJECT2\_W_sexlinked/$PROJECT2\_W_sexlinked_LOF.tsv;

# Remove temporary files
rm -r $OUTDIR/shared_genes \
$OUTDIR/merge_files \
$OUTDIR/final_input;

### Concatenate impact information
cat $(find $OUTDIR/ -name "*_warnings.tsv") | sort -u \
> $OUTDIR/metadata/annotation_warnings.tsv;
cat $(find $OUTDIR/ -name "*_high_impacts.tsv") | sort \
> $OUTDIR/metadata/annotation_high_impacts.tsv;
cat $(find $OUTDIR/ -name "*_LOF.tsv") | sort \
> $OUTDIR/metadata/annotation_LOF.tsv;

