#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -C mem256GB
#SBATCH -J align_shared_genes
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT1=Skylark_2021;
PROJECT2=Rasolark_2021;
WORKDIR1=$MAINDIR/data/$PROJECT1;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/Simon/$PROJECT1;
GENES1=$WORKDIR1/D3_codon_aware_alignment/metadata/genes_info.tsv
GENES2=$WORKDIR2/D3_codon_aware_alignment/metadata/genes_info.tsv;
FEATURES1=$WORKDIR1/D1_filter_genes/$PROJECT1\_features.tsv;
FEATURES2=$WORKDIR2/D1_filter_genes/$PROJECT2\_features.tsv;
DB_NAME=Tgut;
METADATA1=$WORKDIR1/metadata;
METADATA2=$WORKDIR2/metadata;
OUTGROUPS=$(echo -e \
"Pmaj \
${MAINDIR}/data/reference/Parus_major/GCF_001522545.3_Parus_major1.1_genomic.fasta \
${WORKDIR1}/D1_filter_genes/Pmaj_uniqID.gtf
:\
Falb \
${MAINDIR}/data/reference/Ficedula_albicollis/GCA_000247815.2_FicAlb1.5_genomic.fasta \
${WORKDIR1}/D1_filter_genes/Falb_uniqID.gtf");

### Load modules
module load bioinfo-tools SeqKit/0.15.0 BEDTools/2.29.2 gnuparallel/20180822;

### Activate conda environment
conda activate gene_alignment;

### Create folders
mkdir $OUTDIR/D4_align_shared_genes;
OUTDIR=$OUTDIR/D4_align_shared_genes;
OUTDIR=$(echo /$(echo $OUTDIR | cut -d'/' -f 3-));
mkdir $OUTDIR/metadata \
$OUTDIR/shared_genes \
$OUTDIR/shared_genes/autosomal_autosomal \
$OUTDIR/shared_genes/sex_dropout_sex_dropout \
$OUTDIR/shared_genes/sex_phase_sex_phase \
$OUTDIR/shared_genes/sex_phase_sex_dropout \
$OUTDIR/shared_genes/sex_dropout_sex_phase \
$OUTDIR/shared_genes/sex_phase_autosomal;

# Set up metadata files
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_A_sampleID.txt \
$WORKDIR2/D3_codon_aware_alignment/metadata/female_A_sampleID.txt | sort -u \
> $OUTDIR/metadata/female_A_A_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_Z_sampleID.txt \
$WORKDIR2/D3_codon_aware_alignment/metadata/female_Z_sampleID.txt | sort -u \
> $OUTDIR/metadata/female_Z_Z_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_sexchrom_sampleID.txt \
$WORKDIR2/D3_codon_aware_alignment/metadata/female_sexchrom_sampleID.txt | sort -u \
> $OUTDIR/metadata/female_ZW_ZW_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_sexchrom_sampleID.txt \
$WORKDIR2/D3_codon_aware_alignment/metadata/female_Z_sampleID.txt | sort -u \
> $OUTDIR/metadata/female_ZW_Z_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_Z_sampleID.txt \
$WORKDIR2/D3_codon_aware_alignment/metadata/female_sexchrom_sampleID.txt | sort -u \
> $OUTDIR/metadata/female_Z_ZW_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_sexchrom_sampleID.txt \
$WORKDIR2/D3_codon_aware_alignment/metadata/female_A_sampleID.txt | sort -u \
> $OUTDIR/metadata/female_ZW_A_sampleID.txt;

# All sequences
NB_SEQ1=$(echo "($(cat $METADATA1/sample_info.txt $METADATA2/sample_info.txt | wc -l) * 2) \
- (2*2) + 1 + $(echo $OUTGROUPS | tr ':' '\n' | wc -l)" | bc);

# Second species dropout
NB_SEQ2=$(echo "$NB_SEQ1 - $(cat $METADATA2/sample_info.txt | \
awk -F'\t' '{if($3=="Female") print $1}' | wc -l)" | bc);

# First species dropout
NB_SEQ3=$(echo "$NB_SEQ1 - $(cat $METADATA1/sample_info.txt | \
awk -F'\t' '{if($3=="Female") print $1}' | wc -l)" | bc);

# Both species dropout
NB_SEQ4=$(echo "$NB_SEQ3 - $(cat $METADATA2/sample_info.txt | \
awk -F'\t' '{if($3=="Female") print $1}' | wc -l)" | bc);

### Find shared genes and transcripts
cat $GENES1 | cut -f2 > $OUTDIR/metadata/$PROJECT1\_transcripts.txt;
cat $GENES2 | cut -f2 > $OUTDIR/metadata/$PROJECT2\_transcripts.txt;
cat $OUTDIR/metadata/$PROJECT1\_transcripts.txt | \
grep -Ff $OUTDIR/metadata/$PROJECT2\_transcripts.txt \
> $OUTDIR/metadata/$PROJECT1\_$PROJECT2\_transcripts.txt;
cat $GENES1 | grep -Ff $OUTDIR/metadata/$PROJECT2\_transcripts.txt | sort -k2 | cut -f1-4 \
> $OUTDIR/metadata/$PROJECT1\_shared_genes.tsv;
cat $GENES2 | grep -Ff $OUTDIR/metadata/$PROJECT2\_transcripts.txt | sort -k2 | cut -f1-4,6,9 \
> $OUTDIR/metadata/$PROJECT2\_shared_genes.tsv;

join -j2 -o1.1,1.2,1.3,2.3,1.4,2.4,2.5,2.6 \
$OUTDIR/metadata/$PROJECT1\_shared_genes.tsv \
$OUTDIR/metadata/$PROJECT2\_shared_genes.tsv | \
tr ' ' '\t' | sort -k1 \
> $OUTDIR/metadata/$PROJECT1\_$PROJECT2\_shared_genes.tsv;

# Genes with different transcripts but overlapping exons
cat $GENES1 | grep -vFf $OUTDIR/metadata/$PROJECT1\_$PROJECT2\_transcripts.txt \
> $OUTDIR/metadata/$PROJECT1\_not_shared_genes.tsv;
cat $GENES2 | grep -vFf $OUTDIR/metadata/$PROJECT1\_$PROJECT2\_transcripts.txt \
> $OUTDIR/metadata/$PROJECT2\_not_shared_genes.tsv;

# Input for function
GENES=$OUTDIR/metadata/$PROJECT1\_$PROJECT2\_shared_genes.tsv;
touch $OUTDIR/metadata/$PROJECT1\_$PROJECT2\_genes_not_work.tsv;

# Remove temporary files
rm $OUTDIR/metadata/$PROJECT1\_transcripts.txt \
$OUTDIR/metadata/$PROJECT2\_transcripts.txt \
$OUTDIR/metadata/$PROJECT1\_$PROJECT2\_transcripts.txt \
$OUTDIR/metadata/$PROJECT1\_shared_genes.tsv \
$OUTDIR/metadata/$PROJECT2\_shared_genes.tsv;

export WORKDIR1 WORKDIR2 OUTDIR FEATURES1 FEATURES2 PROJECT1 PROJECT2 DB_NAME NB_SEQ1 NB_SEQ2 NB_SEQ3 NB_SEQ4;

### Align shared genes
align_shared_genes() {

  GENE=$1;
  TRANS=$2;
  REGION0_1=$3;
  REGION0_2=$4;
  STRATA1=$5;
  STRATA2=$6;

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

  echo $GENE;
  mkdir $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE;
  cd $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE;

  # Align shared genes
  java -Xmx12g -jar \
  /proj/snic2020-2-25/nobackup/simon/conda/envs/gene_alignment/share/macse-2.07-0/macse_v2.07.jar \
  -prog alignTwoProfiles \
  -p1 $WORKDIR1/D3_codon_aware_alignment/genes/$REGION1/$GENE/$PROJECT1\_$GENE\_alignment_NT_original_all.fasta \
  -p2 $WORKDIR2/D3_codon_aware_alignment/genes/$REGION2/$GENE/$PROJECT2\_$GENE\_alignment_NT_original_all.fasta \
  -out_NT $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp.fasta \
  -out_AA $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_AA_temp.fasta;

  # If alignment did not work, report and abort
  if ! test -f $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp.fasta; then
    echo -e "${GENE}\t${TRANS}\t${REGION0_1}\t${REGION0_2}" \
    >> $OUTDIR/metadata/genes_not_work.tsv;
    rm -r $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE;
    exit;
  fi;

  # Make sure the DB sequence is first in the alignment
  cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp.fasta | \
  awk -v DB_NAME=$DB_NAME 'BEGIN {notprint="T"} {if($0 ~ DB_NAME) notprint="F"; \
  if(notprint == "F" && $0 ~ /^>/ && $0 !~ DB_NAME) notprint="T"; if(notprint=="F") print}' \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original.fasta;
  cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp.fasta | \
  awk -v DB_NAME=$DB_NAME 'BEGIN {notprint="T"} {if($0 ~ DB_NAME) notprint="T"; \
  if(notprint == "T" && $0 ~ /^>/ && $0 !~ DB_NAME) notprint="F"; if(notprint=="F") print}' \
  >> $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original.fasta;

  # Remove leading gaps "-"" and frameshift insertions "!" relative to the DB species
  LEADING1=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original.fasta | \
  head -n2 | tail -n1 | grep -ob "-" | cut -d":" -f1 | \
  awk 'BEGIN {count=0} {if($0==count) {count+=1; print}}' | wc -l);
  LEADING2=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original.fasta | \
  head -n2 | tail -n1 | grep -ob "\!" | cut -d":" -f1 | \
  awk -v LEADING1=$LEADING1 'BEGIN {count=LEADING1} {if($0==count) {count+=1; print}}' | wc -l);
  cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original.fasta | \
  sed -r "/^>/! s/.{$LEADING1}//" | sed -r "/^>/! s/.{$LEADING2}//"\
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff.fasta;

  # Convert missing data "N" and frameshift insertions "!" to gaps "-"
  cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff.fasta | \
  awk '{if($0 ~ ">") {print} else {gsub(/[!N]/, "-"); print}}' \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp1.fasta;

  # Filter alignment
  if [ $REGION1 != sex_dropout ] && [ $REGION2 != sex_dropout ]; then
      Gblocks $PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp1.fasta \
      -t=c -p=n -b2=$NB_SEQ1 -b1=$NB_SEQ1 -b3=10 -b4=5 -b5=n;
  elif [ $REGION1 != sex_dropout ] && [ $REGION2 == sex_dropout ]; then
      Gblocks $PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp1.fasta \
      -t=c -p=n -b2=$NB_SEQ2 -b1=$NB_SEQ2 -b3=10 -b4=5 -b5=n;
  elif [ $REGION1 == sex_dropout ] && [ $REGION2 != sex_dropout ]; then
      Gblocks $PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp1.fasta \
      -t=c -p=n -b2=$NB_SEQ3 -b1=$NB_SEQ3 -b3=10 -b4=5 -b5=n;
  elif [ $REGION1 == sex_dropout ] && [ $REGION2 == sex_dropout ]; then
      Gblocks $PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp1.fasta \
      -t=c -p=n -b2=$NB_SEQ4 -b1=$NB_SEQ4 -b3=10 -b4=5 -b5=n;
  fi;

  # Modify to single line fasta
  cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp1.fasta-gb | \
  tr -d ' ' | awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp2.fasta;

  if [ $(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp2.fasta | head -n2 | tail -n1 | tr -d '\n' | wc -c) -gt 0 ]; then

    # Convert stop codons to gaps "---"
    java -Xmx12g -jar \
    /proj/snic2020-2-25/nobackup/simon/conda/envs/gene_alignment/share/macse-2.07-0/macse_v2.07.jar \
    -prog exportAlignment \
    -align $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp2.fasta \
    -codonForFinalStop --- -codonForInternalStop --- \
    -out_NT $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp3.fasta \
    -out_AA $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_AA_temp3.fasta;

    # Filter alignment from stop codon gaps
    if [ $REGION1 != sex_dropout ] && [ $REGION2 != sex_dropout ]; then
        Gblocks $PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp3.fasta \
        -t=c -p=n -b2=$NB_SEQ1 -b1=$NB_SEQ1 -b3=10 -b4=5 -b5=n;
    elif [ $REGION1 != sex_dropout ] && [ $REGION2 == sex_dropout ]; then
        Gblocks $PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp3.fasta \
        -t=c -p=n -b2=$NB_SEQ2 -b1=$NB_SEQ2 -b3=10 -b4=5 -b5=n;
    elif [ $REGION1 == sex_dropout ] && [ $REGION2 != sex_dropout ]; then
        Gblocks $PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp3.fasta \
        -t=c -p=n -b2=$NB_SEQ3 -b1=$NB_SEQ3 -b3=10 -b4=5 -b5=n;
    elif [ $REGION1 == sex_dropout ] && [ $REGION2 == sex_dropout ]; then
        Gblocks $PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp3.fasta \
        -t=c -p=n -b2=$NB_SEQ4 -b1=$NB_SEQ4 -b3=10 -b4=5 -b5=n;
    fi;

    # Modify to single line fasta
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp3.fasta-gb | \
    tr -d ' ' | awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT.fasta;

  else

    # Create dummy files
    cp $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp2.fasta \
    $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT.fasta;
    touch $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp3.fasta-gb;
    touch $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_AA_temp3.fasta;

  fi;

  # Count transcript lengths, number of exons, and lift over success and callable start codons
  NB_SITES1=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original.fasta | \
  head -n2 | tail -n1 | tr -d '\n' | wc -c);
  NB_SITES2=$7;
  NB_SITES3=$(cat $FEATURES1 $FEATURES2 | \
  awk -F'\t' -v GENE=$GENE '{if($13 ~ "^"GENE"_" && $9 != "missing") print $6"\t"$7"\t"$13}' | sort -k3 | \
  awk -F'\t' -v GENE=$GENE 'BEGIN {sum=0} {if(NR!=1 && $3==exon) {sum+=$2-$1+1; exon=$3} else {exon=$3}} END {print sum}');
  NB_SITES4=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT.fasta | \
  head -n2 | tail -n1 | tr -d '\n' | wc -c);
  NB_EXONS1=$8;
  NB_EXONS2=$(cat $FEATURES1 $FEATURES2 | \
  awk -F'\t' -v GENE=$GENE '{if($13 ~ "^"GENE"_" && $9 == "callable") print $6"\t"$7"\t"$13}' | sort -k3 | \
  awk -F'\t' -v GENE=$GENE 'BEGIN {sum=0} {if(NR!=1 && $3==exon) {sum+=1; exon=$3} else {exon=$3}} END {print sum}');
  if [ $(cat $WORKDIR1/D1_filter_genes/transcripts_nonmissing.txt | grep $GENE\_$TRANS | wc -l) -gt 0 ]; then
    INFO=complete_liftover;
  elif [ $(cat $WORKDIR1/D1_filter_genes/transcripts_start_codon_trans_missing.txt $WORKDIR2/D1_filter_genes/transcripts_start_codon_trans_missing.txt | grep $GENE\_$TRANS | wc -l) == 2 ]; then
    INFO=start_codon_incomplete_liftover;
  else
    INFO=no_start_incomplete_liftover;
  fi;
  echo -e "${GENE}\t${TRANS}\t${REGION0_1}\t${REGION0_2}\t${STRATA1}\t${STRATA2}\t$NB_SITES1\t$NB_SITES2\t$NB_SITES3\t$NB_SITES4\t$NB_EXONS1\t$NB_EXONS2\t$INFO\t${PROJECT1}_${PROJECT2}" \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_info.tsv;

  # Remove temporary files
  rm $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_temp* \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_AA_temp*;

  # Subset alignments
  mv $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original.fasta \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_all.fasta;
  mv $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT.fasta \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_all.fasta;
  mv $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff.fasta \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta;

  if [ $REGION1 == "autosomal" ] && [ $REGION2 == "autosomal" ]; then

    # One sequence per female
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_A_A_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_females_A_A.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_A_A_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females_A_A.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_A_A_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females_A_A.fasta;

  elif [ $REGION1 == "sex_dropout" ] && [ $REGION2 == "sex_dropout" ]; then

    # Female Z sequences only
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_females_Z_Z.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females_Z_Z.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females_Z_Z.fasta;

  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "sex_phase" ]; then

    # Female sequences only
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_ZW_ZW_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_females_ZW_ZW.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_ZW_ZW_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females_ZW_ZW.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_ZW_ZW_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females_ZW_ZW.fasta;

  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "sex_dropout" ]; then

    # Female sequences only
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_ZW_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_females_ZW_Z.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_ZW_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females_ZW_Z.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_ZW_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females_ZW_Z.fasta;

  elif [ $REGION1 == "sex_dropout" ] && [ $REGION2 == "sex_phase" ]; then

    # Female sequences only
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_ZW_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_females_Z_ZW.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_ZW_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females_Z_ZW.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_ZW_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females_Z_ZW.fasta;

  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "autosomal" ]; then

    # Female sequences only
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_ZW_A_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_original_females_ZW_A.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_ZW_A_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females_ZW_A.fasta;
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_ZW_A_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_females_ZW_A.fasta;

  fi;
  cd -;

};

## Excecute function in parallell
export -f align_shared_genes;
parallel --colsep '\t' 'align_shared_genes {}' :::: $GENES;

# Make file with info of each gene alignment
cat $(find $OUTDIR/shared_genes -name "*_info.tsv") | sort -k1 \
> $OUTDIR/metadata/genes_info.tsv;
cat $OUTDIR/metadata/genes_info.tsv | cut -f1,2 | tr '\t' '_' \
> $OUTDIR/metadata/genes_work.tsv;
cat $GENES | cut -f1,2 | tr '\t' '_' | grep -vFf $OUTDIR/metadata/genes_work.tsv \
> $OUTDIR/metadata/genes_not_work.tsv;

# Genes with different transcripts but overlapping exons
cat $OUTDIR/metadata/$PROJECT1\_not_shared_genes.tsv | \
awk -F'\t' -v PROJECT1=$PROJECT1 '{print $1"\t"$2"\t"$3"\tNA\t"$4"\tNA\tNA\t"$6"\t"$7"\tNA\t"$9"\t"$10"\t"$11"\t"PROJECT1}' \
>> $OUTDIR/metadata/genes_info.tsv;
cat $OUTDIR/metadata/$PROJECT2\_not_shared_genes.tsv | \
awk -F'\t' -v PROJECT2=$PROJECT2 '{print $1"\t"$2"\t"$3"\tNA\t"$4"\tNA\tNA\t"$6"\t"$7"\tNA\t"$9"\t"$10"\t"$11"\t"PROJECT2}' \
>> $OUTDIR/metadata/genes_info.tsv;
rm $(find $OUTDIR/shared_genes -name "*_info.tsv");

