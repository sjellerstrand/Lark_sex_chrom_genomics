#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -C mem256GB
#SBATCH -J codon_aware_alignment
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
GENES=$WORKDIR/D2_extract_sequences/metadata/genes.tsv;
FEATURES=$WORKDIR/D1_filter_genes/$PROJECT\_features.tsv;
DB_NAME=Tgut;
METADATA=$WORKDIR/metadata;
OUTGROUPS=$(echo -e \
"Pmaj \
${MAINDIR}/data/reference/Parus_major/GCF_001522545.3_Parus_major1.1_genomic.fasta \
${WORKDIR}/D1_filter_genes/Pmaj_uniqID.gtf
:\
Falb \
${MAINDIR}/data/reference/Ficedula_albicollis/GCA_000247815.2_FicAlb1.5_genomic.fasta \
${WORKDIR}/D1_filter_genes/Falb_uniqID.gtf");

### Load modules
module load bioinfo-tools SeqKit/0.15.0 gnuparallel/20180822;

### Activate conda environment
conda activate gene_alignment;

### Create folders
mkdir $OUTDIR/D3_codon_aware_alignment;
OUTDIR=$OUTDIR/D3_codon_aware_alignment;
mkdir $OUTDIR/genes \
$OUTDIR/genes/autosomal \
$OUTDIR/genes/sex_dropout \
$OUTDIR/genes/sex_phase \
$OUTDIR/metadata;

# Set up metadata files
echo $DB_NAME \
> $OUTDIR/metadata/female_A_sampleID.txt;
echo $OUTGROUPS | tr ":" "\n" | cut -d' ' -f1 \
>> $OUTDIR/metadata/female_A_sampleID.txt;
cp $OUTDIR/metadata/female_A_sampleID.txt \
$OUTDIR/metadata/female_W_sampleID.txt;
cp $OUTDIR/metadata/female_A_sampleID.txt \
$OUTDIR/metadata/female_Z_sampleID.txt;
for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do
  echo $IND\_A$(shuf -i 1-2 -n 1);
done >> $OUTDIR/metadata/female_A_sampleID.txt;
for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do
  echo $IND\_W;
done >> $OUTDIR/metadata/female_W_sampleID.txt;
for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do
  echo $IND\_Z;
done >> $OUTDIR/metadata/female_Z_sampleID.txt;
cp $OUTDIR/metadata/female_Z_sampleID.txt $OUTDIR/metadata/all_Z_sampleID.txt;
for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Male") print $1}'); do
  echo $IND\_Z1;
  echo $IND\_Z2;
done >> $OUTDIR/metadata/all_Z_sampleID.txt;
cat $OUTDIR/metadata/female_W_sampleID.txt $OUTDIR/metadata/female_Z_sampleID.txt | sort -u \
> $OUTDIR/metadata/female_sexchrom_sampleID.txt;
NB_SEQ2=$(cat $OUTDIR/metadata/all_Z_sampleID.txt | wc -l);
NB_SEQ1=$(echo "$(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}' | wc -l) + $NB_SEQ2" | bc);

export WORKDIR OUTDIR FEATURES PROJECT DB_NAME METADATA NB_SEQ1 NB_SEQ2;

codon_aware_alignment() {

  GENE=$1;
  TRANS=$2;
  REGION0=$3;
  STRATA=$4;

  if [ $REGION0 == "sex_phase" ] || [ $REGION0 == "partial_dropout" ]; then
    REGION="sex_phase";
  else
    REGION=$REGION0;
  fi;

  echo $GENE;
  mkdir $OUTDIR/genes/$REGION/$GENE;

  if [ $REGION != "sex_phase" ]; then
    java -Xmx12g -jar \
    /proj/snic2020-2-25/nobackup/simon/conda/envs/gene_alignment/share/macse-2.07-0/macse_v2.07.jar \
    -prog alignSequences \
    -seq $WORKDIR/D2_extract_sequences/genes/$REGION/$GENE/$PROJECT\_$GENE\_all_sequences.fasta \
    -out_NT $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp.fasta \
    -out_AA $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_AA_temp.fasta;
  else
    java -Xmx12g -jar \
    /proj/snic2020-2-25/nobackup/simon/conda/envs/gene_alignment/share/macse-2.07-0/macse_v2.07.jar \
    -prog alignSequences \
    -seq $WORKDIR/D2_extract_sequences/genes/$REGION0/$GENE/$PROJECT\_$GENE\_nonW_sequences.fasta \
    -seq_lr $WORKDIR/D2_extract_sequences/genes/$REGION0/$GENE/$PROJECT\_$GENE\_W_sequences.fasta \
    -out_NT $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp.fasta \
    -out_AA $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_AA_temp.fasta;
  fi;

  # Make sure the DB sequence is first in the alignment
  cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp.fasta | \
  awk -v DB_NAME=$DB_NAME 'BEGIN {notprint="T"} {if($0 ~ DB_NAME) notprint="F"; \
  if(notprint == "F" && $0 ~ /^>/ && $0 !~ DB_NAME) notprint="T"; if(notprint=="F") print}' \
  > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original.fasta;
  cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp.fasta | \
  awk -v DB_NAME=$DB_NAME 'BEGIN {notprint="T"} {if($0 ~ DB_NAME) notprint="T"; \
  if(notprint == "T" && $0 ~ /^>/ && $0 !~ DB_NAME) notprint="F"; if(notprint=="F") print}' \
  >> $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original.fasta;

  # Remove leading gaps "-"" and frameshift insertions "!" relative to the DB species
  LEADING1=$(cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original.fasta | \
  head -n2 | tail -n1 | grep -ob "-" | cut -d":" -f1 | \
  awk 'BEGIN {count=0} {if($0==count) {count+=1; print}}' | wc -l);
  LEADING2=$(cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original.fasta | \
  head -n2 | tail -n1 | grep -ob "\!" | cut -d":" -f1 | \
  awk -v LEADING1=$LEADING1 'BEGIN {count=LEADING1} {if($0==count) {count+=1; print}}' | wc -l);
  cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original.fasta | \
  sed -r "/^>/! s/.{$LEADING1}//" | sed -r "/^>/! s/.{$LEADING2}//"\
  > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff.fasta;

  # Convert missing data "N" and frameshift insertions "!" to gaps "-"
  cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff.fasta | \
  awk '{if($0 ~ ">") {print} else {gsub(/[!N]/, "-"); print}}' \
  > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp.fasta;

  # Filter alignment
  if [ $REGION != sex_dropout ]; then
      Gblocks $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp.fasta \
      -t=c -p=n -b2=$NB_SEQ1 -b1=$NB_SEQ1 -b3=10 -b4=5 -b5=n;
  else
    Gblocks $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp.fasta \
    -t=c -p=n -b2=$NB_SEQ2 -b1=$NB_SEQ2 -b3=10 -b4=5 -b5=n;
  fi;

  # Modify to single line fasta
  cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp.fasta-gb | \
  tr -d ' ' | awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
  > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT.fasta;

  # Count transcript lengths, number of exons, and lift over success and callable start codons
  NB_SITES1=$(cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp.fasta | \
  head -n2 | tail -n1 | tr -d '\n' | wc -c);
  NB_SITES2=$(cat $FEATURES | \
  awk -F'\t' -v GENE=$GENE 'BEGIN {sum=0} {if($13 ~ "^"GENE"_") {sum+=$3-$2+1}} END {print sum}');
  NB_SITES3=$(cat $FEATURES | \
  awk -F'\t' -v GENE=$GENE 'BEGIN {sum=0} {if($13 ~ "^"GENE"_" && $5 != "missing") {sum+=$7-$6+1}} END {print sum}');
  NB_SITES4=$(cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT.fasta | \
  head -n2 | tail -n1 | tr -d '\n' | wc -c);
  NB_EXONS1=$(cat $FEATURES | \
  awk -F'\t' -v GENE=$GENE 'BEGIN {sum=0} {if($13 ~ "^"GENE"_") {sum+=1}} END {print sum}');
  NB_EXONS2=$(cat $FEATURES | \
  awk -F'\t' -v GENE=$GENE 'BEGIN {sum=0} {if($13 ~ "^"GENE"_" && $5 != "missing") {sum+=1}} END {print sum}');
  if [ $(cat $WORKDIR/D1_filter_genes/transcripts_nonmissing.txt | grep $GENE\_$TRANS | wc -l) -gt 0 ]; then
    INFO=complete_liftover;
  elif [ $(cat $WORKDIR/D1_filter_genes/transcripts_start_codon_trans_missing.txt | grep $GENE\_$TRANS | wc -l) -gt 0 ]; then
    INFO=start_codon_incomplete_liftover;
  else
    INFO=no_start_incomplete_liftover;
  fi;
  echo -e "${GENE}\t${TRANS}\t${REGION0}\t${STRATA}\t$NB_SITES1\t$NB_SITES2\t$NB_SITES3\t$NB_SITES4\t$NB_EXONS1\t$NB_EXONS2\t$INFO" \
  > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_info.tsv;

  # Remove temporary files
  rm $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_temp* \
  $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_AA_temp.fasta;

  # Subset alignments
  mv $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original.fasta \
  $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_all.fasta;
  mv $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT.fasta \
  $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_all.fasta;
  mv $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff.fasta \
  $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_all.fasta;

  if [ $REGION == "autosomal" ]; then

    # One sequence per female
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_A_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_females_A.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_A_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_females_A.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_A_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_females_A.fasta;

  elif [ $REGION == "sex_dropout" ]; then

    # Female Z sequences only
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_females_Z.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_females_Z.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_females_Z.fasta;

  elif [ $REGION == "sex_phase" ]; then

    # Female sequences only
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_sexchrom_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_females_ZW.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_sexchrom_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_females_ZW.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_sexchrom_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_females_ZW.fasta;

    # Female Z sequences only
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_females_Z.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_females_Z.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_Z_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_females_Z.fasta;

    # Female W sequences only
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_W_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_original_females_W.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_W_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_females_W.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_all.fasta | \
    seqkit grep -f $OUTDIR/metadata/female_W_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_alignment_NT_snpeff_females_W.fasta;

  fi;

};

## Excecute function in parallell
export -f codon_aware_alignment;
parallel --colsep '\t' 'codon_aware_alignment {}' :::: $GENES;

# Make file with info of each gene alignment
cat $(find $OUTDIR/genes -name "*_info.tsv") | sort -k1 | awk -F'\t' 'BEGIN {OFS="\t"} {if($5!=0) print}' \
> $OUTDIR/metadata/genes_info.tsv;
cat $OUTDIR/metadata/genes_info.tsv | cut -f1,2 | tr '\t' '_' \
> $OUTDIR/metadata/genes_work.tsv;
cat $GENES | cut -f1,2 | tr '\t' '_' | grep -vFf  $OUTDIR/metadata/genes_work.tsv \
> $OUTDIR/metadata/genes_not_work.tsv;
rm $(find $OUTDIR/genes -name "*_info.tsv");

