#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -J haploinsufficiency_scores
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
GENES=$WORKDIR/D4_align_shared_genes/metadata/genes_info.tsv;
HUMAN_PROT=$MAINDIR/data/reference/Homo_sapiens/Homo_sapiens.GRCh37.75.pep.all.fa;
HUMAN_GTF=$MAINDIR/data/reference/Homo_sapiens/Homo_sapiens.GRCh37.75.gtf;
DB_NAME=Tgut;
DB_REF=$MAINDIR/data/Skylark_2021/B3_annotation_lift_over/$DB_NAME.fa;
DB_TRANS=$MAINDIR/data/reference/Taeniopygia_guttata/GCF_003957565.2_bTaeGut1.4.pri_rna.fna;
DB_GTF=$MAINDIR/data/reference/Taeniopygia_guttata/GCF_003957565.2_bTaeGut1.4.pri_genomic.no_W.gtf;
METADATA=$WORKDIR/metadata;
RESOURCES=$METADATA/resources;

### Load modules
module load bioinfo-tools blast/2.15.0+ gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/haploinsufficiency_scores;
OUTDIR=$OUTDIR/haploinsufficiency_scores;
mkdir $OUTDIR/database \
$OUTDIR/matches \
$OUTDIR/metadata;

echo -e "Gene\tTranscript\tHuman_gene_name\tHuman_gene_ID\tMethod\t%HI\tpHaplo" \
> $OUTDIR/haploinsufficiency_scores.tsv;

# Create Database
cp $HUMAN_PROT $OUTDIR/database/human_protein.faa;
cp $DB_TRANS $OUTDIR/database/$DB_NAME\_rna.fna;
makeblastdb -in $OUTDIR/database/human_protein.faa -dbtype prot;

export WORKDIR OUTDIR PROJECT HUMAN_GTF DB_NAME METADATA RESOURCES;

haploinsufficiency_scores() {

  GENE=$1;
  TRANS=$2;

  echo $GENE;
  mkdir $OUTDIR/matches/$GENE\_$TRANS;

  # Match gene to human
  GENE_INFO=$(cat $RESOURCES/biomart_human_homology.tsv | awk -F'\t' -v GENE=$GENE '{if(tolower($1)==tolower(GENE) && $5=="ortholog_one2one" && $6==1) print}');
  HUMAN_GENE=$(echo $GENE_INFO | cut -d' ' -f3);
  HUMAN_GENE_ID=$(echo $GENE_INFO | cut -d' ' -f2);

  if [ $(echo $HUMAN_GENE | tr -d '\n' | wc -c) -gt 0 ]; then

    # Find haploinsufficieny scores
    HI=$(cat $RESOURCES/HI_Predictions_Version3.bed | awk -F'[\t|]' -v HUMAN_GENE=$HUMAN_GENE '{if(tolower($4) == tolower(HUMAN_GENE)) print $5}' | head -n1);
    pHaplo=$(cat $RESOURCES/pHaplo_collins.txt | awk -F'[\t]' -v HUMAN_GENE=$HUMAN_GENE '{if(tolower($1) == tolower(HUMAN_GENE)) print $2}' | tr ',' '.' | head -n1);

  fi;
  if  [ $(echo $HI | tr -d '\n' | wc -c) -gt 0 ] || [ $(echo $pHaplo | tr -d '\n' | wc -c) -gt 0 ]; then

    METHOD=1;

    if  [ $(echo $HI | tr -d '\n' | wc -c) -eq 0 ]; then
      HI="NA";
    fi;

    if  [ $(echo $pHaplo | tr -d '\n' | wc -c) -eq 0 ]; then
      pHaplo="NA";
    fi;

  else

    # Extract database sequence
    cat $OUTDIR/database/$DB_NAME\_rna.fna | awk -v TRANS=$TRANS 'BEGIN {printline="F"} \
    {if($0 ~ "^" ">"TRANS) {printline="T"} else if($0 ~ "^" ">") {printline="F"}; if(printline=="T") {print}}' \
    > $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS.fasta;

    # Blastx against Human
    blastx -query $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS.fasta -db $OUTDIR/database/human_protein.faa \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle" \
    -out $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_human_match -evalue 1e-6 -max_target_seqs 10;

    # Evaluate top match. Try at least 10 attempts
    if [ $(cat $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_human_match | wc -l) -gt 10 ]; then
      ATTEMPTS=10;
    else
      ATTEMPTS=$(cat $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_human_match | wc -l);
    fi;

    if [ $(echo $ATTEMPTS) -gt 0 ]; then

      for ATTEMPT in $(seq 1 1 $ATTEMPTS); do

        HUMAN_GENE_ID=$(cat $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_human_match | head -n$ATTEMPT | tail -n1 | cut -f9 | \
        awk -F'[ :]' '{for(i=1; i<=NF; i++) {if($i=="gene") {FIELD=i+1; break}}; {print $FIELD; exit}}' | cut -d'.' -f1);

        if [ $(echo $ATTEMPT) -gt 1 ] && [ $HUMAN_GENE_ID0 == $HUMAN_GENE_ID ]; then
          continue;
        fi;

        HUMAN_GENE=$(cat $HUMAN_GTF | grep $HUMAN_GENE_ID | awk -F'\t' '{if($3 == "CDS") print $9}' | tr -d '"' | \
        awk -F'[ ;]' -v HUMAN_GENE_ID=$HUMAN_GENE_ID '{for(i=1; i<=NF; i++) {if($i=="gene_id") {FIELD=i+1; break}}; if($FIELD==HUMAN_GENE_ID) {print; exit}}' | \
        awk -F'[ ;]' '{for(i=1; i<=NF; i++) {if($i=="gene_name") {GENE=i+1}}; print $GENE}');

        if [ $(echo $HUMAN_GENE | tr -d '\n' | wc -c) -gt 0 ]; then

          # Find haploinsufficieny scores
          HI=$(cat $RESOURCES/HI_Predictions_Version3.bed | awk -F'[\t|]' -v HUMAN_GENE=$HUMAN_GENE '{if(tolower($4) == tolower(HUMAN_GENE)) print $5}' | head -n1);
          pHaplo=$(cat $RESOURCES/pHaplo_collins.txt | awk -F'[\t]' -v HUMAN_GENE=$HUMAN_GENE '{if(tolower($1) == tolower(HUMAN_GENE)) print $2}' | tr ',' '.' | head -n1);

        fi;
        if  [ $(echo $HI | tr -d '\n' | wc -c) -gt 0 ] || [ $(echo $pHaplo | tr -d '\n' | wc -c) -gt 0 ]; then

          METHOD=2.$ATTEMPT;

          if  [ $(echo $HI | tr -d '\n' | wc -c) -eq 0 ]; then
            HI="NA";
          fi;

          if  [ $(echo $pHaplo | tr -d '\n' | wc -c) -eq 0 ]; then
            pHaplo="NA";
          fi;

          break;

        else

          HUMAN_GENE_ID0=$HUMAN_GENE_ID;

        fi;

      done;

    fi;
    if  [ $(echo $HI | tr -d '\n' | wc -c) -eq 0 ] && [ $(echo $pHaplo | tr -d '\n' | wc -c) -eq 0 ]; then

      HUMAN_GENE="NA";
      HUMAN_GENE_ID="NA";
      METHOD="NA";
      HI="NA";
      pHaplo="NA";

    fi;

  fi;

  # Output info
  echo -e "$GENE\t$TRANS\t$HUMAN_GENE\t$HUMAN_GENE_ID\t$METHOD\t$HI\t$pHaplo" \
  > $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_info.tsv;

};

## Excecute function in parallell
export -f haploinsufficiency_scores;
parallel --colsep '\t' 'haploinsufficiency_scores {}' :::: $GENES;

# Concatenate data
cat $(find $OUTDIR/matches/ -name "*_info.tsv") | sort -k1 \
>> $OUTDIR/haploinsufficiency_scores.tsv;

# Remove temporary files
rm -r $OUTDIR/matches;

