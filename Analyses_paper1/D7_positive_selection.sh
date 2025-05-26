#!/bin/bash -l

#SBATCH -A naiss2024-5-340
#SBATCH -p shared
#SBATCH -n 2
#SBATCH -t 20:00:00
#SBATCH -J positive_selection
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/cfs/klemming/projects/supr/snic2020-2-25/user_data/simon/Sylvioidea;
PROJECT1=Skylark_2021;
PROJECT2=Rasolark_2021;
WORKDIR1=$MAINDIR/data/$PROJECT1;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/$PROJECT1;
GENES=$WORKDIR1/organise_data1/$PROJECT1\_$PROJECT2\_organised_data1.tsv;
METADATA1=$WORKDIR1/metadata;
METADATA2=$WORKDIR2/metadata;
TREE=$METADATA1/tree_topologies/Skylark_2021_Rasolark_2021_sex_phase_sex_phase_females.nwk;

DB_NAME=Tgut;
DB_NAME=Tgut;
OUTGROUPS=$(echo -e \
"Pmaj \
${MAINDIR}/data/reference/Parus_major/GCF_001522545.3_Parus_major1.1_genomic.fasta \
${WORKDIR}/D1_filter_genes/Pmaj_uniqID.gtf
:\
Falb \
${MAINDIR}/data/reference/Ficedula_albicollis/GCA_000247815.2_FicAlb1.5_genomic.fasta \
${WORKDIR}/D1_filter_genes/Falb_uniqID.gtf");

### Activate conda environment
conda activate hyphy;

### Create folders
mkdir $OUTDIR/D7_positive_selection3;
OUTDIR=$OUTDIR/D7_positive_selection3;
mkdir $OUTDIR/metadata \
$OUTDIR/genes;

# Subset genes
cat $GENES | awk -F'\t' '{if($98 == "OK" && $100=="OK" && $101=="OK" && $102=="OK" && \
$4!="S0" && $4!="S1" && $4!="S2" && $4!="S3" && $4!="Ancestral unknown" && $3=="sex_phase") print $1}' | \
sort | uniq -D | uniq \
> $OUTDIR/metadata/genes_in.tsv;
GENES2=$OUTDIR/metadata/genes_in.tsv;

# Annotate branches in tree
cat $TREE | sed -E 's/\b(BranchZ|BranchW|BranchA)\b/&{Foreground}/g' \
> $OUTDIR/metadata/tree.nwk;
TREE=$OUTDIR/metadata/tree.nwk;

for GENE in $(cat $GENES2); do

  echo $GENE;
  mkdir $OUTDIR/genes/$GENE;
  cd $OUTDIR/genes/$GENE;
  
  hyphy absrel \
  --alignment $WORKDIR1/D4_align_shared_genes/shared_genes/sex_phase_sex_phase/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females_ZW_ZW.fasta \
  --code Universal --tree $TREE --branches Foreground \
  --output $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_ABSREL.json \
  --save-fit $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_ABSREL.txt;

done;


# Retrieve output data
echo -e "geneID\tp-valueZ\tpvalue_W\tpvalue_A\tcorrected_p-valueZ1\tcorrected_p-valueW1\tcorrected_p-valueA1\tcorrected_p-valueZ2\tcorrected_p-valueW2\tcorrected_p-valueA2" \
> $OUTDIR/metadata/ABSREL_results.tsv;
NGENES=$(cat $GENES2 | wc -l);

for GENE in $(cat $GENES2); do

  PVAL_Z=$(echo $(jq '.["branch attributes"]["0"]["BranchZ"]["Uncorrected P-value"]' \
  $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_ABSREL.json));
  
  if [[ -n "$PVAL_Z" && "$PVAL_Z" != "null" ]]; then
    
    PVALCORR_Z1=$(echo $(jq '.["branch attributes"]["0"]["BranchZ"]["Corrected P-value"]' \
    $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_ABSREL.json));
    PVALCORR_Z2=$(echo $PVALCORR_Z1 | awk -v NGENES=$NGENES '{pcorr=$1*NGENES; if(pcorr > 1) {pcorr=1}; print pcorr}');
    
  else 
    PVAL_Z="NA";
    PVALCORR_Z1="NA";
    PVALCORR_Z2="NA";
  
  fi;
   
  PVAL_W=$(echo $(jq '.["branch attributes"]["0"]["BranchW"]["Uncorrected P-value"]' \
  $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_ABSREL.json));
  
 if [[ -n "$PVAL_W" && "$PVAL_W" != "null" ]]; then
   
    PVALCORR_W1=$(echo $(jq '.["branch attributes"]["0"]["BranchW"]["Corrected P-value"]' \
    $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_ABSREL.json));
    PVALCORR_W2=$(echo $PVALCORR_W1 | awk -v NGENES=$NGENES '{pcorr=$1*NGENES; if(pcorr > 1) {pcorr=1}; print pcorr}');
    
  else
    PVAL_W="NA";
    PVALCORR_W1="NA";
    PVALCORR_W2="NA";
  
  fi;

  
  PVAL_A=$(echo $(jq '.["branch attributes"]["0"]["BranchA"]["Uncorrected P-value"]' \
  $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_ABSREL.json));
  
  if [[ -n "$PVAL_A" && "$PVAL_A" != "null" ]]; then
   
    PVALCORR_A1=$(echo $(jq '.["branch attributes"]["0"]["BranchA"]["Corrected P-value"]' \
    $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_ABSREL.json));
    PVALCORR_A2=$(echo $PVALCORR_A1 | awk -v NGENES=$NGENES '{pcorr=$1*NGENES; if(pcorr > 1) {pcorr=1}; print pcorr}');
    
  else
    PVAL_A="NA";
    PVALCORR_A1="NA";
    PVALCORR_A2="NA";
  
  fi;
  echo -e "${GENE}\t${PVAL_Z}\t${PVAL_W}\t${PVAL_A}\t${PVALCORR_Z1}\t${PVALCORR_W1}\t${PVALCORR_A1}\t${PVALCORR_Z2}\t${PVALCORR_W2}\t${PVALCORR_A2}";

done >> $OUTDIR/metadata/ABSREL_results.tsv;