#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J Skylark_HWE
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## General
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
PROJECT2=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
VCFS=$WORKDIR/A8_PhaseWY/final_output/vcfs;
BEDS=$WORKDIR/A8_PhaseWY/final_output/beds;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
METADATA=$WORKDIR/metadata;
FUNCTIONS=$MAINDIR/scripts/$PROJECT/analyses;
PLOT=$MAINDIR/scripts/$PROJECT/plotting;

### Load modules
module load bioinfo-tools bcftools/1.17 vcftools/0.1.16 vcflib/1.0.1 ADMIXTURE/1.3.0 plink/1.90b4.9 R_packages/4.0.0;

### Create folders
mkdir $OUTDIR/Skylark_HWE;
OUTDIR=$OUTDIR/Skylark_HWE;
mkdir $OUTDIR/metadata \
$OUTDIR/data \
$OUTDIR/HWE \
$OUTDIR/Fst \
$OUTDIR/PCA \
$OUTDIR/admixture;

# Set up data
cat $METADATA/sample_info.txt | tail -n+2 | \
awk -F'\t' 'BEGIN {print FID"\t"IID"\t"Population} {print $1"\t"$1"\t"$4}' \
> $OUTDIR/metadata/pop_file.txt;
cat $METADATA/sample_info.txt | tail -n+2 | awk -F'\t' '{print $1"\t"$1}' \
> $OUTDIR/metadata/samples_autosomal.txt;
cp $OUTDIR/metadata/samples_autosomal.txt \
$OUTDIR/metadata/samples_homogametic.txt;
cat $METADATA/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3 == "Female") print $1"\t"$1}' \
> $OUTDIR/metadata/samples_heterogametic.txt;


DATASETS=$(echo autosomal);
for DATA in $DATASETS; do

  ### Linkage prune data
  if [ $DATA == "autosomal" ]; then
    GENOMES=$(cat $METADATA/sample_info.txt | tail -n+2 | \
    awk -F'\t' 'BEGIN {SUM=0} {SUM+=2} END {print SUM}');
  fi;

  bcftools view $VCFS/$PROJECT\_$DATA\_snps.vcf.gz | \
  vcffilter -f  "AC > 1 & AC < ( $GENOMES - 1 )" | \
  bgzip -c > $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.vcf.gz;
  tabix $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.vcf.gz;

  plink --vcf $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.vcf.gz \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --vcf-half-call missing \
  --indep-pairwise 50 10 0.1 --out $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2;


  ### Caluclate HWE
  mkdir $OUTDIR/HWE/$DATA;

  ## Create bed files
  plink --vcf $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.vcf.gz \
  --keep $OUTDIR/metadata/samples_$DATA.txt \
  --extract $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.prune.in \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --vcf-half-call missing \
  --make-bed --out $OUTDIR/HWE/$DATA/$PROJECT\_$DATA;

  ## Calculate HWE
  plink --bfile $OUTDIR/HWE/$DATA/$PROJECT\_$DATA \
  --keep $OUTDIR/metadata/samples_$DATA.txt \
  --extract $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.prune.in \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --hardy --out $OUTDIR/HWE/$DATA/$PROJECT\_$DATA;

  # Perform Fisher's Combined Probability Test
  Rscript -e  \
  'args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
  for(i in 1:length(args)) {
    assign(args[[i]][1], args[[i]][2])
  }
  pval <-  read.table(paste(OUTDIR, "/HWE/", DATA, "/", PROJECT, "_", DATA, ".hwe", sep=""), header=TRUE)$P
  chisq_stat <- -2 * sum(log(pval))
  df <- 2 * length(pval)
  combined_pval <- pchisq(chisq_stat, df, lower.tail = FALSE)
  write(combined_pval, paste(OUTDIR, "/HWE/", DATA, "/", PROJECT, "_", DATA, "_global_pval.txt", sep=""))' \
  --args OUTDIR=$OUTDIR DATA=$DATA PROJECT=$PROJECT;

  ### Calculate Fst

  ## Create bed files
  mkdir $OUTDIR/Fst/$DATA;
  plink --vcf $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.vcf.gz \
  --keep $OUTDIR/metadata/samples_$DATA.txt \
  --extract $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.prune.in \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --vcf-half-call missing \
  --make-bed --out $OUTDIR/Fst/$DATA/$PROJECT\_$DATA;

  ## Calculte Fst
  plink --bfile $OUTDIR/Fst/$DATA/$PROJECT\_$DATA \
  --keep $OUTDIR/metadata/samples_$DATA.txt \
  --extract $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.prune.in \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --within $OUTDIR/metadata/pop_file.txt \
  --fst --out $OUTDIR/Fst/$DATA/$PROJECT\_$DATA;

  ### Calculate principal components
  mkdir $OUTDIR/PCA/$DATA;
  plink --vcf $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.vcf.gz \
  --keep $OUTDIR/metadata/samples_$DATA.txt \
  --extract $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.prune.in \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --vcf-half-call missing \
  --pca 20 --out $OUTDIR/PCA/$DATA/$PROJECT\_$DATA;


  ### Perform admixture analyses
  mkdir $OUTDIR/admixture/$DATA;
  VCF_IN=$VCFS/$PROJECT\_$DATA\_neutral_pruned.vcf.gz

  ## Create bed files
  plink --vcf $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.vcf.gz \
  --keep $OUTDIR/metadata/samples_$DATA.txt \
  --extract $OUTDIR/data/$PROJECT\_$DATA\_snps_mac2.prune.in \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --vcf-half-call missing \
  --make-bed --out $OUTDIR/admixture/$DATA/$PROJECT\_$DATA;

  awk '{$1=0;print $0}' $OUTDIR/admixture/$DATA/$PROJECT\_$DATA.bim \
  > $OUTDIR/admixture/$DATA/$PROJECT\_$DATA.bim.tmp;
  mv $OUTDIR/admixture/$DATA/$PROJECT\_$DATA.bim.tmp \
  $OUTDIR/admixture/$DATA/$PROJECT\_$DATA.bim;

  ## Run admixture
  cd $OUTDIR/admixture/$DATA;
  K=$(seq 1 1 4  | sort -n -r);
  for k in $K; do
    admixture $OUTDIR/admixture/$DATA/$PROJECT\_$DATA.bed $k --cv=20 -j20 -s=12345 \
    > $OUTDIR/admixture/$DATA/$PROJECT\_$DATA\_K$k;
  done;
  cat $OUTDIR/admixture/$DATA/$PROJECT\_$DATA\_K* | grep "CV error" | cut -d' ' -f 3,4 | \
  sed -e 's/(K=\(.*\)):/\1/' | sort -n > $OUTDIR/admixture/$DATA/$PROJECT\_$DATA\_CVerror;
  cd -;

done;

# Remove temporary files
rm -r $OUTDIR/data $(find $OUTDIR -name "*.bim");
