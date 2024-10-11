#!/bin/bash -l

#SBATCH -A snic2022-5-484
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J define_regions
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
METADATA=$WORKDIR/metadata;
BEDS=$WORKDIR/A8_PhaseWY/final_output/beds;

### Load modules
module load bioinfo-tools BEDTools/2.29.2;

### Create folders
mkdir $OUTDIR/A9_define_regions;
OUTDIR=$OUTDIR/A9_define_regions;

### Target regions
bedtools subtract -a $BEDS/$PROJECT\_target_region.bed \
-b $METADATA/region_unknown.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_target_region.bed;

### Autosomal and PAR
cat $METADATA/region_auto.bed $METADATA/region_PAR_3.bed $METADATA/region_PAR_3_unk.bed \
$METADATA/region_PAR_5.bed $METADATA/region_PAR_5_unk.bed \
> $OUTDIR/region_autosomal_and_PAR.bed;
bedtools intersect -a $OUTDIR/region_autosomal_and_PAR.bed \
-b $BEDS/$PROJECT\_autosomal.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_autosomal_and_PAR.bed;
rm $OUTDIR/region_autosomal_and_PAR.bed;

### Autosomal
bedtools intersect -a $METADATA/region_auto.bed \
-b $OUTDIR/$PROJECT\_autosomal_and_PAR.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_autosomal.bed;

### Heterogametic drop out
bedtools subtract -a $BEDS/$PROJECT\_hetgam_dropout.bed \
-b $METADATA/region_unknown.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_hetgam_dropout.bed;

### Homogametic
bedtools subtract -a $BEDS/$PROJECT\_homogametic.bed \
-b $METADATA/region_unknown.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_homogametic.bed;

### Heterogametic
bedtools subtract -a $BEDS/$PROJECT\_heterogametic.bed \
-b $METADATA/region_unknown.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_heterogametic.bed;

### PAR 3
bedtools intersect -a $METADATA/region_PAR_3.bed \
-b $OUTDIR/$PROJECT\_autosomal_and_PAR.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_PAR_3.bed;

### PAR 3 unknown
bedtools intersect -a $METADATA/region_PAR_3_unk.bed \
-b $OUTDIR/$PROJECT\_autosomal_and_PAR.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_PAR_3_unk.bed;

### PAR 5
bedtools intersect -a $METADATA/region_PAR_5.bed \
-b $OUTDIR/$PROJECT\_autosomal_and_PAR.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_PAR_5.bed;

### PAR 5 unknown
bedtools intersect -a $METADATA/region_PAR_5_unk.bed \
-b $OUTDIR/$PROJECT\_autosomal_and_PAR.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_PAR_5_unk.bed;

### Ancestral Z
bedtools intersect -a $METADATA/region_Z.bed \
-b $OUTDIR/$PROJECT\_homogametic.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_Z.bed;

### 4A
bedtools intersect -a $METADATA/region_4A.bed \
-b $OUTDIR/$PROJECT\_homogametic.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_4A.bed;

### 3a
bedtools intersect -a $METADATA/region_3a.bed \
-b $OUTDIR/$PROJECT\_homogametic.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_3a.bed;

### 3b
bedtools intersect -a $METADATA/region_3b.bed \
-b $OUTDIR/$PROJECT\_homogametic.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_3b.bed;

### 3c
bedtools intersect -a $METADATA/region_3c.bed \
-b $OUTDIR/$PROJECT\_homogametic.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_3c.bed;


### 5
bedtools intersect -a $METADATA/region_5.bed \
-b $OUTDIR/$PROJECT\_homogametic.bed | \
bedtools sort | bedtools merge \
> $OUTDIR/$PROJECT\_5.bed;

### Unknown
cp $METADATA/region_unknown.bed $OUTDIR/$PROJECT\_unknown.bed;
