#!/bin/bash
#$ -pe smp 6
#$ -l h_vmem=16G
#$ -l h_rt=72:0:0
#$ -cwd
#$ -j y
#$ -N  filter_stats
#$ -o /data/scratch/bt24076
#$ -m bea
#$ -l rocky

set -e 

#load modules
module load bcftools
module load vcftools

#loop for all chromsomes
for i in {1..22} X;
do

path to vcfs
path_to_vcf="/data/scratch/bt24076/03-02-2025-T2D_pjct_Data/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

bangladesh
bcftools view --threads $NSLOTS -S /data/scratch/bt24076/03-02-2025-T2D_pjct_Data/benbang.txt -o /data/scratch/bt24076/03-02-2025-T2D_pjct_Data/beb_filt/benbang_CHR${i}_filt.vcf $path_to_vcf

indian in hou
bcftools view --threads $NSLOTS -S /data/scratch/bt24076/03-02-2025-T2D_pjct_Data/gujIndHou.txt -o /data/scratch/bt24076/03-02-2025-T2D_pjct_Data/gih_filt/gujIndHou_CHR${i}_filt.vcf $path_to_vcf

pakistan
bcftools view --threads $NSLOTS -S /data/scratch/bt24076/03-02-2025-T2D_pjct_Data/pjlPak.txt -o /data/scratch/bt24076/03-02-2025-T2D_pjct_Data/pjl_filt/pjlPak_CHR${i}_filt.vcf $path_to_vcf

indian Tel UK
bcftools view --threads $NSLOTS -S /data/scratch/bt24076/03-02-2025-T2D_pjct_Data/indTelUK.txt -o /data/scratch/bt24076/03-02-2025-T2D_pjct_Data/itu_filt/indTelUK_CHR${i}_filt.vcf $path_to_vcf
done


#Tajim's D on all Chrom and Populations
for i in {1..22} X;
do

#paths to filtered south asain popultaions
beb_vcf="/data/scratch/bt24076/03-02-2025-T2D_pjct_Data/beb_filt/benbang_CHR${i}_filt.vcf"
gih_vcf="/data/scratch/bt24076/03-02-2025-T2D_pjct_Data/gih_filt/gujIndHou_CHR${i}_filt.vcf"
itu_vcf="/data/scratch/bt24076/03-02-2025-T2D_pjct_Data/itu_filt/indTelUK_CHR${i}_filt.vcf"
pjl_vcf="/data/scratch/bt24076/03-02-2025-T2D_pjct_Data/pjl_filt/pjlPak_CHR${i}_filt.vcf"

#Tajimas D calculation
vcftools --vcf $beb_vcf --TajimaD 10000 --out /data/scratch/bt24076/11-02-2025_psResults/beb/beb_CHR${i}

vcftools --vcf $gih_vcf --TajimaD 10000 --out /data/scratch/bt24076/11-02-2025_psResults/gih/gih_CHR${i}

vcftools --vcf $itu_vcf --TajimaD 10000 --out /data/scratch/bt24076/11-02-2025_psResults/itu/itu_CHR${i}

vcftools --vcf $pjl_vcf --TajimaD 10000 --out /data/scratch/bt24076/11-02-2025_psResults/pjl/pjl_CHR${i}

done
