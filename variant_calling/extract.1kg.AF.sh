#!/bin/sh
#SBATCH --time=20:00:00 --mem=5GB

set -o nounset

# set number of processes

KG_DIR="$1"
outdir="$2"
wgs_samples_dir="$3"

# get AF of the 5 superpopulations from the 1kg vcfs

# for each chromosome, extract the AF information from the vcf
# only do this for the sites in GTEx


gtexsnp=${wgs_samples_dir}GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller_446EAonly_SNPs.frq



export gtexsnp
export outdir


for i in chr{1..22}; do 
    echo $i
    vcftools --gzvcf ${KG_DIR}ALL.$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --out ${outdir}$i.SNPs --remove-indels --positions $gtexsnp --get-INFO EAS_AF --get-INFO AMR_AF --get-INFO AFR_AF --get-INFO EUR_AF --get-INFO SAS_AF
done


python process.1kg.AF.py $outdir | sort --temporary-directory=$outdir -k1,1 -k2,2n > ${outdir}SNPs.1kg.AF.bed

if false; then

cat ${outdir}*SNPs.INFO | ./process.1kg.AF.py | sort --temporary-directory=$outdir -k1,1 -k2,2n > ${outdir}SNPs.1kg.AF.bed


rm ${outdir}chr*.INFO
gzip ${outdir}*.bed
fi