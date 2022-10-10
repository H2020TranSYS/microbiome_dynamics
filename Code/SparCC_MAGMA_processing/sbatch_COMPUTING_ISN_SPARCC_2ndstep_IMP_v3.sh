#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=2 # each task uses 1 cpu
#SBATCH --partition=urtgen_unlimit
#SBATCH --mem-per-cpu=10000 #10GB


source /home/u/f058977/miniconda3/etc/profile.d/conda.sh
conda activate cpp_sparcc_env

# directory
specific_file_name=$1
directori=/massstorage/URT/GEN/BIO3/Federico/LucKi_project/SparCC_MAGMA_processing
data=DATA/OTU_TABLE_input_sparCC_${specific_file_name}.tsv

mkdir -p ${directori}/LooData${specific_file_name}
mkdir -p ${directori}/LooNet${specific_file_name}
mkdir -p ${directori}/GlobNet${specific_file_name}

a="$( head -2 ${directori}/${data} |tail -1 |tr '\t' '\n' |wc -l)"
for((i=1;i<$((a));i++)) 
do
echo $i
Ind=$((i))
if [ "$(($i+2))" -le "$a" ]; then ## to have also the case in which we have the last observation coverd

cut -f1-${i},$((i+2))-$a ${directori}/${data} > ${directori}/LooData${specific_file_name}/DataInd_${Ind}.tsv
## to have the mapping idnividual new name

cut -f$((i+1)) ${directori}/${data} | head -n1 >> ${directori}/GlobNet${specific_file_name}/name.txt
name=$(cut -f$((i+1))  ${directori}/${data} | head -n1 | tr -d '\r')

echo $((i)) >> ${directori}/GlobNet${specific_file_name}/individuals.txt 
else
cut -f1-${i} ${directori}/${data} > ${directori}/LooData${specific_file_name}/DataInd_${Ind}.tsv

#mapping individual new name 
cut -f$((i+1))  ${directori}/${data} | head -n1 >> ${directori}/GlobNet${specific_file_name}/name.txt
name=$(cut -f$((i+1))  ${directori}/${data} | head -n1 | tr -d '\r')
echo $((i)) >> ${directori}/GlobNet${specific_file_name}/individuals.txt

fi
fastspar --otu_table ${directori}/LooData${specific_file_name}/DataInd_${Ind}.tsv --correlation  ${directori}/LooNet${specific_file_name}/median_correlation_${name}.tsv --covariance ${directori}/LooNet${specific_file_name}/median_covariance_${name}.tsv
#tr -d '\r' < ${directori}/LooNet${specific_file_name}/median_correlation_${name}.tsv > ${directori}/LooNet${specific_file_name}/median_correlation_${name}.tsv
# tr -d '\r' < ${directori}/LooNet${specific_file_name}/median_covariance_${name}.tsv > ${directori}/LooNet${specific_file_name}/median_covariance_${name}.tsv


done

## COMPUTE GLOBAL NET
fastspar --otu_table ${directori}/${data} --correlation ${directori}/GlobNet${specific_file_name}/median_correlation.tsv --covariance ${directori}/GlobNet${specific_file_name}/median_covariance.tsv
