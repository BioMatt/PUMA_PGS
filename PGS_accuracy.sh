#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=40G
#SBATCH --array=1-100
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0
module load gemma/0.98.3
module load plink/1.9b_6.21-x86_64

cd /home/biomatt/scratch/PGS_accuracy

gemma -bfile /home/biomatt/scratch/PGS_accuracy/gwas_data/spring74_revision${SLURM_ARRAY_TASK_ID} -gk 1 -miss 1 -maf 0 -notsnp -o relatedness_matrix_${SLURM_ARRAY_TASK_ID}

gemma -bfile /home/biomatt/scratch/PGS_accuracy/gwas_data/spring74_revision${SLURM_ARRAY_TASK_ID} -k /home/biomatt/scratch/PGS_accuracy/output/relatedness_matrix_${SLURM_ARRAY_TASK_ID}.cXX.txt -miss 1 -maf 0 -lmm 1 -c /home/biomatt/scratch/PGS_accuracy/gwas_data/spring74_${SLURM_ARRAY_TASK_ID}_cov7.colony.sex.year.age.pc1.spdur.spdist.txt -o spring74_${SLURM_ARRAY_TASK_ID}

plink --bfile /home/biomatt/scratch/PGS_accuracy/gwas_data/spring74_revision${SLURM_ARRAY_TASK_ID} --allow-extra-chr --clump-p1 1 --clump-r2 0.2 --clump-kb 250 --clump /home/biomatt/scratch/PGS_accuracy/output/spring74_${SLURM_ARRAY_TASK_ID}.assoc.txt --clump-snp-field rs --clump-field p_wald --allow-no-sex --out /home/biomatt/scratch/PGS_accuracy/output/spring74_clumped_${SLURM_ARRAY_TASK_ID}

awk 'NR!=1{print $3}' /home/biomatt/scratch/PGS_accuracy/output/spring74_clumped_${SLURM_ARRAY_TASK_ID}.clumped > /home/biomatt/scratch/PGS_accuracy/output/spring74_clumped_${SLURM_ARRAY_TASK_ID}.valid.snp

plink --bfile /home/biomatt/scratch/PGS_accuracy/gwas_data/spring74_revision${SLURM_ARRAY_TASK_ID} --allow-extra-chr --allow-no-sex --score /home/biomatt/scratch/PGS_accuracy/output/spring74_${SLURM_ARRAY_TASK_ID}.assoc.txt 2 5 8 --q-score-range /home/biomatt/scratch/PGS_accuracy/range_list /home/biomatt/scratch/PGS_accuracy/output/spring74_${SLURM_ARRAY_TASK_ID}.assoc.txt 2 12 header --extract /home/biomatt/scratch/PGS_accuracy/output/spring74_clumped_${SLURM_ARRAY_TASK_ID}.valid.snp --out /home/biomatt/scratch/PGS_accuracy/output/spring74_PRS_${SLURM_ARRAY_TASK_ID}

plink --bfile /home/biomatt/scratch/PGS_accuracy/gwas_data/spring74_test_${SLURM_ARRAY_TASK_ID} --allow-extra-chr --allow-no-sex --score /home/biomatt/scratch/PGS_accuracy/output/spring74_${SLURM_ARRAY_TASK_ID}.assoc.txt 2 5 8 --q-score-range /home/biomatt/scratch/PGS_accuracy/range_list /home/biomatt/scratch/PGS_accuracy/output/spring74_${SLURM_ARRAY_TASK_ID}.assoc.txt 2 12 header --extract /home/biomatt/scratch/PGS_accuracy/output/spring74_clumped_${SLURM_ARRAY_TASK_ID}.valid.snp --out /home/biomatt/scratch/PGS_accuracy/test_out/spring74_PRS_test_${SLURM_ARRAY_TASK_ID}
