# Purple martin polygenic scores
Scripts used to calculate polygenic scores and estimate polygenic score accuracy in spring migration timing for the purple martin. Preprint here: https://www.biorxiv.org/content/10.1101/2022.09.14.508039v1

Published paper here:
https://doi.org/10.1038/s41598-023-29470-7

Also, the rest of the code used for the paper is in these repos:
https://github.com/edegreef/PUMA-resequencing-data

https://github.com/edegreef/PUMA-reference-genome

The script revision_PRS.R was used on a series of .profile files from plink to calculate polygenic scores at different p-value thresholds, and assess R^2 among them. 

The script PGS_revision_random_sampling.R took genotype and covariate files in bed format and randomly created 100 training and test datasets with a 85/15 split. It calculated new principal component scores for the covariate files based on the training data to simulate a smaller GWAS

The script PGS_accuracy.sh was run on the Graham cluster of the Digital Research Alliance of Canada, and used gemma for GWAS and plink for calculating polygenic scores. Gemma first needs a relatedness matrix, then to run the GWAS. Plink then takes GWAS results, clumps SNPs based on linkage, and calculates polygenic scores at a range of adjusted p-values (0.001, 0.05, 0.1, 0.2, 0.3, 0.4, and 0.5). Plink was then used to calculate polygenic scores for the test individuals based on results from the training individuals.

The script PGS_revision_accuracy_test.R takes the .profile polygenic score output files from the training and test datasets, the covariate files, and the .fam phenotype files to assess whether or not there is predictive utility in polygenic scores calculated in the 100 test datasets. 

### The steps to get polygenic scores were run in Ubuntu, following the guide at https://choishingwan.github.io/PRS-Tutorial/  
Here is the specific code used:  

Clump the SNPs from the GWAS for calculating polygenic scores, using the unique SNP ID files and the new association results  
`plink --bfile spring87_uniqueSNPIDs --clump-p1 1 --clump-r2 0.2 --clump-kb 250 --clump /mnt/e/puma_lmm/sci_reports_revision/spring_87_reformat.lmm.assoc.txt --clump-snp-field rs --clump-field p_wald --allow-no-sex --allow-extra-chr --out spring_revision_clumped`

Extract clumped SNP IDs  
`awk 'NR!=1{print $3}' spring_revision_clumped.clumped > clumped.valid.snp`

Create a range list of p values for significance thresholds  
`echo "0.001 0 0.001" > range_list`  
`echo "0.05 0 0.05" >> range_list`  
`echo "0.1 0 0.1" >> range_list`  
`echo "0.2 0 0.2" >> range_list`  
`echo "0.3 0 0.3" >> range_list`  
`echo "0.4 0 0.4" >> range_list`  
`echo "0.5 0 0.5" >> range_list`  

Get the polygenic scores. A warning about variant ID mismatches will appear because of the clumping process
`plink --bfile spring87_uniqueSNPIDs --score /mnt/e/puma_lmm/sci_reports_revision/spring_87_reformat.lmm.assoc.txt 2 5 8 header --q-score-range range_list /mnt/e/puma_lmm/sci_reports_revision/spring_87_reformat.lmm.assoc.txt 2 13 header --extract clumped.valid.snp --allow-no-sex --allow-extra-chr --out spring_revision_PRS`  

The above line makes the .profile files used with linear models to calculate R^2 in revision_PRS.R
