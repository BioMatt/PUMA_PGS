# PUMA_PGS
Scripts used to calculate polygenic scores and estimate polygenic score accuracy in spring migration timing for the purple martin. Preprint here: https://www.biorxiv.org/content/10.1101/2022.09.14.508039v1

The script PGS_revision_random_sampling.R took genotype and covariate files in bed format and randomly created 100 training and test datasets with a 85/15 split. It calculated new principal component scores for the covariate files based on the training data to simulate a smaller GWAS

The script PGS_accuracy.sh was run on the Graham cluster of the Digital Research Alliance of Canada, and used gemma for GWAS and plink for calculating polygenic scores. Gemma first needs a relatedness matrix, then to run the GWAS. Plink then takes GWAS results, clumps SNPs based on linkage, and calculates polygenic scores at a range of adjusted p-values (0.001, 0.05, 0.1, 0.2, 0.3, 0.4, and 0.5). Plink was then used to calculate polygenic scores for the test individuals based on results from the training individuals.

The script PGS_revision_accuracy_test.R takes the .profile polygenic score output files from the training and test datasets, the covariate files, and the .fam phenotype files to assess whether or not there is predictive utility in polygenic scores calculated in the 100 test datasets. 
