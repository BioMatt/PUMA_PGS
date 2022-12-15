# This script is used to look at the 87 individuals 

# Use genio to handle the plink bed, bim, and fam files because it has nice read and write functions.
library(genio)
library(SNPRelate)
library(tidyverse)

# Read the fam, bim, and bed files. Use individual names from the fam and SNP IDs from the bim to add context to the bed file
fam <- read_fam("spring87_uniqueSNPIDs.fam")
bim <- read_bim("spring87_uniqueSNPIDs.bim")
bed <- read_bed("spring87_uniqueSNPIDs.bed", m_loci = 4484670, names_loci = bim$id, n_ind = 87, names_ind = fam$fam)

# Set up a function to write out random subsets of individuals based on which plink fileset you need, how many subsamples you need, and how many individuals in each subsample
write_random_beds <- function(bed_file, bim_file, fam_file, iterations, sample_size) {
  for (i in 1:iterations) { # Loop through for as many file set iterations as you want to make.
    
    temp_fam <- fam[sample(nrow(fam), size = sample_size, replace = FALSE),] # Pull a certain number of individuals randomly from the fam file, with their phenotypes
    
    temp_bed <- bed[,temp_fam$fam] # Take those sample individuals and pull their genotypes from the bed file
    
    file_name <- paste0("spring74_revision", i)
    write_plink(file = file_name, X = temp_bed, bim = bim, fam = temp_fam) # Write out the subsampled bed and fam files, with the normal bim file because that does not rely on individual information.
  }
}


# Use the function just created to write out 100 subsampled files with 74 individuals. Choosing 74 because there are 87 individuals in the dataset, so we're using 85.06% of individuals as training data
write_random_beds(bed_file = bed, bim_file = bim, fam_file = fam, iterations = 100, sample_size = 74)


# Set up a function to write out test files, complementary to the subsetted files created with "write_random_beds". The files from that function will be used as training data, and the files from this function will be used as test data.
write_test_data <- function(bed_file, bim_file, fam_file, iterations) {
  for (i in 1:iterations) { # Loop through for as many file set iterations as you want to make.
    # Read the fam file already created
    training_filename <- paste0("spring74_revision", i, ".fam")
    training_fam <- read_fam(training_filename)
    
    # Pull individual names *not* in the training fam from the bed file. These are the test individuals.
    test_individuals <- setdiff(colnames(bed), unlist(c(training_fam)))
    
    temp_bed <- bed[,test_individuals] # Take those sample individuals and pull their genotypes from the bed file
    temp_fam <- fam[fam$fam %in% test_individuals,] # Do the same thing with the fam file, but use the individual IDs in the first row instead
    
    file_name <- paste0("spring74_test_", i)
    write_plink(file = file_name, X = temp_bed, bim = bim, fam = temp_fam) # Write out the subsampled bed and fam files, with the normal bim file because that does not rely on individual information.
  }
}

# Use the function just created to write out the test files for each of the 100 iterations created previously.
write_test_data(bed_file = bed, bim_file = bim, fam_file = fam, iterations = 100)

# This function is to read in the subsetted plink files, LD prune the subsetted SNP data, run a PCA, and update the covariate file with the new PCA information. Then write the new covariate file with consistent filenames.

write_covariates <- function(prefix, iterations){
  # Read in the overall phenotype file. This file does not change between iterations.
  fam_file <- read_table("spring87_uniqueSNPIDs.fam", col_names = c("ind", "fam", "father", "mother", "sex", "phenotype"))
  # Read in the covariate file. Set column names of the covariate file to the individual names. This file doesn't change between iterations.
  cov_file <- read_tsv("cov7.colony.sex.year.age.pc1.spdur.spdist.txt", col_names = c("intercept", "colony", "sex", "year", "age", "pc1", "spdur", "spdist")) %>% 
    mutate(ind = as.character(fam_file$ind))
  for (i in 1:iterations){
    # Set filenames for the three plink files worked on this iteration.
    temp_fam_name <- paste0(prefix, i, ".fam")
    temp_bed_name <- paste0(prefix, i, ".bed")
    temp_bim_name <- paste0(prefix, i, ".bim")
    temp_gds_name <- paste0(prefix, i, ".gds")
    
    # Read in the temporary fam file to be worked on for this iteration
    temp_fam <- read_table(temp_fam_name, col_names = c("ind", "fam", "father", "mother", "sex", "phenotype")) %>% 
      select(ind)
    
    # Filter the covariate file for just the individuals in the temporary fam file. Use left_join to also re-order the covariate file
    temp_cov <- left_join(temp_fam, cov_file, by = "ind") %>% 
      select(-ind)
    
    # Read in the subsetted plink fileset
    snpgdsBED2GDS(bed.fn = temp_bed_name, fam.fn = temp_fam_name, bim.fn = temp_bim_name, out.gdsfn = temp_gds_name, cvt.chr = "char")
    
    genofile <- snpgdsOpen(temp_gds_name, readonly = FALSE, allow.duplicate = TRUE)
    
    # Prune the subsetted group of SNPs
    train_snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only = FALSE, verbose = FALSE)
    train_snpset.id <- unlist(train_snpset)
    # Run the PCA
    train_pca <- snpgdsPCA(genofile, snp.id=train_snpset.id, num.thread=8,need.genmat=TRUE, autosome.only = FALSE, sample.id = as.factor(temp_fam$ind))
    # Replace the principal component 1 column with the new PC1
    temp_cov$pc1 <- train_pca$eigenvect[,1]
    # Close the genofile
    snpgdsClose(genofile)
    # Write out the new subsetted covariate file
    write_delim(x = temp_cov, file = paste0("spring74_", i, "_cov7.colony.sex.year.age.pc1.spdur.spdist.txt"), delim = "\t", col_names = FALSE)
  }
}

write_covariates(prefix = "spring74_revision", 100)



# Write a for loop in bash to run gemma. Copy and paste this code into Ubuntu to run a GWAS with a new relatedness matrix for each GWAS
for i in {1..100}
do
gemma -bfile /mnt/e/puma_lmm/sci_reports_revision/revision_accuracy/spring74_revision${i} -gk 1 -miss 1 -maf 0 -notsnp -o relatedness_matrix_${i}
gemma -bfile /mnt/e/puma_lmm/sci_reports_revision/revision_accuracy/spring74_revision${i} -k output/relatedness_matrix_${i}.cXX.txt -miss 1 -maf 0 -lmm 1 -c /mnt/e/puma_lmm/sci_reports_revision/revision_accuracy/spring74_${i}_cov7.colony.sex.year.age.pc1.spdur.spdist.txt -o spring74_${i}
plink --bfile /mnt/e/puma_lmm/sci_reports_revision/revision_accuracy/spring74_revision${i} --allow-extra-chr --clump-p1 1 --clump-r2 0.2 --clump-kb 250 --clump output/spring74_${i}.assoc.txt --clump-snp-field rs --clump-field p_wald --allow-no-sex --out spring74_clumped_${i}
awk 'NR!=1{print $3}' spring74_clumped_${i}.clumped > spring74_${i}.valid.snp
plink --bfile /mnt/e/puma_lmm/sci_reports_revision/revision_accuracy/spring74_revision${i} --allow-extra-chr --allow-no-sex --score output/spring74_${i}.assoc.txt 2 5 8 --q-score-range range_list output/spring74_${i}.assoc.txt 2 12 header --extract spring74_${i}.valid.snp --out spring74_PRS_${i}
done


# A loop to clump SNPs. 
for i in {1..100}
do
plink --bfile /mnt/e/puma_lmm/sci_reports_revision/revision_accuracy/spring74_revision${i} --allow-extra-chr --clump-p1 1 --clump-r2 0.2 --clump-kb 250 --clump output/spring74_${i}.assoc.txt --clump-snp-field rs --clump-field p_wald --allow-no-sex --out spring74_clumped_${i}
done

# Extract valid SNPs from the clumped SNP files
for i in {1..100}
do
awk 'NR!=1{print $3}' spring74_clumped_${i}.clumped > spring74_${i}.valid.snp
done

# Calculate polygenic scores
for i in {1..100}
do
plink --bfile /mnt/e/puma_lmm/sci_reports_revision/revision_accuracy/spring74_revision${i} --allow-extra-chr --allow-no-sex --score output/spring74_${i}.assoc.txt 2 5 8 --q-score-range range_list output/spring74_${i}.assoc.txt 2 12 header --extract spring74_${i}.valid.snp --out spring74_PRS_${i}
done
