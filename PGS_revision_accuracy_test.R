# R script based on guide at https://choishingwan.github.io/PRS-Tutorial/plink/
# However, the script is also meant to test the accuracy of PGS data with PUMA migration by looking at r^2 across different training models, and associating predicted phenotype deciles with real spring migration phenotypes.
# It draws ideas from Fuller 2020, at https://github.com/zfuller5280/CoralGenomes/blob/master/pgs_accuracy.R
library(tidyverse)


p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
# Initialize an empty dataframe
prs_results <- data.frame(Threshold = as.character(rep(NA, 100)), R2 = as.character(rep(NA, 100)), P = as.character(rep(NA, 100)), BETA = as.character(rep(NA, 100)), SE = as.character(rep(NA, 100)))

# Writing a loop to read in multiple .profile files, find the best p threshold for each of 100 iterations, and write the r^2 and p threshold information for each run to a data frame
for(run in 1:100){
  # Read in the phenotype file 
  fam_file <- paste0("../spring74_revision", run, ".fam")
  phenotype <- read.table(fam_file, col.names = c("FID", "IID", "foo1", "foo2", "foo3", "spring_migration"))
  # Remove everything except basin residency and FID since this was the phenotype tested for
  phenotype$foo1 <- phenotype$foo2 <- phenotype$foo3 <- NULL
  
  # Read in the covariates
  cov_file <- paste0("../spring74_", run, "_cov7.colony.sex.year.age.pc1.spdur.spdist.txt")
  covariate <- read.table("../spring74_1_cov7.colony.sex.year.age.pc1.spdur.spdist.txt", header=F, col.names = c("intercept", "colony", "sex", "year", "age", "PC1", "spdur", "spdist"))
  # Remove the intercept column
  covariate$intercept  <- NULL
  # Add in FID columns to match the other tables
  covariate$FID <- phenotype$FID
  covariate$IID <- phenotype$IID
  
  # Now merge the files
  pheno <- merge(phenotype, covariate, by=c("FID", "IID"))
  
  
  # We can then calculate the null model (model with PRS) using a linear regression 
  # (as height is quantitative)
  null.model <- lm(spring_migration~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
  # And the R2 of the null model is 
  null.r2 <- summary(null.model)$r.squared
  null.r2
  prs.result <- NULL
  
  for(i in p.threshold){
    # Go through each p-value threshold
    prs <- read.table(paste0("./training_profile_files/spring74_PRS_", run, ".",i,".profile"), header=T)
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
    # relevant columns
    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    # Now perform a linear regression on Height with PRS and the covariates
    # ignoring the FID and IID from our model
    model <- lm(spring_migration~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
    # model R2 is obtained as 
    model.r2 <- summary(model)$r.squared
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    prs.coef <- summary(model)$coeff["SCORE",]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    # We can then store the results
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
  }
  
  # Best result is:
  print(prs.result[which.max(prs.result$R2),])
  #prs_results <<- rbind(prs_results, prs.result[which.max(prs.result$R2),])
  prs_results[run, ] <- prs.result[which.max(prs.result$R2),]
  
  
  #q() # exit R
}

# Clean up the environment
rm(cov_file, fam_file, model.r2, null.r2, p.threshold, prs.beta, prs.coef, prs.p, prs.r2, prs.se, i, covariate, model, null.model, pheno, pheno.prs, phenotype, prs, prs.result, run)

# Add in a column to use as the X position for plotting r^2 values to show how the r2 results come from the same model
prs_results <- data.frame(model = rep("seven_covariates", nrow(prs_results)), prs_results)
prs_results$R2 <- as.numeric(prs_results$R2) # Make R2 numeric
# Plot accuracy of the PGS across different training datasets, as measured by r^2
ggplot(prs_results, aes(x = model, y = R2)) +
  geom_boxplot(notch=TRUE, width = 0.25) +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  ylim(0.1, 0.35) +
  theme_bw()

# Read in phenotype data for all individuals, training and test
# Remove everything except basin residency and FID since this was the phenotype tested for
phenotypes <- read_table("../spring87_uniqueSNPIDs.fam", col_names = c("FID", "IID", "foo1", "foo2", "foo3", "spring_migration")) %>% 
  select(-foo1, -foo2, -foo3)

# Create an empty dataframe for just test individuals from each iteration of the GWAS
test_individuals <- data.frame(FID = as.character(), IID = as.character(), PHENO = as.numeric(), CNT = as.numeric(), CNT2 = as.numeric(), SCORE = as.numeric(), data = as.character(), decile = as.integer(), spring_migration = as.integer(), iteration = as.integer())

# This loop is intended to find predicted phenotype deciles for each individual in the test data sets (with respect to the PGS values in the respective training dataset), and then attach that information to the real observed phenotype for the test individuals.
for(run in 1:100){
  # Check which threshold p value led to the lowest model p value for each run
  p_threshold <- prs_results[run, 2]
  # With the p_threshold and iteration information, read in the scores for the best p threshold for that run
  training_data <- read.table(paste0("spring74_PRS_test_", run, ".", p_threshold, ".profile"), header=T)
  # Add in a column describing the training data as created for training
  training_data$data <- rep("training", nrow(training_data))
  
  # Add a column to the training data to see the decile for each individual's PGS
  training_data$decile <- ntile(training_data$SCORE, 10)
  
  # Group the training data by decile
  training_decile <- training_data %>% 
    group_by(decile) %>% 
    summarise(min_PGS = min(SCORE, na.rm = TRUE), max_PGS = max(SCORE, na.rm = TRUE))
  
  # Also read in the scores for the test data, using the best p_threshold from the training data
  test_data <- read.table(paste0("spring74_PRS_test_", run, ".", p_threshold, ".profile"), header = TRUE)
  # Add in a column describing the test dataset as created for testing
  test_data$data <- rep("test", nrow(test_data))
  
  # Join all combinations of the data together
  all_join <- full_join(training_decile, test_data, by = character())
  # Filter for only the rows where the polygenic score is correctly placed between in a decile
  all_join <- dplyr::filter(all_join, SCORE >= min_PGS, SCORE <= max_PGS) %>% 
    select(FID, IID, PHENO, CNT, CNT2, SCORE, data, decile)
  
  
  # Bind the two datasets together
  combined_data <- rbind(training_data, all_join)
  
  # Add a column to the combined data to see the decile for each individual's PGS
  combined_data$decile <- ntile(combined_data$SCORE, 10)
  # Turn the combined_data dataframe into a tibble
  combined_data <- as_tibble(combined_data)
  
  # Attach the real phenotype to the combined_data 
  combined_data <- left_join(combined_data, phenotypes)
  
  # Filter for just the test individuals
  combined_data <- dplyr::filter(combined_data, data == "test")
  
  # Add a column showing which iteration of the program these test individuals are from
  combined_data$iteration <- rep(run, nrow(combined_data))
  
  # Add these test individuals to a dataframe of all test individuals to plot for prediction accuracy
  test_individuals <- rbind(test_individuals, combined_data)
}

# Clean up the environment
rm(run, p_threshold, training_data, test_data, combined_data, all_join)

# Plot prediction accuracy for PGS
PGS_test_plot <- ggplot(test_individuals, aes(x = decile, y = spring_migration, group= decile)) + 
  geom_boxplot(notch = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.1), alpha = 0.6) +
  ylab("Spring Migration Date") +
  xlab("Decile") +
  scale_x_continuous(breaks = 1:10, labels = c(as.character(1:10))) +
  theme_bw(base_size = 12)

# Write out the plot
ggsave(filename = "PGS_test_plot.pdf", plot = PGS_test_plot, dpi = 900)

