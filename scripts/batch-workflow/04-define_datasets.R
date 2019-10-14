
source("init.R", chdir = TRUE)

# Read German QC3 fam file
fam <- data.table::fread(sprintf("%s.fam", german_qc3_file_prefix),
                         col.names = c("FID", "IID", "PID", "MID", "SEX", "STATUS"))

# Define validation and training samples
validation_ids <- sample(1:fam[, .N], validation_size, replace = FALSE)
training_ids <- setdiff(1:fam[, .N], validation_ids)

# Save ids
saveRDS(validation_ids, validation_ids_rds_file)
saveRDS(training_ids, training_ids_rds_file)

# Create fam subsets
validation_fam <- fam[validation_ids]
training_fam <- fam[training_ids]

# Save fam subsets
saveRDS(validation_fam, validation_fam_rds_file)
saveRDS(training_fam, training_fam_rds_file)
