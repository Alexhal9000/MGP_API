library(plumber)
library(Morpho)
library(ddsPLS)
library(Jovid)#devtools::install_github("https://github.com/J0vid/Jovid")
library(GenomicFeatures)#BiocManager::install("GenomicFeatures")
library(org.Mm.eg.db)#BiocManager::install("org.Mm.eg.db")
library(dplyr)
library(dbplyr)
library(future)
library(promises)
library(abind)
future::plan("multicore")

library(devtools)
library(ggplot2)
library(cmocean)
library(pals)
library(RColorBrewer)

library(Rvcg)
library(Morpho)
library(geomorph)
library(shapes)

library(foreach)
library(parallel)
library(doParallel)
library(org.Hs.eg.db)#BiocManager::install("org.Hs.eg.db")
library(GO.db)
library(biomaRt)
library(biomartr)
library(GOSemSim)#BiocManager::install("GOSemSim")
library(plinkFile)
library(LDlinkR)
library(paran)
library(gmodels)
library(jsonlite)

library(DBI)
library(RSQLite)

DEBUG <- FALSE

if(DEBUG){
  library(plot3D)
  library(rgl)
}

options("plumber.port" = 6234)

#save(combined.markers, DO.go, giga.pca, mutant.db, mutant.lms, shape.mean, Y, DO.probs, file = "/data/MGP_data/offline_data.Rdata")
#save(combined.markers, giga.pca, mutant.db, mutant.lms, shape.mean, Y, file = "~/shiny/shinyapps/MGP/shiny_data2.Rdata")


#local dirs
mmusculusEnsembl <- loadDb(file="~/Documents/MGP/jovid/MGP_API/ensemble.sqlite")
load("~/Documents/MGP/jovid/MGP_API/offline_data_no_biggie.Rdata")
load("~/Documents/MGP/jovid/MGP_API/cached.results.Rdata")
load("~/Documents/MGP/jovid/MGP_API/Mouse_Data/MGP_landmarks_mus.Rdata")
load("~/Documents/MGP/jovid/MGP_API/Mouse_Data/MGP_landmarks_mus_mutant.Rdata")

# New code using DBI/RSQLite:
DO_probs_DB <- dbConnect(RSQLite::SQLite(), "~/Documents/MGP/jovid/MGP_API/MGP_genotypes.sqlite")
DO_probs_DB_cleaner <- dbConnect(RSQLite::SQLite(), "~/Documents/MGP/jovid/MGP_API/Mouse_Data/MGP_genotypes_mus_A.sqlite")
DO_probs_DB_cleaner_M <- dbConnect(RSQLite::SQLite(), "~/Documents/MGP/jovid/MGP_API/Mouse_Data/MGP_genotypes_mus_M.sqlite")
DO_probs_DB_cleaner_F <- dbConnect(RSQLite::SQLite(), "~/Documents/MGP/jovid/MGP_API/Mouse_Data/MGP_genotypes_mus_F.sqlite")

# Count Markers:
# List all marker tables
marker_tables <- dbListTables(DO_probs_DB_cleaner)
# Print the number of markers (tables)
cat("Number of markers (tables):", length(marker_tables), "\n")
# For the first marker, print the number of subjects (rows)
first_marker <- marker_tables[1]
num_subjects <- dbGetQuery(DO_probs_DB_cleaner, paste0("SELECT COUNT(*) FROM '", first_marker, "'"))[[1]]
cat("Number of subjects for marker", first_marker, ":", num_subjects, "\n")


## ── offline annotation tables ─────────────────────────────────────────
base_path <- "~/Documents/MGP/jovid/MGP_API/humanGeneDBs"

load(file.path(base_path, "ensembl_gene_coords_GRCh37.RData"))   # object coords
ENSEMBL_COORDS_HUMAN <- coords ; rm(coords)

load(file.path(base_path, "symbol2ensembl_hs.RData"))            # SYMBOL2ENSEMBL
SYMBOL2ENSEMBL_HUMAN <- SYMBOL2ENSEMBL ; rm(SYMBOL2ENSEMBL)

load(file.path(base_path, "ensembl2symbol_hs.RData"))            # ENSEMBL2SYMBOL
ENSEMBL2SYMBOL_HUMAN <- ENSEMBL2SYMBOL ; rm(ENSEMBL2SYMBOL)

load(file.path(base_path, "go2symbol_hs.RData"))                 # GO2SYMBOL
GO2SYMBOL_HUMAN <- GO2SYMBOL ; rm(GO2SYMBOL)

## ── helper mappers (purely local look-ups) ────────────────────────────
map_symbol_to_ensembl <- function(sym)
  subset(SYMBOL2ENSEMBL_HUMAN, toupper(SYMBOL) %in% toupper(sym))

map_ensembl_to_symbol <- function(ens)
  subset(ENSEMBL2SYMBOL_HUMAN, ENSEMBL %in% ens)

map_go_to_symbol <- function(go)
  subset(GO2SYMBOL_HUMAN, GO %in% go, select = c("ENSEMBL","SYMBOL"))



load("~/Documents/MGP/jovid/MGP_API/combined_markers_cleaner.Rdata")

load("~/Documents/MGP/jovid/MGP_API/humanMGP-master/R/YHuman.Rdata")

# # ---------- process covariate data ---------------
# load("~/Documents/MGP/jovid/MGP_API/humanMGP-master/R/DENSEDATA_TANZ.RData")
# lm = "dense"
# covs = c("age", "sex", "height", "weight", "csize")
# # Number of landmarks
# nlandmarks <- 5629
# if (lm == "sparse") {
#   nsplandmarks <- 65
#   nlandmarks <- 2565
# }
# 
# # Covariate adjustment
# cat("\033[33m", "== Standardize shapes by covariates", "\033[0m", "\n")
# 
# # Only keep individuals that have covariate info
# sort_idx1 <- pheno.id %in% meta.cov$IID
# pheno.id <- pheno.id[sort_idx1]
# pheno.coeff <- pheno.coeff[sort_idx1,]
# 
# # Make sure shape data and meta data are in the same order
# sort_idx2 <- match(pheno.id, meta.cov$IID)
# meta.cov <- meta.cov[sort_idx2,]
# 
# # Standardization
# if ("none" %in% covs) {
#   # No standardization
#   pheno.coeff.adj <- pheno.coeff
#   pheno.avg.adj <- pheno.avg
# } else {
#   cov_idx <- match(tolower(covs), tolower(colnames(meta.cov)))
#   cov_keep <- meta.cov[, cov_idx, drop = FALSE]
# 
#   # Remove individuals with missing data
#   keep_id <- rowSums(is.na(cov_keep)) < 1
#   cov_keep <- cov_keep[keep_id, , drop = FALSE]
#   pheno.coeff <- pheno.coeff[keep_id, , drop = FALSE]
#   pheno.id <- pheno.id[keep_id]
# 
#   # Make geomorph dataframe
#   df <- geomorph.data.frame(cov_keep)
#   df$Shape <- pheno.coeff
# 
#   # Fit model
#   f1 <- as.formula(paste(names(df)[length(df)], paste(names(df)[1:(length(df) - 1)], collapse = " + "), sep = " ~ "))
#   mod <- procD.lm(f1, data = df, RRPP = TRUE, iter = 99, Parallel = 1)
# 
#   # Get residuals
#   pheno.coeff.adj <- mod$residuals + matrix(data = colMeans(df$Shape), nrow = nrow(pheno.coeff), ncol = ncol(pheno.coeff), byrow = TRUE)
#   pheno.avg.adj <- colMeans(pheno.coeff.adj)
#   print("Done with regressing out covariates")
# }
# 
# # Store the processed data in YHuman
# YHumanTANZDense <- pheno.coeff.adj
# pheno.id.human = pheno.id
# save(YHumanTANZDense, pheno.id.human, pca.eigstd, pca.eigvec, pheno.avg, file = "~/Documents/MGP/jovid/MGP_API/humanMGP-master/R/YHumanTANZDense.Rdata")
# # ---------------------------------------------------------------


load("~/Documents/MGP/jovid/MGP_API/humanMGP-master/R/YHumanDense.Rdata")
nlandmarks <- 5629
load("~/Documents/MGP/jovid/MGP_API/humanMGP-master/R/YHumanTANZ.Rdata")



#docker dirs
# setwd("/srv/shiny-server/")
# mmusculusEnsembl <- loadDb(file="ensemble.sqlite")
# load("shiny_data.Rdata")
# load("cached.results.Rdata")
# DO_probs_DB <- src_sqlite("MGP_genotypes.sqlite")

#genopheno deployment dirs
# mmusculusEnsembl <- loadDb(file="/data/MGP_data/ensemble.sqlite")
# load("/data/MGP_data/shiny_data.Rdata")
# load("/data/MGP_data/cached.results.Rdata")
# DO_probs_DB <- src_sqlite("/data/MGP_data/MGP_genotypes.sqlite")

# Create a named list to map pheno values to corresponding DO_probs_DB tables
pheno_tables <- list(
  "Y" = DO_probs_DB,
  "A_lm_raw" = DO_probs_DB_cleaner,
  "A_lm_gen" = DO_probs_DB_cleaner,
  "F_lm_gen" = DO_probs_DB_cleaner_F,
  "M_lm_gen" = DO_probs_DB_cleaner_M,
  "A_lm_gen_sex" = DO_probs_DB_cleaner,
  "A_lm_gen_sex_size" = DO_probs_DB_cleaner,
  "A_lm_gen_size" = DO_probs_DB_cleaner,
  "F_lm_gen_size" = DO_probs_DB_cleaner_F,
  "M_lm_gen_size" = DO_probs_DB_cleaner_M,
  "A_lm" = DO_probs_DB_cleaner,
  "F_lm" = DO_probs_DB_cleaner_F,
  "M_lm" = DO_probs_DB_cleaner_M, 
  "YHuman" = "NA",
  "YHumanTANZ" = "NA",
  "YHumanDense" = "NA",
  "YHumanTANZDense" = "NA"
  # Add more mappings as needed
)

# Create a named list to map marker values to corresponding DO_probs_DB tables
combined_markers <- list(
  "Y" = combined.markers,
  "A_lm_raw" = combined.markers.cleaner,
  "A_lm_gen" = combined.markers.cleaner,
  "F_lm_gen" = combined.markers.cleaner,
  "M_lm_gen" = combined.markers.cleaner,
  "A_lm_gen_sex" = combined.markers.cleaner,
  "A_lm_gen_sex_size" = combined.markers.cleaner,
  "A_lm_gen_size" = combined.markers.cleaner,
  "F_lm_gen_size" = combined.markers.cleaner,
  "M_lm_gen_size" = combined.markers.cleaner,
  "A_lm" = combined.markers.cleaner,
  "F_lm" = combined.markers.cleaner,
  "M_lm" = combined.markers.cleaner,
  "YHuman" = "NA",
  "YHumanTANZ" = "NA",
  "YHumanDense" = "NA",
  "YHumanTANZDense" = "NA"
  # Add more mappings as needed
)

# Create a named list to map marker values to corresponding scaling factors
corrective_scales <- list(
  "Y" = 67.626370731,
  "A_lm_raw" = 1.0,
  "A_lm_gen" = 232.0139,
  "F_lm_gen" = 232.0139,
  "M_lm_gen" = 232.0139,
  "A_lm_gen_sex" = 232.0139,
  "A_lm_gen_sex_size" = 232.0139,
  "A_lm_gen_size" = 232.0139,
  "F_lm_gen_size" = 232.0139,
  "M_lm_gen_size" = 232.0139,
  "A_lm" = 232.0139,
  "F_lm" = 232.0139,
  "M_lm" = 232.0139,
  "YHuman" = 1.0,
  "YHumanTANZ" = 1.0,
  "YHumanDense" = 1.0,
  "YHumanTANZDense" = 1.0
  # Add more mappings as needed
)


# Create a list with the names of the phenotypes
pheno_names_mouse <- list(
  "Y",
  "A_lm_raw",
  "A_lm_gen",
  "F_lm_gen",
  "M_lm_gen",
  "A_lm_gen_sex",
  "A_lm_gen_sex_size",
  "A_lm_gen_size",
  "F_lm_gen_size",
  "M_lm_gen_size",
  "A_lm",
  "F_lm",
  "M_lm"
  # Add more elements as needed
)

# Create a list with the names of the phenotype groups
pheno_names_mouse_groups <- list(
  "Manual" = c("Y"),
  "Automated" = c(
    "A_lm_raw",
    "A_lm_gen",
    "F_lm_gen",
    "M_lm_gen",
    "A_lm_gen_sex",
    "A_lm_gen_sex_size",
    "A_lm_gen_size",
    "F_lm_gen_size",
    "M_lm_gen_size",
    "A_lm",
    "F_lm",
    "M_lm")
  # Add more elements as needed
)

pheno_names_mouse_mutants <- list(
  "Manual" = c("mutant.db", "mutant.lms"),
  "Automated" = c("auto.mutant.db", "auto.mutant.lms")
)

pheno_names_human <- list(
  "YHuman",
  "YHumanTANZ",
  "YHumanDense",
  "YHumanTANZDense"
)


# Function to reverse PCA and approximate dense landmark data
reverse_pca_dense <- function(coeff_row) {
  # Ensure coeff_row is a matrix
  coeff_row <- as.matrix(coeff_row)
  
  # Initialize a matrix to store the approximated data
  approx_data <- matrix(0, nrow = nrow(coeff_row), ncol = nrow(pca.eigvec))
  
  # Loop through each row of coeff_row
  for (i in 1:nrow(coeff_row)) {
    approx_data[i, ] <- t(pca.eigvec %*% coeff_row[i, ]) + as.vector(pheno.avg)
  }
  
  return(approx_data)
}



# Function to remove rows from all_genes based on the first column element
remove_rows <- function(mat, bad_mat) {
  # Extract the first column from both matrices
  all_genes_first_col <- mat[, 1]
  bad_genes_first_col <- bad_mat[, 1]
  
  # Identify rows in all_genes to be removed
  rows_to_remove <- which(all_genes_first_col %in% bad_genes_first_col)
  
  # Remove the identified rows
  cleaned_mat <- mat[-rows_to_remove, ]
  
  return(cleaned_mat)
}

get_curated_gene_list <- function() {
  # Check if the curated gene list is already in memory
  if (!exists("curated_gene_list")) {
    # If not, compile the gene list
    keys <- keys(org.Mm.eg.db)
    all_genes <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys, columns = c("ENSEMBL", "SYMBOL"), keytype = "ENTREZID")))
    
    # Create rangenes.gene2symbol
    all_genes <- all_genes[, c("SYMBOL", "ENSEMBL")]
    
    # Remove bad genes
    bad_genes <- read.csv("~/Documents/MGP/jovid/MGP_API/bad_genes_mm.csv", stringsAsFactors = FALSE)
    all_genes <- remove_rows(all_genes, bad_genes)
    
    # Remove chromosome 2 genes
    bad_genes <- read.csv("~/Documents/MGP/jovid/MGP_API/chrom_two_genes_mm.csv", stringsAsFactors = FALSE)
    all_genes <- remove_rows(all_genes, bad_genes)
    
    # Save the curated gene list in memory
    curated_gene_list <- all_genes
    
    # Group by SYMBOL and take the first occurrence of each group
    curated_gene_list <- curated_gene_list %>%
      group_by(SYMBOL) %>%
      slice(1) %>%
      ungroup()
  }
  
  # Return the curated gene list
  return(curated_gene_list)
}



get_curated_gene_list_human <- function(pheno,
                                        batch_size    = 500,
                                        force_rebuild = FALSE,
                                        cache_env     = .GlobalEnv) {

  ## ----- helper ---------------------------------------------------- ##
  clr           <- function(code, txt) paste0("\033[", code, "m", txt, "\033[0m")
  json_filename <- switch(pheno,
                          "YHuman"          = "~/Documents/MGP/jovid/MGP_API/pitt_whitelist.json",
                          "YHumanDense"     = "~/Documents/MGP/jovid/MGP_API/pitt_whitelist.json",
                          "YHumanTANZ"      = "~/Documents/MGP/jovid/MGP_API/tanz_whitelist.json",
                          "YHumanTANZDense" = "~/Documents/MGP/jovid/MGP_API/tanz_whitelist.json",
                          NULL)
  cache_var     <- paste0("curated_gene_list_human_", ifelse(is.null(pheno),"default",pheno))


  ## ----- quick-load ------------------------------------------------ ##
  if (!force_rebuild && exists(cache_var, envir = cache_env))
    return(get(cache_var, envir = cache_env))

  if (!force_rebuild && !is.null(json_filename) && file.exists(json_filename)) {
    cat(clr(32,paste("Loading curated list from",basename(json_filename))),"\n")
    js  <- jsonlite::fromJSON(json_filename)
    res <- data.frame(SYMBOL=js$SYMBOL, ENSEMBL=js$ENSEMBL, stringsAsFactors=FALSE)
    assign(cache_var, res, envir = cache_env)
    return(res)
  }


  ## ----- build from scratch --------------------------------------- ##
  cat(clr(33,"Rebuilding curated gene list – chunked mode"),"\n")

  # dummy shape matrix (big enough → sub-setting never fails)
  dummy_Y <- matrix(0, nrow = 20000, ncol = 2)

  all_syms   <- unique(SYMBOL2ENSEMBL_HUMAN$SYMBOL)
  good_syms  <- character(0)

  batches <- split(all_syms, ceiling(seq_along(all_syms)/batch_size))
  pb <- txtProgressBar(min = 0, max = length(batches), style = 3)

  start_time <- Sys.time()  # Record the start time

  for (i in seq_along(batches)) {
    batch <- batches[[i]]

    batch_ok <- tryCatch({
      res <- getHumanGenomeX_curate(GOterm = batch, pheno = pheno)

      survived <- unlist(res)
      # survived must be present in batch to be considered survived
      survived <- survived[survived %in% batch]
      good_syms <- c(good_syms, survived)
      TRUE
    }, error = function(e) { FALSE })

    setTxtProgressBar(pb, i)

    # Estimate time remaining
    elapsed_time <- Sys.time() - start_time
    avg_time_per_batch <- elapsed_time / i
    batches_left <- length(batches) - i
    estimated_remaining_time <- avg_time_per_batch * batches_left
    cat("\nProcessed", i*batch_size, "out of", length(all_syms), "genes - (survived:", length(good_syms), ")\n")
    cat("\nEstimated time remaining:", round(estimated_remaining_time, 2), "seconds\n")
  }
  close(pb)

  good_syms <- unique(good_syms)
  cat(clr(32, paste(length(good_syms),"symbols passed the full pipeline")),"\n")


  ## ----- SYMBOL → ENSEMBL mapping --------------------------------- ##
  curated <- SYMBOL2ENSEMBL_HUMAN[SYMBOL2ENSEMBL_HUMAN$SYMBOL %in% good_syms, ]
  #curated <- unique(na.omit(AnnotationDbi::select(
  #             org.Hs.eg.db,
  #              keys    = good_syms,
  #              columns = c("SYMBOL","ENSEMBL"),
  #              keytype = "SYMBOL")))


  ## ----- save + cache --------------------------------------------- ##
  if (!is.null(json_filename)) {
    jsonlite::write_json(
      list(creation_date = Sys.time(),
           pheno         = pheno,
           total_genes   = nrow(curated),
           SYMBOL        = curated$SYMBOL,
           ENSEMBL       = curated$ENSEMBL),
      json_filename, pretty = TRUE)
    cat(clr(32,paste("Saved curated list to",basename(json_filename))),"\n")
  }

  assign(cache_var, curated, envir = cache_env)
  curated    # final value returned
}



filter_bad_genes <- function(original_genes) {
  
  # Remove bad genes
  bad_genes <- read.csv("~/Documents/MGP/jovid/MGP_API/bad_genes_mm.csv", stringsAsFactors = FALSE)
  # Extract the first column from both matrices
  all_genes_third_col <- original_genes[, "SYMBOL"]
  bad_genes_first_col <- bad_genes[, 1]
  
  # Identify rows in all_genes to be removed
  rows_to_remove <- all_genes_third_col %in% bad_genes_first_col
  
  # Remove the identified rows
  clean_genes <- original_genes[!rows_to_remove, ]
  
  # Return the curated gene list
  return(clean_genes)
}

filter_chrom_two <- function(original_genes) {
  
  # Remove bad genes
  bad_genes <- read.csv("~/Documents/MGP/jovid/MGP_API/chrom_two_genes_mm.csv", stringsAsFactors = FALSE)
  # Extract the first column from both matrices
  all_genes_third_col <- original_genes[, "SYMBOL"]
  bad_genes_first_col <- bad_genes[, 1]
  
  # Identify rows in all_genes to be removed
  rows_to_remove <- all_genes_third_col %in% bad_genes_first_col
  
  # Remove the identified rows
  clean_genes <- original_genes[!rows_to_remove, ]
  
  # Return the curated gene list
  return(clean_genes)
}


cumulative_percent <- function(arr, percentage) {
  total_sum <- sum(arr)
  sorted_arr <- sort(arr, decreasing = TRUE)
  cumulative_sum <- 0
  for (i in 1:length(sorted_arr)) {
    cumulative_sum <- cumulative_sum + sorted_arr[i]
    if (cumulative_sum >= percentage * total_sum) {
      return(sorted_arr[i])
    }
  }
}

# Function to compute cosine distance
cos_dist <- function(v1, v2, t) {
  magnitudes_v1 <- sqrt((v1[,1]^2) + (v1[,2]^2) + (v1[,3]^2))
  magnitudes_v2 <- sqrt((v2[,1]^2) + (v2[,2]^2) + (v2[,3]^2))
  threshold_v1 = 0.003
  threshold_v2 = 0.003 
  
  # Cumulative percent rank 90% cutoff
  percent_cutoff <- cumulative_percent(magnitudes_v1, 0.8)
  mask <- magnitudes_v1 >= percent_cutoff
  
  print(percent_cutoff)
  
  # make the small noisy v2 magnitudes zero
  percent_cutoff_b <- cumulative_percent(magnitudes_v2, 0.8)
  zero_v2_mask <- magnitudes_v2 < percent_cutoff_b
  v2[zero_v2_mask, ] = c(0.0, 0.0, 0.0)
  
  #print(percent_cutoff_b)
  
  #mask <- (magnitudes_v1 >= (max(magnitudes_v1) - 3*sd(magnitudes_v1))) #& magnitudes_v2 >= threshold_v1)
  #print(max(magnitudes_v1)*0.3)
  #print(max(magnitudes_v2)*0.3)
  if(sum(mask, na.rm = TRUE)==0){
    return(0.0)
  }else{
    proc_dist <- sqrt(sum((magnitudes_v1[mask] - magnitudes_v2[mask])^2))
    dot_prod <- sum(v1[mask,] * v2[mask,])
    mag_v1 <- sqrt(sum(v1[mask,]^2))
    mag_v2 <- sqrt(sum(v2[mask,]^2))
    #print(paste("We keep ",sum(mask, na.rm=TRUE), " out of ",nrow(v1)," landmarks."))
    #print(paste("We zero out ",sum(zero_v2_mask, na.rm=TRUE), " out of ",nrow(v2)," V2 landmarks."))
    #print(proc_dist)
    cos_dis <- abs(dot_prod / (mag_v1 * mag_v2))
    return(cos_dis) #(1-(proc_dist))*cos_dis
  }
}

reversePCA <- function(Yshapes, pca_result, use_standardized_PCA) {
  if(use_standardized_PCA){
    # Reverse PCA standardization
    for(i in 1:length(pca_result$sdev)) {
      Yshapes[, i] <- Yshapes[, i] * pca_result$sdev[i]
    }
    
    # Reconstruct the original data
    reconstructed <- Yshapes %*% t(pca_result$rotation)
    reconstructed <- reconstructed + matrix(pca_result$center, nrow = nrow(reconstructed), ncol = ncol(reconstructed), byrow = TRUE)
    
    return(reconstructed)
  }else {
    return(Yshapes)
  }
}

standardizePCA <- function(Yshapes, use_standardized_PCA, x) {
  if(use_standardized_PCA){
    # Use PCA as phenotype
    Y_pca_result <- prcomp(Yshapes, center = TRUE, scale = FALSE)
    scores <- Y_pca_result$x
    
    # Calculate the proportion of variance explained by each component
    prop_var <- Y_pca_result$sdev^2
    prop_var <- prop_var / sum(prop_var)
    
    # Calculate the cumulative proportion of variance explained
    cum_var <- cumsum(prop_var)
    
    # Determine the threshold for x% of accumulated variability
    threshold <- 0.01 * x # Replace x with your desired percentage
    
    # Find the components that fall below the threshold
    components_to_zero_out <- which(cum_var > threshold)
    
    # Zero out the variability in the components below the threshold
    for(i in components_to_zero_out) {
      scores[, i] <- 0
    }
    
    # Normalize the remaining scores by dividing by the standard deviation
    for(i in 1:length(Y_pca_result$sdev)) {
      if(!(i %in% components_to_zero_out)) {
        scores[, i] <- scores[, i] / Y_pca_result$sdev[i]
      }
    }
    
    return(list(Yshapes = scores, Y_pca_result = Y_pca_result))
  }else {
    return(list(Yshapes = Yshapes, Y_pca_result = "NA"))
  }
}


# ---------------------------------------------------------------------
# Fast gene-filter for the curation step only
# ---------------------------------------------------------------------
getHumanGenomeX_curate <- function(GOterm,
                                   pheno,
                                   window    = 0,
                                   cache_env = .GlobalEnv) {

  ## ---------- genotype BIM (cached) ---------------------------------
  mntpath <- "~/Documents/MGP/jovid/MGP_API/humanMGP-master/R/"
  loadpath2 <- paste0(mntpath, "Genetics/")
  bfile <- if (pheno %in% c("YHuman", "YHumanDense"))
             paste0(loadpath2, "PITT/Marazita_imputed_qc_prune_rmrel")
           else
             paste0(loadpath2, "TANZ/Spritz_imputed_qc_prune_rmrel")

  bim_var <- paste0("bim_cache_", pheno)
  if (!exists(bim_var, envir = cache_env))
    assign(bim_var, readBIM(bfile), envir = cache_env)
  geno.bim <- get(bim_var, envir = cache_env)



  ## ---------- SYMBOL ⇢ ENSEMBL --------------------------------------
  map <- SYMBOL2ENSEMBL_HUMAN
  sel <- map[map$SYMBOL %in% GOterm | map$ENSEMBL %in% GOterm, ]
  if (nrow(sel) == 0) stop("No SYMBOL/ENSEMBL match found")



  ## ---------- coordinates -------------------------------------------
  coords <- ENSEMBL_COORDS_HUMAN[
             ENSEMBL_COORDS_HUMAN$ensembl_gene_id %in% sel$ENSEMBL, ]



  # clean & unify chromosome names
  coords$chromosome_name <- gsub("HSCHR", "", coords$chromosome_name)
  coords$chromosome_name <- sub("_.*", "", coords$chromosome_name)
  coords$chromosome_name <- gsub("X", "23", coords$chromosome_name)
  coords$chromosome_name <- gsub("Y", "24", coords$chromosome_name)
  coords <- coords[!grepl("HG", coords$chromosome_name), ]



  ## ---------- SNP presence test -------------------------------------
  survived <- character(0)
  for (g in unique(sel$SYMBOL)) {
    ens <- sel$ENSEMBL[match(g, sel$SYMBOL)]
    cc  <- coords[coords$ensembl_gene_id == ens, ]
    if (nrow(cc) == 0) next
    
    # --- handle chromosome conversion (consistent with existing code style) ---
    chr <- suppressWarnings(as.numeric(cc$chromosome_name[1]))
    if (is.na(chr)) next  # Skip genes on non-standard chromosomes

    # --- collapse all transcripts of the gene to a single range ----------
    gene_start <- min(cc$transcript_start, na.rm = TRUE)
    gene_end   <- max(cc$transcript_end,   na.rm = TRUE)

    # --- SNP-presence test (no vector recycling warnings) ----------------
    hit <- any(with(geno.bim[geno.bim$chr == chr, ],
                    bps > (gene_start - window) &
                    bps < (gene_end   + window)))
    if (hit) survived <- c(survived, g)
  }
  survived <- unique(survived)


  as.list(survived)
}



# Extract the human genes of interest and reduce them with PCA, returns X

getHumanGenomeX <- function(GOterm, gene_type, pheno, Yshapes, single_match=FALSE, cache_env = .GlobalEnv){
  
  # manual mode
  # GOterm = selection.vector
  # gene_type = "custom"
  # Yshapes = YHuman
  # ---
  
  
  mntpath = "~/Documents/MGP/jovid/MGP_API/humanMGP-master/R/"
  loadpath1 = paste(mntpath,"HumanMGP/Data/",sep="")
  loadpath2 = paste(mntpath,"Genetics/",sep="")
  pls = "cov"
  window = 0
  ncomp = 2
  npc = 1
  signif = FALSE
  nperm = 99
  ncores = 10

  # Genotypes
  if (pheno == "YHuman" || pheno == "YHumanDense"){ 
    bfile = paste(loadpath2,"PITT/Marazita_imputed_qc_prune_rmrel",sep="") 
    pheno.id = pheno.id.human
  } else if (pheno == "YHumanTANZ" || pheno == "YHumanTANZDense"){ 
    bfile = paste(loadpath2,"TANZ/Spritz_imputed_qc_prune_rmrel",sep="") 
    pheno.id = pheno.id.human.TANZ
  }
  
  ## ---------- genotype BIM (cached) ---------------------------------
  bim_var <- paste0("bim_cache_", pheno)
  if (!exists(bim_var, envir = cache_env))
    assign(bim_var, readBIM(bfile), envir = cache_env)
  geno.full.bim <- get(bim_var, envir = cache_env)


  #geno.full.fam =  readFAM(bfile)
  GRCh = 37
  
  
  ## GENE/GO SEARCH TERM ####
  cat("\033[33m", "=== Match Genes to Search Term", "\033[0m", "\n")
  if (sum(tolower(GOterm) %in% tolower(GO2SYMBOL_HUMAN$GO_TERM) | tolower(GOterm) %in% tolower(GO2SYMBOL_HUMAN$GO))) {
    # Use offline GO database
    go_matches <- GO2SYMBOL_HUMAN[tolower(GO2SYMBOL_HUMAN$GO_TERM) %in% tolower(GOterm) | 
                                  tolower(GO2SYMBOL_HUMAN$GO) %in% tolower(GOterm), ]
    go2symbol <- data.frame(
      ENSEMBL = go_matches$ENSEMBL,
      SYMBOL = go_matches$SYMBOL,
      stringsAsFactors = FALSE
    )
    go2symbol <- unique(na.omit(go2symbol))
  }

  # # List all IDs of GO terms available in org.Hs object, and reduce to biological processes
  # GO_ID = toTable(org.Hs.egGO)
  # GO_ID = unique(GO_ID$go_id[GO_ID$Ontology == "BP"])
  # # Match GO IDs with process names
  # GO_all = toTable(GOTERM)
  # GO_BP = GO_all[GO_all$Ontology == "BP",]
  # GO_BP = GO_BP[match(GO_ID,GO_BP$go_id),]
  
  # # ---------------------save CSV of human ensembl and GO dbs---------------------------
  # # Initialize a new column for the number of genes
  # GO_BP$numGenes <- 0
  # for(nn in 1:nrow(GO_BP)){
  #     
  #     GO = GO_BP$go_id[nn]  
  #   
  #     suppressMessages({go2symbol = unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = GO, columns = c("ENSEMBL", "SYMBOL"), keytype = "GO")[,-2:-3]))})
  #     
  #     numGenes = nrow(go2symbol)
  #     
  #     # Update the numGenes column
  #     GO_BP$numGenes[nn] <- numGenes
  #     
  #     print(paste0(nn, ", ", numGenes))
  #     
  # }
  # 
  # write.csv(GO_BP, file = "~/Documents/MGP/jovid/MGP_API/DO_go_human.csv", row.names = FALSE)
  # 
  # # Load the necessary libraries
  # library(AnnotationDbi)
  # library(org.Hs.eg.db)
  # 
  # # Define the columns you want to retrieve
  # desired_columns <- c("ENSEMBL", "GO", "EVIDENCE", "ONTOLOGY", "SYMBOL", "ENTREZID")
  # 
  # # Retrieve data for Biological Process (BP) ontology
  # bp_data <- AnnotationDbi::select(org.Hs.eg.db,
  #                                  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  #                                  columns = desired_columns,
  #                                  keytype = "ENTREZID")
  # 
  # # Filter for Biological Process ontology
  # bp_data_bp_only <- bp_data[bp_data$ONTOLOGY == "BP", ]
  # 
  # # Remove NA values and get unique entries
  # bp_data_unique <- unique(na.omit(bp_data_bp_only))
  # 
  # # Order the columns as per your request
  # bp_data_ordered <- bp_data_unique[, desired_columns]
  # write.csv(bp_data_ordered, file = "~/Documents/MGP/jovid/MGP_API/human_all_GO_genes.csv", row.names = FALSE)
  # 
  # 
  # # Retrieve data with the specified columns
  # data <- AnnotationDbi::select(org.Hs.eg.db,
  #                               keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
  #                               columns = c("ENTREZID", "SYMBOL", "ENSEMBL", ),
  #                               keytype = "SYMBOL")
  # 
  # 
  # 
  # # Remove duplicates to get unique SYMBOL entries
  # unique_data <- data[!duplicated(data$SYMBOL), ]
  # 
  # write.csv(unique_data, file = "~/Documents/MGP/jovid/MGP_API/humanEnsembl.csv", row.names = FALSE)
  # ------------------------------------------------
  

  # Biological process (NAME / ID) - direct offline lookup
  if (sum(tolower(GOterm) %in% tolower(GO2SYMBOL_HUMAN$GO_TERM) | tolower(GOterm) %in% tolower(GO2SYMBOL_HUMAN$GO))) {
    # Use offline GO database directly
    go_matches <- GO2SYMBOL_HUMAN[tolower(GO2SYMBOL_HUMAN$GO_TERM) %in% tolower(GOterm) | 
                                  tolower(GO2SYMBOL_HUMAN$GO) %in% tolower(GOterm), ]
    go2symbol <- data.frame(
      ENSEMBL = go_matches$ENSEMBL,
      SYMBOL = go_matches$SYMBOL,
      stringsAsFactors = FALSE
    )
    go2symbol <- unique(na.omit(go2symbol))
  # Gene (NAME / ID) - FAST BATCH VERSION
  } else {
    
    # Pre-categorize genes using offline databases
    symbol_genes = GOterm[toupper(GOterm) %in% toupper(SYMBOL2ENSEMBL_HUMAN$SYMBOL)]
    ensembl_genes = GOterm[toupper(GOterm) %in% toupper(SYMBOL2ENSEMBL_HUMAN$ENSEMBL)]
    
    # Batch queries (only 2 total queries max)
    all_results = data.frame()
    
    # Function to perform batch querying and mapping
    perform_batch_query <- function(genes, keytype, column_names) {
      if (length(genes) > 0) {
        cat("\033[32m", paste("Offline lookup for", length(genes), keytype, "genes"), "\033[0m", "\n")
        
        if (keytype == "SYMBOL") {
          result <- map_symbol_to_ensembl(genes)
        } else if (keytype == "ENSEMBL") {
          result <- SYMBOL2ENSEMBL_HUMAN[SYMBOL2ENSEMBL_HUMAN$ENSEMBL %in% genes, ]
        }
        return(result)
      }
      return(data.frame())
    }
    
    # Batch 1: SYMBOL genes
    result <- perform_batch_query(symbol_genes, "SYMBOL", c("ENSEMBL", "SYMBOL"))
    if (nrow(result) > 0) {
      result$INPUT_NAME <- result$SYMBOL
      all_results = dplyr::bind_rows(all_results, result)
    }
    
    # Batch 2: ALIAS genes
    # result <- perform_batch_query(alias_genes, "ALIAS", c("ENSEMBL", "SYMBOL"))
    # if (nrow(result) > 0) {
    #   alias_map <- setNames(ALIAS_all$alias_symbol, ALIAS_all$gene_id)
    #   gene_to_input <- setNames(alias_genes, alias_genes)
    #   result$INPUT_NAME <- gene_to_input[rownames(result)]
      
    #   for (alias_gene in alias_genes) {
    #     mask <- result$SYMBOL %in% ALIAS_all$alias_symbol[toupper(ALIAS_all$alias_symbol) == toupper(alias_gene)]
    #     if (any(mask)) result$INPUT_NAME[mask] <- alias_gene
    #   }
      
    #   all_results = dplyr::bind_rows(all_results, result)
    # }
    
    # Batch 3: ENSEMBL genes
    result <- perform_batch_query(ensembl_genes, "ENSEMBL", c("ENSEMBL", "SYMBOL"))
    if (nrow(result) > 0) {
      result$INPUT_NAME <- ensembl_genes[match(result$ENSEMBL, ensembl_genes)]
      all_results = dplyr::bind_rows(all_results, result)
    }
    
    # Report missing genes
    found_inputs = unique(all_results$INPUT_NAME)
    missing_genes = GOterm[!GOterm %in% found_inputs]
    if(length(missing_genes) > 0) {
      cat("\033[33m", paste("No match found for:", paste(missing_genes, collapse=", ")), "\033[0m", "\n")
    }
    
    # Check if we found anything
    if(nrow(all_results) == 0) {
      stop("No match was found for any genes")
    }
    
    # Create the final mapping preserving input names
    go2symbol = data.frame(
      ENSEMBL = all_results$ENSEMBL,
      SYMBOL = all_results$INPUT_NAME,  # Use original input names
      stringsAsFactors = FALSE
    )
    go2symbol = unique(na.omit(go2symbol))

    cat("Number of unique genes found in the database: ", length(unique(go2symbol$SYMBOL)), "\n")
    
    # Non coding [to be added] 
    #gene.list2 <- chdr.list[grepl("TCONS_", chdr.list$gene_id2, fixed = TRUE),]
    #db2 <- TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts
    #symbol2info2 <- transcriptsBy(db2, by=c("gene"), use.names=FALSE)[gene.list2]
    #symbol2info2 <- symbol2info2@unlistData
  }
    
  
  
  
  ## GENES TX START-END SITE ####   
  cat("\033[33m", "==== Map SNPs to Genes", "\033[0m", "\n")
  #ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=GRCh)
  #symbol2info = getBM(attributes = c('ensembl_gene_id','chromosome_name','transcript_start', 'transcript_end', 'hgnc_symbol'), 
  #                    mart = ensembl, 
  #                    filters = 'ensembl_gene_id', 
  #                    values = go2symbol$ENSEMBL)

  # Local retrieval of gene marker data
  symbol2info = ENSEMBL_COORDS_HUMAN[ENSEMBL_COORDS_HUMAN$ensembl_gene_id %in% go2symbol$ENSEMBL, ]
  
  # Create mapping from ENSEMBL back to original input symbols
  original_symbol_map = setNames(go2symbol$SYMBOL, go2symbol$ENSEMBL)
  
  # Replace empty gene names
  tmp.ind = symbol2info$hgnc_symbol==""
  tmp.symbol = go2symbol$SYMBOL[match(symbol2info$ensembl_gene_id[tmp.ind],go2symbol$ENSEMBL)]
  symbol2info$hgnc_symbol[tmp.ind] = tmp.symbol 
  
  # Add original input name to each row
  symbol2info$original_input = original_symbol_map[symbol2info$ensembl_gene_id]
  
  # remove weird chr names
  chr = symbol2info$chromosome_name
  chr = gsub("HSCHR","",chr);  chr = sapply(strsplit(chr,"_"),"[[",1)
  symbol2info$chromosome_name = chr
  if (length(grep(symbol2info$chromosome_name, pattern = "HG")) > 0) { symbol2info = symbol2info[-grep(symbol2info$chromosome_name, pattern = "HG"),] }
  
  # Rename sex chromosomes
  symbol2info$chromosome_name = gsub("X","23",symbol2info$chromosome_name)
  symbol2info$chromosome_name = gsub("Y","24",symbol2info$chromosome_name)
  
  # Skip if no match is found
  if (dim(symbol2info)[1] == 0) { stop("No genes could be mapped to GO search term") }
  
  # SMART GENE NAME RESOLUTION
  final_genes = data.frame()

  cat("Number of unique genes: ", length(unique(symbol2info$original_input)), "\n")
  
  for(original_gene in unique(symbol2info$original_input)) {
    # Get all Biomart genes that came from this original input
    biomart_rows = symbol2info[symbol2info$original_input == original_gene, ]
    biomart_genes = unique(biomart_rows$hgnc_symbol)
    
    cat("\033[36m", paste("Original:", original_gene, "-> Biomart found:", paste(biomart_genes, collapse=", ")), "\033[0m", "\n")
    
    if(length(biomart_genes) == 1) {
      # Case 1: One gene returned -> keep original name as proxy
      biomart_rows$final_symbol = original_gene
      cat("\033[32m", paste("  Keeping original name:", original_gene), "\033[0m", "\n")
      
    } else if(length(biomart_genes) > 1 && single_match == FALSE) {
      # Case 2: Multiple genes returned
      if(original_gene %in% biomart_genes) {
        # Case 2a: Original name is among returned genes -> keep ALL with their Biomart names
        biomart_rows$final_symbol = biomart_rows$hgnc_symbol
        cat("\033[32m", paste("  Original found among results -> keeping all:", paste(biomart_genes, collapse=", ")), "\033[0m", "\n")
        
      } else {
        # Case 2b: Original name NOT among returned -> keep original as proxy for first + extras
        for(i in 1:length(biomart_genes)) {
          gene_rows = biomart_rows[biomart_rows$hgnc_symbol == biomart_genes[i], ]
          
          if(i == 1) {
            # First gene: use original name as proxy
            gene_rows$final_symbol = original_gene
            cat("\033[32m", paste("  Using", original_gene, "as proxy for:", biomart_genes[i]), "\033[0m", "\n")
          } else {
            # Additional genes: keep their Biomart names
            gene_rows$final_symbol = biomart_genes[i]
            cat("\033[32m", paste("  Adding extra gene:", biomart_genes[i]), "\033[0m", "\n")
          }
          
          # Add to final results
          if(i == 1) {
            processed_rows = gene_rows
          } else {
            processed_rows = rbind(processed_rows, gene_rows)
          }
        }
        biomart_rows = processed_rows
      }
    } else if(length(biomart_genes) > 1 && single_match == TRUE) {
      if(original_gene %in% biomart_genes) {
        # Case 3: Multiple genes returned and single_match is TRUE
        biomart_rows$final_symbol = original_gene
        cat("\033[32m", paste("  Keeping original name:", original_gene), "\033[0m", "\n")
      } else {
        # Case 4: Multiple genes returned and single_match is TRUE but original gene is not among them
        biomart_rows$final_symbol = biomart_genes[1]
        cat("\033[32m", paste("  Keeping first gene:", biomart_genes[1]), "\033[0m", "\n")
      }
    }
    
    final_genes = rbind(final_genes, biomart_rows)
  }

  
  # Now use final_symbol for the rest of the processing
  symbol2info = final_genes
  
  # Select maximum transcript size per gene using final symbols
  n = length(unique(symbol2info$final_symbol))

  cat("Number of unique genes post multiple matches filtering: ", n, "\n")
  
  genvar.name = rep(NA, n)
  genvar.start = rep(NA,  n)
  genvar.end = rep(NA,  n)
  chr.name = rep(NA,  n)
  for(i in 1:n){
    ind = symbol2info$final_symbol == unique(symbol2info$final_symbol)[i]
    
    genvar.name[i] = unique(symbol2info$final_symbol[ind])
    genvar.start[i] = min(symbol2info$transcript_start[ind])
    genvar.end[i] = max(symbol2info$transcript_end[ind])
    chr.name[i] = unique(symbol2info$chromosome_name[ind])
  }
  seq.info = data.frame(mgi_symbol = genvar.name, chromosome_name = chr.name, start_position = genvar.start, end_position = genvar.end)
  
  ## FIND MARKERS WITHIN WINDOW ####
  cat("\033[33m", "==== Finding markers within window", "\033[0m", "\n")
  seq.indexes = matrix(nrow = 0, ncol = 3)
  for (i in 1 : dim(seq.info)[1]) {
    tmp.markers = geno.full.bim[which(geno.full.bim$chr == as.numeric(seq.info$chromosome_name[i])),]
    tmp.dist1 = (tmp.markers$bps) > (as.numeric(seq.info$start_position[i]) - window)
    tmp.dist2 = (tmp.markers$bps) < (as.numeric(seq.info$end_position[i]) + window)
    ind = which(tmp.dist1 & tmp.dist2)
    
    if (length(ind) > 0) { seq.indexes = rbind(seq.indexes,as.matrix(cbind(seq.info$mgi_symbol[i],tmp.markers$vid[ind],tmp.markers$bps[ind]))) }
  }
  # Make dataframe
  seq.indexes = data.frame(gene = seq.indexes[,1], snp = seq.indexes[,2], pos = seq.indexes[,3])
  
  # Skip if no markers are found within window
  if (dim(seq.indexes)[1] == 0) { stop("No SNPs could be mapped to GO search term") }
  
  
  ## REDUCE GENOTYPES IN PLINK #### 
  plink_path = path.expand("~/Documents/MGP/jovid/MGP_API/humanMGP-master/R/Plink/plink")
  tmp_path = paste(mntpath,"tmp_mgp/",sep = "/") 
  if (!dir.exists(tmp_path)){ dir.create(tmp_path) }
  
  # Temporarily export SNPs to keep
  rnd_i = sample(1:1e6,1)
  tmp = seq.indexes$snp
  write.csv(tmp, file = paste(tmp_path,"tmp_snp_list_",rnd_i,".txt",sep=""),row.names = F, quote = F)
  
  # Run PLINK
  # Extract SNPs + remove related individuals
  system2(plink_path, args = paste(" --bfile ", bfile, " --threads ", ncores, " --extract ", tmp_path, "tmp_snp_list_", rnd_i, ".txt --make-bed --out ", tmp_path, "tmp_geno_", rnd_i, " > NUL", sep=""))
  
  # Import geno files
  geno.bim = readBIM(paste(tmp_path,"tmp_geno_",rnd_i,sep=""))
  geno.fam = readFAM(paste(tmp_path,"tmp_geno_",rnd_i,sep=""))
  geno.bed = readBED(paste(tmp_path,"tmp_geno_",rnd_i,sep="")); geno.bed = as.data.frame(geno.bed)
  
  # Remove tmp files
  unlink(paste(tmp_path,"tmp_snp_list_",rnd_i,".txt",sep=""))
  unlink(paste(tmp_path,"tmp_geno_",rnd_i,".bed",sep=""))
  unlink(paste(tmp_path,"tmp_geno_",rnd_i,".bim",sep=""))
  unlink(paste(tmp_path,"tmp_geno_",rnd_i,".fam",sep=""))
  unlink(paste(tmp_path,"tmp_geno_",rnd_i,".log",sep=""))
  
  
  
  
  ## CHECK ORDER OF GENO AND PHENO ####
  ## ENSURE MATCHING SUBJECTS IN GENOTYPE AND SHAPE DATA ####
  common_subjects = intersect(geno.fam$iid, pheno.id)
  
  # Filter Yshapes and pheno.id
  pheno_indices = match(common_subjects, pheno.id)
  Yshapes = Yshapes[pheno_indices, ]
  pheno.id = pheno.id[pheno_indices]
  
  # Filter geno.fam and geno.bed
  geno_indices = match(common_subjects, geno.fam$iid)
  geno.fam = geno.fam[geno_indices, , drop=F]
  geno.bed = geno.bed[geno_indices, , drop=F]
  

  ## REPLACE MISSING GENO WITH MAJOR GENOTYPE #### 
  geno.0 = colSums(geno.bed==0,na.rm = T)
  geno.1 = colSums(geno.bed==1,na.rm = T)
  geno.2 = colSums(geno.bed==2,na.rm = T)
  tmp.geno = rbind(geno.0,geno.1,geno.2)
  
  max.ind = apply(tmp.geno, 2, function(x) which.max(x))
  max.val = gsub("geno.","",rownames(tmp.geno)[max.ind])
  max.val = matrix(data=as.numeric(max.val),nrow=dim(geno.bed)[1],ncol=dim(geno.bed)[2],byrow = T)
  
  geno.bed[is.na(geno.bed)] = max.val[is.na(geno.bed)]
  
  
  
  
  ## GET GENE COMPOSITE SCORE ####
  cat("\033[33m", "===== Get Gene Composite Score", "\033[0m", "\n")
  
  #  Do PCA per gene
  gene.names <- unique(seq.indexes$gene)
  n_gen <- length(gene.names)
  
  gene.score = matrix(NA,dim(Yshapes)[1],1e4)
  gene.id = matrix(NA,1,1e4)
  gene.n = matrix(NA,1,1e4)
  count = 0
  for (i in 1:n_gen) {
    ind <- which(seq.indexes$gene == gene.names[i])
    if (length(ind) > 0) {
      snp_indices <- match(seq.indexes$snp[ind], colnames(geno.bed))
      snp_indices <- snp_indices[!is.na(snp_indices)]
      
      if (length(snp_indices) > 0) {
        # PCA
        tmp.pca <- fast.prcomp(geno.bed[,snp_indices])
        
        pc_axis = 1
        n = 1
        if (npc == "all"){
          # Parallel analysis (only needed if more than 1 SNP is defined within gene)
          if (length(snp_indices) > 1){
            tmp.pca.pa = paran(geno.bed[,snp_indices], quietly = T, status = F, all = T, iterations = 50, centile = 5)
            pc_axis = which(tmp.pca.pa$Ev > tmp.pca.pa$RndEv)
            n = length(pc_axis)
            if (n == 0){ pc_axis = 1; n = 1 }
          }
        }
        
        gene.score[,(count+1):(count+n)] = as.matrix(tmp.pca$x[,pc_axis])
        tmp.gene.id = as.data.frame(cbind(rep(gene.names[i],1,n),as.character(1:n)))
        gene.id[(count+1):(count+n)] = apply(tmp.gene.id[ , 1:2 ], 1, paste, collapse = "_" )
        gene.n[(count+1):(count+n)] = matrix(n,1,n)
        count = count + n
      } else {
        cat("Skipping gene", gene.names[i], "due to no matching SNPs in geno.bed\n")
      }
    } else {
      cat("Skipping gene", gene.names[i], "due to no matching entries in seq.indexes\n")
    }
  }
  gene.score = gene.score[,1:count]
  gene.id = gene.id[1:count]
  gene.n = gene.n[1:count]

  cat("\033[32m", "DEBUG: Final gene.id list:", paste(gene.id, collapse=", "), "\033[0m", "\n")
  cat("\033[32m", "DEBUG: Returning gene list:", paste(gene.names, collapse=", "), "\033[0m", "\n")
  
  return(list(as.matrix(gene.score), as.list(gene.names), Yshapes))
  
}



#core mgp function####
mgp <- function(GO.term = "chondrocyte differentiation", Yshapes, cv = F, lambda = .06, pls_axis = 1, pheno="Y", remove_pc1=F, use_standardized_PCA=F, permutation = 0, doPermutation = F, rangenes = 0, doRangenes = F){
  
  # GO.term = "chondrocyte differentiation"
  # cv = F
  # lambda = 0
  # pls_axis = 1
  # pheno = "A_lm_gen_sex_size"
  # remove_pc1 = F
  # use_standardized_PCA = F
  # permutation = 0
  # doPermutation = F
  # rangenes = 3
  # doRangenes = F
  # pheno_index = "1:65"
  # coordinate.table <- matrix(1:(ncol(get(pheno))), ncol = 3, byrow = T)
  # selected.pheno <- eval(parse(text = pheno_index))
  # pheno.xyz <- as.numeric(t(coordinate.table[selected.pheno,]))
  # Yshapes = subset(get(pheno), select = pheno.xyz)
  
  # Remove PC1 if True
  if(remove_pc1){
    # Step 1: do PCA
    pca_result <- prcomp(Yshapes, center = TRUE, scale = FALSE)
    
    # Step 2: Remove PC1 from the scores
    scores <- pca_result$x
    scores[, "PC1"] <- 0  # Set PC1 scores to zero
    
    reconstructed <- scores %*% t(pca_result$rotation)
    
    # Add the mean back to the reconstructed data
    reconstructed <- reconstructed + matrix(pca_result$center, nrow = nrow(reconstructed), ncol = ncol(reconstructed), byrow = TRUE)
    
    # Step 4: Replace the original data with the reconstructed data
    Yshapes <- reconstructed
  } 
  

  if(pheno=="YHuman" || pheno=="YHumanTANZ" || pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
    res = getHumanGenomeX(GO.term, "GO", pheno, Yshapes)
    probs.rows = res[[1]]
    gene.names = res[[2]]
    Yshapes = res[[3]]
    
    Y_before = Yshapes
    stddPCA = standardizePCA(Yshapes, use_standardized_PCA, 97)
    Yshapes = stddPCA$Yshapes
    Y_pca_result = stddPCA$Y_pca_result
  }else{
    Y_before = Yshapes
    stddPCA = standardizePCA(Yshapes, use_standardized_PCA, 97)
    Yshapes = stddPCA$Yshapes
    Y_pca_result = stddPCA$Y_pca_result
    
    selection.vector <- GO.term
    # selection.vector <- process.list()[[1]][process.list()[[2]] %in% input$variables2]
    #print(length(selection.vector))
    process.ano <- NULL
    # for(i in 1: length(selection.vector)) process.ano <- c(process.ano, as.character(DO.go[DO.go[,3] == selection.vector[i], 2]))
    process.ano <- as.vector(as.character(DO.go[DO.go[,3] == selection.vector, 2]))
    # CHANGE: No multiple go terms allowed with comma since some of them have commas. 
    coi <- c("ENSEMBL", "SYMBOL")
    go2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = process.ano, columns = coi, keytype = "GO")[,-2:-3]))
    go2symbol <- filter_bad_genes(go2symbol)
    go2symbol <- filter_chrom_two(go2symbol) # Part of chromosome 2 has been deemed as un-specific due to high linkage
    # unique is not operating on the symbol for 'insulin receptor signaling 
    # pathway' resulting in duplicate gene later removed. 
    # Thus additional duplicate index for SYMBOL added.
    go2symbol <- go2symbol[!duplicated(go2symbol$SYMBOL), ]
    coi2 <- c("TXCHROM", "TXSTART", "TXEND")
    symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = go2symbol[,2], columns = coi2, keytype="GENEID")
    
    transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
    
    #symbol, chr, start, end
    chr_name <- rep(NA,  length(unique(symbol2info$GENEID)))
    gene.start <- rep(NA,  length(unique(symbol2info$GENEID)))
    gene.end <- rep(NA,  length(unique(symbol2info$GENEID)))
    
    for(i in 1:length(unique(symbol2info$GENEID))){
      
      tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
      chr_name[i] <- tmp.transcript$TXCHROM
      gene.start[i] <- tmp.transcript$TXSTART
      gene.end[i] <- tmp.transcript$TXEND
      
    }
    
    seq.info <- data.frame(mgi_symbol = unique(go2symbol$SYMBOL), chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
    seq.info[,2] <- as.character(seq.info[,2])
    seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6  
    
    #biomart method for getting gene metadata
    # seq.info <- getBM(attributes = c("mgi_symbol", "chromosome_name", "start_position", "end_position") , filters = "go" , values = process.ano ,mart = mouse)
    # seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6
    
    # print(seq.info)
    # Remove all rows where the chromosome_name is not X or a number.
    seq.info <- seq.info[grepl("^\\d+|^X$", seq.info$chromosome_name), ]
    
    #get rid of weird chromosome names
    if(length(grep(seq.info$chromosome_name, pattern = "CHR")) > 0) seq.info <- seq.info[-grep(seq.info$chromosome_name, pattern = "CHR"),]
    
    seq.indexes <- matrix(NA, ncol = 3, nrow = dim(seq.info)[1])
    #we have seq.info which gives us a gene name and its location on the chromosome
    
    the_markers = combined_markers[[pheno]]
    
    for(j in 1 : dim(seq.info)[1]){
      #seq.indexes <- rbind(seq.indexes, cbind(seq.info[j,1],MM_snps[which(MM_snps$chr == seq.info[j,2] & MM_snps$pos > mean(as.numeric(seq.info[j,3:4])) - .07 & MM_snps$pos < mean(as.numeric(seq.info[j,3:4])) + .07), c(1,3)]))
      tmp.indexes <-  the_markers[which(the_markers$chr == seq.info[j,2] & the_markers$Mbp_mm10 > mean(as.numeric(seq.info[j,3:4])) - 2 & the_markers$Mbp_mm10 < mean(as.numeric(seq.info[j,3:4])) + 2), c(1,3)]
      #for each gene, select the marker closest to the middle of the gene
      seq.indexes[j,] <- as.matrix(cbind(seq.info[j,1],tmp.indexes[which.min(abs(tmp.indexes[,2] - mean(as.numeric(seq.info[j,3:4])))),]))
    }
    
    probs.rows <- matrix(NA, nrow = nrow(Yshapes), ncol = nrow(seq.indexes) * 8)
    probrowseq <- seq(1, ncol(probs.rows) + 8, by = 8)
    
    gene.names <- seq.info[,1]
    
    if (pheno %in% names(pheno_tables)) {
      pheno_table <- pheno_tables[[pheno]]
      
      # Iterate over seq.indexes
      for (i in 1:nrow(seq.indexes)) {
        # Access the appropriate table based on pheno value and update probs.rows
        probs.rows[, probrowseq[i]:(probrowseq[i+1] - 1)] <- as.matrix(collect(tbl(pheno_table, seq.indexes[i, 2])))
      }
    }
  }
  
  
  #fit pls2B, need duv, gene names, seq.indexes
  # if(cv){
  #   pls.svd.cv <- perf_mddsPLS(Xs = probs.rows, Y = Yshapes, lambda_min = .03, lambda_max = .15, n_lambda = 4, kfolds = 10, R = pls_axis, mode = "reg", NCORES = 11)
  #   pls.svd <- mddsPLS(Xs = probs.rows, Y = Yshapes, R = pls_axis, lambda = pls.svd.cv$Optim$optim_para_one[1])
  # } else {pls.svd <- ddsPLS(X = probs.rows, Y = Yshapes, lambdas = rep(lambda, pls_axis),NCORES=6, verbose=TRUE, doBoot = FALSE)
  # }
  
  print("Computing PLS")
  successful_components = TRUE
  assign("pls_axis", pls_axis, env=globalenv())
  while (successful_components) {
    tryCatch({
      pls_axis = get("pls_axis", env=globalenv())
      print(pls_axis)
      pls.svd <- ddsPLS(X = probs.rows, Y = Yshapes, lambdas = rep(lambda, pls_axis), NCORES=6, verbose=TRUE, doBoot = FALSE)
      successful_components = FALSE
    }, error = function(e) {
      # Check if the error is due to a singular matrix
      if (grepl("Lapack routine dgesv: system is exactly singular", e$message)) {
        # Extract the component number from the error message using regex
        singular_component <- as.numeric(gsub(".*U\\[([0-9]+),\\1\\] = 0.*", "\\1", e$message))
        
        print(paste("Error: Singular matrix encountered at component", singular_component))
        
        # Reduce the number of components based on the error (component that failed - 1)
        pls_axis = singular_component - 1
        if (pls_axis==0){
          pls_axis = 1
        }
        assign("pls_axis", pls_axis, env=globalenv())
        
      }else {
        # If another error occurs, throw it
        stop(e)
      }
    })
  }
  
  # ------------------------------
  # 1.   PRE-ALLOCATION & SETTINGS
  # ------------------------------
  approximate.p <- -1        # Final permutation p-value for the full model, indicating the probability that the observed variance explained is due to chance.
  perm.r2 <- -1              # Vector to store both the observed R² and the null R²s from permutations, representing the proportion of variance explained.
  perm.pdist.mean <- -1      # Distance statistic for the “3 SD shape”, measuring how far an extrapolated shape is from the mean shape.
  # Lists to collect component-wise and landmark-wise permutation results
  perm.r2.components <- list()   # R² per component for each permutation, showing variance explained by each component.
  perm.r2.landmark <- list() # R2 per landmark for each permutation, showing variance explained by each landmark.
  component.p.values <- c()      # Stores formatted p-values for each component, reflecting the probability that the component's effect is due to chance.
  component.effect.sizes <- c()  # Component effect sizes, calculated as the absolute difference in variance explained divided by the standard deviation.
  landmark.p.values <- c()       # P-values for landmarks, indicating the probability that the observed landmark error is due to chance.
  landmark.effect.sizes <- c()   # Landmark effect sizes, representing the standardized difference in landmark errors.

  if(doPermutation){

    # ---- Safety: keep number of permutations within reasonable bounds ----
    if(as.numeric(permutation) > 5000){
      permutation <- 5000  # Limit permutations to 5000 to prevent excessive computation.
    }else if(as.numeric(permutation) < 2){
      permutation <- 2     # Ensure at least 2 permutations for meaningful results.
    }    
    permutation <- as.numeric(permutation)
    
    perm.result <- rep(NA, permutation)     # Pre-allocate space for global R² from each permutation, representing variance explained.

    # ---- Reserve space: one vector per PLS component ---------------------
    for(comp in 1:pls_axis){
      perm.r2.components[[comp]] <- rep(NA, permutation)  # Pre-allocate space for R² of each component across permutations.
    }

    # ---- Reserve space: one vector per landmark --------------------------
    # Determine number of landmarks based on phenotype type.
    if(pheno %in% c("YHumanDense","YHumanTANZDense")){
      for(lmk in 1:nlandmarks){
        perm.r2.landmark[[lmk]] <- rep(NA, permutation)  # Pre-allocate space for landmark errors across permutations.
      }
    }else{
      for(lmk in 1:(ncol(Y_before)/3)){
        perm.r2.landmark[[lmk]] <- rep(NA, permutation)  # Pre-allocate space for landmark errors across permutations.
      }
    }
            
    # ------------------------------
    # 2.   PREPARE ORIGINAL SHAPES
    # ------------------------------
    YshapesReversed <- reversePCA(Yshapes, Y_pca_result, use_standardized_PCA)  # Reverse PCA to get original shapes.
    if(pheno %in% c("YHumanDense","YHumanTANZDense")){
      YshapesReversed <- reverse_pca_dense(YshapesReversed)  # Undo dense encoding for dense phenotypes.
      Y_mat <- row2array3d(YshapesReversed, Nlandmarks = nlandmarks)
    } else {
      Y_mat <- row2array3d(YshapesReversed, Nlandmarks = ncol(Y_before)/3)
    }
    data.mean.shape <- apply(Y_mat, c(1,2), mean)  # Calculate the grand mean configuration of the shapes.

    # ==============================
    # 3.   PERMUTATION LOOP
    # ==============================
    for(i in 1:length(perm.result)){
      print(paste("Permutation #",i,"/",length(perm.result),"\n"))
      
      # 3a) Break the X-Y link, fit ddsPLS, predict Y
      process.svd <- ddsPLS(
          X       = probs.rows,
          Y       = as.matrix(Yshapes[ sample(1:nrow(Yshapes)), ]), # Randomly shuffle Y to break the link with X.
          lambdas = rep(lambda, pls_axis),
          NCORES  = 6,
          verbose = TRUE,
          doBoot  = FALSE
      )
      
      # 3b) GLOBAL R² under permutation
      full.pred <- predict(process.svd, probs.rows)$y
      full.pred <- reversePCA(full.pred, Y_pca_result, use_standardized_PCA)
      Y_col_means <- matrix(colMeans(Y_before), nrow = nrow(Y_before), ncol = ncol(Y_before), byrow = TRUE)
      ess <- sum((full.pred - Y_col_means)^2)            # Explained sum of squares (ESS), variance explained by the model.
      tss <- sum((Y_before - Y_col_means)^2)             # Total sum of squares (TSS), total variance in the data.
      perm.result[i] <- sqrt(ess/tss)                    # √R², the square root of the proportion of variance explained. The square root is used to provide a measure of correlation, which is more interpretable in terms of linear relationships.

      # 3c) COMPONENT-WISE R² under permutation
      for(comp in 1:pls_axis){
        # Reconstruct Y using only component 'comp'
        comp.pred <- process.svd$model$t[,comp,drop=FALSE] %*% t(process.svd$model$C[,comp,drop=FALSE])
        comp.pred <- comp.pred * t(replicate(nrow(comp.pred), process.svd$model$sdY)) + t(replicate(nrow(comp.pred), process.svd$model$muY))
        comp.pred <- reversePCA(comp.pred, Y_pca_result, use_standardized_PCA)
        comp.ess  <- sum((comp.pred - Y_col_means)^2)  # ESS for the component.
        perm.r2.components[[comp]][i] <- sqrt(comp.ess / tss)  # √R² for the component, variance explained by this component. The square root here also provides a correlation-like measure for each component.
      }
      
      # 3d) LANDMARK-WISE √R² under permutation (same methodology as global)
      if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
        perm.pred.dense <- reverse_pca_dense(full.pred)
        perm.true.dense <- reverse_pca_dense(Y_before)
      }else{
        perm.pred.dense <- full.pred
        perm.true.dense <- Y_before
      }

      if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
        n_landmarks <- nlandmarks
      }else{
        n_landmarks <- ncol(Y_before)/3
      }
      
      print(paste("DEBUG: In permutation loop, n_landmarks =", n_landmarks))
      print(paste("DEBUG: length(perm.r2.landmark) =", length(perm.r2.landmark)))
      
      for(lmk in 1:n_landmarks){
        # Extract the 3 columns (x,y,z) for this landmark
        lmk_cols <- ((lmk-1)*3 + 1):(lmk*3)
        
        print(paste("DEBUG: Processing landmark", lmk, "of", n_landmarks))
        
        # Same ESS/TSS calculation as observed model
        lmk_pred <- perm.pred.dense[, lmk_cols, drop=FALSE]
        lmk_true <- perm.true.dense[, lmk_cols, drop=FALSE]
        lmk_means <- matrix(colMeans(lmk_true), nrow = nrow(lmk_true), ncol = 3, byrow = TRUE)
        
        ess_lmk <- sum((lmk_pred - lmk_means)^2)
        tss_lmk <- sum((lmk_true - lmk_means)^2)
        
        perm.r2.landmark[[lmk]][i] <- sqrt(ess_lmk / tss_lmk)  # Store √R² instead of Euclidean error
      }
      
      
      # 3e) Pdist: distance of a 3 SD extrapolated shape to the mean shape
      # Calculate projection of 3 standard deviations for each permutation
      # Get the mean and standard deviation of the first PLS axis
      mean_score <- apply(process.svd$model$t, 2, mean)
      sd_score <- apply(process.svd$model$t, 2, sd)
      new_score <- matrix(0, ncol = length(mean_score), nrow = 1)
      new_score[1, ] <- mean_score + 3 * sd_score  # Extrapolate 3 standard deviations along the PLS axis.
      
      # Get the coefficient vector from the first response block
      c <- process.svd$model$C
      
      # Perform the prediction
      y <- new_score %*% t(c)
      y = y * t(process.svd$model$sdY)
      y = (y + t(process.svd$model$muY))
      
      y = reversePCA(y, Y_pca_result, use_standardized_PCA)
      
      if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
        y = reverse_pca_dense(y)
        pred.shape <- row2array3d(y, Nlandmarks = nlandmarks)
      }else{
        pred.shape <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
      }
      
      pred.shape <- pred.shape[,,1]
      
      perm.pdist.mean[i] <- sqrt(sum((pred.shape - data.mean.shape)^2))  # Distance of the extrapolated shape from the mean shape, providing insight into the variability and extremity of shape changes.
    }
    
    # ==============================
    # 4.   OBSERVED (REAL) MODEL
    # ==============================
    full.pred <- predict(pls.svd, probs.rows)$y
    full.pred = reversePCA(full.pred, Y_pca_result, use_standardized_PCA)
    Y_col_means = matrix(colMeans(Y_before), nrow = nrow(Y_before), ncol = ncol(Y_before), byrow=TRUE)
    ess <- sum((full.pred - Y_col_means)^2)  # ESS for the observed model.
    tss <- sum((Y_before - Y_col_means)^2)   # TSS for the observed data.
    perm.r2 <- c(sqrt(ess/tss), perm.result)  # Prepend observed √R² to permutation results. The square root provides a correlation measure, connecting variance explained to model fit.
    approximate.p <- sum(sqrt(ess/tss) < perm.result)/length(perm.result)  # Calculate p-value for the observed model, probability that observed R² is due to chance.
    approximate.effect.size <- round(abs((sqrt(ess/tss) - mean(perm.result))/sd(perm.result)), 3)

    # ------------------------------
    # 5.   COMPONENT-LEVEL INFERENCE
    # ------------------------------
    for(comp in 1:pls_axis){
      # Calculate actual R² for this component
      comp.pred <- pls.svd$model$t[, comp, drop=FALSE] %*% t(pls.svd$model$C[, comp, drop=FALSE])
      comp.pred <- comp.pred * t(replicate(nrow(comp.pred), pls.svd$model$sdY))
      comp.pred <- comp.pred + t(replicate(nrow(comp.pred), pls.svd$model$muY))
      comp.pred <- reversePCA(comp.pred, Y_pca_result, use_standardized_PCA)
      
      comp.ess <- sum((comp.pred - Y_col_means)^2)  # ESS for the component.
      actual.comp.r2 <- sqrt(comp.ess/tss)  # √R² for the component, variance explained by this component. The square root here provides a correlation-like measure for each component.
      
      # P-value for component using permutation method
      component.p.values[comp] <- sum(actual.comp.r2 < perm.r2.components[[comp]]) / length(perm.r2.components[[comp]])  # Probability that the component's effect is due to chance.

      
      # Effect size for component
      component.effect.sizes[comp] <- round(abs((actual.comp.r2 - mean(perm.r2.components[[comp]]))/sd(perm.r2.components[[comp]])), 3)  # Standardized effect size for the component, showing the magnitude of the component's effect relative to variability.
    }
    
    # ------------------------------
    # 6.   LANDMARK-LEVEL INFERENCE
    # ------------------------------
    # Convert back to coordinate space for landmark-wise analysis
    if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      full.pred.dense <- reverse_pca_dense(full.pred)
      full.pred.coords <- full.pred.dense
      Y_before.coords <- reverse_pca_dense(Y_before)
    }else{
      full.pred.coords <- full.pred
      Y_before.coords <- Y_before
    }
    
    landmark.actual.R2 <- numeric(nlandmarks)
    for(lmk in 1:nlandmarks){
      # Extract the 3 columns (x,y,z) for this landmark
      lmk_cols <- ((lmk-1)*3 + 1):(lmk*3)
      
      # Apply same ESS/TSS logic as global model, but only for this landmark's coordinates
      lmk_pred <- full.pred.coords[, lmk_cols, drop=FALSE]
      lmk_true <- Y_before.coords[, lmk_cols, drop=FALSE]
      lmk_means <- matrix(colMeans(lmk_true), nrow = nrow(lmk_true), ncol = 3, byrow = TRUE)
      
      ess_lmk <- sum((lmk_pred - lmk_means)^2)  # ESS for this landmark
      tss_lmk <- sum((lmk_true - lmk_means)^2)  # TSS for this landmark
      
      landmark.actual.R2[lmk] <- sqrt(ess_lmk / tss_lmk)  # √R² for this landmark
      
      # P-value: upper tail (higher √R² is better, same as global and component tests)
      lmk.p <- sum(landmark.actual.R2[lmk] < perm.r2.landmark[[lmk]]) / 
                length(perm.r2.landmark[[lmk]])
      landmark.p.values[lmk] <- lmk.p
      
      # Effect size (standardized difference)
      landmark.effect.sizes[lmk] <- round(
          abs((landmark.actual.R2[lmk] - mean(perm.r2.landmark[[lmk]])) / 
              sd(perm.r2.landmark[[lmk]])), 3)
    }
    
    # Debug visualization: plot landmarks colored by p-value gradient
    if(DEBUG){ # Set to TRUE to enable debug visualization
      
      # Get mean shape for visualization
      if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
        Y_before.array <- row2array3d(reverse_pca_dense(Y_before), Nlandmarks = nlandmarks)
      }else{
        Y_before.array <- row2array3d(Y_before, Nlandmarks = ncol(Y_before)/3)
      }
      mean.shape <- apply(Y_before.array, c(1,2), mean)
      
      # Create gradient based on actual p-value range (0 to max observed)
      max_p <- max(landmark.p.values, na.rm = TRUE)
      min_p <- 0  # Always start from 0 for p-values
      
      # Color gradient: green (p=0, most significant) to red (p=max, least significant)
      colors <- colorRampPalette(c("darkgreen", "green", "yellow", "orange", "red"))(100)
      
      # Scale p-values to [0,1] based on actual range
      if(max_p > min_p){
        p.values.scaled <- (landmark.p.values - min_p) / (max_p - min_p)
      } else {
        p.values.scaled <- rep(0, length(landmark.p.values))  # All p-values are the same
      }
      
      # Map to color indices (1-100)
      color.indices <- pmax(1, pmin(100, ceiling(p.values.scaled * 99) + 1))
      
      # Open 3D window
      open3d()
      
      # Plot landmarks with gradient colors
      for(lmk in 1:nlandmarks){
        points3d(mean.shape[lmk,1], mean.shape[lmk,2], mean.shape[lmk,3], 
                  col=colors[color.indices[lmk]], size=10)
      }
      
      # Create informative legend with actual p-value ranges
      legend_p_values <- c(min_p, max_p/4, max_p/2, 3*max_p/4, max_p)
      legend_labels <- sprintf("p = %.3f", legend_p_values)
      legend_colors <- colors[c(1, 25, 50, 75, 100)]
      
      legend3d("topright", 
                legend = legend_labels,
                col = legend_colors, 
                pch = 16, 
                cex = 1.2)
      
      title3d(main = sprintf("Landmarks colored by p-value (0 to %.3f)", max_p))
      
      # Print summary statistics
      cat(sprintf("P-value range: %.4f to %.4f\n", min(landmark.p.values), max(landmark.p.values)))
      cat(sprintf("Mean p-value: %.4f\n", mean(landmark.p.values)))
      cat(sprintf("Significant landmarks (p < 0.05): %d/%d\n", 
                  sum(landmark.p.values < 0.05), length(landmark.p.values)))
    }
  }

  
  # Random gene analysis, same number of genes but random, compared to the original 
  # list result in terms of the coherence, which is defined as the salient vectors 
  # pointing in the same direction
  approximate.p.rangenes <- -1
  rangenes.coherence <- -1
  rangenes.coherence.shapes <- -1
  if(doRangenes){
    if(as.numeric(rangenes)>5000){
      rangenes = 5000
    }else if(as.numeric(rangenes)<2){
      rangenes = 2
    }
    rangenes.result <- rep(NA, rangenes)
    
    rangenes.coherence.shapes <- array(NA, dim = c(length(rangenes.result), ncol(Yshapes) / 3, 3)) 
    
    YshapesReversed = reversePCA(Yshapes, Y_pca_result,  use_standardized_PCA)      
      
    # Mean shape of data
    if (pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      YshapesReversed = reverse_pca_dense(YshapesReversed)
      Y_mat = row2array3d(YshapesReversed, Nlandmarks = nlandmarks)
    }else{
      Y_mat = row2array3d(YshapesReversed, Nlandmarks = ncol(Y_before)/3)
    }
    data.mean.shape <- apply(Y_mat, c(1,2), mean)
    
    # #FIND ALL FAULTY GENES
    # curated_genes <- get_curated_gene_list()
    # 
    # coi2 <- c("TXCHROM", "TXSTART", "TXEND")
    # faultyGenes <- list()
    # # Record the start time
    # start_time <- Sys.time()
    # elapsed_times <- c()
    # convert_seconds <- function(seconds) {
    #   hours <- floor(seconds / 3600)
    #   minutes <- floor((seconds %% 3600) / 60)
    #   return(sprintf("%02d:%02d", hours, minutes))
    # }
    # for(jj in 1:nrow(curated_genes)){
    #   tryCatch(
    #     {
    #       rangenes.gene2symbol <- curated_genes[jj, ]
    #       suppressMessages(symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = curated_genes[jj,2]$ENSEMBL, columns = coi2, keytype="GENEID"))
    #       transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
    #       
    #       chr_name <- rep(NA,  length(unique(symbol2info$GENEID)))
    #       gene.start <- rep(NA,  length(unique(symbol2info$GENEID)))
    #       gene.end <- rep(NA,  length(unique(symbol2info$GENEID)))
    #       
    #       for(i in 1:length(unique(symbol2info$GENEID))){
    #         tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
    #         chr_name[i] <- tmp.transcript$TXCHROM
    #         gene.start[i] <- tmp.transcript$TXSTART
    #         gene.end[i] <- tmp.transcript$TXEND
    #       }
    #       
    #       rangenes.seq.info <- data.frame(mgi_symbol = rangenes.gene2symbol[match(unique(symbol2info$GENEID), rangenes.gene2symbol$ENSEMBL),1], chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
    #       rangenes.seq.info[,2] <- as.character(rangenes.seq.info[,2])
    #       rangenes.seq.info[,3:4] <- as.matrix(rangenes.seq.info[,3:4])/1e6
    #       gene.names <- rangenes.seq.info[,1]
    #       
    #       # print(rangenes.seq.info)
    #       # Remove all rows where the chromosome_name is not X or a number.
    #       rangenes.seq.info <- rangenes.seq.info[grepl("^\\d+|^X$", rangenes.seq.info$chromosome_name), ]
    #       
    #       if(length(grep(rangenes.seq.info$chromosome_name, pattern = "CHR")) > 0) rangenes.seq.info <- rangenes.seq.info[-grep(rangenes.seq.info$chromosome_name, pattern = "CHR"),]
    #       
    #       seq.indexes <- matrix(NA, ncol = 3, nrow = dim(rangenes.seq.info)[1])
    #       #we have rangenes.seq.info which gives us a gene name and its location on the chromosome
    #       
    #       the_markers = combined_markers[[pheno]]
    #       
    #       for(j in 1 : dim(rangenes.seq.info)[1]){
    #         tmp.indexes <-  the_markers[which(the_markers$chr == rangenes.seq.info[j,2] & the_markers$Mbp_mm10 > mean(as.numeric(rangenes.seq.info[j,3:4])) - 2 & the_markers$Mbp_mm10 < mean(as.numeric(rangenes.seq.info[j,3:4])) + 2), c(1,3)]
    #         # for each gene, select the marker closest to the middle of the gene
    #         seq.indexes[j,] <- as.matrix(cbind(rangenes.seq.info[j,1],tmp.indexes[which.min(abs(tmp.indexes[,2] - mean(as.numeric(rangenes.seq.info[j,3:4])))),]))
    #       }
    #       
    #       #put together selected genotype data
    #       random.probs.rows <- matrix(NA, nrow = nrow(Yshapes), ncol = nrow(seq.indexes) * 8)
    #       probrowseq <- seq(1, ncol(random.probs.rows) + 8, by = 8)
    #       if (pheno %in% names(pheno_tables)) {
    #         pheno_table <- pheno_tables[[pheno]]
    #         
    #         # Iterate over seq.indexes
    #         for (i in 1:nrow(seq.indexes)) {
    #           # Access the appropriate table based on pheno value and update probs.rows
    #           random.probs.rows[, probrowseq[i]:(probrowseq[i+1] - 1)] <- as.matrix(collect(tbl(pheno_table, seq.indexes[i, 2])))
    #         }
    #       }
    #       
    # 
    #       if (jj %% 100 == 0) {
    #         end_time <- Sys.time()
    #         elapsed_time <- end_time - start_time
    #         elapsed_times <- c(elapsed_times, elapsed_time)
    # 
    #         # Calculate the average time per 100 iterations
    #         avg_elapsed_time <- mean(elapsed_times)
    # 
    #         # Calculate and print estimated time for completion
    #         total_iterations <- nrow(curated_genes)
    #         iterations_left <- total_iterations - jj
    #         estimated_remaining_time <- (iterations_left / 100) * avg_elapsed_time
    # 
    #         cat("Iteration", jj, "completed in", elapsed_time, "seconds\n")
    #         cat("Estimated time remaining:", convert_seconds(estimated_remaining_time), " hrs\n")
    #         start_time <- Sys.time()  # Reset the start time
    #       }
    #     },error = function(cond) {
    #       #message(conditionMessage(cond))
    #       message(paste("Gene at fault: ", curated_genes[jj, 1], ", ENSEMBL: ", curated_genes[jj, 2]))
    #       faultyGenes <<- append(faultyGenes, curated_genes[jj, 1]$SYMBOL)
    #       cat(paste0(curated_genes[jj, 1]$SYMBOL, ","), file = "~/Documents/MGP/jovid/MGP_API/faultyGenes.txt", sep = "", append = TRUE)
    # 
    #       if (jj %% 100 == 0) {
    #         end_time <- Sys.time()
    #         elapsed_time <- end_time - start_time
    #         elapsed_times <- c(elapsed_times, elapsed_time)
    # 
    #         # Calculate the average time per 100 iterations
    #         avg_elapsed_time <- mean(elapsed_times)
    # 
    #         # Calculate and print estimated time for completion
    #         total_iterations <- nrow(curated_genes)
    #         iterations_left <- total_iterations - jj
    #         estimated_remaining_time <- (iterations_left / 100) * avg_elapsed_time
    # 
    #         cat("Iteration", jj, "completed in", elapsed_time, "seconds\n")
    #         cat("Estimated time remaining:", convert_seconds(estimated_remaining_time), " hrs\n")
    #         start_time <- Sys.time()  # Reset the start time
    #       }
    #     }
    #   )
    # }
    # 
    k=1
    while(k <= length(rangenes.result)){
      tryCatch(
        {
          # Start randomizing the gene selection
          # select the ENSEMBL and SYMBOL columns from the org.Mm.eg.db object
          if(pheno %in% pheno_names_human){
            curated_genes <- get_curated_gene_list_human(pheno)
          }else{
            curated_genes <- get_curated_gene_list()
          }

          rangenes.gene2symbol <- curated_genes[sample(nrow(curated_genes), length(selection.vector)), ]
          
          if(pheno %in% pheno_names_human){
            res = getHumanGenomeX(rangenes.gene2symbol$SYMBOL, "custom", pheno, Yshapes, single_match = TRUE)
            random.probs.rows = res[[1]]
            random.gene.names = res[[2]]
          }else{
            coi2 <- c("TXCHROM", "TXSTART", "TXEND")
            suppressMessages(symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = rangenes.gene2symbol[,2]$ENSEMBL, columns = coi2, keytype="GENEID"))
            transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
            
            chr_name <- rep(NA,  length(unique(symbol2info$GENEID)))
            gene.start <- rep(NA,  length(unique(symbol2info$GENEID)))
            gene.end <- rep(NA,  length(unique(symbol2info$GENEID)))
            
            for(i in 1:length(unique(symbol2info$GENEID))){
              tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
              chr_name[i] <- tmp.transcript$TXCHROM
              gene.start[i] <- tmp.transcript$TXSTART
              gene.end[i] <- tmp.transcript$TXEND
            }
            
            rangenes.seq.info <- data.frame(mgi_symbol = rangenes.gene2symbol[match(unique(symbol2info$GENEID), rangenes.gene2symbol$ENSEMBL),1], chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
            rangenes.seq.info[,2] <- as.character(rangenes.seq.info[,2])
            rangenes.seq.info[,3:4] <- as.matrix(rangenes.seq.info[,3:4])/1e6
            gene.names <- rangenes.seq.info[,1]
            
            # print(rangenes.seq.info)
            # Remove all rows where the chromosome_name is not X or a number.
            rangenes.seq.info <- rangenes.seq.info[grepl("^\\d+|^X$", rangenes.seq.info$chromosome_name), ]
            
            if(length(grep(rangenes.seq.info$chromosome_name, pattern = "CHR")) > 0) rangenes.seq.info <- rangenes.seq.info[-grep(rangenes.seq.info$chromosome_name, pattern = "CHR"),]
            
            seq.indexes <- matrix(NA, ncol = 3, nrow = dim(rangenes.seq.info)[1])
            #we have rangenes.seq.info which gives us a gene name and its location on the chromosome
            
            the_markers = combined_markers[[pheno]]
            
            for(j in 1 : dim(rangenes.seq.info)[1]){
              tmp.indexes <-  the_markers[which(the_markers$chr == rangenes.seq.info[j,2] & the_markers$Mbp_mm10 > mean(as.numeric(rangenes.seq.info[j,3:4])) - 2 & the_markers$Mbp_mm10 < mean(as.numeric(rangenes.seq.info[j,3:4])) + 2), c(1,3)]
              # for each gene, select the marker closest to the middle of the gene
              seq.indexes[j,] <- as.matrix(cbind(rangenes.seq.info[j,1],tmp.indexes[which.min(abs(tmp.indexes[,2] - mean(as.numeric(rangenes.seq.info[j,3:4])))),]))
            }
            
            #put together selected genotype data
            random.probs.rows <- matrix(NA, nrow = nrow(Yshapes), ncol = nrow(seq.indexes) * 8)
            probrowseq <- seq(1, ncol(random.probs.rows) + 8, by = 8)
            if (pheno %in% names(pheno_tables)) {
              pheno_table <- pheno_tables[[pheno]]
              
              # Iterate over seq.indexes
              for (i in 1:nrow(seq.indexes)) {
                # Access the appropriate table based on pheno value and update probs.rows
                random.probs.rows[, probrowseq[i]:(probrowseq[i+1] - 1)] <- as.matrix(collect(tbl(pheno_table, seq.indexes[i, 2])))
              }
            }
          }
          process.svd <- ddsPLS(X = random.probs.rows, Y = Yshapes, lambdas = rep(lambda, pls_axis), NCORES=6, verbose=TRUE, doBoot = FALSE)
          
          # Calculate projection of 3 stdev for each randomization
          # Get the mean and standard deviation of the PLS axes
          mean_score <- apply(process.svd$model$t, 2, mean)
          sd_score <- apply(process.svd$model$t, 2, sd)
          new_score <- matrix(0, ncol = length(mean_score), nrow = 1)
          new_score[1, ] <- mean_score + 3 * sd_score
          
          # Get the coefficient vector from the first response block
          c <- process.svd$model$C
          
          # Perform the prediction
          y <- new_score %*% t(c)
          y = y * t(process.svd$model$sdY)
          y = (y + t(process.svd$model$muY))
          
          y = reversePCA(y, Y_pca_result, use_standardized_PCA)
          
          if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
            y = reverse_pca_dense(y)
            pred.shape <- row2array3d(y, Nlandmarks = nlandmarks)
          }else{
            pred.shape <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
          }
          
          pred.shape <- pred.shape[,,1]
          
          rangenes.coherence.shapes[k,,] <- pred.shape
          k = k + 1
        },error = function(cond) {
          if(pheno %in% pheno_names_human){
            message(rangenes.gene2symbol[,1]$SYMBOL)
          }else{
            matched_symbol <- subset(rangenes.gene2symbol, ENSEMBL==unique(symbol2info$GENEID[i]), select = SYMBOL)$SYMBOL
            message(paste("skipping ", k, " - gene at fault: ", matched_symbol))
            message(conditionMessage(cond))
          }
        },
        warning = function(cond) {
          message(conditionMessage(cond))
        },
        finally = {
          message("")
        }
      )
    }
    
    # Calculate projection of 3 stdev for each randomization
    # Get the mean and standard deviation of the PLS axes
    mean_score <- apply(pls.svd$model$t, 2, mean)
    sd_score <- apply(pls.svd$model$t, 2, sd)
    new_score <- matrix(0, ncol = length(mean_score), nrow = 1)
    new_score[1, ] <- mean_score + 3 * sd_score
    
    # Get the coefficient vector from the first response block
    c <- pls.svd$model$C
    
    # Perform the prediction
    y <- new_score %*% t(c)
    y = y * t(pls.svd$model$sdY)
    y = (y + t(pls.svd$model$muY))
    
    y = reversePCA(y, Y_pca_result, use_standardized_PCA)
    
    if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      y = reverse_pca_dense(y)
      pred.shape <- row2array3d(y, Nlandmarks = nlandmarks)
    }else{
      pred.shape <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
    }
    
    pred.shape <- pred.shape[,,1]
    
    rangenes.coherence.cos <- numeric(length = length(rangenes.result))
    
    for (l in 1:length(rangenes.result)) {
      rangenes.coherence.cos[l] <- cos_dist(pred.shape-data.mean.shape, rangenes.coherence.shapes[l,,]-data.mean.shape, 0.007)
    }
    
    # Summary statistics
    mean_cos_distance <- mean(rangenes.coherence.cos)
    std_dev_cos_distance <- sd(rangenes.coherence.cos)
    
    # Visualize the distribution
    #hist(rangenes.coherence.cos, main = "Cosine Distances", xlab = "Cosine Distance")
    
    # Perform hypothesis testing (for example, one-sample t-test against a reference value)
    reference_value <- 1.0
    # Generate an artificial uniform distribution for each iteration with the same number of values
    reference_values <- matrix(runif(length(rangenes.coherence.cos), 0.5, 1.0), ncol = 1)
    
    t_test_result <- t.test(rangenes.coherence.cos, mu = reference_value)
    approximate.p.rangenes <- t_test_result$p.value
  }
  
  
  if(pheno %in% pheno_names_human){
    #cache to pls list for new analyses
    results <- list(pls.svd, gene.names, "NA", probs.rows)
    
    tmp.reactive <- results
    reactive.svd <- tmp.reactive[[1]]$model$P
    gene.names <- tmp.reactive[[2]]
    seq.info <- tmp.reactive[[3]]
    
    # Determine the number of components based on the size of the second axis
    num_of_components <- ncol(reactive.svd)
    
    # --- Hanne's weighted sum of loadings ---
    # Sum of absolute values of gene loading, weighted by number of PCs, per gene
    genes = unlist(unique(gene.names))
    gloadings = matrix(0,length(genes), num_of_components)
    comp = 1
    for (comp in 1:num_of_components) {
      for (i in 1:length(genes)){
        ind = (gene.names %in% genes[i])
        gloadings[i, comp] = (sum(abs(reactive.svd[ind,comp]))) / (sum(ind))
      }
    }
    # ----------------------------------------
    
    # Initialize an empty list to store the data frames
    loadings_list <- list()
    
    # Calculate the combined loading for component 0
    combined_loading <- data.frame(
      gloadings = apply(abs(gloadings), 1, sum), # Multiply the loadings across each row
      gnames = genes
    )
    
    # Add the combined loading to the list as component 0
    loadings_list[["component_0"]] <- combined_loading
    
    # Loop through each component and create a data frame for it
    for (i in 1:num_of_components) {
      # Create a data frame for the i-th component
      component_loading <- data.frame(
        gloadings = gloadings[, i],
        gnames = genes
      )
      
      # Add the data frame to the list
      loadings_list[[paste0("component_", i)]] <- component_loading
    }
    
  }else{
    #cache to pls list for new analyses
    results <- list(pls.svd, gene.names, seq.info, probs.rows)
    
    tmp.reactive <- results
    reactive.svd <- tmp.reactive[[1]]$model$P
    gene.names <- tmp.reactive[[2]]
    seq.info <- tmp.reactive[[3]]
    
    #now we should be able to take pls.svd directly and maybe label them by founder in a new column, then barplot by family, by gene
    do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
    do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
    
    pathway.loadings <- data.frame(gloadings = reactive.svd[,1], gnames = as.character(rep(seq.info[,1], each = 8)), founders = rep(do.names, nrow(seq.info)))
    # above it was reactive.svd[,pls_axis]
    
    # Determine the number of components based on the size of the second axis
    num_of_components <- ncol(reactive.svd)
    
    # Initialize an empty list to store the data frames
    loadings_list <- list()
    
    # Calculate the combined loading for component 0
    combined_loading <- data.frame(
      gloadings = apply(abs(reactive.svd), 1, sum), # Multiply the loadings across each row
      gnames = as.character(rep(seq.info[, 1], each = 8)),
      founders = rep(do.names, nrow(seq.info))
    )
    
    bar_order <- pathway.loadings %>% 
      group_by(gnames) %>%
      summarise(test = diff(range(gloadings))) %>%
      arrange(-test) 
    
    combined_loading$gnames <- factor(combined_loading$gnames, levels = lapply(bar_order, as.character)$gnames)
    
    # Add the combined loading to the list as component 0
    loadings_list[["component_0"]] <- combined_loading
    
    # Loop through each component and create a data frame for it
    for (i in 1:num_of_components) {
      # Create a data frame for the i-th component
      component_loading <- data.frame(
        gloadings = reactive.svd[, i],
        gnames = as.character(rep(seq.info[, 1], each = 8)),
        founders = rep(do.names, nrow(seq.info))
      )
      
      bar_order <- pathway.loadings %>% 
        group_by(gnames) %>%
        summarise(test = diff(range(gloadings))) %>%
        arrange(-test) 
      
      component_loading$gnames <- factor(component_loading$gnames, levels = lapply(bar_order, as.character)$gnames)
      
      # Add the data frame to the list
      loadings_list[[paste0("component_", i)]] <- component_loading
    }
  }
  
  tmp.reactive <- results
  probs.rows <- tmp.reactive[[4]]
  
  snp.dim = pls_axis
  
  # Calculate projection

  # Get the mean and standard deviation of all PLS axes
  mean_score <- apply(pls.svd$model$t, 2, mean)
  sd_score <- apply(pls.svd$model$t, 2, sd)
  new_score <- matrix(0, ncol = length(mean_score), nrow = 2)
  new_score[1, ] <- mean_score - 3 * sd_score
  new_score[2, ] <- mean_score + 3 * sd_score
  
  print(paste("all - means:", mean_score, " - sds:", sd_score))
  
  # Get the coefficient vector from the first response block
  c <- pls.svd$model$C
  
  # Perform the prediction
  y <- new_score %*% t(c)
  y = y * t(replicate(2, pls.svd$model$sdY))
  y = (y + t(replicate(2, pls.svd$model$muY)))
  
  y = reversePCA(y, Y_pca_result, use_standardized_PCA) 
  
  if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
    y = reverse_pca_dense(y)
    proj.coords.a1 <- row2array3d(y, Nlandmarks = nlandmarks)
  }else{
    proj.coords.a1 <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
  }
  
  #print(dim(proj.coords.a1))
  proj.coords.a2 <- proj.coords.a1[,,2]
  proj.coords.a1 <- proj.coords.a1[,,1]
  
  if(pheno=="Y"){
    # Swapping rows 17 and 18
    proj.coords.a1[c(17, 18), ] <- proj.coords.a1[c(18, 17), ]
    proj.coords.a2[c(17, 18), ] <- proj.coords.a2[c(18, 17), ]
    
    # Swapping rows 48 and 49
    proj.coords.a1[c(48, 49), ] <- proj.coords.a1[c(49, 48), ]
    proj.coords.a2[c(48, 49), ] <- proj.coords.a2[c(49, 48), ]
  }
  
  
  # ----------Perform shape prediction for each component ----------------
  
  
  # Assuming tmp.reactive[[1]]$model$t is a matrix of PLS scores
  # and tmp.reactive[[1]]$model$C is a matrix of PLS loadings
  
  # Get the number of PLS components
  num_components <- pls_axis
  
  # Initialize an array to store the shapes for each PLS component
  if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
    proj.coords_components_a1 <- array(0, dim = c(nrow(proj.coords.a1), 3, num_components))
    proj.coords_components_a2 <- array(0, dim = c(nrow(proj.coords.a1), 3, num_components))
  }else{
    proj.coords_components_a1 <- array(0, dim = c(ncol(Y_before)/3, 3, num_components))
    proj.coords_components_a2 <- array(0, dim = c(ncol(Y_before)/3, 3, num_components))
  }

  for(i in 1:num_components) {
    # Calculate mean and standard deviation for the current PLS axis
    mean_score <- mean(pls.svd$model$t[, i])
    sd_score <- sd(pls.svd$model$t[, i])
    
    print(paste(i, " mean:", mean_score, " - sd:", sd_score))
    
    # Create a new score matrix for the current PLS axis
    new_score <- matrix(0, ncol = 1, nrow = 2)
    new_score[1, ] <- mean_score - 3 * sd_score
    new_score[2, ] <- mean_score + 3 * sd_score
    
    # Get the coefficient vector for the current PLS axis
    c <- pls.svd$model$C[, i]
    
    # Perform the prediction for the current PLS axis
    y <- new_score %*% t(c)
    y = y * t(replicate(2, pls.svd$model$sdY))
    y = (y + t(replicate(2, pls.svd$model$muY)))
    
    y = reversePCA(y, Y_pca_result, use_standardized_PCA)
    
    if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      y = reverse_pca_dense(y)
      proj.coords <- row2array3d(y, Nlandmarks = nlandmarks)
    }else{
      # Convert the predicted values into the appropriate shape format
      proj.coords <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
    }
    
    # Store the shapes in the list
    proj.coords_components_a1[,,i] <- proj.coords[,,1]
    proj.coords_components_a2[,,i] <- proj.coords[,,2]
    
    if(pheno=="Y"){
      # Swapping rows 17 and 18
      proj.coords_components_a1[c(17, 18),,i] <- proj.coords_components_a1[c(18, 17),,i]
      proj.coords_components_a2[c(17, 18),,i] <- proj.coords_components_a2[c(18, 17),,i]
      
      # Swapping rows 48 and 49
      proj.coords_components_a1[c(48, 49),,i] <- proj.coords_components_a1[c(49, 48),,i]
      proj.coords_components_a2[c(48, 49),,i] <- proj.coords_components_a2[c(49, 48),,i]
    }
  }
  
  # Now shapes_list contains the shapes for each PLS component
  
  # Prepend the previous calculation of all the combined components
  proj.coords_components_a1 <- abind::abind(array(proj.coords.a1, dim = c(dim(proj.coords.a1), 1)), proj.coords_components_a1, along=3)
  proj.coords_components_a2 <- abind::abind(array(proj.coords.a2, dim = c(dim(proj.coords.a2), 1)), proj.coords_components_a2, along=3)

  
  # ------------------------------------------------------------------------
  
  
  # Perform the prediction mean
  mean_score_t <- matrix(0, ncol = length(mean_score), nrow = 1)
  mean_score_t[1, ] <- mean_score
  ym <- as.matrix(mean_score_t) %*% t(c)
  ym = ym * t(pls.svd$model$sdY)
  ym = ym + t(pls.svd$model$muY)
  
  ym = reversePCA(ym, Y_pca_result, use_standardized_PCA)
  
  if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
    ym = reverse_pca_dense(ym)
    proj.coords.mean <- row2array3d(ym, Nlandmarks = nlandmarks)[,,1]
  }else{
    # Convert the predicted values into the appropriate shape format
    proj.coords.mean <- row2array3d(ym, Nlandmarks = ncol(Y_before)/3)[,,1]
  }
  
  if(pheno=="Y"){
    # Swapping rows 17 and 18
    proj.coords.mean[c(17, 18), ] <- proj.coords.mean[c(18, 17), ]
    
    # Swapping rows 48 and 49
    proj.coords.mean[c(48, 49), ] <- proj.coords.mean[c(49, 48), ]
  }

  
  full.pred <- predict(pls.svd, probs.rows)$y # predicted shapes from the genes
  full.pred = reversePCA(full.pred, Y_pca_result, use_standardized_PCA)
  Y_col_means = matrix(colMeans(Y_before), nrow = nrow(Y_before), ncol = ncol(Y_before), byrow=TRUE)
  ess <- sum((full.pred - Y_col_means)^2)
  tss <- sum((Y_before - Y_col_means)^2)
  variance_explained = sqrt(ess/tss)
  
  # Calculate Bayesian Information Criterion
  
  # Set-up
  X_scores <- tmp.reactive[[1]]$model$t
  Y_scores <- tmp.reactive[[1]]$model$t
  X_loadings <- tmp.reactive[[1]]$model$P
  Y_loadings <- tmp.reactive[[1]]$model$C
  
  X_fitted <- X_scores %*% t(X_loadings)
  Y_fitted <- Y_scores %*% t(Y_loadings)
  X <- probs.rows
  
  # De-standardize
  if(length(tmp.reactive[[1]]$model$sdX)==1){
    X_fitted = X_fitted * (replicate(nrow(X_fitted), tmp.reactive[[1]]$model$sdX))
    X_fitted = X_fitted + (replicate(nrow(X_fitted), tmp.reactive[[1]]$model$muX))
    
  }else{
    X_fitted = X_fitted * t(replicate(nrow(X_fitted), tmp.reactive[[1]]$model$sdX))
    X_fitted = X_fitted + t(replicate(nrow(X_fitted), tmp.reactive[[1]]$model$muX))
  }
  
  Y_fitted = Y_fitted * t(replicate(nrow(Y_fitted), tmp.reactive[[1]]$model$sdY))
  Y_fitted = Y_fitted + t(replicate(nrow(Y_fitted), tmp.reactive[[1]]$model$muY))
  
  Y_fitted = reversePCA(Y_fitted, Y_pca_result, use_standardized_PCA)
  
  # Calculate BIC
  X_residuals <- X - X_fitted
  Y_residuals <- Y_before - Y_fitted
  n <- nrow(X) # number of observations
  w <- rep(1, n) # weights of the observations (equal to 1 for unweighted data)
  log_likelihood <- 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * Y_residuals^2))))
  num_parameters <- sum(X_loadings != 0) + sum(Y_loadings != 0)
  BIC <- abs(-2 * log_likelihood + num_parameters * log(n))
  
  corrective_scale = corrective_scales[[pheno]]
  
  regression_error = (sqrt(mean((as.vector(X_residuals))^2)) + sqrt(mean((as.vector(corrective_scale * Y_residuals))^2)))/2
  pls_mag =  sqrt(sum((new_score[1] - new_score[2])^2))
  procrustes_dist = sqrt(sum(((proj.coords.a2) - (proj.coords.mean))^2))
  displacement = as.vector(sqrt(rowSums(((corrective_scale * proj.coords.a2) - (corrective_scale * proj.coords.mean))^2)))
  # Cumulative percent
  percent_cutoff_disp <- cumulative_percent(displacement, 0.8)
  displacement = displacement[displacement>percent_cutoff_disp]
  displacement_mean = mean(displacement)
  displacement_max = max(displacement)
  
  if(doPermutation){
    # Global permutation results
    p_vs_scrambledsubjects = approximate.p #pnorm(q = perm.r2[1], mean=mean(perm.result), sd=sd(perm.result), lower.tail = FALSE)
    p_vs_scrambledsubjects = ifelse(p_vs_scrambledsubjects < 0.0001, sprintf("%.2e", p_vs_scrambledsubjects), as.character(round(p_vs_scrambledsubjects, digits = 4)))
    esize_vs_scrambledsubjects = approximate.effect.size
    
    # Component-wise permutation results
    p_vs_scrambledsubjects_components = sapply(component.p.values, function(p) {
      ifelse(p < 0.0001, sprintf("%.2e", p), as.character(round(p, digits = 4)))
    })
    esize_vs_scrambledsubjects_components = component.effect.sizes
    
    # Landmark-wise permutation results  
    p_vs_scrambledsubjects_landmarks = landmark.p.values
    esize_vs_scrambledsubjects_landmarks = landmark.effect.sizes
    
  }else{
    p_vs_scrambledsubjects = -1
    esize_vs_scrambledsubjects = -1
    p_vs_scrambledsubjects_components = -1
    esize_vs_scrambledsubjects_components = -1
    p_vs_scrambledsubjects_landmarks = -1
    esize_vs_scrambledsubjects_landmarks = -1
  }
  
  # rangenes.coherence
  if(doRangenes){
    p_vs_rangenes <- ifelse(approximate.p.rangenes < 0.0001, sprintf("%.2e", approximate.p.rangenes), as.character(round(approximate.p.rangenes, digits = 4)))
    esize_vs_rangenes = round(mean(rangenes.coherence.cos), 2)
    
    # Define COS histogram by bins with intervals of 0.1
    bins <- seq(0, 1, by = 0.1)
    # Cut the data into bins and calculate the frequency
    bin_freq <- table(cut(rangenes.coherence.cos, breaks = bins, include.lowest = TRUE, right = FALSE))
    # Convert the table to a dataframe
    histogram_df <- as.data.frame(bin_freq)
    # Rename the columns
    colnames(histogram_df) <- c("Bins", "Frequency")
    # Print the dataframe
    
  }else{
    p_vs_rangenes = -1
    esize_vs_rangenes = -1
    histogram_df = -1
    
  }
    
  return(list(loadings = loadings_list, pheno1 = proj.coords_components_a1, pheno2 = proj.coords_components_a2, 
            pheno_loadings = t(pls.svd$model$C), p_value = approximate.p, 
            perm_r2 = perm.r2, perm_pdist_mean = perm.pdist.mean,
            variance_explained = variance_explained, BIC=BIC, genelist=gene.names, regression_error=regression_error,
            pls_mag=pls_mag, procrustes_dist=procrustes_dist, displacement_mean=displacement_mean, displacement_max=displacement_max,
            p_vs_rangenes=p_vs_rangenes, p_vs_scrambledsubjects=p_vs_scrambledsubjects, esize_vs_rangenes=esize_vs_rangenes,
            esize_vs_scrambledsubjects=esize_vs_scrambledsubjects, 
            p_vs_scrambledsubjects_components=p_vs_scrambledsubjects_components,
            esize_vs_scrambledsubjects_components=esize_vs_scrambledsubjects_components,
            p_vs_scrambledsubjects_landmarks=p_vs_scrambledsubjects_landmarks,
            esize_vs_scrambledsubjects_landmarks=esize_vs_scrambledsubjects_landmarks,
            histogram_df=histogram_df))
}

#custom MGP function####
custom.mgp <- function(genelist = c("Bmp7", "Bmp2", "Bmp4", "Ankrd11"), Yshapes, cv = F, lambda = .06, pls_axis = 1, pheno="Y", remove_pc1=F, use_standardized_PCA=F, permutation = 0, doPermutation = F, rangenes = 0, doRangenes = F){
  
    # genelist = c("PAX3","EYA4","PKHD1","RPS12","SRBD1","TWIST1","KAT7")
    # cv = F
    # lambda = 0
    # pls_axis = 3
    # pheno = "YHumanDense"
    # remove_pc1 = F
    # use_standardized_PCA = F
    # permutation = 40
    # doPermutation = T
    # rangenes = 40
    # doRangenes = F
    # pheno_index = "1:5629"
    # # coordinate.table <- matrix(1:(ncol(get(pheno))), ncol = 3, byrow = T)
    # # selected.pheno <- eval(parse(text = pheno_index))
    # # pheno.xyz <- as.numeric(t(coordinate.table[selected.pheno,]))
    # # Yshapes = subset(get(pheno), select = pheno.xyz)
    # Yshapes = get(pheno)

    # Remove PC1 if True
    if(remove_pc1){
      # Step 1: do PCA
      pca_result <- prcomp(Yshapes)
      
      # Step 2: Calculate scores for PC1
      scores_pc1 <- pca_result$x[, "PC1"]
      
      # Step 3: Reconstruct the data using only PC1
      reconstructed_pc1 <- scores_pc1 %*% t(pca_result$rotation[, "PC1"])
      
      # Step 4: Subtract the PC1 reconstruction from the original data
      Yshapes <- Yshapes - reconstructed_pc1
    }  
  
    Yshapes_original = Yshapes
    if(pheno=="YHuman" || pheno=="YHumanTANZ" || pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      selection.vector <- as.character(genelist)
      print(selection.vector)
      print(dim(Yshapes))
      res = getHumanGenomeX(selection.vector, "custom", pheno, Yshapes)
      probs.rows = res[[1]]
      gene.names = res[[2]]
      Yshapes = res[[3]]

      cat("selection.vector and gene.names: ", selection.vector, " - gene names: ", as.character(gene.names), "\n")
      cat("number of genes: ", length(gene.names), "\n")
      cat("number of genelist: ", length(genelist), "\n")

      cat("Yshapes: ", dim(Yshapes), "\n")
      cat("probs.rows: ", dim(probs.rows), "\n")
            
      Y_before = Yshapes
      stddPCA = standardizePCA(Yshapes, use_standardized_PCA, 97)
      Yshapes = stddPCA$Yshapes
      Y_pca_result = stddPCA$Y_pca_result
    }else{
      Y_before = Yshapes
      stddPCA = standardizePCA(Yshapes, use_standardized_PCA, 97)
      Yshapes = stddPCA$Yshapes
      Y_pca_result = stddPCA$Y_pca_result
      
      selection.vector <- as.character(genelist)
      coi <- c("ENSEMBL", "SYMBOL")
      gene2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = selection.vector, columns = coi, keytype = "SYMBOL")))
      gene2symbol <- filter_bad_genes(gene2symbol)
      gene2symbol <- filter_chrom_two(gene2symbol)
  
      coi2 <- c("TXCHROM", "TXSTART", "TXEND")
      symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = gene2symbol[,2], columns = coi2, keytype="GENEID")
      transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
  
      chr_name <- rep(NA,  length(unique(symbol2info$GENEID)))
      gene.start <- rep(NA,  length(unique(symbol2info$GENEID)))
      gene.end <- rep(NA,  length(unique(symbol2info$GENEID)))
  
      for(i in 1:length(unique(symbol2info$GENEID))){
          tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
          chr_name[i] <- tmp.transcript$TXCHROM
          gene.start[i] <- tmp.transcript$TXSTART
          gene.end[i] <- tmp.transcript$TXEND
      }
      
      seq.info <- data.frame(mgi_symbol = gene2symbol[match(unique(symbol2info$GENEID), gene2symbol$ENSEMBL),1], chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
      seq.info[,2] <- as.character(seq.info[,2])
      seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6
      gene.names <- seq.info[,1]
      
      # print(seq.info)
      # Remove all rows where the chromosome_name is not X or a number.
      seq.info <- seq.info[grepl("^\\d+|^X$", seq.info$chromosome_name), ]
  
      if(length(grep(seq.info$chromosome_name, pattern = "CHR")) > 0) seq.info <- seq.info[-grep(seq.info$chromosome_name, pattern = "CHR"),]
  
      seq.indexes <- matrix(NA, ncol = 3, nrow = dim(seq.info)[1])
      #we have seq.info which gives us a gene name and its location on the chromosome
      
      the_markers = combined_markers[[pheno]]
      
      for(j in 1 : dim(seq.info)[1]){
          tmp.indexes <-  the_markers[which(the_markers$chr == seq.info[j,2] & the_markers$Mbp_mm10 > mean(as.numeric(seq.info[j,3:4])) - 2 & the_markers$Mbp_mm10 < mean(as.numeric(seq.info[j,3:4])) + 2), c(1,3)]
          #for each gene, select the marker closest to the middle of the gene
          seq.indexes[j,] <- as.matrix(cbind(seq.info[j,1],tmp.indexes[which.min(abs(tmp.indexes[,2] - mean(as.numeric(seq.info[j,3:4])))),]))
      }
  
      #put together selected genotype data
      probs.rows <- matrix(NA, nrow = nrow(Yshapes), ncol = nrow(seq.indexes) * 8)
      probrowseq <- seq(1, ncol(probs.rows) + 8, by = 8)
      
      if (pheno %in% names(pheno_tables)) {
        pheno_table <- pheno_tables[[pheno]]
        
        # Iterate over seq.indexes
        for (i in 1:nrow(seq.indexes)) {
          # Access the appropriate table based on pheno value and update probs.rows
          probs.rows[, probrowseq[i]:(probrowseq[i+1] - 1)] <- as.matrix(collect(tbl(pheno_table, seq.indexes[i, 2])))
        }
      }
    }
  

    print("Computing PLS")
    successful_components = TRUE
    assign("pls_axis", pls_axis, env=globalenv())
    while (successful_components) {
      tryCatch({
        pls_axis = get("pls_axis", env=globalenv())
        print(pls_axis)
        pls.svd <- ddsPLS(X = probs.rows, Y = Yshapes, lambdas = rep(lambda, pls_axis), NCORES=6, verbose=TRUE, doBoot = FALSE)
        successful_components = FALSE
      }, error = function(e) {
        # Check if the error is due to a singular matrix
        if (grepl("Lapack routine dgesv: system is exactly singular", e$message)) {
          # Extract the component number from the error message using regex
          singular_component <- as.numeric(gsub(".*U\\[([0-9]+),\\1\\] = 0.*", "\\1", e$message))
          
          print(paste("Error: Singular matrix encountered at component", singular_component))
          
          # Reduce the number of components based on the error (component that failed - 1)
          pls_axis = singular_component - 1
          if (pls_axis==0){
            pls_axis = 1
          }
          assign("pls_axis", pls_axis, env=globalenv())
          
        }else {
          # If another error occurs, throw it
          stop(e)
        }
      })
    }
    
    print("PLS done")
    
    max_components <- min(ncol(probs.rows), nrow(probs.rows), nrow(Yshapes))
    print(paste("Max number of components:", max_components))
    
    
    # ------------------------------
    # 1.   PRE-ALLOCATION & SETTINGS
    # ------------------------------
    approximate.p <- -1        # Final permutation p-value for the full model, indicating the probability that the observed variance explained is due to chance.
    perm.r2 <- -1              # Vector to store both the observed R² and the null R²s from permutations, representing the proportion of variance explained.
    perm.pdist.mean <- -1      # Distance statistic for the “3 SD shape”, measuring how far an extrapolated shape is from the mean shape.
    # Lists to collect component-wise and landmark-wise permutation results
    perm.r2.components <- list()   # R² per component for each permutation, showing variance explained by each component.
    perm.r2.landmark <- list() # R2 per landmark for each permutation, showing variance explained by each landmark.
    component.p.values <- c()      # Stores formatted p-values for each component, reflecting the probability that the component's effect is due to chance.
    component.effect.sizes <- c()  # Component effect sizes, calculated as the absolute difference in variance explained divided by the standard deviation.
    landmark.p.values <- c()       # P-values for landmarks, indicating the probability that the observed landmark error is due to chance.
    landmark.effect.sizes <- c()   # Landmark effect sizes, representing the standardized difference in landmark errors.

    if(doPermutation){

      # ---- Safety: keep number of permutations within reasonable bounds ----
      if(as.numeric(permutation) > 5000){
        permutation <- 5000  # Limit permutations to 5000 to prevent excessive computation.
      }else if(as.numeric(permutation) < 2){
        permutation <- 2     # Ensure at least 2 permutations for meaningful results.
      }    
      permutation <- as.numeric(permutation)
      perm.result <- rep(NA, permutation)     # Pre-allocate space for global R² from each permutation, representing variance explained.

      # ---- Reserve space: one vector per PLS component ---------------------
      for(comp in 1:pls_axis){
        perm.r2.components[[comp]] <- rep(NA, permutation)  # Pre-allocate space for R² of each component across permutations.
      }

      # ---- Reserve space: one vector per landmark --------------------------
      # Determine number of landmarks based on phenotype type.
      if(pheno %in% c("YHumanDense","YHumanTANZDense")){
        for(lmk in 1:nlandmarks){
          perm.r2.landmark[[lmk]] <- rep(NA, permutation)  # Pre-allocate space for landmark errors across permutations.
        }
      }else{
        for(lmk in 1:(ncol(Y_before)/3)){
          perm.r2.landmark[[lmk]] <- rep(NA, permutation)  # Pre-allocate space for landmark errors across permutations.
        }
      }

      # Add this BEFORE the permutation loop starts
      Y_col_means <- matrix(colMeans(Y_before), nrow = nrow(Y_before), ncol = ncol(Y_before), byrow = TRUE)
      tss <- sum((Y_before - Y_col_means)^2)
              
      # ------------------------------
      # 2.   PREPARE ORIGINAL SHAPES
      # ------------------------------
      YshapesReversed <- reversePCA(Yshapes, Y_pca_result, use_standardized_PCA)  # Reverse PCA to get original shapes.
      if(pheno %in% c("YHumanDense","YHumanTANZDense")){
        YshapesReversed <- reverse_pca_dense(YshapesReversed)  # Undo dense encoding for dense phenotypes.
        Y_mat <- row2array3d(YshapesReversed, Nlandmarks = nlandmarks)
      } else {
        Y_mat <- row2array3d(YshapesReversed, Nlandmarks = ncol(Y_before)/3)
      }
      data.mean.shape <- apply(Y_mat, c(1,2), mean)  # Calculate the grand mean configuration of the shapes.

      # ==============================
      # 3.   PERMUTATION LOOP
      # ==============================
      for(i in 1:length(perm.result)){
        print(paste("Permutation #",i,"/",length(perm.result),"\n"))

        if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
          n_landmarks <- nlandmarks
        }else{
          n_landmarks <- ncol(Y_before)/3
        }
        
        # 3a) Break the X-Y link, fit ddsPLS, predict Y
        process.svd <- ddsPLS(
            X       = probs.rows,
            Y       = as.matrix(Yshapes[ sample(1:nrow(Yshapes)), ]), # Randomly shuffle Y to break the link with X.
            lambdas = rep(lambda, pls_axis),
            NCORES  = 6,
            verbose = TRUE,
            doBoot  = FALSE
        )
        
        # 3b) GLOBAL R² under permutation
        full.pred <- predict(process.svd, probs.rows)$y
        full.pred <- reversePCA(full.pred, Y_pca_result, use_standardized_PCA)
        #Y_col_means <- matrix(colMeans(Y_before), nrow = nrow(Y_before), ncol = ncol(Y_before), byrow = TRUE)
        ess <- sum((full.pred - Y_col_means)^2)            # Explained sum of squares (ESS), variance explained by the model.
        #tss <- sum((Y_before - Y_col_means)^2)             # Total sum of squares (TSS), total variance in the data.
        perm.result[i] <- sqrt(ess/tss)                    # √R², the square root of the proportion of variance explained. The square root is used to provide a measure of correlation, which is more interpretable in terms of linear relationships.

        # 3c) COMPONENT-WISE R² under permutation
        for(comp in 1:pls_axis){
          # Reconstruct Y using only component 'comp'
          comp.pred <- process.svd$model$t[,comp,drop=FALSE] %*% t(process.svd$model$C[,comp,drop=FALSE])
          comp.pred <- comp.pred * t(replicate(nrow(comp.pred), process.svd$model$sdY)) + t(replicate(nrow(comp.pred), process.svd$model$muY))
          comp.pred <- reversePCA(comp.pred, Y_pca_result, use_standardized_PCA)
          comp.ess  <- sum((comp.pred - Y_col_means)^2)  # ESS for the component.
          perm.r2.components[[comp]][i] <- sqrt(comp.ess / tss)  # √R² for the component, variance explained by this component. The square root here also provides a correlation-like measure for each component.
        }
        
        # 3d) LANDMARK-WISE √R² under permutation (same methodology as global)
        if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
          perm.pred.dense <- reverse_pca_dense(full.pred)
          perm.true.dense <- reverse_pca_dense(Y_before)
        }else{
          perm.pred.dense <- full.pred
          perm.true.dense <- Y_before
        }
        
        # for(lmk in 1:nlandmarks){
        #   # Extract the 3 columns (x,y,z) for this landmark
        #   lmk_cols <- ((lmk-1)*3 + 1):(lmk*3)
          
        #   # Same ESS/TSS calculation as observed model
        #   lmk_pred <- perm.pred.dense[, lmk_cols, drop=FALSE]
        #   lmk_true <- perm.true.dense[, lmk_cols, drop=FALSE]
        #   lmk_means <- matrix(colMeans(lmk_true), nrow = nrow(lmk_true), ncol = 3, byrow = TRUE)
          
        #   ess_lmk <- sum((lmk_pred - lmk_means)^2)
        #   tss_lmk <- sum((lmk_true - lmk_means)^2)
          
        #   perm.r2.landmark[[lmk]][i] <- sqrt(ess_lmk / tss_lmk)  # Store √R² instead of Euclidean error
        # }
        # Vectorized landmark calculation
        pred_reshaped <- array(perm.pred.dense, dim = c(nrow(perm.pred.dense), n_landmarks, 3))
        true_reshaped <- array(perm.true.dense, dim = c(nrow(perm.true.dense), n_landmarks, 3))

        # Calculate means and expanded means for all landmarks at once
        lmk_means <- apply(true_reshaped, c(2,3), mean)
        lmk_means_expanded <- array(rep(lmk_means, each = nrow(perm.true.dense)), 
                                  dim = c(nrow(perm.true.dense), n_landmarks, 3))

        # Calculate ESS and TSS for all landmarks at once
        ess_per_landmark <- apply((pred_reshaped - lmk_means_expanded)^2, 2, sum)
        tss_per_landmark <- apply((true_reshaped - lmk_means_expanded)^2, 2, sum)

        # Store results
        for(lmk in 1:n_landmarks){
          perm.r2.landmark[[lmk]][i] <- sqrt(ess_per_landmark[lmk] / tss_per_landmark[lmk])
        }

        
        # 3e) Pdist: distance of a 3 SD extrapolated shape to the mean shape
        # Calculate projection of 3 standard deviations for each permutation
        # Get the mean and standard deviation of the first PLS axis
        mean_score <- apply(process.svd$model$t, 2, mean)
        sd_score <- apply(process.svd$model$t, 2, sd)
        new_score <- matrix(0, ncol = length(mean_score), nrow = 1)
        new_score[1, ] <- mean_score + 3 * sd_score  # Extrapolate 3 standard deviations along the PLS axis.
        
        # Get the coefficient vector from the first response block
        c <- process.svd$model$C
        
        # Perform the prediction
        y <- new_score %*% t(c)
        y = y * t(process.svd$model$sdY)
        y = (y + t(process.svd$model$muY))
        
        y = reversePCA(y, Y_pca_result, use_standardized_PCA)
        
        if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
          y = reverse_pca_dense(y)
        }        
        pred.shape <- row2array3d(y, Nlandmarks = n_landmarks)        
        pred.shape <- pred.shape[,,1]
        
        perm.pdist.mean[i] <- sqrt(sum((pred.shape - data.mean.shape)^2))  # Distance of the extrapolated shape from the mean shape, providing insight into the variability and extremity of shape changes.
      }
      
      # ==============================
      # 4.   OBSERVED (REAL) MODEL
      # ==============================
      full.pred <- predict(pls.svd, probs.rows)$y
      full.pred = reversePCA(full.pred, Y_pca_result, use_standardized_PCA)
      #Y_col_means = matrix(colMeans(Y_before), nrow = nrow(Y_before), ncol = ncol(Y_before), byrow=TRUE)
      ess <- sum((full.pred - Y_col_means)^2)  # ESS for the observed model.
      #tss <- sum((Y_before - Y_col_means)^2)   # TSS for the observed data.
      perm.r2 <- c(sqrt(ess/tss), perm.result)  # Prepend observed √R² to permutation results. The square root provides a correlation measure, connecting variance explained to model fit.
      approximate.p <- sum(sqrt(ess/tss) < perm.result)/length(perm.result)  # Calculate p-value for the observed model, probability that observed R² is due to chance.
      approximate.effect.size <- round(abs((sqrt(ess/tss) - mean(perm.result))/sd(perm.result)), 3)

      # ------------------------------
      # 5.   COMPONENT-LEVEL INFERENCE
      # ------------------------------
      for(comp in 1:pls_axis){
        # Calculate actual R² for this component
        comp.pred <- pls.svd$model$t[, comp, drop=FALSE] %*% t(pls.svd$model$C[, comp, drop=FALSE])
        comp.pred <- comp.pred * t(replicate(nrow(comp.pred), pls.svd$model$sdY))
        comp.pred <- comp.pred + t(replicate(nrow(comp.pred), pls.svd$model$muY))
        comp.pred <- reversePCA(comp.pred, Y_pca_result, use_standardized_PCA)
        
        comp.ess <- sum((comp.pred - Y_col_means)^2)  # ESS for the component.
        actual.comp.r2 <- sqrt(comp.ess/tss)  # √R² for the component, variance explained by this component. The square root here provides a correlation-like measure for each component.
        
        # P-value for component using permutation method
        component.p.values[comp] <- sum(actual.comp.r2 < perm.r2.components[[comp]]) / length(perm.r2.components[[comp]])  # Probability that the component's effect is due to chance.

        
        # Effect size for component
        component.effect.sizes[comp] <- round(abs((actual.comp.r2 - mean(perm.r2.components[[comp]]))/sd(perm.r2.components[[comp]])), 3)  # Standardized effect size for the component, showing the magnitude of the component's effect relative to variability.
      }
      
      
      # # ------------------------------
      # # 6.   LANDMARK-LEVEL INFERENCE
      # # ------------------------------
      # # Convert back to coordinate space for landmark-wise analysis
      # if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      #   full.pred.dense <- reverse_pca_dense(full.pred)
      #   full.pred.coords <- full.pred.dense
      #   Y_before.coords <- reverse_pca_dense(Y_before)
      # }else{
      #   full.pred.coords <- full.pred
      #   Y_before.coords <- Y_before
      # }
      
      # landmark.actual.R2 <- numeric(nlandmarks)
      # for(lmk in 1:nlandmarks){
      #   # Extract the 3 columns (x,y,z) for this landmark
      #   lmk_cols <- ((lmk-1)*3 + 1):(lmk*3)
        
      #   # Apply same ESS/TSS logic as global model, but only for this landmark's coordinates
      #   lmk_pred <- full.pred.coords[, lmk_cols, drop=FALSE]
      #   lmk_true <- Y_before.coords[, lmk_cols, drop=FALSE]
      #   lmk_means <- matrix(colMeans(lmk_true), nrow = nrow(lmk_true), ncol = 3, byrow = TRUE)
        
      #   ess_lmk <- sum((lmk_pred - lmk_means)^2)  # ESS for this landmark
      #   tss_lmk <- sum((lmk_true - lmk_means)^2)  # TSS for this landmark
        
      #   landmark.actual.R2[lmk] <- sqrt(ess_lmk / tss_lmk)  # √R² for this landmark
        
      #   # P-value: upper tail (higher √R² is better, same as global and component tests)
      #   lmk.p <- sum(landmark.actual.R2[lmk] < perm.r2.landmark[[lmk]]) / 
      #            length(perm.r2.landmark[[lmk]])
      #   landmark.p.values[lmk] <- lmk.p
        
      #   # Effect size (standardized difference)
      #   landmark.effect.sizes[lmk] <- round(
      #       abs((landmark.actual.R2[lmk] - mean(perm.r2.landmark[[lmk]])) / 
      #           sd(perm.r2.landmark[[lmk]])), 3)
      # }
      
      # -------------------------------------------------------------------
      # 6.   LANDMARK-LEVEL INFERENCE - BENJAMINI-HOCHBERG FDR ADJUSTMENT
      # -------------------------------------------------------------------
      # Convert back to coordinate space for landmark-wise analysis
      if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
        full.pred.dense <- reverse_pca_dense(full.pred)
        full.pred.coords <- full.pred.dense
        Y_before.coords <- reverse_pca_dense(Y_before)
      }else{
        full.pred.coords <- full.pred
        Y_before.coords <- Y_before
      }


      landmark.actual.R2 <- numeric(n_landmarks)
      for(lmk in 1:n_landmarks){
        # Extract the 3 columns (x,y,z) for this landmark
        lmk_cols <- ((lmk-1)*3 + 1):(lmk*3)

        # Apply same ESS/TSS logic as global model, but only for this landmark's coordinates
        lmk_pred <- full.pred.coords[, lmk_cols, drop=FALSE]
        lmk_true <- Y_before.coords[, lmk_cols, drop=FALSE]
        lmk_means <- matrix(colMeans(lmk_true), nrow = nrow(lmk_true), ncol = 3, byrow = TRUE)

        ess_lmk <- sum((lmk_pred - lmk_means)^2)  # ESS for this landmark
        tss_lmk <- sum((lmk_true - lmk_means)^2)  # TSS for this landmark

        landmark.actual.R2[lmk] <- sqrt(ess_lmk / tss_lmk)  # √R² for this landmark

        # P-value: upper tail (higher √R² is better, same as global and component tests)
        lmk.p <- sum(landmark.actual.R2[lmk] < perm.r2.landmark[[lmk]]) / 
                 length(perm.r2.landmark[[lmk]])
        landmark.p.values[lmk] <- lmk.p

        # Effect size (standardized difference)
        landmark.effect.sizes[lmk] <- round(
            abs((landmark.actual.R2[lmk] - mean(perm.r2.landmark[[lmk]])) / 
                sd(perm.r2.landmark[[lmk]])), 3)
      }
      
      # Apply Benjamini-Hochberg FDR adjustment to landmark p-values
      landmark.p.values.adjusted <- numeric(n_landmarks)
      
      # Sort p-values and keep track of original order
      m <- length(landmark.p.values)
      p_sorted_indices <- order(landmark.p.values)
      p_sorted <- landmark.p.values[p_sorted_indices]
      
      # Calculate adjusted p-values using BH procedure
      q_prime <- p_sorted * m / (1:m)  # p(i) * m/i for each sorted p-value
      
      # Apply step-down procedure: work backwards to ensure monotonicity
      adjusted_sorted <- numeric(m)
      adjusted_sorted[m] <- q_prime[m]  # Last (largest) p-value
      
      for(i in (m-1):1){
        adjusted_sorted[i] <- min(adjusted_sorted[i+1], q_prime[i])
      }
      
      # Ensure adjusted p-values don't exceed 1
      adjusted_sorted <- pmin(adjusted_sorted, 1)
      
      # Map back to original order
      landmark.p.values[p_sorted_indices] <- adjusted_sorted

      
      # Debug visualization: plot landmarks colored by p-value gradient
      if(DEBUG){ # Set to TRUE to enable debug visualization
        
        # Get mean shape for visualization
        if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
          Y_before.array <- row2array3d(reverse_pca_dense(Y_before), Nlandmarks = nlandmarks)
        }else{
          Y_before.array <- row2array3d(Y_before, Nlandmarks = ncol(Y_before)/3)
        }
        mean.shape <- apply(Y_before.array, c(1,2), mean)
        
        # Create gradient based on actual p-value range (0 to max observed)
        max_p <- max(landmark.p.values, na.rm = TRUE)
        min_p <- 0  # Always start from 0 for p-values
        
        # Color gradient: green (p=0, most significant) to red (p=max, least significant)
        colors <- colorRampPalette(c("darkgreen", "green", "yellow", "orange", "red"))(100)
        
        # Scale p-values to [0,1] based on actual range
        if(max_p > min_p){
          p.values.scaled <- (landmark.p.values - min_p) / (max_p - min_p)
        } else {
          p.values.scaled <- rep(0, length(landmark.p.values))  # All p-values are the same
        }
        
        # Map to color indices (1-100)
        color.indices <- pmax(1, pmin(100, ceiling(p.values.scaled * 99) + 1))
        
        # Open 3D window
        open3d()
        
        # Plot landmarks with gradient colors
        for(lmk in 1:nlandmarks){
          points3d(mean.shape[lmk,1], mean.shape[lmk,2], mean.shape[lmk,3], 
                   col=colors[color.indices[lmk]], size=10)
        }
        
        # Create informative legend with actual p-value ranges
        legend_p_values <- c(min_p, max_p/4, max_p/2, 3*max_p/4, max_p)
        legend_labels <- sprintf("p = %.3f", legend_p_values)
        legend_colors <- colors[c(1, 25, 50, 75, 100)]
        
        legend3d("topright", 
                 legend = legend_labels,
                 col = legend_colors, 
                 pch = 16, 
                 cex = 1.2)
        
        title3d(main = sprintf("Landmarks colored by p-value (0 to %.3f)", max_p))
        
        # Print summary statistics
        cat(sprintf("P-value range: %.4f to %.4f\n", min(landmark.p.values), max(landmark.p.values)))
        cat(sprintf("Mean p-value: %.4f\n", mean(landmark.p.values)))
        cat(sprintf("Significant landmarks (p < 0.05): %d/%d\n", 
                    sum(landmark.p.values < 0.05), length(landmark.p.values)))
      }
    }
    
    
    # Random gene analysis, same number of genes but random, compared to the original 
    # list result in terms of the coherence, which is defined as the salient vectors 
    # pointing in the same direction
    approximate.p.rangenes <- -1
    rangenes.coherence <- -1
    rangenes.coherence.shapes <- -1
    if(doRangenes){
      if(as.numeric(rangenes)>5000){
        rangenes = 5000
      }else if(as.numeric(rangenes)<2){
        rangenes = 2
      }
      rangenes.result <- rep(NA, rangenes)
      
      if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
        rangenes.coherence.shapes <- array(NA, dim = c(length(rangenes.result), nlandmarks, 3)) 
      }else{
        rangenes.coherence.shapes <- array(NA, dim = c(length(rangenes.result), ncol(Y_before) / 3, 3)) 
      }
      
      # Mean shape of data
      YshapesReversed = reversePCA(Yshapes, Y_pca_result,  use_standardized_PCA)      
      
      # If human dense landmarks, reverse pca dense
      if (pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
        YshapesReversed = reverse_pca_dense(YshapesReversed)
        Y_mat = row2array3d(YshapesReversed, Nlandmarks = nlandmarks)
      }else{
        Y_mat = row2array3d(YshapesReversed, Nlandmarks = ncol(Y_before)/3)
      }

      data.mean.shape <- apply(Y_mat, c(1,2), mean)
      k=1
      while(k <= length(rangenes.result)){
        tryCatch(
          {
            # Start randomizing the gene selection
            # select the ENSEMBL and SYMBOL columns from the org.Mm.eg.db object
            if(pheno %in% pheno_names_human){
              curated_genes <- get_curated_gene_list_human(pheno)
            }else{
              curated_genes <- get_curated_gene_list()
            }
            cat(length(gene.names), "\n")
            cat(length(selection.vector), "\n")
            rangenes.gene2symbol <- curated_genes[sample(nrow(curated_genes), length(gene.names)), ]
            cat(length(rangenes.gene2symbol$SYMBOL), "\n")
          
            # #FIND ALL FAULTY GENES
            # coi2 <- c("TXCHROM", "TXSTART", "TXEND")
            # faultyGenes <- list()
            # # Record the start time
            # start_time <- Sys.time()ss
            # elapsed_times <- c()
            # convert_seconds <- function(seconds) {
            #   hours <- floor(seconds / 3600)
            #   minutes <- floor((seconds %% 3600) / 60)
            #   return(sprintf("%02d:%02d", hours, minutes))
            # }
            # for(jj in 1:nrow(curated_genes)){
            #   tryCatch(
            #     {
            #       suppressMessages(symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = curated_genes[jj,2]$ENSEMBL, columns = coi2, keytype="GENEID"))
            #       
            #       if (jj %% 100 == 0) {
            #         end_time <- Sys.time()
            #         elapsed_time <- end_time - start_time
            #         elapsed_times <- c(elapsed_times, elapsed_time)
            #         
            #         # Calculate the average time per 100 iterations
            #         avg_elapsed_time <- mean(elapsed_times)
            #         
            #         # Calculate and print estimated time for completion
            #         total_iterations <- nrow(curated_genes)
            #         iterations_left <- total_iterations - jj
            #         estimated_remaining_time <- (iterations_left / 100) * avg_elapsed_time
            #         
            #         cat("Iteration", jj, "completed in", elapsed_time, "seconds\n")
            #         cat("Estimated time remaining:", convert_seconds(estimated_remaining_time), " hrs\n")
            #         start_time <- Sys.time()  # Reset the start time
            #       }
            #     },error = function(cond) {
            #       #message(conditionMessage(cond))
            #       message(paste("Gene at fault: ", curated_genes[jj, 1], ", ENSEMBL: ", curated_genes[jj, 2]))
            #       faultyGenes <<- append(faultyGenes, curated_genes[jj, 1]$SYMBOL)
            #       cat(paste0(curated_genes[jj, 1]$SYMBOL, ","), file = "~/Documents/MGP/jovid/MGP_API/faultyGenes.txt", sep = "", append = TRUE)
            #       
            #       if (jj %% 100 == 0) {
            #         end_time <- Sys.time()
            #         elapsed_time <- end_time - start_time
            #         elapsed_times <- c(elapsed_times, elapsed_time)
            #         
            #         # Calculate the average time per 100 iterations
            #         avg_elapsed_time <- mean(elapsed_times)
            #         
            #         # Calculate and print estimated time for completion
            #         total_iterations <- nrow(curated_genes)
            #         iterations_left <- total_iterations - jj
            #         estimated_remaining_time <- (iterations_left / 100) * avg_elapsed_time
            #         
            #         cat("Iteration", jj, "completed in", elapsed_time, "seconds\n")
            #         cat("Estimated time remaining:", convert_seconds(estimated_remaining_time), " hrs\n")
            #         start_time <- Sys.time()  # Reset the start time
            #       }
            #     }
            #   )
            # }
            
            coi2 <- c("TXCHROM", "TXSTART", "TXEND")
            if(pheno %in% pheno_names_human){

              # Get genotypes + already-synchronised phenotypes for the
              # random gene set;  getHumanGenomeX() guarantees that
              # random.probs.rows and the returned Y matrix have identical
              # rows / subjects.
              res = getHumanGenomeX(rangenes.gene2symbol$SYMBOL,"custom", pheno, Yshapes_original, single_match = TRUE)
              random.probs.rows  <- res[[1]]
              Yshapes_for_random <- res[[3]]

              cat("Yshapes_for_random: ", dim(Yshapes_for_random), "\n")
              cat("random.probs.rows: ", dim(random.probs.rows), "\n")
             
              stddPCA_rand = standardizePCA(Yshapes_for_random, use_standardized_PCA, 97)
              Yshapes_rand_pca = stddPCA_rand$Yshapes
              Y_pca_result_rand = stddPCA_rand$Y_pca_result

              process.svd <- ddsPLS(X = random.probs.rows, Y = Yshapes_rand_pca, lambdas = rep(lambda, pls_axis), NCORES=6, verbose=TRUE, doBoot = FALSE)

              mean_score <- apply(process.svd$model$t, 2, mean)
              sd_score <- apply(process.svd$model$t, 2, sd)
              new_score <- matrix(0, ncol = length(mean_score), nrow = 1)
              new_score[1, ] <- mean_score + 3 * sd_score
              
              c <- process.svd$model$C
              y <- new_score %*% t(c)
              y = y * t(process.svd$model$sdY)
              y = (y + t(process.svd$model$muY))
              
              y = reversePCA(y, Y_pca_result_rand, use_standardized_PCA)

              if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
                y = reverse_pca_dense(y)
                cat("y: ", dim(y), "\n")
                cat("nlandmarks: ", nlandmarks, "\n")
                cat("y before: ", dim(Y_before), "\n")
                pred.shape <- row2array3d(y, Nlandmarks = nlandmarks)
              }else{
                cat("y: ", dim(y), "\n")
                cat("ncol(Y_before)/3: ", ncol(Y_before)/3, "\n")
                pred.shape <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
              }
            }else{
              suppressMessages(symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = rangenes.gene2symbol[,2]$ENSEMBL, columns = coi2, keytype="GENEID"))
              transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
              
              chr_name <- rep(NA,  length(unique(symbol2info$GENEID)))
              gene.start <- rep(NA,  length(unique(symbol2info$GENEID)))
              gene.end <- rep(NA,  length(unique(symbol2info$GENEID)))
              
              for(i in 1:length(unique(symbol2info$GENEID))){
                tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
                chr_name[i] <- tmp.transcript$TXCHROM
                gene.start[i] <- tmp.transcript$TXSTART
                gene.end[i] <- tmp.transcript$TXEND
              }
              
              rangenes.seq.info <- data.frame(mgi_symbol = rangenes.gene2symbol[match(unique(symbol2info$GENEID), rangenes.gene2symbol$ENSEMBL),1], chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
              rangenes.seq.info[,2] <- as.character(rangenes.seq.info[,2])
              rangenes.seq.info[,3:4] <- as.matrix(rangenes.seq.info[,3:4])/1e6
              gene.names <- rangenes.seq.info[,1]
              
              # print(rangenes.seq.info)
              # Remove all rows where the chromosome_name is not X or a number.
              rangenes.seq.info <- rangenes.seq.info[grepl("^\\d+|^X$", rangenes.seq.info$chromosome_name), ]
              
              if(length(grep(rangenes.seq.info$chromosome_name, pattern = "CHR")) > 0) rangenes.seq.info <- rangenes.seq.info[-grep(rangenes.seq.info$chromosome_name, pattern = "CHR"),]
              
              seq.indexes <- matrix(NA, ncol = 3, nrow = dim(rangenes.seq.info)[1])
              #we have rangenes.seq.info which gives us a gene name and its location on the chromosome
              
              the_markers = combined_markers[[pheno]]
              
              for(j in 1 : dim(rangenes.seq.info)[1]){
                tmp.indexes <-  the_markers[which(the_markers$chr == rangenes.seq.info[j,2] & the_markers$Mbp_mm10 > mean(as.numeric(rangenes.seq.info[j,3:4])) - 2 & the_markers$Mbp_mm10 < mean(as.numeric(rangenes.seq.info[j,3:4])) + 2), c(1,3)]
                # for each gene, select the marker closest to the middle of the gene
                seq.indexes[j,] <- as.matrix(cbind(rangenes.seq.info[j,1],tmp.indexes[which.min(abs(tmp.indexes[,2] - mean(as.numeric(rangenes.seq.info[j,3:4])))),]))
              }
              
              #put together selected genotype data
              random.probs.rows <- matrix(NA, nrow = nrow(Yshapes), ncol = nrow(seq.indexes) * 8)
              probrowseq <- seq(1, ncol(random.probs.rows) + 8, by = 8)
              if (pheno %in% names(pheno_tables)) {
                pheno_table <- pheno_tables[[pheno]]
                
                # Iterate over seq.indexes
                for (i in 1:nrow(seq.indexes)) {
                  # Access the appropriate table based on pheno value and update probs.rows
                  random.probs.rows[, probrowseq[i]:(probrowseq[i+1] - 1)] <- as.matrix(collect(tbl(pheno_table, seq.indexes[i, 2])))
                }
              }

              process.svd <- ddsPLS(X = random.probs.rows, Y = Yshapes, lambdas = rep(lambda, pls_axis), NCORES=6, verbose=TRUE, doBoot = FALSE)

              mean_score <- apply(process.svd$model$t, 2, mean)
              sd_score <- apply(process.svd$model$t, 2, sd)
              new_score <- matrix(0, ncol = length(mean_score), nrow = 1)
              new_score[1, ] <- mean_score + 3 * sd_score
              
              c <- process.svd$model$C
              y <- new_score %*% t(c)
              y = y * t(process.svd$model$sdY)
              y = (y + t(process.svd$model$muY))
              
              y = reversePCA(y, Y_pca_result, use_standardized_PCA)
              
              pred.shape <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
            }
            
            pred.shape <- pred.shape[,,1]

            cat("pred.shape: ", dim(pred.shape), "\n")
            cat("rangenes.coherence.shapes: ", dim(rangenes.coherence.shapes), "\n")
            
            rangenes.coherence.shapes[k,,] <- pred.shape
            k = k + 1
          },error = function(cond) {
            if(pheno %in% pheno_names_human){
              message(paste(rangenes.gene2symbol$SYMBOL, collapse=", "))
            }else{
              matched_symbol <- subset(rangenes.gene2symbol, ENSEMBL==unique(symbol2info$GENEID[i]), select = SYMBOL)$SYMBOL
              message(paste("skipping ", k, " - gene at fault: ", matched_symbol))
            }
            message(conditionMessage(cond))
          },
          warning = function(cond) {
            message(conditionMessage(cond))
          },
          finally = {
            message("")
          }
        )
      }
      
      # Calculate projection of 3 stdev for each randomization
      # Get the mean and standard deviation of the PLS axes
      mean_score <- apply(pls.svd$model$t, 2, mean)
      sd_score <- apply(pls.svd$model$t, 2, sd)
      new_score <- matrix(0, ncol = length(mean_score), nrow = 1)
      new_score[1, ] <- mean_score + 3 * sd_score
      
      # Get the coefficient vector from the first response block
      c <- pls.svd$model$C
      
      # Perform the prediction
      y <- new_score %*% t(c)
      y = y * t(pls.svd$model$sdY)
      y = (y + t(pls.svd$model$muY))
      
      y = reversePCA(y, Y_pca_result, use_standardized_PCA)
      
      if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
        y = reverse_pca_dense(y)
        pred.shape <- row2array3d(y, Nlandmarks = nlandmarks)
      }else{
        pred.shape <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
      }
      
      pred.shape <- pred.shape[,,1]
      
      rangenes.coherence.cos <- numeric(length = length(rangenes.result))
      
      for (l in 1:length(rangenes.result)) {
        rangenes.coherence.cos[l] <- cos_dist(pred.shape-data.mean.shape, rangenes.coherence.shapes[l,,]-data.mean.shape, 0.007)
      }
      
      # Summary statistics
      mean_cos_distance <- mean(rangenes.coherence.cos)
      std_dev_cos_distance <- sd(rangenes.coherence.cos)
      
      # Visualize the distribution
      #hist(rangenes.coherence.cos, main = "Cosine Distances", xlab = "Cosine Distance")
      
      # Perform hypothesis testing (for example, one-sample t-test against a reference value)
      reference_value <- 1.0
      # Generate an artificial uniform distribution for each iteration with the same number of values
      reference_values <- matrix(runif(length(rangenes.coherence.cos), 0.5, 1.0), ncol = 1)
      
      t_test_result <- t.test(rangenes.coherence.cos, mu = reference_value)
      approximate.p.rangenes <- t_test_result$p.value
    }
    
    
    if(pheno == "YHuman" || pheno == "YHumanTANZ" || pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      #cache to pls list for new analyses
      results <- list(pls.svd, gene.names, "NA", probs.rows)
      
      tmp.reactive <- results
      reactive.svd <- tmp.reactive[[1]]$model$P
      gene.names <- tmp.reactive[[2]]
      seq.info <- tmp.reactive[[3]]
      
      str(reactive.svd)
      str(seq.info)
      
      # Determine the number of components based on the size of the second axis
      num_of_components <- ncol(reactive.svd)
      
      # --- Hanne's weighted sum of loadings ---
      # Sum of absolute values of gene loading, weighted by number of PCs, per gene
      genes = unlist(unique(gene.names))
      gloadings = matrix(0,length(genes), num_of_components)
      comp = 1
      for (comp in 1:num_of_components) {
        for (i in 1:length(genes)){
          ind = (gene.names %in% genes[i])
          gloadings[i, comp] = (sum(abs(reactive.svd[ind,comp]))) / (sum(ind))
        }
      }
      # ----------------------------------------
      
      # Initialize an empty list to store the data frames
      loadings_list <- list()
      
      # Calculate the combined loading for component 0
      combined_loading <- data.frame(
        gloadings = apply(abs(gloadings), 1, sum), # Multiply the loadings across each row
        gnames = genes
      )
      
      # Add the combined loading to the list as component 0
      loadings_list[["component_0"]] <- combined_loading
      
      # Loop through each component and create a data frame for it
      for (i in 1:num_of_components) {
        # Create a data frame for the i-th component
        component_loading <- data.frame(
          gloadings = gloadings[, i],
          gnames = genes
        )
        
        # Add the data frame to the list
        loadings_list[[paste0("component_", i)]] <- component_loading
      }
      
    }else{
    
      #cache to pls list for new analyses
      results <- list(pls.svd, gene.names, seq.info, probs.rows)
      
      tmp.reactive <- results
      reactive.svd <- tmp.reactive[[1]]$model$P
      gene.names <- tmp.reactive[[2]]
      seq.info <- tmp.reactive[[3]]
      
      str(reactive.svd)
      str(seq.info)
      
      #now we should be able to take pls.svd directly and maybe label them by founder in a new column, then barplot by family, by gene
      do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
      do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
      
      pathway.loadings <- data.frame(gloadings = reactive.svd[,pls_axis], gnames = as.character(rep(seq.info[,1], each = 8)), founders = rep(do.names, nrow(seq.info)))
      
      # Determine the number of components based on the size of the second axis
      num_of_components <- ncol(reactive.svd)
      
      # Initialize an empty list to store the data frames
      loadings_list <- list()
      
      # Calculate the combined loading for component 0
      combined_loading <- data.frame(
        gloadings = apply(abs(reactive.svd), 1, sum), # Multiply the loadings across each row
        gnames = as.character(rep(seq.info[, 1], each = 8)),
        founders = rep(do.names, nrow(seq.info))
      )
      
      bar_order <- pathway.loadings %>% 
        group_by(gnames) %>%
        summarise(test = diff(range(gloadings))) %>%
        arrange(-test) 
      
      combined_loading$gnames <- factor(combined_loading$gnames, levels = lapply(bar_order, as.character)$gnames)
      
      # Add the combined loading to the list as component 0
      loadings_list[["component_0"]] <- combined_loading
    
      # Loop through each component and create a data frame for it
      for (i in 1:num_of_components) {
        # Create a data frame for the i-th component
        component_loading <- data.frame(
          gloadings = reactive.svd[, i],
          gnames = as.character(rep(seq.info[, 1], each = 8)),
          founders = rep(do.names, nrow(seq.info))
        )
        
        bar_order <- pathway.loadings %>% 
          group_by(gnames) %>%
          summarise(test = diff(range(gloadings))) %>%
          arrange(-test) 
        
        component_loading$gnames <- factor(component_loading$gnames, levels = lapply(bar_order, as.character)$gnames)
        
        # Add the data frame to the list
        loadings_list[[paste0("component_", i)]] <- component_loading
      }
      
    }
    
    
    tmp.reactive <- results
    probs.rows <- tmp.reactive[[4]]
    
    snp.dim = pls_axis
    
    # Calculate projection
    
    mean_score <- apply(tmp.reactive[[1]]$model$t, 2, mean)
    sd_score <- apply(tmp.reactive[[1]]$model$t, 2, sd)
    new_score <- matrix(0, ncol = length(mean_score), nrow = 2)
    new_score[1, ] <- mean_score - 3 * sd_score
    new_score[2, ] <- mean_score + 3 * sd_score
    
    #print(new_score)
    
    # Get the coefficient vector from the first response block
    c <- tmp.reactive[[1]]$model$C
    
    # Perform the prediction
    y <- new_score %*% t(c)
    y = y * t(replicate(2, tmp.reactive[[1]]$model$sdY))
    y = (y + t(replicate(2, tmp.reactive[[1]]$model$muY)))
    
    y = reversePCA(y, Y_pca_result, use_standardized_PCA)
    
    if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      y = reverse_pca_dense(y)
      proj.coords.a1 <- row2array3d(y, Nlandmarks = nlandmarks)
    }else{
      proj.coords.a1 <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
    }

    
    
    proj.coords.a2 <- proj.coords.a1[,,2]
    proj.coords.a1 <- proj.coords.a1[,,1]
    
    
    # # Test if shapes look good
    # 
    # mean_score <- apply(tmp.reactive[[1]]$model$t, 2, mean)
    # sd_score <- apply(tmp.reactive[[1]]$model$t, 2, sd)
    # new_score <- matrix(0, ncol = length(mean_score), nrow = 2)
    # new_score[1, ] <- mean_score - 3 * 10 * sd_score
    # new_score[2, ] <- mean_score + 3 * 10 * sd_score
    # c <- tmp.reactive[[1]]$model$C
    # y <- new_score %*% t(c)
    # y = y * t(replicate(2, tmp.reactive[[1]]$model$sdY))
    # y = (y + t(replicate(2, tmp.reactive[[1]]$model$muY)))
    # y = reversePCA(y, Y_pca_result, use_standardized_PCA)
    # proj.coords.a1 <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
    # proj.coords.a2 <- proj.coords.a1[,,2]
    # proj.coords.a1 <- proj.coords.a1[,,1]
    # points3d(proj.coords.a1)
    # 
    # atlas_2k <- file2mesh("~/Documents/MGP/jovid/MGP_API/dense_2k_ears_atlas.ply")
    # sparse_reference <- as.matrix(read.csv("~/Documents/MGP/jovid/MGP_API/mgpMeshLMsHumanSparse.csv", row.names=NULL))
    # 
    # # Close existing windows
    # while (rgl.cur() > 0) { close3d() }
    # 
    # # Open window
    # open3d(zoom = 0.7)
    # shade3d(atlas_2k, col = "grey", specular = 1)
    # points3d(sparse_reference, size=10, col = "green", add=TRUE)
    # lmk_scheme = "65_lmk"
    # 
    # # Open window
    # open3d(zoom = 0.7)
    # # Define the size of the window
    # par3d("windowRect"= c(0, 100, 600, 800))
    # # Select subject
    # if (lmk_scheme == "65_lmk") {
    #   subject1 <- Morpho::tps3d(atlas_2k, as.matrix(sparse_reference), proj.coords.a1)
    # } else {
    #   subject1 <- atlas_5k
    #   subject$vb[-4,] <- t(converted[,,to_plot[i]])
    # }
    # # Plot
    # shade3d(subject1, col = "grey", specular = 1)
    # points3d(proj.coords.a1, size=10, col = "red", add=TRUE)
    # 
    # # Open window
    # open3d(zoom = 0.7)
    # # Define the size of the window
    # par3d("windowRect"= c(0, 100, 600, 800))
    # # Select subject
    # if (lmk_scheme == "65_lmk") {
    #   subject2 <- Morpho::tps3d(atlas_2k, as.matrix(sparse_reference), proj.coords.a2)
    # } else {
    #   subject2 <- atlas_5k
    #   subject$vb[-4,] <- t(converted[,,to_plot[i]])
    # }
    # # Plot
    # shade3d(subject2, col = "grey", specular = 1)
    # points3d(proj.coords.a2, size=10, col = "blue", add=TRUE)
    # 
    # 
    # meshDist(subject1, subject2)
    # 
    # # Save the image
    # rgl.snapshot(paste0("Morphs/", lmk_scheme, "_", tolower(to_plot[i]), "_front.png"), top = TRUE)
    # # Close
    # close3d()
    #-------------------------------------------------------
    
    if(pheno=="Y"){
      # Swapping rows 17 and 18
      proj.coords.a1[c(17, 18), ] <- proj.coords.a1[c(18, 17), ]
      proj.coords.a2[c(17, 18), ] <- proj.coords.a2[c(18, 17), ]
      
      # Swapping rows 48 and 49
      proj.coords.a1[c(48, 49), ] <- proj.coords.a1[c(49, 48), ]
      proj.coords.a2[c(48, 49), ] <- proj.coords.a2[c(49, 48), ]
    }
    
    
    # ----------Perform shape prediction for each component ----------------
    
    
    # Assuming tmp.reactive[[1]]$model$t is a matrix of PLS scores
    # and tmp.reactive[[1]]$model$C is a matrix of PLS loadings
    
    # Get the number of PLS components
    num_components <- pls_axis
    
    # Initialize an array to store the shapes for each PLS component
    if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      proj.coords_components_a1 <- array(0, dim = c(nrow(proj.coords.a1), 3, num_components))
      proj.coords_components_a2 <- array(0, dim = c(nrow(proj.coords.a1), 3, num_components))
    }else{
      proj.coords_components_a1 <- array(0, dim = c(ncol(Y_before)/3, 3, num_components))
      proj.coords_components_a2 <- array(0, dim = c(ncol(Y_before)/3, 3, num_components))
    }
    
    
    
    for(i in 1:num_components) {
      # Calculate mean and standard deviation for the current PLS axis
      mean_score <- mean(tmp.reactive[[1]]$model$t[, i])
      sd_score <- sd(tmp.reactive[[1]]$model$t[, i])
      
      # Create a new score matrix for the current PLS axis
      new_score <- matrix(0, ncol = 1, nrow = 2)
      new_score[1, ] <- mean_score - 3 * sd_score
      new_score[2, ] <- mean_score + 3 * sd_score
      
      # Get the coefficient vector for the current PLS axis
      c <- tmp.reactive[[1]]$model$C[, i]
      
      # Perform the prediction for the current PLS axis
      y <- new_score %*% t(c)
      y = y * t(replicate(2, tmp.reactive[[1]]$model$sdY))
      y = (y + t(replicate(2, tmp.reactive[[1]]$model$muY)))
      
      y = reversePCA(y, Y_pca_result, use_standardized_PCA)
      
      if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
        y = reverse_pca_dense(y)
        proj.coords <- row2array3d(y, Nlandmarks = nlandmarks)
      }else{
        # Convert the predicted values into the appropriate shape format
        proj.coords <- row2array3d(y, Nlandmarks = ncol(Y_before)/3)
      }
      
      
      # Store the shapes in the list
      proj.coords_components_a1[,,i] <- proj.coords[,,1]
      proj.coords_components_a2[,,i] <- proj.coords[,,2]
      
      if(pheno=="Y"){
        # Swapping rows 17 and 18
        proj.coords_components_a1[c(17, 18),,i] <- proj.coords_components_a1[c(18, 17),,i]
        proj.coords_components_a2[c(17, 18),,i] <- proj.coords_components_a2[c(18, 17),,i]
        
        # Swapping rows 48 and 49
        proj.coords_components_a1[c(48, 49),,i] <- proj.coords_components_a1[c(49, 48),,i]
        proj.coords_components_a2[c(48, 49),,i] <- proj.coords_components_a2[c(49, 48),,i]
      }

    }
    
    # Now shapes_list contains the shapes for each PLS component

    # Prepend the previous calculation of all the combined components
    proj.coords_components_a1 <- abind::abind(proj.coords.a1, proj.coords_components_a1, along=3)
    proj.coords_components_a2 <- abind::abind(proj.coords.a2, proj.coords_components_a2, along=3)
    

    
    # ------------------------------------------------------------------------

    # Perform the prediction mean
    mean_score_t <- matrix(0, ncol = length(mean_score), nrow = 1)
    mean_score_t[1, ] <- mean_score
    ym <- as.matrix(mean_score_t) %*% t(c)
    ym = ym * t(tmp.reactive[[1]]$model$sdY)
    ym = ym + t(tmp.reactive[[1]]$model$muY)
    
    ym = reversePCA(ym, Y_pca_result, use_standardized_PCA)
    
    if(pheno == "YHumanDense" || pheno == "YHumanTANZDense"){
      ym = reverse_pca_dense(ym)
      proj.coords.mean <- row2array3d(ym, Nlandmarks = nlandmarks)[,,1]
    }else{
      # Convert the predicted values into the appropriate shape format
      proj.coords.mean <- row2array3d(ym, Nlandmarks = ncol(Y_before)/3)[,,1]
    }
    
    
    
    if(pheno=="Y"){
      # Swapping rows 17 and 18
      proj.coords.mean[c(17, 18), ] <- proj.coords.mean[c(18, 17), ]
      
      # Swapping rows 48 and 49
      proj.coords.mean[c(48, 49), ] <- proj.coords.mean[c(49, 48), ]
    }
    
    full.pred <- predict(pls.svd, probs.rows)$y
    full.pred = reversePCA(full.pred, Y_pca_result, use_standardized_PCA)
    Y_col_means = matrix(colMeans(Y_before), nrow = nrow(Y_before), ncol = ncol(Y_before), byrow=TRUE)
    ess <- sum((full.pred - Y_col_means)^2)
    tss <- sum((Y_before - Y_col_means)^2)
    variance_explained = sqrt(ess/tss)

    
    # Calculate Bayesian Information Criterion
    
    # Set-up
    X_scores <- tmp.reactive[[1]]$model$t
    Y_scores <- tmp.reactive[[1]]$model$t
    X_loadings <- tmp.reactive[[1]]$model$P
    Y_loadings <- tmp.reactive[[1]]$model$C
    X_fitted <- X_scores %*% t(X_loadings)
    Y_fitted <- Y_scores %*% t(Y_loadings)
    X <- probs.rows
    
    # De-standardize
    if(length(tmp.reactive[[1]]$model$sdX)==1){
      X_fitted = X_fitted * (replicate(nrow(X_fitted), tmp.reactive[[1]]$model$sdX))
      X_fitted = X_fitted + (replicate(nrow(X_fitted), tmp.reactive[[1]]$model$muX))
      
    }else{
      X_fitted = X_fitted * t(replicate(nrow(X_fitted), tmp.reactive[[1]]$model$sdX))
      X_fitted = X_fitted + t(replicate(nrow(X_fitted), tmp.reactive[[1]]$model$muX))
    }
    
    Y_fitted = Y_fitted * t(replicate(nrow(Y_fitted), tmp.reactive[[1]]$model$sdY))
    Y_fitted = Y_fitted + t(replicate(nrow(Y_fitted), tmp.reactive[[1]]$model$muY))

    Y_fitted = reversePCA(Y_fitted, Y_pca_result, use_standardized_PCA)
    
    # Calculate BIC
    X_residuals <- X - X_fitted
    Y_residuals <- Y_before - Y_fitted
    n <- nrow(X) # number of observations
    w <- rep(1, n) # weights of the observations (equal to 1 for unweighted data)
    log_likelihood <- 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * Y_residuals^2))))
    num_parameters <- sum(X_loadings != 0) + sum(Y_loadings != 0)
    BIC <- abs(-2 * log_likelihood + num_parameters * log(n))
    
    corrective_scale = corrective_scales[[pheno]]
    
    regression_error = (sqrt(mean((as.vector(X_residuals))^2)) + sqrt(mean((as.vector(corrective_scale * Y_residuals))^2)))/2
    pls_mag =  sqrt(sum((new_score[1] - new_score[2])^2))
    procrustes_dist = sqrt(sum(((proj.coords.a2) - (proj.coords.mean))^2))
    displacement = as.vector(sqrt(rowSums(((corrective_scale * proj.coords.a2) - (corrective_scale * proj.coords.mean))^2)))
    # Cumulative percent
    percent_cutoff_disp <- cumulative_percent(displacement, 0.8)
    displacement = displacement[displacement>percent_cutoff_disp]
    displacement_mean = mean(displacement)
    displacement_max = max(displacement)
    
    if(doPermutation){
      # Global permutation results (using pre-calculated empirical p-value)
      p_vs_scrambledsubjects = ifelse(approximate.p < 0.0001, sprintf("%.2e", approximate.p), as.character(round(approximate.p, digits = 4)))
      esize_vs_scrambledsubjects = approximate.effect.size
      
      # Component-wise permutation results
      p_vs_scrambledsubjects_components = sapply(component.p.values, function(p) {
        ifelse(p < 0.0001, sprintf("%.2e", p), as.character(round(p, digits = 4)))
      })
      esize_vs_scrambledsubjects_components = component.effect.sizes
      
      # Landmark-wise permutation results  
      p_vs_scrambledsubjects_landmarks = landmark.p.values
      esize_vs_scrambledsubjects_landmarks = landmark.effect.sizes
      
    }else{
      p_vs_scrambledsubjects = -1
      esize_vs_scrambledsubjects = -1
      p_vs_scrambledsubjects_components = -1
      esize_vs_scrambledsubjects_components = -1
      p_vs_scrambledsubjects_landmarks = -1
      esize_vs_scrambledsubjects_landmarks = -1
    }


    
    # rangenes.coherence
    if(doRangenes){
      p_vs_rangenes <- ifelse(approximate.p.rangenes < 0.0001, sprintf("%.2e", approximate.p.rangenes), as.character(round(approximate.p.rangenes, digits = 4)))
      esize_vs_rangenes = round(mean(rangenes.coherence.cos), 2)
      
      # Define COS histogram by bins with intervals of 0.1
      bins <- seq(0, 1, by = 0.1)
      # Cut the data into bins and calculate the frequency
      bin_freq <- table(cut(rangenes.coherence.cos, breaks = bins, include.lowest = TRUE, right = FALSE))
      # Convert the table to a dataframe
      histogram_df <- as.data.frame(bin_freq)
      # Rename the columns
      colnames(histogram_df) <- c("Bins", "Frequency")
      # Print the dataframe
    }else{
      p_vs_rangenes = -1
      esize_vs_rangenes = -1
      histogram_df = -1
    }
    
    return(list(loadings = loadings_list, pheno1 = proj.coords_components_a1, pheno2 = proj.coords_components_a2, 
            pheno_loadings = t(pls.svd$model$C), p_value = approximate.p, 
            perm_r2 = perm.r2, perm_pdist_mean = perm.pdist.mean,
            variance_explained = variance_explained, BIC=BIC, genelist=gene.names, regression_error=regression_error,
            pls_mag=pls_mag, procrustes_dist=procrustes_dist, displacement_mean=displacement_mean, displacement_max=displacement_max,
            p_vs_rangenes=p_vs_rangenes, p_vs_scrambledsubjects=p_vs_scrambledsubjects, esize_vs_rangenes=esize_vs_rangenes,
            esize_vs_scrambledsubjects=esize_vs_scrambledsubjects, 
            p_vs_scrambledsubjects_components=p_vs_scrambledsubjects_components,
            esize_vs_scrambledsubjects_components=esize_vs_scrambledsubjects_components,
            p_vs_scrambledsubjects_landmarks=p_vs_scrambledsubjects_landmarks,
            esize_vs_scrambledsubjects_landmarks=esize_vs_scrambledsubjects_landmarks,
            histogram_df=histogram_df))
}

#set CORS parameters####
#* @filter cors
cors <- function(res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "Origin, X-Requested-With, Content-Type, Accept, Authorization")
  res$setHeader("Access-Control-Allow-Credentials", "true")
  plumber::forward()
}

#API definition####
#* @apiTitle MGP API

#* Run an MGP model
#* @param GO.term GO term to run
#* @param lambda Regularization strength
#* @param pls_axis how many axes to return
#* @param pheno what phenotype to use
#* @param remove_pc1 should pc1 of pheno be removed
#* @param permutation how many permutations should we use for testing? Max 5000.
#* @get /mgp

function(GO.term = "chondrocyte differentiation", lambda = .06, pls_axis = 1, pheno = "Y", pheno_index = "1:54", permutation = 0, remove_pc1 = F, use_standardized_PCA=F, doPermutation = F, rangenes = 0, doRangenes = F) {
  # future::future({
    print(GO.term)
  
    if(pheno=="YHuman" && pheno_index=="1:5629"){
      pheno = "YHumanDense"
    }
  
    if(pheno=="YHumanTANZ" && pheno_index=="1:5629"){
      pheno = "YHumanTANZDense"
    }
    
    # No string splitting due to GO terms usually containing commas
    # GO.term <- strsplit(GO.term, split = ",")[[1]]
    if(pheno=="YHumanDense" || pheno=="YHumanTANZDense"){
      Y = get(pheno)
    }else{
      coordinate.table <- matrix(1:(ncol(get(pheno))), ncol = 3, byrow = T)
      selected.pheno <- eval(parse(text = pheno_index))
      pheno.xyz <- as.numeric(t(coordinate.table[selected.pheno,]))
      Y = subset(get(pheno), select = pheno.xyz)
    }
    
    mgp(
      GO.term = GO.term, 
      lambda = as.numeric(lambda), 
      Y = Y, 
      pls_axis = as.numeric(pls_axis), 
      permutation = permutation,
      pheno = pheno,
      remove_pc1 = remove_pc1,
      use_standardized_PCA=use_standardized_PCA,
      doPermutation = doPermutation, 
      rangenes = rangenes, 
      doRangenes = doRangenes
    )
    
  # })
}

#* Run an automatic Lambda exploration
#* @param GO.term GO term to run
#* @param pls_axis how many axes to return
#* @param pheno what phenotype to use
#* @post /autoLambda

function(GO.term = "chondrocyte differentiation", genelist = vector(mode = "character", length = 0), pls_axis = 1, pheno = "Y", pheno_index = "1:54", remove_pc1 = F, use_standardized_PCA=F) {
  # future::future({
  
  # GO.term <- strsplit(GO.term, split = ",")[[1]]
  coordinate.table <- matrix(1:(ncol(get(pheno))), ncol = 3, byrow = T)
  selected.pheno <- eval(parse(text = pheno_index))
  pheno.xyz <- as.numeric(t(coordinate.table[selected.pheno,]))
  
  best_lambda = 0.0001
  largest_magnitude = 9999999
  
  # Generating sequence from 0.005 to 0.1 in steps of 0.005
  sequence <- seq(0.005, 0.1, by = 0.005)
  
  # Using a for loop to assign each member of the sequence to lambda
  average_time = 0
  for (lambda in sequence) {
    start.time <- Sys.time()
    tryCatch({
      if (length(genelist) > 0) {
        print(paste("Custom list", genelist))
        mymgp <- custom.mgp(genelist = genelist, 
                            lambda = as.numeric(lambda), 
                            pls_axis = as.numeric(pls_axis),
                            pheno = pheno,
                            remove_pc1 = remove_pc1,
                            use_standardized_PCA=use_standardized_PCA,
                            Y = subset(get(pheno), select = pheno.xyz))
      } else {
        mymgp <- mgp(GO.term = GO.term, 
                     lambda = as.numeric(lambda), 
                     Y = subset(get(pheno), select = pheno.xyz), 
                     pheno = pheno,
                     remove_pc1 = remove_pc1,
                     use_standardized_PCA=use_standardized_PCA,
                     pls_axis = as.numeric(pls_axis))
      }
    
      #max(abs(mymgp$pheno2 - mymgp$pheno1))
      magnitude_var = mymgp$variance_explained 
      magnitude_proc = sum((mymgp$pheno2 - mymgp$pheno1)^2)
      magnitude_bic = mymgp$BIC
      
      print(magnitude_var)
      print(magnitude_proc)
      print(magnitude_bic)
      
      if (magnitude_bic < largest_magnitude) {
        best_lambda = lambda
        largest_magnitude = magnitude_bic
      }
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      average_time = average_time + time.taken
      print(time.taken)
    }, error = function(e) {
      message("Error in iteration with lambda = ", lambda, ": ", e$message)
    })
  }
  
  
  print(best_lambda)
  print(largest_magnitude)
  print((average_time/20))
  print(0.0121705*length(mymgp$genelist)+1.5)
  return(best_lambda)
  # })
}

#* Run a custom MGP model
#* @param genelist comma-separated list of gene names
#* @param lambda Regularization strength
#* @param pls_axis how many axes to return
#* @param pheno what phenotype to use
#* @param permutation should we run a permutation test?
#* @get /custom_mgp


function(genelist = "Bmp7, Bmp2, Bmp4, Ankrd11", lambda = .06, pls_axis = 1,  pheno = "Y", pheno_index = "1:54", remove_pc1 = F, use_standardized_PCA=F, permutation = 0, doPermutation = F, rangenes = 0, doRangenes = F) {
  # future::future({
  
    if (is.character(genelist) && length(genelist) == 1 && grepl(",", genelist)) {
      genelist <- strsplit(genelist, "\\s*,\\s*")[[1]]
    }
  
    print(genelist)
    
    if(pheno=="YHuman" && pheno_index=="1:5629"){
      pheno = "YHumanDense"
    }
    
    if(pheno=="YHumanTANZ" && pheno_index=="1:5629"){
      pheno = "YHumanTANZDense"
    }
  
    if(pheno=="YHumanDense" || pheno=="YHumanTANZDense"){
      Y = get(pheno)
    }else{
      coordinate.table <- matrix(1:(ncol(get(pheno))), ncol = 3, byrow = T)
      selected.pheno <- eval(parse(text = pheno_index))
      pheno.xyz <- as.numeric(t(coordinate.table[selected.pheno,]))
      Y = subset(get(pheno), select = pheno.xyz)
    }
    
    custom.mgp(
      genelist = genelist, 
      lambda = as.numeric(lambda), 
      Y = Y, 
      pls_axis = as.numeric(pls_axis), 
      permutation = permutation,
      doPermutation = doPermutation, 
      pheno = pheno,
      remove_pc1 = remove_pc1,
      use_standardized_PCA=use_standardized_PCA,
      rangenes = rangenes, 
      doRangenes = doRangenes
    )
  # })

}


#* Draw a random correlation matrix for MGP effects in the cache
#* @post /mgp_cor

function(num.processes = 5){
  future::future({
    chosen.processes <- sample(1:nrow(cached.params), as.numeric(num.processes))
    tmp.proc.names <- rep(NA, num.processes)
    tmp.cor <- matrix(NA, nrow = 162, ncol = as.numeric(num.processes))
    for(i in 1:num.processes){
      tmp.cor[,i] <- cached.results[[chosen.processes[i]]][[1]][[1]]$model$C[,1]
      tmp.split <- strsplit(cached.params[chosen.processes[i],1], split = "\\|")[[1]]
      name.buffer <- NULL
      for(j in 1 : length(tmp.split)) {
        name.buffer <- c(name.buffer, as.character(DO.go[DO.go[,2] == tmp.split[j],3]))
      }
      tmp.proc.names[i] <- paste(name.buffer, collapse = " + ") 
    }
    return(list(cormat = cor(tmp.cor), process.names = tmp.proc.names))
  }, seed = NULL)
}

#* Query the Ensembl database
#* @param process
#* @get /mgi

function(process = "chondrocyte differentiation") {
  future::future({
    process.call <- DO.go[DO.go[,3] == process, -1]
    colnames(process.call) <- c("go_id", "go_term", "ngenes")
    
    selection.vector <- c(process)
    process.ano <- NULL
    for(i in 1: length(selection.vector)) process.ano <- c(process.ano, as.character(DO.go[DO.go[,3] == selection.vector[i], 2]))
    
    coi <- c("ENSEMBL", "SYMBOL")
    go2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = process.ano, columns = coi, keytype = "GO")[,-2:-3]))
    coi2 <- c("TXCHROM", "TXSTART", "TXEND")
    symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = go2symbol[,2], columns = coi2, keytype="GENEID")
    
    print(symbol2info)
    
    transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
    
    #symbol, chr, start, end
    chr_name <- rep(NA,  length(unique(symbol2info$GENEID)))
    gene.start <- rep(NA,  length(unique(symbol2info$GENEID)))
    gene.end <- rep(NA,  length(unique(symbol2info$GENEID)))
    
    for(i in 1:length(unique(symbol2info$GENEID))){
      
      tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
      
      chr_name[i] <- tmp.transcript$TXCHROM
      gene.start[i] <- tmp.transcript$TXSTART
      gene.end[i] <- tmp.transcript$TXEND
      
    }
    
    seq.info <- data.frame(mgi_symbol = go2symbol$SYMBOL, chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
    seq.info[,2] <- as.character(seq.info[,2])
    seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6  
    
    #biomart method for getting gene metadata
    # seq.info <- getBM(attributes = c("mgi_symbol", "chromosome_name", "start_position", "end_position") , filters = "go" , values = process.ano ,mart = mouse)
    # seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6
    #get rid of weird chromosome names
    if(length(grep(seq.info$chromosome_name, pattern = "CHR")) > 0) seq.info <- seq.info[-grep(seq.info$chromosome_name, pattern = "CHR"),]
    
    return(list(process_call = process.call, seq_info = seq.info))
    
  })
}


#* Get the entire GO network for a process or gene list
#* @param genelist custom comma-separated list of genes
#* @param process process name
#* @get /GO_network
function(genelist = NULL, process = NULL){
future_promise({
  if(is.null(genelist) == F){
selection.vector <- as.character(strsplit(genelist, ", ")[[1]])
coi <- c("ENSEMBL", "SYMBOL")
gene2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = selection.vector, columns = coi, keytype = "SYMBOL")))

coi2 <- c("GO", "SYMBOL")
ensembl2go <- AnnotationDbi::select(org.Mm.eg.db, keys = gene2symbol[,2], columns = coi2, keytype="ENSEMBL")

GO_ensembl_join <- right_join(DO.go, ensembl2go, by = c("V2" = "GO"))
}

  if(is.null(process) == F){
    
    selection.vector <- c(process)
  
    process.ano <- NULL
    for(i in 1: length(selection.vector)) process.ano <- c(process.ano, as.character(DO.go[DO.go[,3] == selection.vector[i], 2]))
    coi <- c("ENSEMBL", "SYMBOL")
    go2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = process.ano, columns = coi, keytype = "GO")[,-2:-3]))
    coi2 <- c("GO", "SYMBOL")
    ensembl2go <- AnnotationDbi::select(org.Mm.eg.db, keys = go2symbol[,2], columns = coi2, keytype="ENSEMBL")
    GO_ensembl_join <- right_join(DO.go, ensembl2go, by = c("V2" = "GO"))
  }
  
  return(GO_ensembl_join)
  
})
}


#### KRISTEN ####

# Add a conditional here: give a parameter pheno; and then return list of unique mutants

# Alejandro, here are the actual names
# Note that generation is not a variable in the mutants, only DO
pheno_names_mouse_groups_mutants <- list(
  "Manual" = c("Y"),
  "Automated" = c(
    "A_lm_raw",
    "A_lm_age",
    "A_lm_sex",
    "F_lm_age",
    "M_lm_age",
    "A_lm_age_sex",
    "A_lm_age_sex_size",
    "A_lm_age_size",
    "F_lm_age_size",
    "M_lm_age_size",
    "A_lm",
    "F_lm",
    "M_lm")
)

# mappings for equivalent pheno to mutant pheno groups
pheno_names_mouse_groups_mutants_equivalences <- c(
  "A_lm_raw" = "A_lm_raw",
  "A_lm_gen" = "A_lm_age",
  "A_lm_gen_sex" = "A_lm_age_sex",
  "F_lm_gen" = "F_lm_age",
  "M_lm_gen" = "M_lm_age",
  "A_lm_gen_sex_size" = "A_lm_age_sex_size",
  "A_lm_gen_size" = "A_lm_age_size",
  "F_lm_gen_size" = "F_lm_age_size",
  "M_lm_gen_size" = "M_lm_age_size",
  "A_lm" = "A_lm",
  "F_lm" = "F_lm",
  "M_lm" = "M_lm",
  "Y" = "Y"
)

pheno_names_mouse_mutants <- list(
  "Manual" = c("mutant.db", "mutant.lms"),
  "Automated" = c("auto.mutant.db", "auto.mutant.lms")
)


#* Get the list of mutants to make comparisons with
#* @param pheno Phenotype name. Determines which mutant DB to use.
#* @get /mutant_list
function(pheno = NULL){
  print(paste("pheno: ", pheno))
  
  # Check that pheno is provided
  if(is.null(pheno)) {
    stop("Parameter 'pheno' must be specified.")
  }
  
  # Determine group (Manual or Automated) based on phenotype
  group <- NULL
  for(grp in names(pheno_names_mouse_groups)) {
    if(pheno %in% pheno_names_mouse_groups[[grp]]) {
      group <- grp
      break
    }
  }
  if(is.null(group)){
    # Return a descriptive list if pheno is not found but specified
    return(c("No mutants found for phenotype: " = pheno))
  }
  
  # Get the appropriate mutant database variable name
  mutant_db_name <- pheno_names_mouse_mutants[[group]][1]
  
  if(!exists(mutant_db_name, inherits = TRUE)){
    stop(sprintf("Mutant database '%s' not found in environment.", mutant_db_name))
  }
  mutant_db <- get(mutant_db_name, inherits = TRUE)
  
  # Ensure mutant_db has a Genotype column
  if(!"Genotype" %in% colnames(mutant_db)){
    stop("The mutant database does not contain a 'Genotype' column.")
  }
  
  # Return unique genotypes as character vector
  as.character(unique(mutant_db$Gene))
}


#* Calculate correlation between mutant vector and MGP vector
#* @param mutant Selected mutant from /mutant_list
#* @param pheno Selected phenotype
#* @param pheno_index Range of landmarks to use (e.g., "1:54")
#* @param pheno_diff Shape difference between -3 and +3 standard deviations from the population mean on each component
#* @post /mutant_comparison
function(mutant = NULL, pheno="A_lm_gen_sex_size", pheno_index="1:93", pheno_diff=NULL){

  if (is.null(pheno_diff)) {
    pheno_diff <- matrix(c(
      -0.00040, 0.00040, 0.00020, -0.00040, -0.00040, 0.00020, -0.00040, 0.00030, 0.00010, -0.00040, -0.00030, 0.00010, 0.00040, 0.00010, 0.00050, 0.00040, -0.00010, 0.00050, -0.00050, 0.00000, 0.00050, 0.00010, -0.00000, -0.00050, 0.00010, -0.00000, -0.00020, -0.00030, -0.00000, -0.00060, -0.00030, 0.00000, -0.00060, -0.00020, 0.00000, 0.00020, -0.00020, 0.00010, 0.00030, -0.00020, -0.00010, 0.00030, 0.00050, 0.00020, 0.00020, 0.00050, -0.00020, 0.00020, 0.00030, 0.00030, 0.00030, 0.00030, -0.00030, 0.00030, 0.00030, 0.00080, 0.00000, 0.00030, -0.00080, 0.00000, -0.00040, 0.00060, 0.00010, 0.00000, -0.00060, 0.00010, 0.00110, 0.00000, -0.00090, 0.00110, 0.00000, -0.00090, 0.00070, 0.00000, -0.00020, 0.00070, 0.00000, -0.00020, 0.00020, -0.00040, 0.00000, 0.00020, 0.00040, 0.00000, 0.00030, -0.00050, 0.00000, 0.00030, 0.00050, 0.00000, 0.00030, -0.00020, 0.00020, 0.00030, 0.00020, 0.00020, 0.00060, -0.00020, 0.00010, 0.00060, 0.00020, 0.00010, 0.00040, -0.00020, 0.00040, 0.00040, 0.00020, 0.00040, 0.00020, -0.00020, 0.00060, 0.00020, 0.00020, 0.00060, 0.00040, -0.00050, 0.00060, 0.00040, 0.00050, 0.00060, 0.00030, -0.00070, 0.00060, 0.00030, 0.00070, 0.00060, -0.00010, -0.00010, 0.00050, -0.00010, 0.00010, 0.00050, 0.00000, -0.00050, 0.00050, 0.00000, 0.00050, 0.00050, 0.00030, -0.00090, 0.00050, 0.00030, 0.00090, 0.00050, 0.00030, -0.00030, 0.00040, 0.00030, 0.00030, 0.00040, 0.00000, -0.00070, 0.00020, 0.00000, 0.00070, 0.00020, -0.00010, 0.00000, 0.00020, -0.00010, 0.00000, 0.00020, -0.00020, 0.00010, 0.00010, -0.00020, -0.00010, 0.00010, -0.00060, 0.00010, -0.00030, -0.00060, -0.00010, -0.00030, 0.00090, -0.00010, -0.00080, 0.00090, 0.00010, -0.00080, 0.00030, -0.00030, 0.00030, 0.00010, -0.00100, 0.00020, 0.00000, -0.00110, 0.00040, -0.00030, -0.00050, 0.00000, -0.00050, -0.00040, -0.00010, 0.00110, 0.00010, -0.00070, 0.00110, -0.00010, -0.00070, 0.00030, 0.00030, 0.00030, 0.00010, 0.00100, 0.00020, 0.00000, 0.00110, 0.00040, -0.00030, 0.00050, 0.00000, -0.00050, 0.00040, -0.00010, -0.00070, 0.00000, -0.00020, -0.00030, -0.00000, 0.00030, -0.00040, 0.00000, 0.00070, -0.00040, 0.00000, 0.00070, -0.00050, 0.00000, 0.00020, 0.00090, 0.00000, -0.00120, -0.00030, 0.00040, 0.00000, 0.00040, -0.00020, 0.00000, 0.00040, 0.00020, 0.00000, 0.00030, -0.00010, 0.00040, 0.00060, -0.00020, 0.00050, 0.00060, -0.00020, 0.00060, 0.00030, 0.00010, 0.00040, 0.00060, 0.00020, 0.00050, 0.00060, 0.00020, 0.00060, 0.00090, 0.00000, -0.00080, -0.00030, -0.00030, 0.00010, -0.00030, 0.00030, 0.00010, -0.00030, -0.00040, 0.00000, 0.00000, -0.00030, -0.00020, 0.00000, 0.00030, -0.00020,
      -0.00030, 0.00040, 0.00020, -0.00030, -0.00040, 0.00020, -0.00020, 0.00040, 0.00040, -0.00020, -0.00040, 0.00040, 0.00000, 0.00000, 0.00020, 0.00000, 0.00000, 0.00020, -0.00010, 0.00000, 0.00000, -0.00010, -0.00000, -0.00030, 0.00000, -0.00000, -0.00010, -0.00030, -0.00000, -0.00040, -0.00030, -0.00000, -0.00050, -0.00010, 0.00000, -0.00020, -0.00010, 0.00000, -0.00010, -0.00010, 0.00000, -0.00010, 0.00030, 0.00020, 0.00000, 0.00030, -0.00020, 0.00000, 0.00030, 0.00030, 0.00010, 0.00030, -0.00030, 0.00010, 0.00050, 0.00040, 0.00000, 0.00050, -0.00040, 0.00000, -0.00020, 0.00040, -0.00030, -0.00020, -0.00040, -0.00030, 0.00090, 0.00000, -0.00010, 0.00090, 0.00000, -0.00010, 0.00070, -0.00010, 0.00000, 0.00070, 0.00010, 0.00000, 0.00020, -0.00040, 0.00020, 0.00020, 0.00040, 0.00020, 0.00030, -0.00030, 0.00000, 0.00030, 0.00030, 0.00000, 0.00050, -0.00020, 0.00010, 0.00050, 0.00020, 0.00010, 0.00070, -0.00010, -0.00010, 0.00070, 0.00010, -0.00010, 0.00060, -0.00010, 0.00000, 0.00060, 0.00010, 0.00000, 0.00040, -0.00020, 0.00000, 0.00040, 0.00020, 0.00000, 0.00060, -0.00040, 0.00010, 0.00060, 0.00040, 0.00010, 0.00040, -0.00050, 0.00010, 0.00040, 0.00050, 0.00010, 0.00010, -0.00010, 0.00010, 0.00010, 0.00010, 0.00010, 0.00000, -0.00050, 0.00010, 0.00000, 0.00050, 0.00010, 0.00000, -0.00080, 0.00000, 0.00000, 0.00080, 0.00000, 0.00000, -0.00020, 0.00030, 0.00000, 0.00020, 0.00030, -0.00020, -0.00070, 0.00040, -0.00020, 0.00070, 0.00040, -0.00010, 0.00000, 0.00020, -0.00010, 0.00000, 0.00020, -0.00010, -0.00010, 0.00030, -0.00010, 0.00010, 0.00030, -0.00030, 0.00000, 0.00020, -0.00030, 0.00000, 0.00020, 0.00090, -0.00010, -0.00020, 0.00090, 0.00010, -0.00020, 0.00010, -0.00020, 0.00010, -0.00010, -0.00070, -0.00020, -0.00020, -0.00090, -0.00010, -0.00030, -0.00050, 0.00010, -0.00030, -0.00040, 0.00010, 0.00090, 0.00000, -0.00010, 0.00090, 0.00000, -0.00010, 0.00010, 0.00020, 0.00010, -0.00010, 0.00070, -0.00020, -0.00020, 0.00090, -0.00010, -0.00030, 0.00050, 0.00010, -0.00030, 0.00040, 0.00010, -0.00020, 0.00000, 0.00020, -0.00010, -0.00000, 0.00030, -0.00010, 0.00000, 0.00020, -0.00020, 0.00000, 0.00010, -0.00020, 0.00000, -0.00020, 0.00090, 0.00000, -0.00020, -0.00030, 0.00050, 0.00020, 0.00040, -0.00030, 0.00020, 0.00040, 0.00030, 0.00020, 0.00050, -0.00010, 0.00000, 0.00070, -0.00010, 0.00000, 0.00070, -0.00020, 0.00010, 0.00050, 0.00010, 0.00000, 0.00070, 0.00010, 0.00000, 0.00070, 0.00020, 0.00010, 0.00090, 0.00000, -0.00020, -0.00020, -0.00050, 0.00020, -0.00020, 0.00050, 0.00020, -0.00030, -0.00050, 0.00020, -0.00020, -0.00040, 0.00000, -0.00020, 0.00040, 0.00000,
      -0.00010, 0.00000, 0.00000, -0.00010, 0.00000, 0.00000, -0.00020, -0.00010, -0.00020, -0.00020, 0.00010, -0.00020, 0.00000, 0.00000, 0.00020, 0.00000, 0.00000, 0.00020, -0.00030, 0.00000, 0.00030, 0.00010, -0.00000, -0.00010, 0.00020, 0.00000, 0.00000, 0.00010, -0.00000, 0.00000, 0.00030, 0.00000, 0.00000, 0.00010, 0.00000, 0.00020, 0.00010, 0.00000, 0.00010, 0.00010, 0.00000, 0.00010, 0.00020, 0.00000, 0.00010, 0.00020, 0.00000, 0.00010, 0.00010, 0.00000, 0.00010, 0.00010, 0.00000, 0.00010, -0.00010, 0.00040, 0.00000, -0.00010, -0.00040, 0.00000, 0.00020, 0.00000, 0.00010, 0.00020, 0.00000, 0.00010, -0.00010, 0.00000, -0.00010, -0.00010, 0.00000, -0.00010, -0.00030, 0.00000, 0.00020, -0.00030, 0.00000, 0.00020, -0.00020, -0.00010, -0.00010, -0.00020, 0.00010, -0.00010, -0.00010, -0.00020, 0.00010, -0.00010, 0.00020, 0.00010, -0.00030, 0.00000, 0.00030, -0.00030, 0.00000, 0.00030, -0.00010, -0.00010, 0.00030, -0.00010, 0.00010, 0.00030, -0.00020, 0.00000, 0.00040, -0.00020, 0.00000, 0.00040, -0.00020, 0.00000, 0.00040, -0.00020, 0.00000, 0.00040, -0.00020, -0.00010, 0.00060, -0.00020, 0.00010, 0.00060, -0.00020, -0.00010, 0.00040, -0.00020, 0.00010, 0.00040, -0.00010, 0.00000, 0.00010, -0.00010, 0.00000, 0.00010, 0.00000, 0.00000, 0.00030, 0.00000, 0.00000, 0.00030, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, -0.00010, -0.00010, 0.00000, 0.00010, -0.00010, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, -0.00010, 0.00000, 0.00000, -0.00010, 0.00000, 0.00000, -0.00020, 0.00020, -0.00020, -0.00020, -0.00020, -0.00020, -0.00040, 0.00010, -0.00050, -0.00040, -0.00010, -0.00050, -0.00020, 0.00000, -0.00010, -0.00020, 0.00000, -0.00010, 0.00010, 0.00000, 0.00010, -0.00020, -0.00010, 0.00010, -0.00010, -0.00010, 0.00020, -0.00010, 0.00010, -0.00010, -0.00010, 0.00000, -0.00020, -0.00030, 0.00000, -0.00010, -0.00030, 0.00000, -0.00010, 0.00010, 0.00000, 0.00010, -0.00020, 0.00010, 0.00010, -0.00010, 0.00010, 0.00020, -0.00010, -0.00010, -0.00010, -0.00010, 0.00000, -0.00020, -0.00030, 0.00000, -0.00040, -0.00020, 0.00000, 0.00010, 0.00000, 0.00000, 0.00020, 0.00000, 0.00000, 0.00030, 0.00000, -0.00000, 0.00020, -0.00030, 0.00000, -0.00020, -0.00010, -0.00010, -0.00010, -0.00020, 0.00000, 0.00020, -0.00020, 0.00000, 0.00020, -0.00030, -0.00010, 0.00040, -0.00020, -0.00010, 0.00050, -0.00020, 0.00000, 0.00060, -0.00030, 0.00010, 0.00040, -0.00020, 0.00010, 0.00050, -0.00020, 0.00000, 0.00060, -0.00030, -0.00000, 0.00000, -0.00020, 0.00010, -0.00010, -0.00020, -0.00010, -0.00010, -0.00010, 0.00010, -0.00010, 0.00010, 0.00010, -0.00020, 0.00010, -0.00010, -0.00020,
      -0.00010, 0.00000, 0.00000, -0.00010, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00040, 0.00000, 0.00010, 0.00040, 0.00000, 0.00010, -0.00010, -0.00000, 0.00010, 0.00010, -0.00000, 0.00000, 0.00000, -0.00000, -0.00010, -0.00010, -0.00000, -0.00020, -0.00030, -0.00000, 0.00000, -0.00020, -0.00000, 0.00020, -0.00020, 0.00000, 0.00010, -0.00020, 0.00000, 0.00010, 0.00010, 0.00000, 0.00010, 0.00010, 0.00000, 0.00010, -0.00010, 0.00010, 0.00010, -0.00010, -0.00010, 0.00010, 0.00010, 0.00010, 0.00000, 0.00010, -0.00010, 0.00000, -0.00040, 0.00020, 0.00030, 0.00000, -0.00020, 0.00030, 0.00030, 0.00000, -0.00070, 0.00030, 0.00000, -0.00070, 0.00030, 0.00000, -0.00040, 0.00030, 0.00000, -0.00040, 0.00020, 0.00010, -0.00010, 0.00020, -0.00010, -0.00010, 0.00000, 0.00010, -0.00020, 0.00000, -0.00010, -0.00020, 0.00010, 0.00000, -0.00020, 0.00010, 0.00000, -0.00020, 0.00010, 0.00000, -0.00010, 0.00010, 0.00000, -0.00010, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00010, 0.00000, 0.00000, 0.00010, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00010, 0.00000, 0.00000, 0.00010, 0.00000, 0.00000, -0.00010, 0.00000, 0.00010, -0.00010, 0.00000, 0.00010, 0.00010, -0.00010, 0.00010, 0.00010, 0.00010, 0.00010, 0.00030, -0.00010, 0.00050, 0.00030, 0.00010, 0.00050, 0.00030, 0.00000, 0.00010, 0.00030, 0.00000, 0.00010, 0.00020, 0.00000, 0.00000, 0.00020, 0.00000, 0.00000, 0.00010, 0.00000, 0.00000, 0.00010, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00010, 0.00000, 0.00000, 0.00010, 0.00000, 0.00000, 0.00030, 0.00000, -0.00050, 0.00030, 0.00000, -0.00050, -0.00010, -0.00010, 0.00020, 0.00030, -0.00010, 0.00020, 0.00030, -0.00010, 0.00020, 0.00010, -0.00010, 0.00000, -0.00010, 0.00000, 0.00000, 0.00040, 0.00010, -0.00050, 0.00040, -0.00010, -0.00050, -0.00010, 0.00010, 0.00020, 0.00030, 0.00010, 0.00020, 0.00030, 0.00010, 0.00020, 0.00010, 0.00010, 0.00000, -0.00010, 0.00000, 0.00000, -0.00010, 0.00000, 0.00000, 0.00000, -0.00000, 0.00010, -0.00020, 0.00000, 0.00020, -0.00040, 0.00000, 0.00030, -0.00030, -0.00000, 0.00020, 0.00030, 0.00000, -0.00080, 0.00010, -0.00010, 0.00000, 0.00020, 0.00010, -0.00020, 0.00020, -0.00010, -0.00020, 0.00000, 0.00000, 0.00000, 0.00020, 0.00000, -0.00010, 0.00010, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00020, 0.00000, -0.00010, 0.00010, 0.00000, 0.00000, 0.00030, 0.00000, -0.00060, 0.00000, 0.00010, -0.00010, 0.00000, -0.00010, -0.00010, 0.00010, 0.00010, 0.00000, 0.00010, -0.00010, -0.00010, 0.00010, 0.00010, -0.00010
    ), nrow = 4, byrow = TRUE)
  }

  # Get number of analyses (all combined + each component) from pheno_diff which contains the per-landmark shape difference between -3 and +3 standard deviations from the population mean on each component
  num_of_analyses <- nrow(pheno_diff)

  # create empty list to store results
  all_results <- vector("list", length = num_of_analyses)

  # Get the equivalent pheno group
  if (!pheno %in% names(pheno_names_mouse_groups_mutants_equivalences)) {
    stop(paste0("Unknown phenotype '", pheno, "'. Not found in equivalences."))
  }
  pheno <- pheno_names_mouse_groups_mutants_equivalences[[pheno]]

  # iterate over components
  for(k in 1:num_of_analyses){ 
    
    # Parse the pheno_index parameter to determine landmark range
    pheno_index_range <- eval(parse(text = pheno_index))
    num_landmarks <- length(pheno_index_range)
    
    # Convert landmark indices to column indices (each landmark has x,y,z coordinates)
    # For landmarks 1:54, we need columns 1:162 (54 landmarks * 3 coordinates each)
    coord_indices <- as.vector(sapply(pheno_index_range, function(lm) c(lm*3-2, lm*3-1, lm*3)))
    
    # Determine phenotype group (Manual vs Automated) using the lookup table
    pheno_group <- NULL
    if(!is.null(pheno_names_mouse_groups_mutants$Manual) && pheno %in% pheno_names_mouse_groups_mutants$Manual) {
      pheno_group <- "Manual"
    } else if(!is.null(pheno_names_mouse_groups_mutants$Automated) && pheno %in% pheno_names_mouse_groups_mutants$Automated) {
      pheno_group <- "Automated"
    } else {
      stop(paste0("Unknown phenotype '", pheno, "'. Not found in `pheno_names_mouse_groups`."))
    }

    # Phenotype variable name is the same as the string
    if(pheno %in% pheno_names_mouse_groups_mutants[pheno_group]) stop(paste0("Phenotype variable '", pheno, "' not found in environment"))
    # pheno_data <- get(pheno, inherits = TRUE)
    
    # Determine which mutant DB / mutant LMS variables to use for this group
    mutant_vars <- pheno_names_mouse_mutants[[pheno_group]]
    if(length(mutant_vars) < 2) stop("pheno_names_mouse_mutants must contain two entries per group (db, lms)")
    mutant_db_name  <- mutant_vars[1]
    mutant_lms_name <- mutant_vars[2]
    
    if(!exists(mutant_db_name, inherits = TRUE)) stop(paste0("Mutant DB '", mutant_db_name, "' not found"))
    if(!exists(mutant_lms_name, inherits = TRUE)) stop(paste0("Mutant LMS '", mutant_lms_name, "' not found"))
    
    mutant_db_all  <- get(mutant_db_name, inherits = TRUE)
    mutant_lms_all <- get(mutant_lms_name, inherits = TRUE)
    
    # select specific landmark set
    if(pheno_group == "Automated"){
      mutant_lms_all <- mutant_lms_all[[pheno]]
      alignment = auto.alignment
    } else {alignment = manual.alignment}
    
    # Convert matrix to vector
    alignment <- as.numeric(t(alignment))[coord_indices]
    
    # Subset by coordinate indices
    mutant_data <- as.matrix(mutant_lms_all[, coord_indices])
    
    # If group is manual then swap certain landmarks
    if(pheno_group == "Manual"){
      # Creating new copy of mutant.lms called mutant.lms.fixed
      mutant_data_fixed <- mutant_data
      
      # Swapping relative columns 17 and 18 (each landmark has x,y,z coordinates)
      relative_column_17 <- c(17*3-2, 17*3-1, 17*3)  # Indices for landmark 17
      relative_column_18 <- c(18*3-2, 18*3-1, 18*3)  # Indices for landmark 18
      # Swap the values using vector indexing
      temp <- mutant_data_fixed[relative_column_17]
      mutant_data_fixed[relative_column_17] <- mutant_data_fixed[relative_column_18]
      mutant_data_fixed[relative_column_18] <- temp
      
      # Swapping relative columns 48 and 49
      relative_column_48 <- c(48*3-2, 48*3-1, 48*3)  # Indices for landmark 48
      relative_column_49 <- c(49*3-2, 49*3-1, 49*3)  # Indices for landmark 49
      # Swap the values using vector indexing
      temp <- mutant_data_fixed[relative_column_48]
      mutant_data_fixed[relative_column_48] <- mutant_data_fixed[relative_column_49]
      mutant_data_fixed[relative_column_49] <- temp
      
      # Replace mutant_data with this fixed stuff
      mutant_data <- mutant_data_fixed; rm(mutant_data_fixed, temp)
    }
    
    #### START: CREATE AXIS BETWEEN CONTROL AND AXIS ####
    
    # Set up table to return
    mutant_list <- unique(na.omit(mutant_db_all$Gene))
    
    # generate labels
    genotype <- c()
    for (i in 1:length(mutant_list)){
      genotype <- c(genotype, unique(na.omit(mutant_db_all$Genotype[mutant_db_all$Gene == mutant_list[i] & mutant_db_all$Experimental_Group == "mutant"]))[1])}
    
    mutant_correlations <- data.frame(
      Gene = mutant_list,
      Genotype = genotype,
      Correlation = NA,
      Cosine_similarity = NA
    ); rownames(mutant_correlations) <- mutant_correlations$Gene
    mutant_shapes_fixed <- list()
    mutant_val <- 0
    
    for (j in 1:length(mutant_list)){
      remove_mutant <- F
      # specific group meta
      jmutant <- mutant_db_all[!is.na(mutant_db_all$Gene) & mutant_db_all$Gene == mutant_list[j],]
      
      # Check controls and mutants sample sizes
      control_shape <- na.omit(jmutant$Biosample[jmutant$Experimental_Group == "control"])
      mutant_shape <- na.omit(jmutant$Biosample[jmutant$Experimental_Group == "mutant"])
      
      # Controls first
      if(length(control_shape) < 3){
        # match based on strain, experimental group
        possible_controls <- mutant_db_all[mutant_db_all$Strain %in% jmutant$Strain &
                                            mutant_db_all$Zygosity == "wild" & 
                                            !is.na(mutant_db_all$Zygosity),]
        control_shape <- unique(c(control_shape, possible_controls$Biosample))
        
        # Check if there are enough controls now
        if(length(control_shape) < 3){
          cat("\nNo suitable controls for",mutant_list[j])
          remove_mutant <- T
        }
      }
      
      # Check if mutants exist
      if(length(mutant_shape) == 0){
        cat("\nNo mutants for ", mutant_list[j])
        remove_mutant <- T
      }
      
      # Continue analysis, only if both mutant and controls exist  
      if(length(control_shape) >= 3 && length(mutant_shape) >= 3){
        
        # Remove any NAs
        control_shape <- na.omit(control_shape)
        mutant_shape <- na.omit(mutant_shape)
        
        # If groups are wildly unbalanced, then sample the controls
        if (length(mutant_shape)/length(control_shape) < 0.5){
          control_shape <- sample(control_shape, size = length(mutant_shape))}
        
        # Test for group differences
        test_gdf <- geomorph::geomorph.data.frame(
          shape = mutant_data[c(control_shape, mutant_shape),],
          group = c(rep(0, length(control_shape)),rep(1, length(mutant_shape)))  # set 0 to control, 1 to mutant
        )
        
        test_means <- geomorph::procD.lm(shape ~ group, data = test_gdf, iter = 1000, RRPP = TRUE)
        if(test_means$aov.table$`Pr(>F)`[1] > 0.2){
          remove_mutant <- T
        } else {
          # Calculate axes
          axis <- test_means$coefficients["group",]
          
          # Calculate phenotypic similarity via cosine similarity or correlation
          mutant_correlations[mutant_list[j],"Cosine_similarity"] <- cos(angle.calc(pheno_diff[k,], axis))
          mutant_correlations[mutant_list[j],"Correlation"] <- cor(pheno_diff[k,], axis)
          
          # save landmarks if max cor shape
          if(is.null(mutant)){
            if (abs(mutant_correlations[mutant_list[j],"Correlation"]) > abs(mutant_val)){
              mutant_val <- mutant_correlations[mutant_list[j],"Correlation"]
              mutant_shapes_fixed[["min"]] <- alignment - axis
              mutant_shapes_fixed[["max"]] <- alignment + axis
            }
          } else if(mutant == mutant_list[j]){ # or if provided mutant
            mutant_shapes_fixed[["min"]] <- alignment - axis
            mutant_shapes_fixed[["max"]] <- alignment + axis
          }

        }
      }
      if(remove_mutant){mutant_correlations <- mutant_correlations[setdiff(rownames(mutant_correlations), mutant_list[j]),]}
    }; rm(jmutant, test_means, test_gdf, possible_controls, control_shape, mutant_shape)

    
    # Convert mutant_shapes_fixed vectors to 3D arrays before returning
    if(length(mutant_shapes_fixed) > 0){
      if(!is.null(mutant_shapes_fixed[["min"]])){
        mutant_shapes_fixed[["min"]] <- row2array3d(mutant_shapes_fixed[["min"]], Nlandmarks = num_landmarks)
      }
      if(!is.null(mutant_shapes_fixed[["max"]])){
        mutant_shapes_fixed[["max"]] <- row2array3d(mutant_shapes_fixed[["max"]], Nlandmarks = num_landmarks)
      }
    }

    # Return correlations and landmarks for shapes
    all_results[[k]] <- list(correlation = mutant_correlations, mutant_shapes = mutant_shapes_fixed)
  }

  print(all_results)

  return(all_results)
}
