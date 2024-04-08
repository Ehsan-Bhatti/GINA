#Install required packages
list_of_packages <-c('GEOquery', 'limma', 'openxlsx', 'edgeR')

for (i in list_of_packages){
     if(! i %in% installed.packages()){

       if(i == 'opeopenxlsx') {
         install.packages(i, repos='http://cran.us.r-project.org', dependencies = TRUE)
       }
       else {
             if (!require("BiocManager", quietly = TRUE)) {
              install.packages("BiocManager", repos='http://cran.us.r-project.org', dependencies = TRUE)
           }
         BiocManager::install(i)
       }
     }
}

library(GEOquery)
library(limma)
library(openxlsx)
library(edgeR)

# Method to get data from GEO
get_data_from_geo <- function (data_id, directory_to_save) {
  return (GEOquery::getGEO(data_id, destdir = directory_to_save))

}

# Method to read a saved GEO file
read_geo_file <- function (file_name) {
  print("Reading File")
  out <- (GEOquery::getGEO(filename = file_name))
  print("File Read")
  return(out)
}

# Method to get phenotype data
pheno_data <- function (df, file_name=NULL) {
  pdf <- pData(df)
  pdf <- as.data.frame(pdf, row.names=row.names(pdf), col.names = colnames(pdf))
  if(!is.null(file_name)) {
    write.csv(pdf, file_name, row.names=TRUE)
  }
  return(pdf)
}

# Method to get expression data
expression_data <- function (df, file_name=NULL) {
  edf <- exprs(df)
  edf <- as.data.frame(edf, row.names = row.names(expression_data),
                                   col.names = colnames(edf))
  if(!is.null(file_name)) {
    write.csv(edf, file_name, row.names=TRUE)
  }
  return(edf)
}

# Method to get feature data
feature_data <- function (df, file_name=NULL) {
  fdf <- fData(df)
  fdf <- as.data.frame(fdf, row.names = row.names(fdf), col.names = colnames(fdf))

  if(!is.null(file_name)) {
    write.csv(fdf, file_name, row.names=TRUE)
  }
  return(fdf)
}


# Method to see columns in dataframe.
see_columns <- function (df) {
  return(colnames(df))
}

unique_values_in_col <- function (df, col) {

  print(sprintf("Unique values in: %s", col))

  print(unique(df[[col]]))
}



# Method to run differenetial expression analysis
differential_expression_analysis <-function(geo_data, phenotype_col_name,
                                            list_of_groups,
                                            e_bayes_proportion = 0.01, decide_test_adj_method='BH',
                                            decide_test_p_val=0.05,
                                            file_name=NULL) {

  print("Adding Phenotype Group Column")
  geo_data$GINA_pheno_group <- NA

  group_val <- 1
  for (g in list_of_groups) {
    string_group_val <- as.character(group_val)
    pg <- paste('group', string_group_val, sep = "_")
    geo_data$GINA_pheno_group[geo_data[[phenotype_col_name]] %in% g] <- pg

    group_val <- group_val + 1
  }

  omit_count <- sum(is.na(geo_data$GINA_pheno_group))

  if (omit_count > 0) {
    fil <- colnames(geo_data)[!is.na(geo_data$GINA_pheno_group)]
    test_data <- geo_data[, fil]
  }
  else {
    test_data <- geo_data
  }

  num_for_table <- length(test_data$GINA_pheno_group)

  print(sprintf("Number of samples ommited because phenotype did not appear in any group:   %s.", omit_count ))

  print("Starting Differntial Analysis")

  design <- model.matrix(~GINA_pheno_group + 0, test_data)
  colnames(design) <- sub('^GINA_pheno_group', '', colnames(design))


  l_fit <- limma::lmFit(test_data, design)

  contrast_name <- combn(colnames(design), 2, function(x) paste0(x[1], "-", x[2]))

  con <- makeContrasts(contrasts=contrast_name, levels=design)
  con_fit <- contrasts.fit(l_fit, con)

  e_fit <- eBayes(con_fit, proportion=e_bayes_proportion)
  e_fit_top_table <- topTable(e_fit, number = num_for_table)
  diff_genes <- decideTests(e_fit, adjust.method=decide_test_adj_method, p.value=decide_test_p_val)

  pheno_data <- pData(geo_data)
  expression_data <- exprs(geo_data)
  feature_data <- fData(geo_data)

  pheno_data <- as.data.frame(pheno_data, row.names=row.names(pheno_data), col.names = colnames(pheno_data))

  feature_data <- as.data.frame(feature_data, row.names = row.names(feature_data), col.names = colnames(feature_data))

  con_fit <- as.data.frame(con_fit, row.names =rownames(con_fit$coefficients, prefix = ""))

  e_fit <- as.data.frame(e_fit, row.names = rownames(e_fit$coefficients, prefix = ""))

  e_fit_top_table <- as.data.frame(e_fit_top_table, row.names = rownames(e_fit_top_table$logFC, prefix = ""))

  expression_data <- as.data.frame(expression_data, row.names = row.names(expression_data),
                                   col.names = colnames(expression_data))

  diff_genes <- as.matrix(diff_genes, row.names=rownames(diff_genes, prefix=""))

  feature_data <- merge(feature_data, diff_genes, by.x =0, by.y=0, all = TRUE)
  rownames(feature_data) <- feature_data$Row.names
  feature_data$Row.names <- NULL

  empirical_bayes_proportion <- e_bayes_proportion
  decide_test_adjust_method <- decide_test_adj_method
  decide_test_p_value <- decide_test_p_val
  num_samples_omitted <- omit_count


  experiment_df <- data.frame(empirical_bayes_proportion, decide_test_adjust_method, decide_test_p_value,
                              num_samples_omitted)

  output <- list(
    "contrast_fit" = con_fit, "empirical_bayes_fit" = e_fit, "table_empirical_bayes" = e_fit_top_table,
    "differentially_expressed_genes" = diff_genes, "calculation_info" = experiment_df,
    "phenotype_data" = pheno_data, "expression_data" = expression_data, "feature_data" = feature_data
  )

  if (!is.null(file_name)) {
    print("Writing To File")
    openxlsx::write.xlsx(output, file=file_name, rowNames=TRUE, colNames=TRUE)
  }

  return(output)
}

# Method to read a single RNA sequence file
read_single_rna_seq_file <- function (path_to_file, sep, header, file_name=NULL) {
  f <- read.delim(file=path_to_file, header = header, sep = sep)
  rna_seq <- edgeR::DGEList(f)

  counts <- rna_seq$counts
  counts <- as.data.frame(counts, row.names=row.names(counts), col.names = colnames(counts))

  samples <- rna_seq$samples
  samples <- as.data.frame(samples, row.names = row.names(samples), col.names = colnames(samples))

  if (exists('genes', where = rna_seq)) {
    genes <- rna_seq$genes
    genes <- as.data.frame(genes, row.names = row.names(genes), col.names = colnames(genes))


  } else {
    genes <- data.frame(name = row.names(counts))
   }

  out <- list('expression_data'=counts, 'phenotype_data'=samples, 'feature_data'=genes)

  if (!is.null(file_name)) {
    print("Writing To File")
    openxlsx::write.xlsx(out, file=file_name, rowNames=TRUE, colNames=TRUE)
  }

  return(out)
}

# Method to read multiple RNA sequence files
read_multiple_rna_seq_files <- function(files, path_to_files, cols, sep, header, file_name=NULL){
  rna_seq <- edgeR::readDGE(files = files, path = path_to_files, columns = cols, sep = sep, header = header)
  counts <- rna_seq$counts
  counts <- as.data.frame(counts, row.names=row.names(counts), col.names = colnames(counts))

  samples <- rna_seq$samples
  samples <- as.data.frame(samples, row.names = row.names(samples), col.names = colnames(samples))

  if (exists('genes', where = rna_seq)) {
    genes <- rna_seq$genes
    genes <- as.data.frame(genes, row.names = row.names(genes), col.names = colnames(genes))


  } else {
    genes <- data.frame(name = row.names(counts))
   }

  out <- list('expression_data'=counts, 'phenotype_data'=samples, 'feature_data'=genes)

  if (!is.null(file_name)) {
    print("Writing To File")
    openxlsx::write.xlsx(out, file=file_name, rowNames=TRUE, colNames=TRUE)
  }

  return(out)
}

# Run differential expression analysis of RNA sequences.
rna_seq_differential_expression_analysis <- function(phenotype_data, expression_data, feature_data,
                                                     phenotype_col_name, list_of_groups, normalising_method,
                                                     e_bayes_proportion = 0.01,
                                                     decide_test_adj_method='BH', decide_test_p_val=0.05,
                                                     file_name=NULL) {

  print("Adding Phenotype Group Column")
  phenotype_data$GINA_pheno_group <- NA

  group_val <- 1
  for (g in list_of_groups) {
    string_group_val <- as.character(group_val)
    pg <- paste('group', string_group_val, sep = "_")
    phenotype_data$GINA_pheno_group[phenotype_data[[phenotype_col_name]] %in% g] <- pg

    group_val <- group_val + 1
  }

  dg <- edgeR::DGEList(counts=expression_data, samples = phenotype_data, genes = feature_data)


  omit_count <- sum(is.na(phenotype_data$GINA_pheno_group))

  if (omit_count > 0) {
    test_data <- dg[,which(!is.na(dg$samples$GINA_pheno_group))]
  }
  else {
    test_data <- dg
  }

  num_for_table <- length(test_data$samples$GINA_pheno_group)

  test_data <- edgeR::calcNormFactors(test_data, method=normalising_method)


  print(sprintf("Number of samples ommited because phenotype did not appear in any group:   %s.", omit_count ))

  print("Starting Differntial Analysis")

  design <- model.matrix(~GINA_pheno_group + 0, test_data$samples)
  colnames(design) <- sub('^GINA_pheno_group', '', colnames(design))

  v_data <- limma::voom(test_data, design)

  l_fit <- limma::lmFit(v_data, design)

  contrast_name <- combn(colnames(design), 2, function(x) paste0(x[1], "-", x[2]))

  con <- makeContrasts(contrasts=contrast_name, levels=design)
  con_fit <- contrasts.fit(l_fit, con)

  e_fit <- eBayes(con_fit, proportion=e_bayes_proportion)
  e_fit_top_table <- topTable(e_fit, number = num_for_table)
  diff_genes <- decideTests(e_fit, adjust.method=decide_test_adj_method, p.value=decide_test_p_val)


  phenotype_data <- as.data.frame(phenotype_data, row.names=row.names(phenotype_data),
                                  col.names = colnames(phenotype_data))

  feature_data <- as.data.frame(feature_data, row.names = row.names(feature_data), col.names = colnames(feature_data))

  con_fit <- as.data.frame(con_fit, row.names =rownames(con_fit$coefficients, prefix = ""))

  e_fit <- as.data.frame(e_fit, row.names = rownames(e_fit$coefficients, prefix = ""))

  e_fit_top_table <- as.data.frame(e_fit_top_table, row.names = rownames(e_fit_top_table$logFC, prefix = ""))

  expression_data <- as.data.frame(expression_data, row.names = row.names(expression_data),
                                   col.names = colnames(expression_data))

  norm_expression <- v_data$E
  norm_expression <- as.data.frame(norm_expression, row.names = row.names(norm_expression),
                                   col.names = colnames(norm_expression))

  norm_weights <- v_data$weights
  norm_weights <- as.data.frame(norm_weights, row.names = row.names(norm_weights), col.names = colnames(norm_weights))

  diff_genes <- as.matrix(diff_genes, row.names=rownames(diff_genes, prefix=""))

  feature_data <- merge(feature_data, diff_genes, by.x =0, by.y=0, all = TRUE)
  rownames(feature_data) <- feature_data$Row.names
  feature_data$Row.names <- NULL

  empirical_bayes_proportion <- e_bayes_proportion
  decide_test_adjust_method <- decide_test_adj_method
  decide_test_p_value <- decide_test_p_val
  num_samples_omitted <- omit_count
  normalisation_method <- normalising_method


  experiment_df <- data.frame(empirical_bayes_proportion, decide_test_adjust_method, decide_test_p_value,
                              num_samples_omitted, normalisation_method)

  output <- list(
    "contrast_fit" = con_fit, "empirical_bayes_fit" = e_fit, "table_empirical_bayes" = e_fit_top_table,
    "differentially_expressed_genes" = diff_genes, "calculation_info" = experiment_df,
    "phenotype_data" = phenotype_data, "expression_data" = expression_data, "feature_data" = feature_data,
    'normalised_expression_data' = norm_expression, 'normalising_expression_weights' = norm_weights)


  if (!is.null(file_name)) {
    print("Writing To File")
    openxlsx::write.xlsx(output, file=file_name, rowNames=TRUE, colNames=TRUE)
  }

  return(output)


}



