#' This is the data for examples
#'
#' \itemize{
#'   \item All_header. ID for all samples that are available in at least one type of data
#'   \item pheno. phenotype file. 1st column is ID, 2nd column is disease status, 3rd column is age, 4th column is gender
#'   \item G. genotypes for 3 genes. G[[i]] is for the ith genes. In each G[[i]], the 1st column is ID and the rest columns are genotypes
#'   \item M. methylated sites 3 genes. M[[i]] is for the ith genes. In each M[[i]], the 1st column is ID and the rest columns are methylated sites
#'   \item R. RNA expression for 3 genes. R[[i]] is for the ith genes. In each R[[i]], the 1st column is ID and the rest columns (usually one column) are RNA expression
#' }
#'
#' @name example_data
#' @aliases All_header pheno G M R
#' @docType data
#' @usage data(example_data)
NULL
