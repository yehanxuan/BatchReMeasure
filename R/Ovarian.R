#' @title Ovarian cancer data measured from Agilent micro-array platform
#' @description The matrix of gene expression of samples
#' @format A data frame with 11861 rows and 306 columns:
#' \describe{
#'   The expression matrix with the gene in rows and samples in column
#'}
#' @source Provided by Mayo clinic
"Agilent.dat"

#' @title Ovarian cancer data measured from Agilent micro-array platform
#' @description The recurrence status of samples: yes or no
#' @format A data frame with 306 rows and 1 variables:
#' \describe{
#'   \item{\code{recur}}{character Whether the samples from Agilent platform are recurrence cases}
#'}
#' @source Provided by Mayo clinic
"Agilent.recur"

#' @title Ovarian cancer data measured from Agilent micro-array platform
#' @description The cancer subtype variables: C1-MES, C2-IMM, C4-DIF, C5-PRO
#' @format A data frame with 306 rows and 1 variables:
#' \describe{
#'   \item{\code{ordered.subtype}}{character Cancer subtype variables}
#'}
#' @source Provided by Mayo clinic
"Agilent.subtype"

#' @title Ovarian cancer data measured from Agilent RNA-Seq platform
#' @description The matrix of gene expression of samples
#' @format A data frame with 11861 rows and 97 variables:
#' \describe{
#'   The expression matrix with the gene in rows and samples in column
#'}
#' @source Provided by Mayo clinic
"RNAseq.dat"

#' @title Ovarian cancer data measured from Agilent RNA-Seq platform
#' @description The recurrence status of samples: yes or no
#' @format A data frame with 97 rows and 1 variables:
#' \describe{
#'   \item{\code{recur}}{character Whether the samples from RNASeq platform are recurrence cases}
#'}
#' @source Provided by Mayo clinic
"RNAseq.recur"

#' @title Ovarian cancer data measured from Agilent RNA-Seq platform
#' @description The cancer subtype variables: C1-MES, C2-IMM, C4-DIF, C5-PRO
#' @format A data frame with 97 rows and 1 variables:
#' \describe{
#'   \item{\code{ordered.subtype}}{character  Cancer subtype variables}
#'}
#' @source Provided by Mayo clinic
"RNAseq.subtype"

#' @title Common remeasured samples between the two data set
#' @description The recurrence status of samples: yes or no
#' @format A data frame with 47 rows and 1 variables:
#' \describe{
#'   \item{\code{recur}}{character Whether the samples are recurrence cases}
#'}
#' @source Provided by Mayo clinic
"remeasure.recur"

#' @title Common remeasured samples between the two data set
#' @description The cancer subtype variables: C1-MES, C2-IMM, C4-DIF, C5-PRO
#' @format A data frame with 47 rows and 1 variables:
#' \describe{
#'   \item{\code{ordered.subtype}}{character Cancer subtype variables}
#'}
#' @source Provided by Mayo clinic
"remeasure.subtype"




