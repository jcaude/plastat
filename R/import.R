library(tibble)
library(dplyr)
library(readxl)

THRESHOLD_LogFC = 1
THRESHOLD_PVALUE = 5e-2

#' Read Shotgun Excel Result Files
#'
#' @param file name of the excel file
#' @param sheet name of the excel worksheet
#'
#' @return a tibble
#' @export
#'
#' @examples
#' SC <- readXL()
readXL <- function(file = "Data/Shotgun_Plastic.xlsx",
                   sheet = "compar spectra") {

  dataset <- read_xlsx(path = file,sheet = sheet,col_names = TRUE)
  colnames(dataset) <-
    c("ProteinID", "P3", "P2", "P1", "C3", "C2", "C1")
  return(dataset)
}

my_ttest <- function(x) {

  # init.
  n <- length(x)
  m <- n/2
  xr <- 1:m
  yr <- (m+1):n

  # t.test
  if (length(unique(x[xr])) == 1 && length(unique(x[yr])) == 1) {
    return(1)
  } else {
    t.results <- t.test(x[xr], x[yr])
    return(t.results$p.value)
  }
}

evalTtest <- function(dataset) {

  # compute FC
  dataset <- dataset %>%
    mutate(P.AVG = (P1+P2+P3)/3) %>%
    mutate(C.AVG = (C1+C2+C3)/3) %>%
    mutate(LogFC = log2(P.AVG) - log2(C.AVG))

  # compute PVAL
  dataset$P.value <- apply(dataset[2:7],1,my_ttest)

  # compute Selection Flag
  dataset <- dataset %>%
    mutate(Selected = if_else(P.value < THRESHOLD_PVALUE &
                                abs(LogFC) > THRESHOLD_LogFC, TRUE,FALSE))

  # end.
  return(dataset)
}

calcDataset <- function() {

  # Spectral Count
  dataset.sc <- readXL(file="Data/Shotgun_Plastic.xlsx",
                       sheet = "compar spectra")
  dataset.sc <- evalTtest(dataset = dataset.sc)

  #PAI
  dataset.PAI <- readXL(file = "Data/Shotgun_Plastic.xlsx",
                        sheet = "compar PAI")
  dataset.PAI <- evalTtest(dataset = dataset.PAI)

  #emPAI
  dataset.emPAI <- readXL(file = "Data/Shotgun_Plastic.xlsx",
                          sheet = "compar emPAI")
  dataset.emPAI <- evalTtest(dataset = dataset.emPAI)

  #NSAF
  dataset.NSAF <- readXL(file = "Data/Shotgun_Plastic.xlsx",
                         sheet = "compar NSAF")
  dataset.NSAF <- evalTtest(dataset = dataset.NSAF)

  # end
  return(list(dataset.sc=dataset.sc,
              dataset.PAI=dataset.PAI,
              dataset.emPAI=dataset.emPAI,
              dataset.NSAF=dataset.NSAF))
}


# Compare selected adsorbed proteins..
CmpSelected <- function() {

  # read data & calc t.test
  data <- calcDataset()
  data.names <- names(data)

  # here is the trick.. we are using the function 'combn' which gives all
  # possible combination of a vector with m elements
  perms <- combn(data.names, m = 2)

  # then we use apply to iterate over the combination matrix
  # - MARGIN = 2, because we are iterating by columns
  t <- apply(X = perms, MARGIN = 2, function(p) {
    # - p[1] == name of the dataset
    # - data[ p[1] ] === list size=1 of data named p[1]
    # - data[[ p[1] ]] the second [ ] unlist the results
    #                  (thus pop out the first item of the list)
    ds1 <- data[[p[1]]]  # first dataset
    ds2 <- data[[p[2]]]  # second dataset

    # print what we compare
    name1 <- gsub("dataset\\.","",p[1])  # -- find & replace, \\ because '.' is a special character
    name2 <- gsub("dataset\\.","",p[2])
    cat("Comparison of",name1,"vs",name2,"\n")
    cat("------------------------------------\n")

    # print the comparison results
    # - FALSE == not selected
    # - TRUE == selected
    print(table(ds1 %>% pull(Selected),
                ds2 %>% pull(Selected)))
    cat("\n")
  })
}
