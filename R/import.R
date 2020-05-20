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

cmpDataset <- function(data) {

  dataset.selected <- data.frame(data$dataset.sc$ProteinID,
                                 data$dataset.sc$Selected,
                                 data$dataset.PAI$Selected,
                                 data$dataset.emPAI$Selected,
                                 data$dataset.NSAF$Selected)

  dataset.selected <- dataset.selected %>%
    rename(SC = dataset.sc.Selected) %>%
    rename(PAI = dataset.PAI.Selected) %>%
    rename(emPAI = dataset.emPAI.Selected) %>%
    rename(NSAF = dataset.NSAF.Selected) %>%

    dataset.selected %>%
    filter(SC == TRUE | PAI == TRUE | emPAI == TRUE | NSAF == TRUE)

  dataset.selected %>%
    filter(SC == TRUE & PAI == TRUE & emPAI == TRUE & NSAF == TRUE)
}