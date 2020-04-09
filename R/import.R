library(readxl)


readSC <- function(file = "Data/Shotgun_Plastic.xlsx") {
  dataset.sc <- read_xlsx(path = file,
                          sheet = "compar spectra",
                          col_names = TRUE)
  colnames(dataset.sc) <- c("ProteinID", "P3","P2","P1","C3","C2","C1")
  return(dataset.sc)
}



