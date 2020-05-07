library(readxl)


readSC <- function(file = "Data/Shotgun_Plastic.xlsx") {
  dataset.sc <- read_xlsx(path = file,
                          sheet = "compar spectra",
                          col_names = TRUE)
  colnames(dataset.sc) <- c("ProteinID", "P3","P2","P1","C3","C2","C1")
  return(dataset.sc)
}

# Spectral Count

dataset.sc <- read_xlsx(path = "Data/Shotgun_Plastic.xlsx",
                        sheet = "compar spectra",
                        col_names = TRUE)

colnames(dataset.sc) <- c("ProteinID", "P3","P2","P1","C3","C2","C1")


dataset.sc$p.value <- apply(dataset.sc[2:7], 1, function(x){
  if(x[1]==x[2] && x[2]==x[3] && x[3]==x[4] && x[4]==x[5] && x[5]==x[6]){return(1)}
  else if ((x[1]==x[2] && x[2]==x[3]) && (x[4]==x[5] && x[5]==x[6])) {return(0)} else {
    t.test(x[1:3], x[4:6])->t.results
    return(t.results$p.value)
  }
})

#PAI

dataset.PAI <- read_xlsx(path = "Data/Shotgun_Plastic.xlsx",
                        sheet = "compar PAI",
                        col_names = TRUE)

colnames(dataset.sc) <- c("ProteinID", "P3","P2","P1","C3","C2","C1")


dataset.PAI$p.value <- apply(dataset.PAI[2:7], 1, function(x){
  if(x[1]==x[2] && x[2]==x[3] && x[3]==x[4] && x[4]==x[5] && x[5]==x[6]){return(1)}
  else if ((x[1]==x[2] && x[2]==x[3]) && (x[4]==x[5] && x[5]==x[6])) {return(0)} else {
    t.test(x[1:3], x[4:6])->t.results
    return(t.results$p.value)
  }
})


#emPAI

dataset.emPAI <- read_xlsx(path = "Data/Shotgun_Plastic.xlsx",
                        sheet = "compar emPAI",
                        col_names = TRUE)

colnames(dataset.emPAI) <- c("ProteinID", "P3","P2","P1","C3","C2","C1")


dataset.emPAI$p.value <- apply(dataset.emPAI[2:7], 1, function(x){
  if(x[1]==x[2] && x[2]==x[3] && x[3]==x[4] && x[4]==x[5] && x[5]==x[6]){return(1)}
  else if ((x[1]==x[2] && x[2]==x[3]) && (x[4]==x[5] && x[5]==x[6])) {return(0)} else {
    t.test(x[1:3], x[4:6])->t.results
    return(t.results$p.value)
  }
})



#NSAF

dataset.NSAF <- read_xlsx(path = "Data/Shotgun_Plastic.xlsx",
                        sheet = "compar NSAF",
                        col_names = TRUE)

colnames(dataset.NSAF) <- c("ProteinID", "P3","P2","P1","C3","C2","C1")


dataset.NSAF$p.value <- apply(dataset.NSAF[2:7], 1, function(x){
  if(x[1]==x[2] && x[2]==x[3] && x[3]==x[4] && x[4]==x[5] && x[5]==x[6]){return(1)}
  else if ((x[1]==x[2] && x[2]==x[3]) && (x[4]==x[5] && x[5]==x[6])) {return(0)} else {
    t.test(x[1:3], x[4:6])->t.results
    return(t.results$p.value)
  }
})










