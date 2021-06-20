#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Updated 2021-06
#HERE: Data aggregation and phenometric estimation


#load  packages
library(tidyverse) #to be tidy

#load functions
source("Rcode/phenometric.functions.R")

############################################
## BEGIN DATA AGGREGATION, PHENOMETRIC ESTIMATION

load("data/curations.RData")

# cured.data2)
occurrences<-split(cured.data2, cured.data2$name2)
occurrences.a<-occurrences[c(1:10)]
#################################################################3
###    WEIBULL PHENOMETRICS


#functions to call in lapply
est.onset.alt<-function(X) {
  we.metric(X,"onset",annual=T,elev.strat=T,7, mindoy = 60) 
}
est.term.alt<-function(X) {
  we.metric(X,"term",annual=T,elev.strat=T,7, maxdoy = 334) 
}

#estimate onset
we.onset.yr.alt<-lapply(occurrences.a, est.onset.alt)

save(we.onset.yr.alt,file="data/phenometrics/we.onset2a.yr.alt.RData")

#estimate termination
we.term.yr.alt<-lapply(occurrences.a, est.term.alt)

save(we.term.yr.alt,file="data/phenometrics/we.term2a.yr.alt.RData")



# the rest of the species
occurrences.b<-occurrences[c(11:length(occurrences))]


#################################################################3

#estimate onset
we.onset.yr.alt<-lapply(occurrences.b, est.onset.alt)

save(we.onset.yr.alt,file="data/phenometrics/we.onset2b.yr.alt.RData")

#estimate termination
we.term.yr.alt<-lapply(occurrences.b, est.term.alt)

save(we.term.yr.alt,file="data/phenometrics/we.term2b.yr.alt.RData")

