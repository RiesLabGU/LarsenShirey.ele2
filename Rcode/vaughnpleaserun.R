#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Updated 2021-06
#HERE: Data aggregation and phenometric estimation


#load  packages
library(tidyverse) #to be tidy
library(phenesse)

#load functions
#source("Rcode/phenometric.functions.R")
#################################################################3
###    WEIBULL PHENOMETRICS FUNCTIONS
we.est.o<-function(doys) {
  y<-suppressWarnings(weib_percentile(doys, 0.01))
  return(y[[1]])
}
we.est.t<-function(doys) {
  y<-suppressWarnings(weib_percentile(doys, 0.99))
  return(y[[1]])
}

#functions to call in lapply
est.onset.alt<-function(X) {
  weib<-X %>%
    dplyr::select(name, region, name2, rndLat, year, alt, doy, elev.yr) %>%
    group_by(name,region,name2, rndLat, year, elev.yr) %>%
    dplyr::mutate(n.occ=length(unique(doy)),alt=min(alt, na.rm=T)) %>%
    dplyr::mutate(sampleEffort=(n.occ/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,7, na.rm=T)))) %>%
    filter(n.occ>=5)
  if(nrow(weib)>=5) {
    #Calculate phenometrics!
    weib<- weib%>% 
      group_by(name,region,name2, rndLat, year, elev.yr) %>%
      summarize(ref=paste("onset.we"),
                n.occ=mean(n.occ,na.rm=T),
                sampleEffort=mean(sampleEffort, na.rm=T), 
                metric=max(1,we.est.o(doy),na.rm=T),
                alt=min(alt, na.rm=T)) 
  } else { 
    #Don't calculate metrics if data are insufficient
    weib<-X %>%
      group_by(name,region,name2, rndLat, year,elev.yr) %>%
      summarize(metric=NA)
  } 
}
#functions to call in lapply
est.term.alt<-function(X) {
  weib<-X %>%
    dplyr::select(name, region, name2, rndLat, year, alt, doy, elev.yr) %>%
    group_by(name,region,name2, rndLat, year, elev.yr) %>%
    dplyr::mutate(n.occ=length(unique(doy)),alt=min(alt, na.rm=T)) %>%
    dplyr::mutate(sampleEffort=(n.occ/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,7, na.rm=T)))) %>%
    filter(n.occ>=5)
  if(nrow(weib)>=5) {
    #Calculate phenometrics!
    weib<- weib%>% 
      group_by(name,region,name2, rndLat, year, elev.yr) %>%
      summarize(ref="term.we",
                n.occ=mean(n.occ,na.rm=T),
                sampleEffort=mean(sampleEffort, na.rm=T), 
                metric=min(365,we.est.t(doy),na.rm=T),
                alt=min(alt, na.rm=T))
  } else { 
    #Don't calculate metrics if data are insufficient
    weib<-X %>%
      group_by(name,region,name2, rndLat, year, elev.yr) %>%
      summarize(metric=NA)
  } 
}

############################################
## BEGIN DATA AGGREGATION, PHENOMETRIC ESTIMATION

load("../data/curations.RData")

# cured.data2)
occurrences<-split(cured.data2, cured.data2$name2)
occurrences.a<-occurrences[c(1:10)]

#################################################################3
###    WEIBULL PHENOMETRICS

#estimate onset
we.onset.yr.alt<-lapply(occurrences.a, est.onset.alt)

save(we.onset.yr.alt,file="../data/phenometrics/we.onset2a.yr.alt.RData")

#estimate termination
we.term.yr.alt<-lapply(occurrences.a, est.term.alt)

save(we.term.yr.alt,file="../data/phenometrics/we.term2a.yr.alt.RData")



# the rest of the species
occurrences.b<-occurrences[c(11:length(occurrences))]


#################################################################3

#estimate onset
we.onset.yr.alt<-lapply(occurrences.b, est.onset.alt)

save(we.onset.yr.alt,file="../data/phenometrics/we.onset2b.yr.alt.RData")

#estimate termination
we.term.yr.alt<-lapply(occurrences.b, est.term.alt)

save(we.term.yr.alt,file="../data/phenometrics/we.term2b.yr.alt.RData")

