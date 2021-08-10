#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Updated 2021-08
#HERE: Data aggregation and phenometric estimation


#load  packages
library(tidyverse) #to be tidy
library(ggplot2) #to make plots
library(viridis) #colorblind friendly color palette for plots

#load functions
source("Rcode/phenometric.functions.R")

############################################
## BEGIN DATA AGGREGATION, PHENOMETRIC ESTIMATION

load("data/curations.RData")

all.datasets<-list(cured.data1, cured.data2)
datasets<-c("cur1","cur2")

###########################################3
##  SET PARAMETERS


#Set parameters for onset analyses
minDOY<- c(1,60) #Weibull distributions are not limited to reasonable day of year values, so any onset estimates prior to day of year 60 are replaced with day of year 60. 
maxDOY<-c(365,334) #Weibull distributions are not limited to reasonable day of year values, so any termination estimates after day of year 334 are replaced with day of year 334. 

#lists to store Extreme Value phenometrics
ev.onset<-list()
ev.term<-list()

for(li in c(1:length(all.datasets))) {
  phenodata<-all.datasets[[li]]
  ## DATA AGGREGATIONS:
  #1. F2021: ONLY BY LATITUDE (RNDLAT) (closest to ELE13419, ELE13739 metrics)
  ev.onset.L<-ev.metric(phenodata, "onset", annual=F, elev.strat=F, mindoy = minDOY[li]) 
  ev.term.L<-ev.metric(phenodata, "term", annual=F, elev.strat=F, maxdoy = maxDOY[li]) 

  #2. LS1: BY LATITUDE AND YEAR 
  ev.onset.LY<-ev.metric(phenodata, "onset", annual=T, elev.strat=F, mindoy = minDOY[li]) 
  ev.term.LY<-ev.metric(phenodata, "term", annual=T, elev.strat=F, maxdoy = maxDOY[li]) 

  #3. NEW: BY LATITUDE, YEAR, ELEVATION
  ev.onset.LYE<-ev.metric(phenodata, "onset", annual=T, elev.strat=T, mindoy = minDOY[li]) 
  ev.term.LYE<-ev.metric(phenodata, "term", annual=T, elev.strat=T, maxdoy = maxDOY[li]) 

  #These collect extreme value metrics for 2 data curations, 3 aggregations
  ev.onset[[li]]<-list(ev.onset.L, ev.onset.LY, ev.onset.LYE)
  ev.term[[li]]<-list(ev.term.L, ev.term.LY, ev.term.LYE)

}
onset.file=paste("data/phenometrics/ev.onset.RData", sep="")
save(ev.onset, file=onset.file)
term.file=paste("data/phenometrics/ev.term.RData", sep="")
save(ev.term, file=term.file)

rm(ev.onset, ev.term, onset.file, term.file)



#################################################################3
###    WEIBULL PHENOMETRICS

## THIS CODE TAKES A LONG TIME TO RUN

##When runall = T, the only loop is between data curations, so if R is interrupted during calculations for a data curation, you have to start over
##When runall = F, there is a 2nd loop through species which saves metrics after each species, in case R is interrupted
runall<-T #run everything w/o breaks?
          # if False, will loop through species and save each 

if(runall==T) {

  #lists to store Weibull phenometrics
  we.onset<-list()
  we.term<-list()

  for(li in c(1:length(all.datasets))) {
    phenodata<-all.datasets[[li]]
    ## DATA AGGREGATIONS:
    #1. F2021: ONLY BY LATITUDE (RNDLAT) 
    we.onset.L<-we.metric(phenodata, "onset", annual=F, elev.strat=F, mindoy = minDOY[li]) 
    we.term.L<-we.metric(phenodata, "term", annual=F, elev.strat=F, maxdoy = maxDOY[li]) 
  
    #2. LS1: BY LATITUDE AND YEAR (closest to ELE13731 metrics)
    we.onset.LY<-we.metric(phenodata, "onset", annual=T, elev.strat=F, mindoy = minDOY[li]) 
    we.term.LY<-we.metric(phenodata, "term", annual=T, elev.strat=F, maxdoy = maxDOY[li]) 
  
    #3. NEW: BY LATITUDE, YEAR, ELEVATION
    we.onset.LYE<-we.metric(phenodata, "onset", annual=T, elev.strat=T, mindoy = minDOY[li]) 
    we.term.LYE<-we.metric(phenodata, "term", annual=T, elev.strat=T, maxdoy = maxDOY[li]) 
  
    #These collect weibull metrics for 2 data curations, 3 aggregations
    we.onset[[li]]<-list(we.onset.L, we.onset.LY, we.onset.LYE)
    we.term[[li]]<-list(we.term.L, we.term.LY, we.term.LYE)
  
  }
  onset.file=paste("data/phenometrics/we.onset.RData", sep="")
  save(we.onset, file=onset.file)
  term.file=paste("data/phenometrics/we.term.RData", sep="")
  save(we.term, file=term.file)

rm(we.onset, we.term, onset.file, term.file)

} elseif(runall==F) {   #if runall=F

  #lists to store Weibull phenometrics
  we.onset.1<-list()
  we.term.1<-list()
  we.onset.2<-list()
  we.term.2<-list()
  
  for(li in c(1:length(all.datasets))) {
    phenodata<-all.datasets[[li]]
    ## DATA AGGREGATIONS:
    #1. F2021: ONLY BY LATITUDE (RNDLAT)
    for(ji in 1:length(sort(unique(phenodata$name2)))) {
      sp.ji<-sort(unique(phenodata$name2))[ji]
      we.onset.lat<-we.metric(filter(phenodata, name2==sp.ji), "onset", annual=F, elev.strat=F, mindoy = minDOY[li]) 
      we.term.lat<-we.metric(filter(phenodata, name2==sp.ji), "term", annual=F, elev.strat=F, maxdoy = maxDOY[li]) #F2020/F2021 metric
  
      #2. LS1: BY LATITUDE AND YEAR
      we.onset.yr<-we.metric(filter(phenodata, name2==sp.ji),"onset",annual=T,elev.strat=F,minFP, mindoy = minDOY[li]) 
      we.term.yr<-we.metric(filter(phenodata, name2==sp.ji),"term",annual=T,elev.strat=F,minFP, maxdoy = maxDOY[li]) 
  
      #3. NEW: BY LATITUDE, YEAR, ELEVATION
      we.onset.yr.alt<-we.metric(filter(phenodata, name2==sp.ji),"onset",annual=T,elev.strat=T,minFP, mindoy = minDOY[li]) 
      we.term.yr.alt<-we.metric(filter(phenodata, name2==sp.ji),"term",annual=T,elev.strat=T,minFP, maxdoy = maxDOY[li]) 
  
      if(ji==1) {
        we.onset.L<-we.onset.lat
        we.onset.LY<-we.onset.yr
        we.onset.LYE<-we.onset.yr.alt
        we.term.L<-we.term.lat
        we.term.LY<-we.term.yr
        we.term.LYE<-we.term.yr.alt
      } else {
        we.onset.L<-bind_rows(we.onset.L, we.onset.lat)
        we.onset.LY<-bind_rows(we.onset.LY, we.onset.yr)
        we.onset.LYE<-bind_rows(we.onset.LYE, we.onset.yr.alt)
        we.term.L<-bind_rows(we.term.L, we.term.lat)
        we.term.LY<-bind_rows(we.term.LY, we.term.yr)
        we.term.LYE<-bind_rows(we.term.LYE, we.term.yr.alt)
      }
      rm(we.onset.lat, we.onset.yr, we.onset.yr.alt)  
      rm(we.term.lat, we.term.yr, we.term.yr.alt)  
    
      onset.tempfile=paste("data/phenometrics/we.onset.temp.RData", sep="")
      save(we.onset.L, we.onset.LY, we.onset.LYE, file=onset.tempfile)
      term.tempfile=paste("data/phenometrics/we.term.temp.RData", sep="")
      save(we.term.L, we.term.LY, we.term.LYE, file=term.tempfile)
    
    }
    #These collect weibull metrics for 2 data curations, 3 aggregations
    we.onset[[li]]<-list(we.onset.L, we.onset.LY, we.onset.LYE)
    we.term[[li]]<-list(we.term.L, we.term.LY, we.term.LYE)
    
  }  
  onset.file=paste("data/phenometrics/we.onset.RData", sep="")
  save(we.onset, file=onset.file)
  term.file=paste("data/phenometrics/we.term.RData", sep="")
  save(we.term, file=term.file)
  rm(we.onset, we.term, onset.file, term.file)
}

## END PHENOMETRIC ESITMATION