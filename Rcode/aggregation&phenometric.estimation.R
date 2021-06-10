#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Updated 2021-06
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
datasets<-c("few.filters","many.filters")

###########################################3
##  SET PARAMETERS


#Set parameters for onset analyses
minFP<-7   #For calculating an effort metric, we correct for the day-of-year span for which data are available. Here we set a minimum day-of-year span at 7, as the effort metric is misleading when that span is less than 1 week.
minDOY<- c(1,60) #Weibull distributions are not limited to reasonable day of year values, so any onset estimates prior to day of year 60 are replaced with day of year 60. 
maxDOY<-c(365,334) #Weibull distributions are not limited to reasonable day of year values, so any termination estimates after day of year 334 are replaced with day of year 334. 

#lists to store Extreme Value phenometrics
ev.onset<-list()
ev.term<-list()

for(li in c(1:length(all.datasets))) {
  phenodata<-all.datasets[[li]]
  ## DATA AGGREGATIONS:
  #1. F2021: ONLY BY LATITUDE (RNDLAT)
  ev.onset.1<-ev.metric(phenodata,"onset",annual=F,elev.strat=F,minFP, mindoy = minDOY[li]) 
  ev.term.1<-ev.metric(phenodata,"term",annual=F,elev.strat=F,minFP, maxdoy = maxDOY[li]) #F2020/F2021 metric

  #2. LS1: BY LATITUDE AND YEAR
  ev.onset.yr<-ev.metric(phenodata,"onset",annual=T,elev.strat=F,minFP, mindoy = minDOY[li]) 
  ev.term.yr<-ev.metric(phenodata,"term",annual=T,elev.strat=F,minFP, maxdoy = maxDOY[li]) 

  #3. NEW: BY LATITUDE, YEAR, ELEVATION
  ev.onset.yr.alt<-ev.metric(phenodata,"onset",annual=T,elev.strat=T,minFP, mindoy = minDOY[li]) 
  ev.term.yr.alt<-ev.metric(phenodata,"term",annual=T,elev.strat=T,minFP, maxdoy = maxDOY[li]) 

  #These collect extreme value metrics for 2 data curations, 3 aggregations
  ev.onset[[li]]<-list(ev.onset.1, ev.onset.yr, ev.onset.yr.alt)
  ev.term[[li]]<-list(ev.term.1, ev.term.yr, ev.term.yr.alt)

  filename=paste("data/phenometrics/ev.metrics.RData", sep="")
  save(ev.onset, ev.term, file=filename)
}
filename=paste("data/phenometrics/all.ev.metrics.RData", sep="")
save(ev.onset, ev.term, file=filename)


#################################################################3
###    WEIBULL PHENOMETRICS

#lists to store Weibull phenometrics
we.onset<-list()
we.term<-list()

for(li in c(1:length(all.datasets))) {
  phenodata<-all.datasets[[li]]
  ## DATA AGGREGATIONS:
  #1. F2021: ONLY BY LATITUDE (RNDLAT)
  we.onset.1<-we.metric(phenodata,"onset",annual=F,elev.strat=F,minFP, mindoy = minDOY[li]) 
  we.term.1<-we.metric(phenodata,"term",annual=F,elev.strat=F,minFP, maxdoy = maxDOY[li]) #F2020/F2021 metric
  
  #2. LS1: BY LATITUDE AND YEAR
  we.onset.yr<-we.metric(phenodata,"onset",annual=T,elwe.strat=F,minFP, mindoy = minDOY[li]) 
  we.term.yr<-we.metric(phenodata,"term",annual=T,elwe.strat=F,minFP, maxdoy = maxDOY[li]) 
  
  #3. NEW: BY LATITUDE, YEAR, ELEVATION
  we.onset.yr.alt<-we.metric(phenodata,"onset",annual=T,elev.strat=T,minFP, mindoy = minDOY[li]) 
  we.term.yr.alt<-we.metric(phenodata,"term",annual=T,elev.strat=T,minFP, maxdoy = maxDOY[li]) 
  
  #These collect extreme value metrics for 2 data curations, 3 aggregations
  we.onset[[li]]<-list(we.onset.1, we.onset.yr, we.onset.yr.alt)
  we.term[[li]]<-list(we.term.1, we.term.yr, we.term.yr.alt)
  
  filename=paste("data/phenometrics/we.metrics.RData", sep="")
  save(we.onset, we.term, file=filename)
}
filename=paste("data/phenometrics/all.we.metrics.RData", sep="")
save(we.onset, we.term, file=filename)
