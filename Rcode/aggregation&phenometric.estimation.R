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


#Weibull phenometric estimation takes time, so these variables allow us to turn the code on & off and select which datasets to estimate weibull phenometrics for
weibull.calc<-T
#select which datasets to estimate weibull phenometrics for
weibull.dc<-2

#lists to store Weibull phenometrics
ev.onset<-list()
ev.term<-list()





```
<h2>Calculate annual phenometrics using weibull distribution, package phenesse</h2>
<p>
Warning: Weibull metric estimation with phenesse takes time. Therefore we have a boolean variable to turn this code on/off. </p>
```{r weibull metrics by year}
weibull.calc<-F  
if(weibull.calc==T) {

  #Loop through datasets, calculate phenometrics
  weibull.datasets<-datasets[[weibull.dc]]
  weibull.phenodata<-all.datasets[[weibull.dc]]
  
for(li in 1:length(weibull.datasets)) {

  dataset<-weibull.datasets[li]
  phenodata<-weibull.phenodata[[li]]
 
  WE.YR.onset<-NULL
  WE.YR.term<-NULL
  
  #PHENOMETRIC CALCULATION: Loop through species, calculate Onset BY YEAR
  for(ni in 1:length(sort(unique(phenodata$name2)))) {
    sp.ni<-sort(unique(phenodata$name2))[ni]
    WE.YR.ni<-NULL
    WE.YR.ni<-we.metric(filter(phenodata, name2==sp.ni), "onset", annual=T, minDOY, maxDOY, minFP) ## THIS TAKES TIME
    if(is.null(WE.YR.onset)) {
      WE.YR.onset<-WE.YR.ni
    } else { 
      WE.YR.onset<-rbind(WE.YR.onset,WE.YR.ni)
    }
  }

  #PHENOMETRIC CALCULATION: Loop through species, calculate Termination BY YEAR
  for(ni in 1:length(sort(unique(phenodata$name2)))) {
    sp.ni<-sort(unique(phenodata$name2))[ni]
    WE.YR.ni<-NULL
    WE.YR.ni<-we.metric(filter(phenodata, name2==sp.ni), "term", annual=T, minDOY, maxDOY, minFP) ## THIS TAKES TIME
    if(is.null(WE.YR.term)) {
      WE.YR.term<-WE.YR.ni
    } else { 
      WE.YR.term<-rbind(WE.YR.term,WE.YR.ni)
    }
  }
  filename=paste("we.yr.",dataset,"metrics.RData")
  save(WE.YR.onset, WE.YR.term, file=filename)
}
}


```



```{r weibull metrics by decade}
  
if(weibull.calc==F) {

  #Loop through datasets, calculate phenometrics
  weibull.datasets<-datasets[[weibull.dc]]
  weibull.phenodata<-all.datasets[[weibull.dc]]
  
for(li in 1:length(weibull.datasets)) {

  dataset<-weibull.datasets[li]
  phenodata<-weibull.phenodata[[li]]
  pheno.decade<-phenodata %>%
    mutate(year=floor(year/nyr.agg)*nyr.agg)
 
  WE.dec.onset<-NULL
  WE.dec.term<-NULL
  
  #PHENOMETRIC CALCULATION: Loop through species, calculate Onset BY DECADE
  for(ni in 1:length(sort(unique(phenodata$name2)))) {
    sp.ni<-sort(unique(phenodata$name2))[ni]
    WE.dec.ni<-NULL
    WE.dec.ni<-we.metric(filter(phenodata, name2==sp.ni), "onset", annual=T, minDOY, maxDOY, minFP) ## THIS TAKES TIME
    if(is.null(WE.dec.onset)) {
      WE.dec.onset<-WE.dec.ni
    } else { 
      WE.dec.onset<-rbind(WE.dec.onset,WE.dec.ni)
    }
  }

  
  #PHENOMETRIC CALCULATION: Loop through species, calculate Termination BY DECADE
  for(ni in 1:length(sort(unique(phenodata$name2)))) {
    sp.ni<-sort(unique(phenodata$name2))[ni]
    WE.dec.ni<-NULL
    WE.dec.ni<-we.metric(filter(phenodata, name2==sp.ni), "term", annual=T, minDOY, maxDOY, minFP) ## THIS TAKES TIME
    if(is.null(WE.dec.term)) {
      WE.dec.term<-WE.dec.ni
    } else { 
      WE.dec.term<-rbind(WE.dec.term,WE.dec.ni)
    }
  }
}
  filename=paste("we.dec.",dataset,"metrics.RData")
  save(WE.dec.onset, WE.dec.term, file=filename)
} 



```

Analysis

```{r analysis: mixed effect models}

estimate.sets<-c("EV.dec","EV.YR","WE.dec","WE.YR")

load("data/phenometrics_TC2/ev.LS2.datametrics.RData")
ev.yr.onset.modelT<-model.lme(EV.YR.onset, effort=T)
(ev.yr.onset.statsT<-model.stats(ev.yr.onset.modelT))
(r.squaredGLMM(ev.yr.onset.modelT))
summary(ev.yr.onset.modelT)$call
(ev.yr.onset.statsT<-model.stats(ev.yr.onset.modelT))
(r.squaredGLMM(ev.yr.onset.modelT))

ev.yr.onset.modelF<-model.lme(EV.YR.onset, effort=F) 
(ev.yr.onset.statsF<-model.stats(ev.yr.onset.modelF))
(r.squaredGLMM(ev.yr.onset.modelF))


## WEIBULL
load("we.F2.datametrics.RData")
we.yr.term<-WE.term[[1]]
for(i in 2:87) {
  we.yr.term<-rbind(we.yr.term, WE.term[[i]])
}


we.yr.term.modelT<-model.lme(we.yr.term, effort=T)
m1<-get_model(step(we.yr.term.modelT))
(we.yrterm.statsT<-model.stats(we.yr.term.modelT))
(r.squaredGLMM(we.yr.term.modelT))
summary(we.yr.term.modelT)$call

we.yr.term.modelF<-model.lme(we.yr.term, effort=F) 
(we.yr.term.statsF<-model.stats(we.yr.term.modelF))
(r.squaredGLMM(we.yr.term.modelF))



```
