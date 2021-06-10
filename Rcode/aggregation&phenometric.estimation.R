#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Updated 2021-06
#HERE: Data aggregation and phenometric estimation


#load  packages
library(tidyverse) #to be tidy
library(lme4) #to run mixed-effects models
library(ggplot2) #to make plots
library(phenesse) #to calculate onset & termination with Weibull distrib.
library(sjPlot) #to plot mixed-effect model results
library(viridis) #colorblind friendly color palette for plots
library(MuMIn) #calculate pseudo r2 values for mixed effects models
library(r2glmm)
library(lmerTest)

##############################
### FUNCTIONS
# %notin% function to remove matches
`%notin%` = function(x,y) !(x %in% y)

##   PHENOMETRIC CALCULATION FUNTIONS
##   No error handling
##   All functions assume a dataframe input of occurrence records with columns:
##   name (Species), region (Continent), name2 (Species-Region dataset), 
##   rndLat (latitudinal band), year, alt, doy (day of year of occurrence)

#Calculation for extreme value  phenometrics per latitude
ev.metric = function(df1,            #dataframe of occurrences
                     pheno="onset",  #which phenometric, default to onset
                     annual=T,      #estimate metrics by year or not
                     minFP=7) {      #minimum # days to calculate effort metric
  
  #for each dataset & latitudinal band, filter to extreme value records & calculate effort metrics (n.occ and Sample Effort)
  if(annual==T){
    ev.m<-df1 %>% 
    select(name, region, name2, rndLat, year, alt, doy) %>%
      group_by(name,region,name2, rndLat,year) %>%
    add_tally(name="n.occ") %>%
    mutate(ref=paste(pheno,".ev",sep=""),
           metric=doy, n.doy=length(unique(doy)),
           sampleEffort=(n.doy/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T))))
  } else {
    ev.m<-df1 %>% 
    select(name, region, name2, rndLat, year, alt, doy) %>%
         group_by(name,region,name2, rndLat)  %>%
    add_tally(name="n.occ") %>%
    mutate(ref=paste(pheno,".ev",sep=""),
           metric=doy, n.doy=length(unique(doy)), 
           sampleEffort=(n.doy/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T))))
  }
  if(pheno=="onset") {
    return(filter(ev.m,metric==min(metric, na.rm=T))) }
  else if(pheno=="term") {
    return(filter(ev.m,metric==max(metric, na.rm=T))) }
}

   
#Calculation for weibull phenometrics
we.est = function(doys,     #vector of occurrence DOYs used in estimation
                  IfTerm) { #IF T, termination gives upper limit. F -> onset
  if(IfTerm) {  suppressWarnings(weib_percentile(doys, 0.99)) } else {
    suppressWarnings(weib_percentile(doys, 0.01))
  }
}

we.metric = function(df1,           #dataframe of occurrences
                     pheno="onset", #which phenometric, defaults to onset
                     annual=T,      #estimate metrics by year or not
                     mindoy=60,     #minimum DOY allowed for onset
                     maxdoy=330,    #maximum DOY allowed for termination
                     minFP=7) {     #minimum # days to calculate effort

      if(annual==T){
        weib<-df1 %>%
        select(name, region, name2, rndLat, year, alt, doy) %>%
        group_by(name,region,name2, rndLat, year) %>%
        mutate(n.occ=length(unique(doy)),alt=min(alt, na.rm=T)) %>%
        mutate(sampleEffort=(n.occ/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T)))) %>%
        filter(n.occ>4)
        if(nrow(weib)>4) {
          weib<- weib%>% 
        group_by(name,region,name2, rndLat, year) %>%
          summarize(ref=paste(pheno,".we",sep=""),
                n.occ=mean(n.occ,na.rm=T),
                sampleEffort=mean(sampleEffort, na.rm=T), 
                metric=if(pheno=="onset"){
                  max(mindoy,we.est(doy,F)[[1]],na.rm=T)
                } else if(pheno=="term"){
                  min(maxdoy,we.est(doy,T)[[1]],na.rm=T)
                  }
                else{NA},
                alt=min(alt, na.rm=T)) 
      } else { 
        weib<-df1 %>%
          group_by(name,region,name2, rndLat, year) %>%
        summarize(metric=NA)
      } 
        } else 
          {
        weib<-df1 %>%
        select(name, region, name2, rndLat, year, alt, doy) %>%
        group_by(name,region,name2, rndLat) %>%
        mutate(n.occ=length(unique(doy)), alt=min(alt, na.rm=T)) %>%
        mutate(sampleEffort=(n.occ/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T)))) %>%
        filter(n.occ>4) 
        if(nrow(weib)>4) {
          weib<-weib%>%
        summarize(ref=paste(pheno,".we",sep=""),
                n.occ=mean(n.occ,na.rm=T),
                sampleEffort=mean(sampleEffort, na.rm=T), 
                metric=if(pheno=="onset"){
                  max(mindoy,we.est(doy,F)[[1]],na.rm=T)
                } else if(pheno=="term"){
                  min(maxdoy,we.est(doy,T)[[1]],na.rm=T)}
                else{NA},
+                alt=min(alt, na.rm=T)) 
        } else {
          weib<-df1 %>%
          group_by(name,region,name2, rndLat, year) %>%
            summarize(metric=NA)

      }
  return(weib)
}
    
      }

############################################
## BEGIN DATA AGGREGATION, PHENOMETRIC ESTIMATION

load("data/curations.RData")

all.datasets<-list(cured.data1, cured.data2)
datasets<-c("few.filters","many.filters")

###########################################3
##  SET PARAMETERS

#Weibull phenometric estimation takes time, so these variables allow us to turn the code on & off and select which datasets to estimate weibull phenometrics for
weibull.calc<-T
#select which datasets to estimate weibull phenometrics for
weibull.dc<-2

#Set parameters for onset analyses
minFP<-7   #For calculating an effort metric, we correct for the day-of-year span for which data are available. Here we set a minimum day-of-year span at 7, as the effort metric is misleading when that span is less than 1 week.
minDOY<-60 #Weibull distributions are not limited to reasonable day of year values, so any onset estimates prior to day of year 60 are replaced with day of year 60. 
maxDOY<-330 #Weibull distributions are not limited to reasonable day of year values, so any termination estimates after day of year 330 are replaced with day of year 330. 




#Loop through datasets, calculate phenometrics
ev.onset<-list()
ev.term<-list()

for(li in 1:length(datasets)) {

  dataset<-datasets[li]
  phenodata<-all.datasets[[li]]
  pheno.decade<-phenodata %>%
    mutate(year=floor(year/nyr.agg)*nyr.agg)
 
  #PHENOMETRIC CALCULATION: Onset
  EV.onset<-ev.metric(phenodata,"onset",annual=F,minFP) #F2020/F2021 metric
  EV.YR.onset<-ev.metric(phenodata,"onset",annual=T,minFP) #annual metrics
  EV.dec.onset<-ev.metric(pheno.decade,"onset",annual=T,minFP) #decadal metrics

  #PHENOMETRIC CALCULATION: Termination
  EV.term<-ev.metric(phenodata,"term",annual=F,minFP) #F2020/F2021 metric
  EV.YR.term<-ev.metric(phenodata,"term",annual=T,minFP) #annual metrics
  EV.dec.term<-ev.metric(pheno.decade,"term",annual=T,minFP) #decadal metrics

  ev.onset[[li]]<-list(EV.onset, EV.YR.onset, EV.dec.onset)
  ev.term[[li]]<-list(EV.term, EV.YR.term, EV.dec.term)
  filename=paste("data/phenometrics_TC2/ev.",dataset,"metrics.RData", sep="")
  save(EV.onset, EV.YR.onset, EV.dec.onset,EV.term, EV.YR.term, EV.dec.term, file=filename)
}

save(ev.onset, ev.term, file="data/phenometrics_TC2/ev.all.RData")

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
