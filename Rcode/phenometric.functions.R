#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Updated 2021-06
#HERE: Functions for phenometric estimation

#load  packages
library(tidyverse) #to be tidy
library(phenesse) #to calculate onset & termination with Weibull distrib.


###### FUNCTIONS

##############################
### FUNCTIONS

##   PHENOMETRIC CALCULATION FUNTIONS
##   No error handling
##   All functions assume a dataframe input of occurrence records with columns:
##   name (Species), region (Continent), name2 (Species-Region dataset), 
##   rndLat (latitudinal band), year, alt, doy (day of year of occurrence)

#Calculation for extreme value  phenometrics per latitude
ev.metric = function(df1,               #dataframe of occurrences
                     pheno="onset",     #which phenometric, default to onset
                     annual=T,          #estimate metrics by year or not
                     elev.strat=F,      #estimate metrics by elevation range or not
                     n.occ.threshold=5, #minimum number of occurrences for calculating phenometrics
                     mindoy=60,         #minimum DOY allowed for onset
                     maxdoy=330,        #maximum DOY allowed for termination
                     minFP=7) {         #minimum # days to calculate effort metric
  
  if(annual==F & elev.strat==T) {message("Elevation stratification is only available for annual metrics. Annual metrics will be calculated.")
    annual<-T}
  
  #for each dataset & latitudinal band, filter to extreme value records & calculate effort metrics (n.occ and Sample Effort)
  if(annual){
    if(elev.strat) {
      #Group & estimate effort by rndLat*year*elevation range
      ev.m<-df1 %>% 
        dplyr::select(name, region, name2, rndLat, year, alt, elev.yr, doy) %>%
        group_by(name, region, name2, rndLat, year, elev.yr, doy) %>%
        add_tally(name="n.m") %>%
        mutate(n.14.onset=ifelse(doy-min(doy, na.rm=T)<15,1,0), n.14.offset=ifelse(max(doy, na.rm=T)-doy<15,1,0)) %>%
        group_by(name,region,name2, rndLat,year, elev.yr) %>%
        add_tally(name="n.occ") %>%
        filter(n.occ>=n.occ.threshold) %>%
        dplyr::mutate(ref=paste(pheno,".ev",sep=""),unit=paste(name2,rndLat,year,elev.yr,sep=""),
               metric=doy, n.doy=length(unique(doy)),
               n.max=max(n.m, na.rm=T), n.14.on=sum(n.14.onset),n.14.term=sum(n.14.offset))
               #sampleEffort=(n.doy/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T))))
    } else {
      #Group & estimate effort by rndLat*year
     ev.m<-df1 %>% 
       dplyr::select(name, region, name2, rndLat, year, alt, doy) %>%
       group_by(name, region, name2, rndLat, year, doy) %>%
       add_tally(name="n.m") %>%
       mutate(n.14.onset=ifelse(doy-min(doy, na.rm=T)<15,1,0), n.14.offset=ifelse(max(doy, na.rm=T)-doy<15,1,0)) %>%
       group_by(name,region,name2, rndLat,year) %>%
        add_tally(name="n.occ") %>%
        filter(n.occ>=n.occ.threshold) %>%
        dplyr::mutate(ref=paste(pheno,".ev",sep=""),unit=paste(name2,rndLat,year,sep=""),
               metric=doy, n.doy=length(unique(doy)),
               n.max=max(n.m, na.rm=T),  n.14.on=sum(n.14.onset),n.14.term=sum(n.14.offset))
              #sampleEffort=(n.doy/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T))))
    }
  } else {
    #Group & estimate effort by rndLat
    ev.m<-df1 %>% 
      dplyr::select(name, region, name2, rndLat, year, alt, doy) %>%
      group_by(name, region, name2, rndLat, year, doy) %>%
      add_tally(name="n.m") %>%
      mutate(n.14.onset=ifelse(doy-min(doy, na.rm=T)<15,1,0), n.14.offset=ifelse(max(doy, na.rm=T)-doy<15,1,0)) %>%
      group_by(name,region,name2, rndLat)  %>%
      add_tally(name="n.occ") %>%
      filter(n.occ>=n.occ.threshold) %>%
      dplyr::mutate(ref=paste(pheno,".ev",sep=""), unit=paste(name2,rndLat,sep=""),
             metric=doy, n.doy=length(unique(doy)), 
             n.max=max(n.m, na.rm=T), 
             n.14.on=sum(n.14.onset),n.14.term=sum(n.14.offset))
             #sampleEffort=(n.doy/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T))))
  }
  if(pheno=="onset") {
    ev.m<-ev.m %>% 
      group_by(unit) %>% filter(metric==min(metric, na.rm=T)) %>%
      dplyr::mutate(metric=ifelse((metric>=mindoy),metric, mindoy), n.14=n.14.on) 
  } else if(pheno=="term") {
    ev.m<-ev.m %>% 
      group_by(unit) %>% filter(metric==max(metric, na.rm=T)) %>%
      dplyr::mutate(metric=ifelse((metric<=maxdoy),metric, maxdoy), n.14=n.14.term) 
  }
  return(ev.m)
}

#Calculation for weibull phenometrics
we.est = function(doys,     #vector of occurrence DOYs used in estimation
                  IfTerm) { #IF T, termination gives upper limit. F -> onset
  if(IfTerm) {  
    try(suppressWarnings(weib_percentile(doys, 0.99)), silent=T)
  } else {
    try(suppressWarnings(weib_percentile(doys, 0.01)), silent=T)
  }
}

we.metric = function(df1,               #dataframe of occurrences
                     pheno="onset",     #which phenometric, defaults to onset
                     annual=T,          #estimate metrics by year or not
                     elev.strat=F,      #estimate metrics by elevation range or not
                     mindoy=60,         #minimum DOY allowed for onset
                     maxdoy=330,        #maximum DOY allowed for termination
                     n.occ.threshold=5, #minimum number of occurrences for calculating phenometrics
                     minFP=7) {         #minimum # days to calculate effort
  
  if(annual==F & elev.strat==T) {
    message("Elevation stratification is only available for annual metrics. Annual metrics will be calculated.")
    annual<-T
  }
  if(annual) {
    if(elev.strat) {
      #Group & estimate effort by rndLat*year*elevation range
      weib<-df1 %>%
        dplyr::select(name, region, name2, rndLat, year, alt, doy, elev.yr) %>%
        group_by(name,region,name2, rndLat, year, elev.yr) %>%
        dplyr::mutate(n.occ=length(unique(doy)),alt=min(alt, na.rm=T)) %>%
        dplyr::mutate(sampleEffort=(n.occ/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T)))) %>%
        filter(n.occ>=n.occ.threshold)
      if(nrow(weib)>=n.occ.threshold) {
        #Calculate phenometrics!
        weib<- weib%>% 
          group_by(name,region,name2, rndLat, year, elev.yr) %>%
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
        #Don't calculate metrics if data are insufficient
        weib<-df1 %>%
          group_by(name,region,name2, rndLat, year) %>%
          summarize(metric=NA)
      } 
    } else if(elev.strat==F) {
      #Group & estimate effort by rndLat*year
      weib<-df1 %>%
        dplyr::select(name, region, name2, rndLat, year, alt, doy) %>%
        group_by(name,region,name2, rndLat, year) %>%
        dplyr::mutate(n.occ=length(unique(doy)),alt=min(alt, na.rm=T)) %>%
        dplyr::mutate(sampleEffort=(n.occ/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T)))) %>%
        filter(n.occ>=n.occ.threshold)
      if(nrow(weib)>=n.occ.threshold) {
        #Calculate phenometrics!
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
        #Don't calculate metrics if data are insufficient
        weib<-df1 %>%
          group_by(name,region,name2, rndLat, year) %>%
          summarize(metric=NA)
      } 
    } #end if elev.strat==F
  } else {
    #Group & estimate effort by rndLat
    weib<-df1 %>%
      dplyr::select(name, region, name2, rndLat, year, alt, doy) %>%
      group_by(name,region,name2, rndLat) %>%
      dplyr::mutate(n.occ=length(unique(doy)), alt=min(alt, na.rm=T)) %>%
      dplyr::mutate(sampleEffort=(n.occ/(max(max(doy, na.rm=T)-min(doy, na.rm=T)+1,minFP, na.rm=T)))) %>%
      filter(n.occ>=n.occ.threshold) 
    if(nrow(weib)>=n.occ.threshold) {
      #Calculate phenometrics!
      weib<-weib%>%
        summarize(ref=paste(pheno,".we",sep=""),
                n.occ=mean(n.occ,na.rm=T),
                sampleEffort=mean(sampleEffort, na.rm=T), 
                metric=if(pheno=="onset") {
                  max(mindoy,we.est(doy,F)[[1]],na.rm=T)
                  } else if(pheno=="term") {
                  min(maxdoy,we.est(doy,T)[[1]],na.rm=T)}
                  else {NA},
                alt=min(alt, na.rm=T)) 
    } else {
      #Don't calculate metrics if data are insufficient
      weib<-df1 %>%
        group_by(name,region,name2, rndLat, year) %>%
        summarize(metric=NA)
    }
  }
  return(weib)
}
  