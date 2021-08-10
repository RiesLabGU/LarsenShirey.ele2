#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Last Update 2021-08
#HERE: Data curation functions

#load packages
library(tidyverse) #to be tidy

###########################
## FUNCTIONS

#Function for removing y from x
`%notin%` = function(x,y) !(x %in% y)

#Function to return summary statistics in format "median (min-max)"
median.range<-function(x) {
  a<-median(x, na.rm=T)
  b<-min(x, na.rm=T)
  c<-max(x, na.rm=T)
  d<-paste(a," (",b,"-",c,")", sep="")
  return(d)
}

#Function to calculate summary statistics for curated datasets
calc.sum<-function(df, nthresh=5) {
  perlatyr<-df %>%
    group_by(name, region, rndLat, year) %>%
    summarize(n.per.AO=n()) %>%
    filter(n.per.AO>=nthresh)
  perlat<-df %>%
    group_by(name, region, rndLat) %>%
    summarize(n.per.lat=n()) 
  nyrs<-df %>%
    group_by(name, region) %>%
    summarize(nyears=length(unique(year)), nlats=length(unique(rndLat)), noccs=n())
  nmetrics<-df %>%
    group_by(name, region, rndLat) %>%
    summarize(n.per.lat=n()) %>%
    filter(n.per.lat>=nthresh)
  nm2<-nmetrics %>% group_by(name, region) %>% tally()
  n.occ.tot<-nrow(df)
  nsp<-length(unique(df$name))
  nsp.reg<-length(unique(paste(df$name,df$region)))
  
  out.fin<-c(nsp, nsp.reg,median.range(nyrs$nlats), median.range(nyrs$nyears),
             median.range(nyrs$noccs),median.range(perlat$n.per.lat),
             median.range(nmetrics$n.per.lat),median.range(nm2$n),
             median.range(perlatyr$n.per.AO))
 names(out.fin)<-c("nsp","nSR","nlat","nyear","nocc","nperlat","metperlat",
                   "nmet","nperlatyrmet")            
             
#             median.range(out3$n.per.dset), 
#             median.range(out2$n.per.lat), median.range(out1$n.per.AO),
#             1, median.range(out2$n.on.we),median.range(out3$nlat),
#             median.range(n.onset$n), median.range(out3$nmetric))
#  names(out.fin)<-c("n.occ","nsp","nsp-reg","n.perSR","n.perLat","n.perAO","1",
#                    "n.latyr.onset","nlat","n.lat.EVonset","n.occ.perSR")
  return(out.fin)
}

#Function to add elevation ranges to data
add.elev<-function(df) {
  df.2<-df %>%
  group_by(name, region, rndLat, year) %>%
  mutate(elev.yr=ifelse(alt>quantile(alt,0.9),10, 
                        ifelse(alt>quantile(alt, 0.8),9,
                          ifelse(alt>quantile(alt,0.7),8, ifelse(alt>quantile(alt, 0.6),7,
                            ifelse(alt>quantile(alt,0.5),6, ifelse(alt>quantile(alt, 0.4),5,
                              ifelse(alt>quantile(alt,0.3),4, ifelse(alt>quantile(alt, 0.2),3,
                                ifelse(alt>quantile(alt,0.1),2,1)))))))))) %>%
  group_by(name, region, rndLat, elev.yr) %>%
  add_tally(name="ny1")

  return(df.2)
}

  #Not super efficient but it does the job ... regroup 
comb.elev<-function(df.2) {
  for(xi in rev(c(1:10))) {
    df.2$elev.yr[which(df.2$elev.yr==xi & df.2$ny1<10)]<-df.2$elev.yr[which(df.2$elev.yr==xi & df.2$ny1<10)]-1
    df.2<-df.2 %>%
      select(row.index:elev.yr) %>%
      group_by(name, region, rndLat, year, elev.yr) %>%
      add_tally(name="ny1")
  }
  return(df.2)
}
