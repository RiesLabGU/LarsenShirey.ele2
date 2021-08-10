#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Last Update 2021-08
#HERE: Data curation functions

#load packages
library(tidyverse) #to be tidy
data.curation

split<-median(f1.effort$rndLat)
summary(lm(nocc~rndLat*name2+year, data=f.effort))
x1b<-summary(lm(nocc~name2+rndLat:name2+year:name2, data=filter(f.effort,rndLat>split)))
x1a<-summary(lm(nocc~name2+rndLat:name2+year:name2, data=filter(f.effort,rndLat<split)))
x2b<-summary(lm(nocc~name2+rndLat:name2, data=filter(f.effort,rndLat>split)))
x2a<-summary(lm(nocc~name2+rndLat:name2, data=filter(f.effort,rndLat<split)))
#x3b<-summary(lm(nocc~name2+rndLat+year, data=filter(f.effort,rndLat>split)))
#x3a<-summary(lm(nocc~name2+rndLat+year, data=filter(f.effort,rndLat<split)))
x4b<-summary(lm(nocc~name2+year:name2, data=filter(f.effort,rndLat>split)))
x4a<-summary(lm(nocc~name2+year:name2, data=filter(f.effort,rndLat<split)))

x1b$adj.r.squared #all
x2b$adj.r.squared #name * lat 
x3b$adj.r.squared #name & year


sig<-x1a$coefficients[x1a$coefficients[,4]<0.05,]
nsplat<-length(grep(":rndLat",rownames(sig)))
nspyr<-length(grep(":year",rownames(sig)))

sig2<-x1b$coefficients[x1b$coefficients[,4]<0.05,]
nsplat2<-length(grep(":rndLat",rownames(sig2)))
nspyr2<-length(grep(":year",rownames(sig2)))

#partial r2 for year
x1a$adj.r.squared - x2a$adj.r.squared
x1b$adj.r.squared - x2b$adj.r.squared

#partial r2 for lat
x1a$adj.r.squared - x4a$adj.r.squared
x1b$adj.r.squared - x4b$adj.r.squared

#partial r2 for genral lat
x1a$adj.r.squared - x4a$adj.r.squared
x1b$adj.r.squared - x4b$adj.r.squared



f1.effort<-comb.elev(add.elev(F1.data)) 

f.effort <- f1.effort %>%
  group_by(name2,rndLat, year,doy) %>%
  add_tally(name='ndoy') %>%
  group_by(name2,rndLat, year) %>%
  mutate(on14=ifelse((doy-min(doy))<15,1,0),off14=ifelse((max(doy)-doy)<15,1,0)) %>%
  summarize(nocc=n(),logn=log(nocc),n14.onset=sum(on14), n14.term=sum(off14), nmax=max(ndoy))

library("corrplot")
cor(f.effort[,c(2:8)])
corrplot(cor(f.effort[,c(2:8)]))



