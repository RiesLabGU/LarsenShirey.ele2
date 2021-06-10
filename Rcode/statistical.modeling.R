#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Updated 2021-06
#HERE: Statistical models


#load  packages
library(tidyverse) #to be tidy
library(lme4) #to run mixed-effects models
library(lmerTest) #to compare mixed effects models
library(MuMIn) #to calculate pseudo r2 values for mixed effects models
library(ggplot2) #to make plots
library(sjPlot) #to plot mixed-effect model results
library(viridis) #colorblind friendly color palette for plots


########################################
##  FUNCTIONS

#MIXED EFFECTS MODEL: ANALYSIS & OUTPUT
model.lme = function(df1, effort=T) {
  df2<-df1 %>% group_by(name2) %>% mutate(nlat=length(unique(rndLat))) %>% 
    #filter(nlat>1) %>%
    group_by(name2, rndLat) %>% mutate(nyear=length(unique(year)), alt2=alt/100)
  if(effort) {
    fullmodel<-lmer(metric~1 + rndLat:name2 +  alt2 +log(nyear) + log(n.occ)*sampleEffort + (1|year)  + (1|name2), data=df2) 
  } else {
    
    fullmodel<-lmer(metric~1 + rndLat:name2  + alt2 + (1|year)  + (1|name2), data=df2)
  }
  step_fm <- step(fullmodel)
  finalmodel <- get_model(step_fm)
  return(finalmodel)
  
}


#OUTPUTS
model.stats<- function(finalmodel, nsp) {
  species.resp<-as.data.frame(summary(finalmodel)$coefficients)
  if(length(c(grep("rndLat:name2",row.names(summary(finalmodel)$coefficients))))> 1) {
      #Calculate # species responses if there is a lat*species interaction
      species.resp<-species.resp[c(grep("rndLat:name2",row.names(summary(finalmodel)$coefficients))),] %>%
        mutate(sign1=ifelse(`Pr(>|t|)`<0.05,sign(Estimate),0))
      sign1<-species.resp$sign1
  } else if(length(c(grep("rndLat",row.names(summary(finalmodel)$coefficients))))==1) {
    #With no species interaction, assign general latitudinal response to all species
    sign.1<-species.resp[c(grep("rndLat",row.names(summary(finalmodel)$coefficients))),]
    species.resp<-species.resp[c(grep("rndLat",row.names(summary(finalmodel)$coefficients))),] 
    sign1<-c(rep(ifelse(sign1$`Pr(>|t|)`<0.05,sign(sign1$Estimate),0),nsp))
  } else {
    #rndLat was not in the final model, so assign 0 to all species
    sign1<-c(rep(0,nsp))
  }
  return(c(length(which(sign1==-1)),length(which(sign1==0)),length(which(sign1==1))))

  }


#
model.compar<- function(data) {
  
  nsp<-length(unique(data$name2))
  #run model without effort
  ev.full.f<-model.lme(data, effort=F)
  ev.final.f<-get_model(step(ev.full.f))
  (ev.stats.f<-model.stats(ev.final.f, nsp))
  (r.squaredGLMM(ev.final.f))[1]
  summary(ev.final.f)$call
  no.effort<-c(paste(formula(summary(ev.final.f)$call)[3]),0, round(extractAIC(ev.final.f)),round(r.squaredGLMM(ev.final.f)[1],2),
           round(r.squaredGLMM(ev.final.f)[2],2),model.stats(ev.final.f))
  
  #run model with effort
  ev.full.t<-model.lme(data, effort=T)
  ev.final.t<-get_model(step(ev.full.t))
  (ev.stats.t<-model.stats(ev.final.t, nsp))
  (r.squaredGLMM(ev.final.t))
  summary(ev.final.t)$call
  effort<-c(paste(formula(summary(ev.final.t)$call)[3]),1, round(extractAIC(ev.final.t)),round(r.squaredGLMM(ev.final.t)[1],2),
           round(r.squaredGLMM(ev.final.t)[2],2),model.stats(ev.final.t))
  names(effort)<-names(no.effort)<-c("Model","Effort","EDF","AIC","r2m","r2c","Nsp.NEG","Nsp.0","Nsp.POS")
  return(rbind(no.effort,effort))
}

########################################
## PARAMETERS

#for visualizations
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = 'A')
v2<-v_colors[c(2,4,6,8)]


########################################
## Load phenometrics


## Extreme Value phenometrics
load("data/phenometrics/all.ev.metrics.RData")

#ev.onset and ev.term each have 2 lists of phenometric datasets, one for each data curation
#each list has 3 datasets, representing different aggregations 

##
aggregs<-c("Lat","Lat*yr","Lat*yr*elev")
model.output<-list()
for(i in 1:length(ev.onset)) {
  pheno<-"onset"
  metric<-"EV"
  curate<-i
  for(j in 1:length(ev.onset[[i]])) {
    aggreg<-  aggregs[j]
    model.input<-ev.onset[[i]][[j]]
    effort<-T
    nsp<-length(unique(model.input$name2))
    ntot<-nrow(model.input)
    params<-as.data.frame(rbind(c(pheno, metric, curate, aggreg, nsp, ntot),c(pheno, metric, curate, aggreg, nsp, ntot)))    
    names(params)<-c("Metric","Estimator","Curation","Aggregation","N.SP","N.REC")
    model.output[[(i-1)*j+j]]<-as.data.frame(cbind(model.compar(model.input), params))
}
}
  
  



