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
  if(effort) {
    df2<-df1 %>% group_by(name2) %>% mutate(nlat=length(unique(rndLat))) %>% 
      filter(nlat>1) %>%
      group_by(name2, rndLat) %>% mutate(nyear=length(unique(year)), alt2=alt/100)
    
    fullmodel<-lmer(metric~1 + rndLat:name2 +  alt2 +log(n.occ) + log(nyear) + sampleEffort + (1|year)  + (1|name2), data=df2) 
  } else {
    fullmodel<-lmer(metric~1 + rndLat:name2  + alt2 + (1|year)  + (1|name2), data=df2)
  }
  step_fm <- step(fullmodel)
  finalmodel <- get_model(step_fm)
  return(finalmodel)
  
}


#OUTPUTS
model.stats<- function(finalmodel) {
  species.resp<-as.data.frame(summary(finalmodel)$coefficients)
  species.resp<-species.resp[c(grep("rndLat",row.names(summary(fullmodel)$coefficients))),] %>%
    mutate(sign1=ifelse(`Pr(>|t|)`<0.05,sign(Estimate),0))
  return(table(species.resp$sign1))
}

########################################
## PARAMETERS

q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = 'A')
v2<-v_colors[c(2,4,6,8)]

