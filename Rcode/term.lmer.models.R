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
library(gridExtra)

########################################
##  FUNCTIONS

#MIXED EFFECTS MODEL: ANALYSIS & OUTPUT
model.lme = function(df1, effort=T) {
  
  df1<-df1 %>%  mutate(alt2=alt/100) %>%
    dplyr::select(metric, name2, rndLat, alt2, n.14, n.max, year)
  
  
  if(effort) {
    fullmodel<-lmer(metric~1 + rndLat:name2 +  alt2  + n.14 + n.max + (1|year)  + (1|name2), data=df1) 
  } else {
    
    fullmodel<-lmer(metric~1 + rndLat:name2  + alt2 + (1|year)  + (1|name2), data=df1)
  }
  try(fullmodel <- get_model(step(fullmodel)))
  return(fullmodel)
  
}


#OUTPUTS
model.coefs<-function(finalmodel) {
  species.resp<-as.data.frame(summary(finalmodel)$coefficients)
  return(species.resp)
}

model.stats.sp<- function(finalmodel, nsp) {
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
  return(sign1)
}

model.stats.sum<-function(sign1) {
  return(c(length(which(sign1==-1)),length(which(sign1==0)),length(which(sign1==1))))
}

#
model.compar1<- function(model1, nsp) {
  model.summary<-c(effort, paste(formula(summary(model1)$call)[3]), round(extractAIC(model1)),round(r.squaredGLMM(model1)[1],2),
                   round(r.squaredGLMM(model1)[2],2),model.stats.sum(model.stats.sp(model1, nsp)))
  return(model.summary)
}


mc.coef<-function(data, eff=F) {
  #run model without effort
  full.f<-model.lme(data, effort=eff)
  final.f<-get_model(step(full.f))
  final.coefs<-summary(final.f)$coefficients
  latrows<-0
  latrows<-c(grep("rndLat:name2",row.names(final.coefs)))
  if(length(latrows)> 1) {
    coef.out<-paste(round(median(final.coefs[latrows,1], na.rm=T),2)," [",
                    round(min(final.coefs[latrows,1], na.rm=T),2),"-",
                    round(max(final.coefs[latrows,1], na.rm=T),2),"]", sep="")
  } else {coef.out<-paste("0 [0,0]")}      
  return(coef.out)
}


########################################
## PARAMETERS & output lists

#for models
aggregs<-c("Lat","Lat*yr","Lat*yr*elev")
ev.term.output<-list()
we.term.output<-list()
ev.sp.term.signs<-list()
we.sp.term.signs<-list()


#for visualizations
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = 'A')
v2<-v_colors[c(2,4,6,8)]



########################################
## Load phenometrics

## Extreme Value phenometrics
load("data/phenometrics/ev.metrics.RData")

## Weibull phenometrics (in process)
load("data/phenometrics/we.metrics.RData")


########################################################
###  RUN EV MIXED EFFECT MODELS
#ev.onset and has 2 lists of phenometric datasets, one for each data curation
#each list has 3 datasets, representing different aggregations 
metric<-"EV"

pheno<-"term"
index<-1
summary.res.ev.term<-list()
sp.resp.ev.term<-list()
ev.term.coefs<-list()
for(i in 1:length(ev.term)) {
  curate<-i
  for(j in 1:length(ev.term[[i]])) {
    aggreg<-  aggregs[j]
    model.input<-ev.term[[i]][[j]] %>% dplyr::mutate(year=year-1850, aggreg=aggreg, curat=curate) %>%
      dplyr::select(name, region, name2, rndLat, year, metric, alt, aggreg, curat, n.14, n.max)
    model.input<-na.omit(model.input)
    if(nrow(model.input)>0) {
      effort<-T
      nsp<-length(unique(model.input$name2[!is.na(model.input$metric)]))
      ntot<-length((model.input$name2[!is.na(model.input$metric)]))
      effort.lme<-model.lme(model.input, effort=T)
      effort.summary<-rbind(c(pheno,metric,curate,aggreg,nsp,ntot,"yes",model.compar1(effort.lme, nsp)))
      effort.sp.resp<-data.frame(spname=sort(unique(model.input$name2[!is.na(model.input$metric)])), metric="ev.off",cur=curate,agg=aggreg,effort="yes", resp=model.stats.sp(effort.lme, nsp))
      #summary.res.ev.on[[index]]<-effort.summary
      #sp.resp.ev.on[[index]]<-effort.sp.resp
      #index<-index+1    
      
      noeffort.lme<-model.lme(model.input, effort=F)
      colnames(effort.summary)<-c("Metric","Estimator","Curation","Aggregation","N.SP","N.REC","effort","blah","best.model","AIC.df","AIC","r2m","r2c","neg","nsig","pos")
      lme.summary<-rbind(effort.summary, c(pheno,metric,curate,aggreg,nsp,ntot,"no",model.compar1(noeffort.lme, nsp)))
      colnames(lme.summary)<-c("Metric","Estimator","Curation","Aggregation","N.SP","N.REC","effort","blah","best.model","AIC.df","AIC","r2m","r2c","neg","nsig","pos")
      sp.resp<-bind_rows(effort.sp.resp,data.frame(spname=sort(unique(model.input$name2[!is.na(model.input$metric)])), metric="ev.off",cur=curate,agg=aggreg,effort="no", resp=model.stats.sp(noeffort.lme, nsp)))
      summary.res.ev.term[[index]]<-as.data.frame(lme.summary)
      sp.resp.ev.term[[index]]<-sp.resp
      temp1<-as.data.frame(summary(noeffort.lme)$coefficients) %>%
        dplyr::mutate(permut=paste(pheno,metric,curate,aggreg,"e0",sep="."), spname=row.names(summary(noeffort.lme)$coefficients))
      temp2<-as.data.frame(summary(effort.lme)$coefficients) %>%
        dplyr::mutate(permut=paste(pheno,metric,curate,aggreg,"e1",sep="."), spname=row.names(summary(effort.lme)$coefficients))
      ev.term.coefs[[index]]<-bind_rows(temp1,temp2)  
      index<-index+1    
    }
  }
}

summary.res.ev.term[[1]]
ev.term.summary<-bind_rows(summary.res.ev.term)

sp.resp.ev.term[[1]]
ev.term.spresp<-bind_rows(sp.resp.ev.term)

ev.term.sprtot<-ev.term.spresp %>%
  group_by(spname) %>%
  summarize(neg=length(which(resp==-1)), ns=length(which(resp==0)),pos=length(which(resp==1)))

write.csv(ev.term.spresp,file="output/ev_term_sp.results.csv")
write.csv(ev.term.summary,file="output/ev_term_model.results.csv")




########################################################
###  RUN WE MIXED EFFECT MODELS
#we.onset has 2 lists of phenometric datasets, one for each data curation
#each list has 3 datasets, representing different aggregations 
load("data/phenometrics/we.term.RData")
we.term[[1]][[1]]<-we.term[[1]][[1]] %>% dplyr::mutate(year=year.y)
metric<-"WE"

pheno<-"term"
index<-1
summary.res.we.term<-list()
sp.resp.we.term<-list()
we.term.coefs<-list()
#onset.coefs<-list()
for(i in 1:length(we.term)) {
  curate<-i
  for(j in 1:length(we.term[[i]])) {
    aggreg<-  aggregs[j]
    model.input<-we.term[[i]][[j]] %>% dplyr::mutate(year=year-1850) %>%
      dplyr::select(name, region, name2, rndLat, year, metric, alt, aggreg, curat, n.14, n.max)
    model.input<-na.omit(model.input) %>% mutate(alt2=alt/100)
    if(nrow(model.input)>0) {
      effort<-T
      nsp<-length(unique(model.input$name2[!is.na(model.input$metric)]))
      ntot<-length((model.input$name2[!is.na(model.input$metric)]))
      df1<-model.input
      effort.lme<-model.lme(model.input, effort=T)
      effort.summary<-rbind(c(pheno,metric,curate,aggreg,nsp,ntot,"yes",model.compar1(effort.lme, nsp)))
      effort.sp.resp<-data.frame(spname=sort(unique(model.input$name2[!is.na(model.input$metric)])), metric="we.off",cur=curate,agg=aggreg,effort="yes", resp=model.stats.sp(effort.lme, nsp))
      #summary.res.we.on[[index]]<-effort.summary
      #sp.resp.we.on[[index]]<-effort.sp.resp
      #index<-index+1    
      
      noeffort.lme<-model.lme(model.input, effort=F)
      colnames(effort.summary)<-c("Metric","Estimator","Curation","Aggregation","N.SP","N.REC","effort","blah","best.model","AIC.df","AIC","r2m","r2c","neg","nsig","pos")
      lme.summary<-rbind(effort.summary, c(pheno,metric,curate,aggreg,nsp,ntot,"no",model.compar1(noeffort.lme, nsp)))
      colnames(lme.summary)<-c("Metric","Estimator","Curation","Aggregation","N.SP","N.REC","effort","blah","best.model","AIC.df","AIC","r2m","r2c","neg","nsig","pos")
      sp.resp<-bind_rows(effort.sp.resp,data.frame(spname=sort(unique(model.input$name2[!is.na(model.input$metric)])), metric="we.off",cur=curate,agg=aggreg,effort="no", resp=model.stats.sp(noeffort.lme, nsp)))
      summary.res.we.term[[index]]<-as.data.frame(lme.summary)
      sp.resp.we.term[[index]]<-sp.resp
      temp1<-as.data.frame(summary(noeffort.lme)$coefficients) %>%
        dplyr::mutate(permut=paste(pheno,metric,curate,aggreg,"e0",sep="."), spname=row.names(summary(noeffort.lme)$coefficients))
      temp2<-as.data.frame(summary(effort.lme)$coefficients) %>%
        dplyr::mutate(permut=paste(pheno,metric,curate,aggreg,"e1",sep="."), spname=row.names(summary(effort.lme)$coefficients))
      we.term.coefs[[index]]<-bind_rows(temp1,temp2)  
      print(index)
      index<-index+1    
    }
  }
}

summary.res.we.term[[1]]
we.term.summary<-bind_rows(summary.res.we.term)

sp.resp.we.term[[1]]
we.term.spresp<-bind_rows(sp.resp.we.term)

we.term.sprtot<-we.term.spresp %>%
  group_by(spname) %>%
  summarize(neg=length(which(resp==-1)), ns=length(which(resp==0)),pos=length(which(resp==1)))

write.csv(we.term.spresp,file="output/we_term_results.csv")
save(we.term.summary, we.term.spresp,we.term.sprtot,we.term.coefs,ev.term.summary, ev.term.spresp,ev.term.sprtot,ev.term.coefs,file="output/term_model_results.RData")
load("output/term_model_results.RData")
write.csv(we.term.summary,file="output/we_term_models.csv")


term.coefs<-bind_rows(bind_rows(ev.term.coefs), bind_rows(we.term.coefs))
term_coefs<-(term.coefs) %>%
  select(permut, spname, Estimate, `Pr(>|t|)`) %>%
  dplyr::mutate(sig=ifelse(`Pr(>|t|)`<0.05,1,0),sign=ifelse(`Pr(>|t|)`<0.05,sign(Estimate),0),
                param=ifelse(str_sub(permut,start=6,end=6)=="E","EV","WE"))

term_coefs<- term_coefs[grep("rndLat:name2", term_coefs$spname),]
n.am.rows<-grep(".N. America", term_coefs$spname)
term_coefs<- term_coefs  %>%
  dplyr::mutate(species=str_sub(spname, start = 13, end = str_locate(spname, "\\.")[,1]-1),
                region=str_sub(spname, start = str_locate(spname, "\\.")[,1]+1),
                aggreg=ifelse(str_length(permut)==24,"LYE",ifelse(str_length(permut)==19,"LY","L")))


shapevec<-c("L"=15,"LY"=17,"LYE"=19)
### TERM RESULTS
(f1.na.evt<-ggplot(data=filter(term_coefs, region=="N. America", param=="EV"), aes(x=species, y=Estimate)) + geom_point(aes(shape=aggreg,color=as.factor(sign)) ) + 
    scale_shape_manual(values=shapevec,guide="none") + ylim(-4,2.5) + 
    scale_color_manual(values=c(c("blue","darkgray", "darkgreen")), guide="none") + 
    scale_x_discrete(limits=rev) +   
    theme_classic() + labs(x="Species (N. America)", title="Termination (Extreme) ~ Latitude") + 
    theme(text = element_text(size=8), axis.title.x = element_blank(), axis.text.x=element_blank()) + 
    geom_hline(yintercept=0, linetype="dashed") + 
    coord_flip()  )

(f1.eur.evt<-ggplot(data=filter(term_coefs, region=="Europe", param=="EV"), aes(x=species, y=Estimate)) + geom_point(aes(shape=aggreg,color=as.factor(sign)) ) + 
    scale_shape_manual(values=shapevec,guide="none") + ylim(-4,2.5) + 
    scale_color_manual(values=c(c("blue","darkgray", "darkgreen")), guide="none") + 
    scale_x_discrete(limits=rev) +   
    theme_classic() + labs(x="Species (Europe)", y="Coefficient estimate") +  
    theme(text = element_text(size=8)) + 
    geom_hline(yintercept=0, linetype="dashed") +
    coord_flip()  )

(f1.na.wet<-ggplot(data=filter(term_coefs, region=="N. America", param=="WE"), aes(x=species, y=Estimate)) + 
    geom_point(aes(shape=aggreg,color=as.factor(sign))) + 
    scale_shape_manual(values=shapevec,guide="none") + ylim(-4,2.5) + 
    scale_color_manual(values=c("blue","darkgray", "darkgreen"), name="Response sign", labels=c("Negative","Non-significant","Positive")) + 
    scale_x_discrete(limits=rev) +   
    theme_classic() + labs(title="Termination (Weibull) ~ Latitude") + 
    geom_hline(yintercept=0, linetype="dashed")  +
    theme(text = element_text(size=12), axis.title = element_blank(),axis.text.y = element_blank(),  axis.text.x=element_blank(), axis.ticks.y = element_blank(), legend.position="bottom") + 
    coord_flip()  )

(f1.eur.wet<-ggplot(data=filter(term_coefs, region=="Europe", param=="WE"), aes(x=species, y=Estimate)) + 
    geom_point(aes(shape=aggreg,color=as.factor(sign)) ) + 
    scale_shape_manual(values=shapevec, name="Aggregation", labels=c("Latitude","Lat*Year","Lat*Year*Elevation")) + 
    ylim(-4,2.5) +
    scale_color_manual(values=c(c("blue","darkgray", "darkgreen")), guide="none") + 
    scale_x_discrete(limits=rev) +   
    theme_classic() + labs( y="Coefficient estimate") +  
    theme(text = element_text(size=12),axis.title.y=element_blank(),  axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position="bottom") + 
    geom_hline(yintercept=0, linetype="dashed")  + labs(x="") + 
    coord_flip()  )

library(ggpubr)
leg1 <- get_legend(f1.na.wet) 
L1<-as_ggplot(leg1)
(f1.na.wet<-f1.na.wet + theme(legend.position = "none", text = element_text(size=8)))


leg2 <- get_legend(f1.eur.wet)
L2<-as_ggplot(leg2)
(f1.eur.wet<-f1.eur.wet + theme(legend.position = "none",text = element_text(size=8)))


library(gridExtra)
fig.termres.s1<-grid.arrange(f1.na.evt,f1.na.wet,f1.eur.evt,f1.eur.wet,L1, L2, 
                              layout_matrix = rbind(c(1, 2), c(3, 4),c(5,5),c(6,6)),
                              widths=c(4.4,3), heights=c(10,8,1,1))
ggsave(filename="output/figs/term.res.summary2.png",fig.termres.s1, width = 6, height = 7)




#### OTHER COEFFICIENTS
term_coefs<-(term.coefs) %>%
  select(permut, spname, Estimate, `Pr(>|t|)`) %>%
  dplyr::mutate(sig=ifelse(`Pr(>|t|)`<0.05,1,0),sign=ifelse(`Pr(>|t|)`<0.05,sign(Estimate),0),
                param=ifelse(str_sub(permut,start=6,end=6)=="E","EV","WE"),
                aggreg=ifelse(str_length(permut)==24,"LYE",ifelse(str_length(permut)==19,"LY","L")))

x1<-c(grep("alt2", term_coefs$spname),grep("n.max", term_coefs$spname),grep("n.14", term_coefs$spname))

other_term_coefs<- term_coefs[x1,]

(fig.texp.coef<-ggplot(data=other_term_coefs, aes(y=Estimate, x=spname)) + geom_boxplot(color="gray") + 
    geom_jitter(width=0.1, height=0,aes(shape=aggreg,color=as.factor(sign)),size=2) + theme_classic() + 
    scale_color_manual(values=c("blue","darkgray", "darkgreen"), name="Termination Response") + 
    scale_shape_manual(values=shapevec, name="Aggregation", labels=c("Lat","Lat*Year","Lat*Year*Elev")) + 
    labs(x="Explanatory variables", y="Coefficient estimate") + 
    scale_x_discrete(labels=c("Elevation","n.14 effort","n.max effort")) + 
    geom_hline(yintercept=0, linetype="dashed") )

ggsave(fig.texp.coef, file="output/figs/fig.supp1.term.expvar.png")



