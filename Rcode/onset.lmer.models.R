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
ev.onset.output<-list()
we.onset.output<-list()
ev.sp.onset.signs<-list()
we.sp.onset.signs<-list()


#for visualizations
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = 'A')
v2<-v_colors[c(2,4,6,8)]



########################################
## Load phenometrics

## Extreme Value phenometrics
load("data/phenometrics/ev.metrics.RData")

## Weibull phenometrics (in process)
load("data/phenometrics/we.onset.RData")


########################################################
###  RUN EV MIXED EFFECT MODELS
#ev.onset and has 2 lists of phenometric datasets, one for each data curation
#each list has 3 datasets, representing different aggregations 
metric<-"EV"

pheno<-"onset"
index<-1
summary.res.ev.on<-list()
sp.resp.ev.on<-list()
ev.onset.coefs<-list()
for(i in 1:length(ev.onset)) {
  curate<-i
  for(j in 1:length(ev.onset[[i]])) {
    aggreg<-  aggregs[j]
    model.input<-ev.onset[[i]][[j]] %>% dplyr::mutate(year=year-1850, aggreg=aggreg, curat=curate) %>%
      dplyr::select(name, region, name2, rndLat, year, metric, alt, aggreg, curat, n.14, n.max)
    model.input<-na.omit(model.input)
    if(nrow(model.input)>0) {
      effort<-T
    nsp<-length(unique(model.input$name2[!is.na(model.input$metric)]))
    ntot<-length((model.input$name2[!is.na(model.input$metric)]))
    effort.lme<-model.lme(model.input, effort=T)
    effort.summary<-rbind(c(pheno,metric,curate,aggreg,nsp,ntot,"yes",model.compar1(effort.lme, nsp)))
    effort.sp.resp<-data.frame(spname=sort(unique(model.input$name2[!is.na(model.input$metric)])), metric="ev.on",cur=curate,agg=aggreg,effort="yes", resp=model.stats.sp(effort.lme, nsp))
    #summary.res.ev.on[[index]]<-effort.summary
    #sp.resp.ev.on[[index]]<-effort.sp.resp
    #index<-index+1    
    
    noeffort.lme<-model.lme(model.input, effort=F)
    colnames(effort.summary)<-c("Metric","Estimator","Curation","Aggregation","N.SP","N.REC","effort","blah","best.model","AIC.df","AIC","r2m","r2c","neg","nsig","pos")
    lme.summary<-rbind(effort.summary, c(pheno,metric,curate,aggreg,nsp,ntot,"no",model.compar1(noeffort.lme, nsp)))
    colnames(lme.summary)<-c("Metric","Estimator","Curation","Aggregation","N.SP","N.REC","effort","blah","best.model","AIC.df","AIC","r2m","r2c","neg","nsig","pos")
    sp.resp<-bind_rows(effort.sp.resp,data.frame(spname=sort(unique(model.input$name2[!is.na(model.input$metric)])), metric="ev.on",cur=curate,agg=aggreg,effort="no", resp=model.stats.sp(noeffort.lme, nsp)))
    summary.res.ev.on[[index]]<-as.data.frame(lme.summary)
    sp.resp.ev.on[[index]]<-sp.resp
    temp1<-as.data.frame(summary(noeffort.lme)$coefficients) %>%
      dplyr::mutate(permut=paste(pheno,metric,curate,aggreg,"e0",sep="."), spname=row.names(summary(noeffort.lme)$coefficients))
    temp2<-as.data.frame(summary(effort.lme)$coefficients) %>%
      dplyr::mutate(permut=paste(pheno,metric,curate,aggreg,"e1",sep="."), spname=row.names(summary(effort.lme)$coefficients))
    ev.onset.coefs[[index]]<-bind_rows(temp1,temp2)  
    index<-index+1    
    }
  }
}

summary.res.ev.on[[1]]
ev.on.sum<-bind_rows(summary.res.ev.on)

sp.resp.ev.on[[1]]
ev.on.spresp<-bind_rows(sp.resp.ev.on)

ev.on.sprtot<-ev.on.spr %>%
  group_by(spname) %>%
  summarize(neg=length(which(resp==-1)), ns=length(which(resp==0)),pos=length(which(resp==1)))

write.csv(ev.on.spresp,file="output/ev_onset_sp.results.csv")
write.csv(ev.on.sum,file="output/ev_onset_model.results.csv")




########################################################
###  RUN WE MIXED EFFECT MODELS
#we.onset has 2 lists of phenometric datasets, one for each data curation
#each list has 3 datasets, representing different aggregations 
load("data/phenometrics/we.onset.RData")
we.onset[[2]][[1]]<-we.onset[[2]][[1]] %>% mutate(curat=2)
we.onset[[2]][[2]]<-we.onset[[2]][[2]] %>% mutate(curat=2)
metric<-"WE"

pheno<-"onset"
index<-1
summary.res.we.on<-list()
sp.resp.we.on<-list()
we.onset.coefs<-list()
#onset.coefs<-list()
for(i in 1:length(we.onset)) {
  curate<-i
  for(j in 1:length(we.onset[[i]])) {
    aggreg<-  aggregs[j]
    model.input<-we.onset[[i]][[j]] %>% dplyr::mutate(year=year-1850) %>%
      dplyr::select(name, region, name2, rndLat, year, metric, alt, aggreg, curat, n.14, n.max)
    model.input<-na.omit(model.input)
    if(nrow(model.input)>0) {
      effort<-T
      nsp<-length(unique(model.input$name2[!is.na(model.input$metric)]))
      ntot<-nrow(filter(model.input, !is.na(metric)))
      effort.lme<-model.lme(model.input, effort=T)
      effort.summary<-rbind(c(pheno,metric,curate,aggreg,nsp,ntot,"yes",model.compar1(effort.lme, nsp)))
      effort.sp.resp<-data.frame(spname=sort(unique(model.input$name2[!is.na(model.input$metric)])), metric="we.on",cur=curate,agg=aggreg,effort="yes", resp=model.stats.sp(effort.lme, nsp))
      #summary.res.we.on[[index]]<-effort.summary
      #sp.resp.we.on[[index]]<-effort.sp.resp
      #index<-index+1    
    
      noeffort.lme<-model.lme(model.input, effort=F)
      colnames(effort.summary)<-c("Metric","Estimator","Curation","Aggregation","N.SP","N.REC","effort","blah","best.model","AIC.df","AIC","r2m","r2c","neg","nsig","pos")
      lme.summary<-rbind(effort.summary, c(pheno,metric,curate,aggreg,nsp,ntot,"no",model.compar1(noeffort.lme, nsp)))
      colnames(lme.summary)<-c("Metric","Estimator","Curation","Aggregation","N.SP","N.REC","effort","blah","best.model","AIC.df","AIC","r2m","r2c","neg","nsig","pos")
      sp.resp<-bind_rows(effort.sp.resp,data.frame(spname=sort(unique(model.input$name2[!is.na(model.input$metric)])), metric="we.on",cur=curate,agg=aggreg,effort="no", resp=model.stats.sp(noeffort.lme, nsp)))
      summary.res.we.on[[index]]<-as.data.frame(lme.summary)
      sp.resp.we.on[[index]]<-sp.resp
      temp1<-as.data.frame(summary(noeffort.lme)$coefficients) %>%
        dplyr::mutate(permut=paste(pheno,metric,curate,aggreg,"e0",sep="."), spname=row.names(summary(noeffort.lme)$coefficients))
      temp2<-as.data.frame(summary(effort.lme)$coefficients) %>%
        dplyr::mutate(permut=paste(pheno,metric,curate,aggreg,"e1",sep="."), spname=row.names(summary(effort.lme)$coefficients))
      we.onset.coefs[[length(onset.coefs)+1]]<-bind_rows(temp1,temp2)  
      print(index)
      index<-index+1    
    }
  }
}

summary.res.we.on[[1]]
we.on.summary<-bind_rows(summary.res.we.on)

sp.resp.we.on[[1]]
we.on.spresp<-bind_rows(sp.resp.we.on)

we.on.sprtot<-we.on.spresp %>%
  group_by(spname) %>%
  summarize(neg=length(which(resp==-1)), ns=length(which(resp==0)),pos=length(which(resp==1)))

write.csv(we.on.spresp,file="output/we_onset_results.csv")
save(we.on.summary, we.on.spresp,we.on.sprtot,we.onset.coefs,ev.on.summary, ev.on.spresp,ev.on.sprtot,ev.onset.coefs,file="output/onset_model_results.RData")
write.csv(we.on.summary,file="output/we_onset_models.csv")


load("output/onset_model_results.RData")
onset.coefs<-bind_rows(bind_rows(ev.onset.coefs), bind_rows(we.onset.coefs))
onset_coefs<-(onset.coefs) %>%
  select(permut, spname, Estimate, `Pr(>|t|)`) %>%
  dplyr::mutate(sig=ifelse(`Pr(>|t|)`<0.05,1,0),sign=ifelse(`Pr(>|t|)`<0.05,sign(Estimate),0),
                param=ifelse(str_sub(permut,start=7,end=7)=="E","EV","WE"))

onset_coefs<- onset_coefs[grep("rndLat:name2", onset_coefs$spname),]



n.am.rows<-grep(".N. America", onset_coefs$spname)
onset_coefs<- onset_coefs  %>%
  dplyr::mutate(species=str_sub(spname, start = 13, end = str_locate(spname, "\\.")[,1]-1),
                region=str_sub(spname, start = str_locate(spname, "\\.")[,1]+1),
                aggreg=ifelse(str_length(permut)==25,"LYE",ifelse(str_length(permut)==20,"LY","L")))



### ONSET RESULTS
(f1.na.ev<-ggplot(data=filter(onset_coefs, region=="N. America", param=="EV"), aes(x=species, y=Estimate)) + geom_point(aes(shape=aggreg,color=as.factor(sign)) ) + 
  scale_shape_manual(values=c(15, 17, 19),guide="none") + ylim(-1,4.5) + 
  scale_color_manual(values=c(c("darkgray", "darkgreen")), guide="none") + 
    scale_x_discrete(limits=rev) +   
    theme_classic() + labs(x="Species (N. America)", title="Onset (Extreme) ~ Latitude coefficients") + 
  theme(text = element_text(size=8), axis.title.x = element_blank(), axis.text.x=element_blank()) + 
  geom_hline(yintercept=0, linetype="dashed") + 
  coord_flip()  )

(f1.eur.ev<-ggplot(data=filter(onset_coefs, region=="Europe", param=="EV"), aes(x=species, y=Estimate)) + geom_point(aes(shape=aggreg,color=as.factor(sign)) ) + 
  scale_shape_manual(values=c(15, 17, 19),guide="none") + ylim(-1,4.5) + 
  scale_color_manual(values=c(c("darkgray", "darkgreen")), guide="none") + 
    scale_x_discrete(limits=rev) +   
    theme_classic() + labs(x="Species (Europe)", y="Coefficient estimate") +  
  theme(text = element_text(size=8)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  coord_flip()  )

(f1.na.we<-ggplot(data=filter(onset_coefs, region=="N. America", param=="WE"), aes(x=species, y=Estimate)) + 
    geom_point(aes(shape=aggreg,color=as.factor(sign))) + 
  scale_shape_manual(values=c(15, 17, 19),guide="none") + ylim(-1,4.5) + 
  scale_color_manual(values=c("darkgray", "darkgreen"), name="Response sign", labels=c("Non-significant","Positive")) + 
    scale_x_discrete(limits=rev) +   
    theme_classic() + labs(title="Onset (Weibull) ~ Latitude") + 
    geom_hline(yintercept=0, linetype="dashed")  +
    theme(text = element_text(size=12), axis.title = element_blank(),axis.text.y = element_blank(),  axis.text.x=element_blank(), axis.ticks.y = element_blank(), legend.position="bottom") + 
  coord_flip()  )

(f1.eur.we<-ggplot(data=filter(onset_coefs, region=="Europe", param=="WE"), aes(x=species, y=Estimate)) + 
    geom_point(aes(shape=aggreg,color=as.factor(sign)) ) + 
  scale_shape_manual(values=c(15, 17, 19), name="Aggregation", labels=c("Latitude","Lat*Year","Lat*Year*Elevation")) + 
    ylim(-1,4.5) + 
  scale_color_manual(values=c(c("darkgray", "darkgreen")), guide="none") + 
  scale_x_discrete(limits=rev) +   
  theme_classic() + labs( y="Coefficient estimate") +  
  theme(text = element_text(size=12),axis.title.y=element_blank(),  axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position="bottom") + 
  geom_hline(yintercept=0, linetype="dashed")  + labs(x="") + 
  coord_flip()  )

library(ggpubr)
leg1 <- get_legend(f1.na.we) 
L1<-as_ggplot(leg1)
(f1.na.we<-f1.na.we + theme(legend.position = "none", text = element_text(size=8)))


leg2 <- get_legend(f1.eur.we)
L2<-as_ggplot(leg2)
(f1.eur.we<-f1.eur.we + theme(legend.position = "none",text = element_text(size=8)))


library(gridExtra)
fig.onsetres.s1<-grid.arrange(f1.na.ev,f1.na.we,f1.eur.ev,f1.eur.we,L1, L2, 
                              layout_matrix = rbind(c(1, 2), c(3, 4),c(5,5),c(6,6)),
                              widths=c(4.4,3), heights=c(10,8,1,1))
ggsave(filename="output/figs/onset.res.summary.png",fig.onsetres.s1, width = 6, height = 8)
ggsave(filename="output/figs/onset.res.summary2.png",fig.onsetres.s1, width = 6, height = 7)


names(ev.on.sprtot)<-c("spname","ev.neg","ev.ns","ev.pos")
names(we.on.sprtot)<-c("spname","we.neg","we.ns","we.pos")
sp.comp<-merge(ev.on.sprtot,we.on.sprtot, by=c("spname"))

(agreement<-sp.comp %>%
  filter((we.pos==12 & ev.pos==12) | (we.ns==12 & ev.ns==12) | (we.neg==12 & ev.neg==12)) )

(comp1<-sp.comp %>%
    filter((we.pos!=12 & ev.pos==12) | (we.ns!=12 & ev.ns==12) | (we.neg!=12 & ev.neg==12)) )
#ag, cc ly and lye SP
#bf always SP
#Gl SP except dc1 lat
#Is usually NS, SP for LY dc2, LYE no effort

(comp2p<-sp.comp %>%
    filter((we.pos==12 & ev.pos!=12) ) )
(comp2ns<-sp.comp %>%
      filter((we.ns==12 & ev.ns!=12) ) )
(comp2n<-sp.comp %>%
      filter((we.neg==12 & ev.neg!=12) ) )


names(ev.on.spr)<-c("spname","metric","agg","effort","resp","cur")
names(we.on.spr)<-c("spname","metric","agg","effort","resp","cur")
on.spr<-bind_rows(ev.on.spr, we.on.spr) %>%
  mutate(m1=paste(spname,metric,effort,cur,sep=".")) %>%
  select(m1, agg, resp) %>%
  pivot_wider(names_from=agg, values_from=resp) %>%
  dplyr::mutate(add.yr=ifelse(`Lat*yr`==Lat,1,0),add.elev=ifelse(`Lat*yr*elev`==`Lat*yr`,1,0))

on.spr.effort<-bind_rows(ev.on.spr, we.on.spr) %>%
  mutate(m1=paste(spname,metric,agg,cur,sep=".")) %>%
  select(m1, effort, resp) %>%
  pivot_wider(names_from=effort, values_from=resp) %>%
  dplyr::mutate(same=ifelse(yes==no,1,0))

on.spr.cur<-bind_rows(ev.on.spr, we.on.spr) %>%
  mutate(m1=paste(spname,metric,agg,effort,sep=".")) %>%
  select(m1, cur, resp) %>%
  pivot_wider(names_from=cur, values_from=resp) %>%
  dplyr::mutate(same=ifelse(`1`==`2`,1,0))
(percent.cur.diff<-round(length(which(on.spr.cur$same==0))/nrow(on.spr.cur),2)*100)

on.spr.agg<-bind_rows(ev.on.spr, we.on.spr) %>%
  mutate(m1=paste(spname,metric,effort,cur,sep=".")) %>%
  select(m1, agg, resp) %>%
  pivot_wider(names_from=agg, values_from=resp) %>%
  dplyr::mutate(same.a=ifelse("Lat"=="Lat*yr",1,0),same.a=ifelse("Lat*yr"=="Lat*yr*elev",1,0))





sp.comp.1<-merge(ev.on.spr[,c(1,3:5)], we.on.spr[,c(1,3:5)], by=c("spname","agg","effort"), all=T) %>%
  dplyr::mutate(same=ifelse(ev.resp==we.resp,1,0))
table(sp.comp.1$same)

#(ev.on.spr,we.on.spr) %>% arrange(spname) %>%
#  pivot_wider(id_cols=c(spname), names_from=metric, values_from=resp)
 



load("output/onset_model_results.RData")
other_coefs<-bind_rows(onset.coefs)
table(other_coefs$spname)
x1<-c(grep("alt2", other_coefs$spname),grep("n.max", other_coefs$spname),grep("n.14", other_coefs$spname))
other_coefs<-other_coefs[x1,] %>%
  dplyr::mutate(sig=ifelse(`Pr(>|t|)`<0.05,1,0),sign=ifelse(`Pr(>|t|)`<0.05,sign(Estimate),0),
                param=ifelse(str_sub(permut,start=7,end=7)=="E","EV","WE"),
                aggreg=ifelse(str_length(permut)==25,"LYE",ifelse(str_length(permut)==20,"LY","L")))


(fig.exp.coef<-ggplot(data=other_coefs, aes(y=Estimate, x=spname)) + geom_boxplot(color="gray") + 
  geom_jitter(width=0.1, height=0,aes(shape=aggreg,color=as.factor(sign)),size=2) + theme_classic() + 
  scale_color_manual(values=c("blue","darkgray", "darkgreen"), name="Onset Response sign") + 
  scale_shape_manual(values=c(15,17,19), name="Aggregation", labels=c("Lat","Lat*Year","Lat*Year*Elev")) + 
  labs(x="Explanatory variables", y="Coefficient estimate") + 
  scale_x_discrete(labels=c("Elevation","n.14 effort","n.max effort")) + 
  geom_hline(yintercept=0, linetype="dashed") )

ggsave(fig.exp.coef, file="output/figs/fig.supp1.onset.expvar.png")



### FIGURE 2??
onset_coefs<-(onset.coefs) %>%
  select(permut, spname, Estimate, `Pr(>|t|)`) %>%
  dplyr::mutate(sig=ifelse(`Pr(>|t|)`<0.05,1,0),sign=ifelse(`Pr(>|t|)`<0.05,sign(Estimate),0),
                param=ifelse(str_sub(permut,start=7,end=7)=="E","EV","WE"),
                aggreg="")

#onset_coefs$aggreg[grep("Lat",onset_coefs$permut)]<-"Lat"
#onset_coefs$aggreg[grep("yr",onset_coefs$permut)]<-"LatYr"
#onset_coefs$aggreg[grep("elev",onset_coefs$permut)]<-"LatYrElev"

onset_coefs<- onset_coefs[grep("rndLat:name2", onset_coefs$spname),]

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

onset_fig<- onset_coefs  %>%
  dplyr::mutate(species=str_sub(spname, start = 13, end = str_locate(spname, "\\.")[,1]-1),
                name2=str_sub(spname, start = 13),
                region=str_sub(spname, start = str_locate(spname, "\\.")[,1]+1),
                aggreg=ifelse(str_length(permut)==25,"Lat*Yr*Elev",ifelse(str_length(permut)==20,"Lat*Yr","Lat")),
                effort=substrRight(permut,1), cur=substr(permut,10,10)) %>%
  group_by(species, region, sign) %>%
    add_tally() %>%
  group_by(species, region) %>%
  dplyr::mutate(nm=max(n),summary=ifelse(nm>20 & median(sign)==1,"1",ifelse(nm>20 & median(sign)==0, "0","mixed")))
  


onset.fig2<-filter(onset_fig,cur==2 & !is.na(summary) & effort==1) %>% 
  dplyr::mutate(agnum=ifelse(aggreg=="Lat",1,ifelse(aggreg=="Lat*Yr",2,3)),
                snum=ifelse(summary=="0",1,ifelse(summary=="1",2,3)))

library(viridis)
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = 'A')
v2<-v_colors[c(10,8,6)]

fig.aggreg<-ggplot(data=onset.fig2, aes(x=summary, y=Estimate), size=2) + 
  geom_boxplot(aes(color=aggreg)) +
  geom_jitter(data=filter(onset.fig2, aggreg=="Lat"),width=0.05,aes(x=snum-0.25, y=Estimate, fill=as.factor(sign), shape=param), size=2) + 
  geom_jitter(data=filter(onset.fig2, aggreg=="Lat*Yr"),width=0.05,aes(x=snum, y=Estimate, fill=as.factor(sign), shape=param), size=2) + 
  geom_jitter(data=filter(onset.fig2, aggreg=="Lat*Yr*Elev"),width=0.05,aes(x=snum+0.25, y=Estimate, fill=as.factor(sign), shape=param), size=2) + 
  scale_fill_manual(values=c("gray","darkgreen"), labels=c(">0.05","<0.05"), name="P-value") + 
  scale_shape_manual(values=c(21,24), labels=c("Extreme","Weibull"), name="Metric") + 
  scale_color_manual(values=v2, labels=c("Lat","Lat*Yr","Lat*Yr*Elev"), name="Aggregation") + #labels=c("Extreme","Weibull"), name="Metric") + 
  geom_hline(yintercept=0, linetype="dashed") + theme_classic() +
  labs(title="Onset~latitude", y="Coefficient estimate") + 
  scale_x_discrete(labels=c("0" = "No trend (n=28)", "1" = "Pos. trend (n=11)","mixed"="variable (n=48)", name=element_blank()))
fig.aggreg



fig.aggreg.1<-ggplot(data=onset.fig2, aes(x=aggreg, y=Estimate, param=summary), size=2) + 
  geom_boxplot(aes(x=aggreg, y=Estimate, fill=summary)) +
  geom_jitter(data=filter(onset.fig2, summary=="0"),width=0.05,aes(x=agnum-0.25, y=Estimate, color=as.factor(sign), shape=param), size=2) + 
  geom_jitter(data=filter(onset.fig2, summary=="1"),width=0.05,aes(x=agnum, y=Estimate, color=as.factor(sign), shape=param), size=2) + 
  geom_jitter(data=filter(onset.fig2, summary=="mixed"),width=0.05,aes(x=agnum+0.25, y=Estimate, color=as.factor(sign), shape=param), size=2) + 
  scale_color_manual(values=c("darkgray","darkgreen"), labels=c(">0.05","<0.05"), name="P-value") + 
  scale_shape_manual(values=c(1,4), labels=c("Extreme","Weibull"), name="Metric") + 
  scale_fill_manual(values=c("white","white","white"), guide="none") + #labels=c("Extreme","Weibull"), name="Metric") + 
  geom_hline(yintercept=0, linetype="dashed") + theme_classic() +
  labs(x="Aggregation",title="Onset~latitude", y="Coefficient estimate")

onset.fig2b<-filter(onset_fig,cur==2) %>% dplyr::mutate(pnum=ifelse(param=="EV",1,2))

(fig.effort<-ggplot(data=filter(onset.fig2b, aggreg=="Lat"), aes(x=param, y=Estimate, color=effort), size=2) + 
  geom_boxplot(aes(x=param, y=Estimate, fill=summary)) +
  geom_jitter(data=filter(onset.fig2b, summary=="0"),width=0.1,aes(x=pnum-0.25, y=Estimate, shape=as.factor(sign), shape=effort)) + 
  geom_jitter(data=filter(onset.fig2b, summary=="1"),width=0.1,aes(x=pnum, y=Estimate, shape=as.factor(sign), shape=effort)) + 
  geom_jitter(data=filter(onset.fig2b, summary=="mixed"),width=0.1,aes(x=pnum+0.25, y=Estimate, shape=as.factor(sign), shape=effort)) 
)
+ 
  scale_color_manual(values=c("darkgray","darkgreen"), labels=c("nonsignificant","positive"), name="Sign") + 
  scale_shape_manual(values=c(15,17,19), labels=c("Extreme","Weibull"), name="Metric") + 
  scale_fill_manual(values=c("white","white","white"), guide="none") + #labels=c("Extreme","Weibull"), name="Metric") + 
  geom_hline(yintercept=0, linetype="dashed") + theme_classic() +
  labs(x="Aggregation",title="Onset~latitude coefficients")


library(tidyverse)
###
load("data/curations.RData")
sp.n<-cured.data2 %>%
  group_by(name2) %>%
  summarize(n.occ=n(),nlat=length(unique(rndLat)))

load("data/phenometrics/ev.metrics.RData")
sp.n2<-ev.onset[[2]][[2]] %>%
  group_by(name2) %>%
  summarize(n.onset=n(),nlat=length(unique(rndLat)))


onset_coefs_f2<-merge(onset_fig, sp.n2, by.x="name2", by.y="name2")
  
(fig.2n<-ggplot(data=filter(onset_coefs_f2, aggreg=="Lat*Yr"), aes(x=nlat, y=Estimate, color=as.factor(sign), shape=param, size=n.onset)) + 
    geom_point() + theme_classic() +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_manual(values=c("darkgray","darkgreen"), labels=c("nonsignificant","positive"), name="Sign") + 
    scale_shape_manual(values=c(1,4), labels=c("Extreme","Weibull"), name="Metric") +
    labs(x="# Latitudes with Lat*Yr phenometrics", y="Onset~Latitude coefficients", title="Latitudinal patterns by # latitudes with data")  
)

sp.n1<-ev.onset[[2]][[1]] %>%
  group_by(name2) %>%
  summarize(n.onset=n(),nlat=length(unique(rndLat)))


onset_coefs_f21<-merge(onset_fig, sp.n1, by.x="name2", by.y="name2")

(fig.2n1<-ggplot(data=filter(onset_coefs_f21, aggreg=="Lat"), aes(x=nlat, y=Estimate, color=as.factor(sign), shape=param), size=3) + 
    geom_point() + theme_classic() +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_manual(values=c("darkgray","darkgreen"), labels=c("nonsignificant","positive"), name="Sign") + 
    scale_shape_manual(values=c(1,4), labels=c("Extreme","Weibull"), name="Metric") +
    labs(x="# Latitudes with Lat phenometrics", y="Onset~Latitude coefficients", title="Latitudinal patterns by # latitudes with data")  
)


sp.n3<-ev.onset[[2]][[3]] %>%
  group_by(name2) %>%
  summarize(n.onset=n(),nlat=length(unique(rndLat)))


onset_coefs_f23<-merge(onset_fig, sp.n3, by.x="name2", by.y="name2")

(fig.2n3<-ggplot(data=filter(onset_coefs_f23, aggreg=="Lat*Yr*Elev"), aes(x=nlat, y=Estimate, color=as.factor(sign), shape=param), size=3) + 
    geom_point() + theme_classic() +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_manual(values=c("darkgray","darkgreen"), labels=c("nonsignificant","positive"), name="Sign") + 
    scale_shape_manual(values=c(1,4), labels=c("Extreme","Weibull"), name="Metric") +
    labs(x="# Latitudes with Lat*Yr*Elev phenometrics", y="Onset~Latitude coefficients", title="Latitudinal patterns by # latitudes with data")  
)







(fig.2p<-ggplot(data=filter(onset_coefs_f2, aggreg=="Lat*Yr"), aes(x=nlat, y=Estimate, color=as.factor(sign), shape=param), size=3) + 
    geom_point() + theme_classic() +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_manual(values=c("darkgray","darkgreen"), labels=c("nonsignificant","positive"), name="Sign") + 
    scale_shape_manual(values=c(1,4), labels=c("Extreme","Weibull"), name="Metric") +
    labs(x="# Latitudes with Lat*Yr phenometrics", y="Onset~Latitude coefficients", title="Latitudinal patterns by # latitudes with data")  
)

onr<-onset_coefs_f2 %>%
  dplyr::mutate(lat.th=ifelse(nlat>5,2,1)) %>%
  group_by(lat.th,name2,sign) %>%
  tally() %>%
  pivot_wider(names_from=sign, values_from=n)


  dplyr::summarize(medsign=median(sign))

summary(lm(sign~nlat, data=onset_coefs_f2))



############################################# EFFORT

aggregs<-c("lat","latyr","latyrelev")

for(i in 1:2) {
  for(j in 1:3) {
    ev.onset[[i]][[j]]<-ev.onset[[i]][[j]] %>%
      dplyr::mutate(curat=i, aggreg=aggregs[j])
    we.onset[[i]][[j]]<-we.onset[[i]][[j]] %>%
      dplyr::mutate(curat=i)
  }
}

ev1<-bind_rows(ev.onset) %>%
  dplyr::select(name2,rndLat,year=year,ev.alt=alt, elev.yr,ev.metric=metric,aggreg, curat)

  
we1<-bind_rows(we.onset) %>%
  dplyr::select(name2,rndLat,year=year,we.alt=alt, elev.yr,we.metric=metric,aggreg, curat)
test.lye<-merge(ev1, we1, by=c("name2","rndLat","year","aggreg","curat","elev.yr"), all=F)
test.ly<-merge(filter(ev1, aggreg=="latyr"), filter(we1, aggreg=="latyr"), by=c("name2","rndLat","year","aggreg","curat"), all=F)
test.l<-merge(filter(ev1, aggreg=="lat"), filter(we1, aggreg=="lat"), by=c("name2","rndLat","aggreg","curat"), all=F)


onetoone<-c(0,50,100,200)

ggplot(data=test.lye, aes(x=ev.metric, y=we.metric)) + geom_point() + 
  geom_line(aes(x=ev.metric, y=ev.metric), linetype="dashed") + facet_wrap(~curat)

ggplot(data=test.ly, aes(x=ev.metric, y=we.metric)) + geom_point() + 
  geom_line(aes(x=ev.metric, y=ev.metric), linetype="dashed") + facet_wrap(~curat)

ggplot(data=test.l, aes(x=ev.metric, y=we.metric)) + geom_point() + 
  geom_line(aes(x=ev.metric, y=ev.metric), linetype="dashed") + facet_wrap(~curat)

ev2<-bind_rows(ev.onset) %>%
  group_by(name2,rndLat,aggreg) %>%
  dplyr::mutate(mednmax=median(n.max, na.rm=T),medn.14=median(on14, na.rm=T))

ggplot(data=ev2, aes(x=as.factor(rndLat), y=mednmax)) + geom_point() + 
  geom_boxplot() + facet_wrap(~aggreg)

we2<-bind_rows(we.onset) %>%
  group_by(name2,rndLat,aggreg) %>%
  dplyr::mutate(mednmax=log(median(n.max, na.rm=T)),medn.14=median(n.14.onset, na.rm=T))

(fignmaxlat<-ggplot(data=filter(we2, aggreg=="lat"), aes(x=as.factor(rndLat), y=mednmax)) + geom_point() + 
  geom_boxplot() + labs(x="Latitudinal band", y="Log median n.max by population") + 
  theme_classic() + 
    annotate(geom="text",x=3,y=5,label="by latitude", hjust="left") + 
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())  )

(fign14lat<-ggplot(data=filter(we2, aggreg=="lat"), aes(x=as.factor(rndLat), y=medn.14)) + geom_point() + 
    geom_boxplot() + labs( y="Median n.14 by population", x="Latitudinal band") + 
    theme_classic() + 
    annotate(geom="text",x=3,y=60,label="by latitude", hjust="left") + 
    theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x = element_text())  )

(fignmaxly<-ggplot(data=filter(we2, aggreg=="latyr"), aes(x=as.factor(rndLat), y=mednmax)) + geom_point() + 
    geom_boxplot() + labs(x="Latitudinal band", y="Log median n.max by population") + 
    theme_classic() + 
    annotate(geom="text",x=3,y=3,label="by latitude*year", hjust="left") + 
    theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())  )

(fign14ly<-ggplot(data=filter(we2, aggreg=="latyr"), aes(x=as.factor(rndLat), y=medn.14)) + geom_point() + 
    geom_boxplot() + labs( y="Median n.14 by population", x="Latitudinal band") + 
    theme_classic() + 
    annotate(geom="text",x=3,y=16,label="by latitude*year", hjust="left") + 
    theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x = element_text())  )

(fignmaxlye<-ggplot(data=filter(we2, aggreg=="latyrelev"), aes(x=as.factor(rndLat), y=mednmax)) + geom_point() + 
    geom_boxplot() + labs(x="Latitudinal band", y="Log median n.max by population") + 
    theme_classic() + 
    annotate(geom="text",x=3,y=2,label="by latitude*year*elevation", hjust="left") + 
    theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())  )

(fign14lye<-ggplot(data=filter(we2, aggreg=="latyrelev"), aes(x=as.factor(rndLat), y=medn.14)) + geom_point() + 
    geom_boxplot() + labs( y="Median n.14 by population", x="Latitudinal band") + 
    theme_classic() + 
    annotate(geom="text",x=3,y=15,label="by latitude*year*elevation", hjust="left") + 
    theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x = element_text())  )

grid.arrange(fignmaxlat,fign14lat,fignmaxly,fign14ly,fignmaxlye,fign14lye)
title="Sampling effort for Latitude aggregation"
