#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Last Update 2021-08
#HERE: Data curation

#load packages
library(tidyverse) #to be tidy
source("Rcode/data.curation.functions.R")
######################################################
### LOAD AND CURATE DATA
#load occurrence data formatted from F2020 raw data
load("data/occurrences.RData")

#initial filter
F1.data<-alldata %>%
  filter(doy>1, name %notin% c("",NA)) %>%
  mutate(name2=paste(name,region,sep="."))

t1<-calc.sum(F1.data)

#Remove species per F2021 & recalculate statistics
f2021.removed.sp<-c("Amblyscirtes vialis","Anthocharis sara","Boloria bellona",
                    "Boloria euphrosyne", "Coenonympha glycerion","Coenonympha tullia",
                    "Lasiommata maera","Limenitis archippus","Limenitis arthemis",
                    "Lycaena phlaeas","Melitaea diamina","Melitaea athalia",
                    "Ochlodes sylvanus","Papilio machaon","Pieris marginalis",
                    "Pyrgus malvae","Thorybes pylades")

Filter.sp<-F1.data %>%
  filter(name%notin%f2021.removed.sp)
t2<-calc.sum(Filter.sp)


#Remove records with altitude<0: definite errors
Filter.alt1<-Filter.sp %>%
  filter(alt>=0)
t3<-calc.sum(Filter.alt1)

#Remove records from January and February
Filter.doy1<-Filter.alt1 %>%
  filter(month>2)
t4<-calc.sum(Filter.doy1)

#Remove records from December
Filter.doy2<-Filter.doy1 %>%
  filter(month<12)
t5<-calc.sum(Filter.doy2)

tab2<-Filter.doy2 %>%
  group_by(name, region, rndLat) %>%
  tally() %>%
  filter(n>=5) 
median.range(tab2$n)

tab2b<-tab2 %>%
  group_by(name, region) %>%
  summarize(nlat=sum(n))
median.range(tab2b$nlat)

tab2c<-tab2 %>% group_by(name, region) %>% summarize(nmet=n(), nocc=sum(n))
median.range(tab2c$nmet)
median.range(tab2c$nocc)

scope<-Filter.doy2 %>%
  group_by(name, region, rndLat, year) %>%
  summarize(nocc=n()) %>%
  filter(nocc>=5) %>%
  group_by(name, region) %>%
  summarize(nlat=length(unique(rndLat)), nyear=length(unique(year)))


#Remove records above alt 500
Filter.alt2<-Filter.doy2 %>%
  filter(alt<=500)
t6<-calc.sum(Filter.alt2)

#Remove records with only 1 observation in a latitudinal band
Filter.data1<-Filter.alt2 %>%
  group_by(name, region, rndLat) %>%
  add_tally() %>%
  filter(n>1)
t7<-calc.sum(Filter.data1)

curations<-list(F1.data, Filter.sp, Filter.alt1, Filter.doy1, Filter.doy2, Filter.alt2, Filter.data1)
#Create table for supplement 1
table.s1<-as.data.frame(cbind(t1,t2,t3,t5,t6,t7))

names(table.s1)<-c('F2020 dataset','Species filter','Alt>=0',
'Month in (2-12)','Alt<=500','N/rndLat > 1')

row.names(table.s1)<-c('Total # obs',
                       'N species',
                       'N species-region datasets',
                       '#obs per dataset',
                       '#obs per latitudinal band / onset AY',
                       '#obs per annual onset',
                       'Number of onsets (all years) per lat band',
                       'Number of annual onsets per lat band',
                       'Number of lat bands per dataset ',
                       'Number of onsets (all years) per dataset',
                       'Number of annual onsets per dataset')

write.csv(table.s1, file="output/table.s1.csv")


#Add elevation ranges to datasets for further comparison
cured.data1<-comb.elev(add.elev(Filter.sp))
cured.data1<-comb.elev(cured.data1)
summary(cured.data1)
table(cured.data1$elev.yr[cured.data1$ny1<5])

cured.data2<-comb.elev(add.elev(Filter.doy2))
cured.data2<-comb.elev(cured.data2)
summary(cured.data2)
table(cured.data2$elev.yr[cured.data2$ny1<5])

cd1<-cured.data2 %>%
  group_by(name2, rndLat, year, elev.yr) %>%
  add_tally() %>%
  filter(n>=5) %>%
  group_by(name2,rndLat,year,elev.yr) %>%
  summarize(noccon=length(which(doy<(min(doy)+15))),noccoff=length(which(doy>(max(doy)-15))))



save(cured.data1, cured.data2, file="data/curations.RData")



thresholds.summary<-Filter.doy2 %>%
  group_by(name2, rndLat, year) %>%
  tally()
length(unique(thresholds.summary$name2[thresholds.summary$n>1]))
length(unique(thresholds.summary$name2[thresholds.summary$n>4]))
length(unique(thresholds.summary$name2[thresholds.summary$n>9]))
summary(filter(thresholds.summary, n>1))

x1<-thresholds.summary %>%
  filter(n>1) %>%
  group_by(name2,rndLat ) %>%
  summarize(nyr=n())
summary(x1)

x2<-thresholds.summary %>%
  filter(n>9) %>%
  group_by(name2) %>%
  summarize(nlat=length(unique(rndLat)), nyr=length(unique(year)))
summary(x2)


##correlating effort and explanatory variables
effort<-Filter.doy2 %>%
  group_by(name2) %>%
  mutate(nlat=length(unique(rndLat)),meanlat=round(mean(rndLat, na.rm=T)), lat.diff=abs(rndLat-meanlat), elev=ceiling(alt/100)) %>%
  filter(nlat>=3) %>%
  group_by(name2, meanlat, rndLat,lat.diff) %>%
  tally() 

eff2<-Filter.doy2 %>%
  group_by(name2, rndLat) %>%
  tally()
  
effort.lm1<-lmer(n~lat.diff:name2+(1|name2),data=effort)
effort.lm1.best<-get_model(step(effort.lm1))
summary(effort.lm1.best)
r.squaredGLMM(effort.lm1.best)


library(ggplot2) #to make plots
library(sjPlot) #to plot mixed-effect model results
library(viridis) #colorblind friendly color palette for plots

q_colors =  90 # for no particular reason
v_colors =  viridis(q_colors, option = 'A')

plot_model(effort.lm1.best, type="pred")

fig1<-ggplot(data=effort, aes(x=rndLat, y=log(n), color=name2)) + 
  geom_smooth(data=filter(effort, rndLat<=meanlat),method = glm, formula = y ~ x, se = F) + 
  geom_smooth(data=filter(effort, rndLat>=meanlat),method = glm, formula = y ~ x, se = F) + 
  scale_color_manual(values=v_colors) + 
  geom_point() + theme(legend.position="none")
  
  #stat_summary(fun=mean, colour="red", geom="line", size = 3) # draw a mean line in the data
fig1

fig1<-ggplot(data=eff2, aes(x=rndLat, y=log(n), color=name2)) + 
  geom_smooth(method = glm, formula = y ~ x, se = F) + 
  scale_color_manual(values=v_colors) + 
  geom_point() + theme(legend.position="none")

#stat_summary(fun=mean, colour="red", geom="line", size = 3) # draw a mean line in the data
fig1



##Use Filter.data1 to create figure panels
library("grid")
library("ggplotify")

pheno.figdata<-Filter.data1 %>%
  group_by(name2, rndLat) %>%
  summarize(n.occ=n(), onset=min(SuccDay), term=max(SuccDay)) %>%
  filter(n.occ>1)
  
phenoyr.figdata<-Filter.data1 %>%
  group_by(name2, rndLat, year) %>%
  summarize(n.occ=n(), onset=min(SuccDay), term=max(SuccDay)) %>%
  filter(n.occ>1)


fd1<-Filter.data1 %>% group_by(name2, rndLat) %>% tally() 
fd2<-Filter.data1 %>% group_by(name2, rndLat) %>% tally() %>% group_by(rndLat) %>% summarize(medn=median(n))
(p1 <- ggplot(fd1, aes(rndLat,name2)) + geom_tile(aes(fill = log(n, base=10)),colour = "black") +
    scale_fill_gradient(low = "white", high="steelblue"))


ggplot(pheno.figdata, aes(x=as.factor(rndLat), y=log(n.occ))) + geom_boxplot() + theme_classic()
hist(Filter.doy2$year)
hist(Filter.doy2$alt)
hist(Filter.doy2$rndLat)

## DATA EXPLORATION
(exp.yr<-ggplot() + 
  geom_boxplot(data=(Filter.doy2%>%group_by(name2,year)%>%tally()), aes(x=as.factor(year), y=log(n,base=10))) + 
  geom_point(data=(Filter.doy2%>%group_by(name2,year)%>%tally()%>%
                    group_by(year)%>%dplyr::summarize(n=median(n))), 
            aes(x=as.factor(year),y=log(n, base=10)), color="red") + 
  annotate("text", y=c(rep(-0.1,4)), x=c(9,43,93,143), label=c("1850","1900","1950","2000")) + 
  theme_classic() + labs(x="Year",y='# Occurrences log10') + 
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) 
)

(exp.lat<-ggplot() + 
  geom_boxplot(data=(Filter.doy2%>%group_by(name2,rndLat)%>%tally()), aes(x=as.factor(rndLat), y=log(n, base=10))) + 
  geom_point(data=(Filter.doy2%>%group_by(name2,rndLat)%>%tally()%>%
                     group_by(rndLat)%>%dplyr::summarize(n=median(n))), 
             aes(x=as.factor(rndLat),y=log(n, base=10)), color="red") + 
  annotate("text", y=c(rep(-0.1,6)), x=c(1,3,11,21,31,41), label=c("6","32","40","50","60","70")) + 
  theme_classic() + labs(x="Latitudinal band",y='# Occurrences (log10)') + 
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) 
)

(exp.elev<-ggplot() + 
  geom_boxplot(data=(Filter.doy2%>%group_by(name2,floor(alt/100))%>%tally()), aes(x=as.factor(`floor(alt/100)`), y=log(n, base=10))) + 
    geom_point(data=(Filter.doy2%>%group_by(name2,floor(alt/100))%>%tally()%>%
                       group_by(`floor(alt/100)`)%>%dplyr::summarize(n=median(n))), 
               aes(x=as.factor(`floor(alt/100)`),y=log(n, base=10)), color="red") + 
    annotate("text", y=c(rep(-0.1,6)), x=c(1,6,11,21,31,40), label=c("0","500","1000","2000","3000","4000")) + 
  theme_classic() + labs(x="Elevation (m)",y='# Occurrences (log10)') + 
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) 
)
##Combine panels into Figure 1:
library(gridExtra)
(figure <- grid.arrange(exp.lat, exp.yr, exp.elev, ncol = 1, nrow = 3) )
ggsave("output/figs/FigS1.1a.png", plot=figure, width = 4, height = 6, dpi = 300, units = "in", device='png')


ggplot(pheno.figdata, aes(x=rndLat, y=log(n.occ), color=name2)) + geom_point() + 
  scale_color_viridis_d(guide=F) + 
  geom_smooth(method="lm", se=F) + theme_classic()

#boxplot of coefficients
modnlat<-lmer(log(n.occ,base=10)~-1+rndLat:name2+(1|name2), data=pheno.figdata)
lm0<-summary(modnlat)
lma1<-as.data.frame(lm0$coefficients)
names(lma1)<-c("est","se","df","tval","pval")
lma1<-lma1 %>% dplyr::mutate(sig=ifelse(pval<0.05,1,0),signn=ifelse(pval<0.05, sign(est),0))
lm0r1<-paste("=",round(r.squaredGLMM(modnlat)[1],2))
pan.nlat<-ggplot(data=lma1, aes(x=1,y=est)) + 
  geom_boxplot() + ylim(-0.015,0.05) + 
  geom_hline(yintercept=0, linetype="dashed") + 
  geom_jitter(aes(color=factor(sig))) + theme_classic() + 
  scale_color_manual(values=c("grey","black"), guide="none") + 
  annotate("text",x=1,y=-0.01,label=bquote(r^2 ~.(lm0r1))) + 
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  labs(x="", y="coefficients", subtitle="log(N) ~ Lat:sp")
inset.nlat<-as.grob(pan.nlat)

ggplot(pheno.figdata, aes(x=rndLat, y=log(n.occ), color=name2)) + geom_point() + 
  scale_color_viridis_d(guide="none") + 
  #geom_smooth(method="lm", se=F) + 
  theme_classic()

#pheno.figdata.ev<-pheno.figdata %>%
#  group_by(name2) %>%
#  filter(rndLat==min(rndLat) | rndLat==max(rndLat))
  
fig1a<-ggplot(pheno.figdata, aes(x=rndLat, y=log(n.occ), color=name2)) + geom_line() + 
  scale_color_viridis_d(guide="none") + 
  #geom_smooth(method="lm", se=F) + 
  #geom_point(data=pheno.figdata.ev)
  labs(x="Latitudinal band", y="# Occurrences (log10)") + 
  #annotation_custom(inset.nlat, xmin=33, xmax=48, ymin=5, ymax=10) + 
  theme_classic() 
fig1a

fig1d<-ggplot(phenoyr.figdata, aes(x=rndLat, y=log(n.occ), color=name2)) + geom_line() + 
  scale_color_viridis_d(guide="none") + 
  #geom_smooth(method="lm", se=F) + 
  #geom_point(data=pheno.figdata.ev)
  labs(x="Latitudinal band", y="# Occurrences (log10)") + 
  #annotation_custom(inset.nlat, xmin=33, xmax=48, ymin=5, ymax=10) + 
  theme_classic() 
fig1d



## ONSET ~ N
modon.n<-lmer(onset~-1+n.occ:name2+(1|name2), data=pheno.figdata)
lm.on<-summary(modon.n)
lm1b<-as.data.frame(lm.on$coefficients)
names(lm1b)<-c("est","se","df","tval","pval")
lm1b<-lm1b %>% dplyr::mutate(sig=ifelse(pval<0.05,1,0),signon=ifelse(pval<0.05,sign(est),0))
lm1br<-paste("=",round(r.squaredGLMM(modon.n)[1],2))
lm1b<-lm1b %>% dplyr::mutate(rd=floor(as.numeric(est))) %>% group_by(signon) %>% tally()
(on.n<-ggplot(data=lm1b, aes(x=1, fill=as.factor(signon), y=n)) + geom_bar(position="stack", stat="identity") + 
  geom_vline(xintercept=(-0.5), linetype="dashed") + theme(x.axis.ticks=element_blank()) + 
  theme_classic() )

fig1b<-ggplot(pheno.figdata, aes(x=log(n.occ), y=onset, color=name2)) + geom_line() + 
  scale_color_viridis_d(guide="none") + 
  #geom_smooth(method="lm", se=F) + 
  #geom_point(data=pheno.figdata.ev)
  labs(x="# occurrences (log10)", y="F2021 Onset") + 
#  annotation_custom(inset.on.n, xmin=6, xmax=10, ymin=200, ymax=300) + 
  #annotate("text",x=7.5, y=200. label=bquote())
  theme_classic() 
fig1b
(fig1c<-ggplot(pheno.figdata, aes(x=log(n.occ), y=term, color=name2)) + geom_line() + 
    scale_color_viridis_d(guide="none") + 
    #geom_smooth(method="lm", se=F) + 
    #geom_point(data=pheno.figdata.ev)
    labs(x="# occurrences (log10)", y="F2021 Termination") + 
    #  annotation_custom(inset.on.n, xmin=6, xmax=10, ymin=200, ymax=300) + 
    theme_classic() )
  

## ONSET ~ N
modoff.n<-lmer(term~-1+n.occ:name2+(1|name2), data=pheno.figdata)
lm.off<-summary(modoff.n)
lm1c<-as.data.frame(lm.off$coefficients)
names(lm1c)<-c("est","se","df","tval","pval")
lm1c<-lm1c %>% dplyr::mutate(sig=ifelse(pval<0.05,1,0),signoff=ifelse(pval<0.05,sign(est),0))
lm1cr<-paste("=",round(r.squaredGLMM(modoff.n)[1],2))
lm1c<-lm1c %>% dplyr::mutate(rd=floor(as.numeric(est))) %>% group_by(signoff) %>% tally()




fig1e<-ggplot(phenoyr.figdata, aes(x=log(n.occ), y=onset, color=name2)) + geom_line() + 
  scale_color_viridis_d(guide="none") + 
  #geom_smooth(method="lm", se=F) + 
  #geom_point(data=pheno.figdata.ev)
  labs(x="# occurrences (log10)", y="Lat-Yr Onset") + 
  #  annotation_custom(inset.on.n, xmin=6, xmax=10, ymin=200, ymax=300) + 
  #annotate("text",x=7.5, y=200. label=bquote())
  theme_classic() 
fig1e
(fig1f<-ggplot(phenoyr.figdata, aes(x=log(n.occ), y=term, color=name2)) + geom_line() + 
    scale_color_viridis_d(guide="none") + 
    #geom_smooth(method="lm", se=F) + 
    #geom_point(data=pheno.figdata.ev)
    labs(x="# occurrences (log10)", y="Lat-Yr Termination") + 
    #  annotation_custom(inset.on.n, xmin=6, xmax=10, ymin=200, ymax=300) + 
    theme_classic() )

##Combine panels into Figure 1:
library(gridExtra)
figure <- grid.arrange(fig1a, fig1b, fig1c, ncol = 3, nrow = 1)
figure
ggsave("output/figs/Fig1.png", plot=figure, width = 10, height = 4, dpi = 300, units = "in", device='png')

  
lma.on<-as.data.frame(summary(lm(onset~log(n.occ,base=10):name2, data=pheno.figdata))$coefficients)[-1,] %>% mutate(model="onset~n")
lma.off<-as.data.frame(summary(lm(term~log(n.occ,base=10):name2, data=pheno.figdata))$coefficients)[-1,] %>% mutate(model="term~n")
lma<-rbind(lma.on, lma.off)
names(lma)<-c("est","se","tval","pval","model")
lma<-lma %>% dplyr::mutate(sig=ifelse(pval<0.05,1,0))

ggplot(data=lma, aes(x=model,y=est,color=as.factor(sig))) + geom_boxplot() + 
  geom_jitter(data=filter(lma, sig==0)) + 
  scale_color_manual(values=c("grey","black"))

fd3<-Filter.data1 %>% group_by(name2, rndLat, year) %>% tally() 
(p2 <- ggplot(fd3, aes(rndLat,name2)) + geom_tile(aes(fill = log(n, base=10)),colour = "black") +
    scale_fill_gradient(low = "white", high="steelblue"))
ggplot(pheno.figdata, aes(x=as.factor(rndLat), y=log(n.occ))) + geom_boxplot()

ggplot(fd3, aes(x=rndLat, y=log(n.occ), color=name2)) + geom_line() + 
  scale_color_viridis_d(guide=F) + 
  geom_smooth(method="lm", se=F) + theme_classic()

ggplot(fd3, aes(x=rndLat, y=log(n.occ), color=name2)) + geom_point() + 
  scale_color_viridis_d(guide=F) + 
  geom_smooth(method="lm", se=F) + theme_classic()

ggplot(data=Filter.doy2, aes(color=name2)) + geom_histogram()

