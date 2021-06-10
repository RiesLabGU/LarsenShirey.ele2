#Larsen & Shirey
#Quantifying impacts of data curation, aggregation, phenometric estimation
#Project started 2021-03, Updated 2021-06
#HERE: Data curation

#load packages
library(tidyverse) #to be tidy

#FUNCTIONS
`%notin%` = function(x,y) !(x %in% y)

median.range<-function(x) {
  a<-median(x, na.rm=T)
  b<-min(x, na.rm=T)
  c<-max(x, na.rm=T)
  d<-paste(a," (",b,"-",c,")", sep="")
  return(d)
}

calc.sum<-function(df) {
  out1<-df %>%
    group_by(name, region, rndLat, year) %>%
    summarize(n.per.AO=n()) 
  out2<-out1 %>%
    group_by(name, region, rndLat) %>%
    summarize(n.per.lat=sum(n.per.AO), n.on.ev=1, n.on.we=n())
  out3<-out2 %>%
    group_by(name, region) %>%
    summarize(n.per.dset=sum(n.per.lat),nlat=n(), nmetric=sum(n.on.we))
    
  n.occ.tot<-nrow(df)
  nsp<-length(unique(df$name))
  nsp.reg<-length(unique(paste(df$name,df$region)))
  
  n.onset<-df %>%
    group_by(name, region, rndLat) %>%
    filter(doy==min(doy, na.rm=T)) %>%
    group_by(name, region) %>%
    tally()
  
  out.fin<-c(n.occ.tot, nsp, nsp.reg,median.range(out3$n.per.dset), 
             median.range(out2$n.per.lat), median.range(out1$n.per.AO),
             1, median.range(out2$n.on.we),median.range(out3$nlat),
             median.range(n.onset$n), median.range(out3$nmetric))
  return(out.fin)
}


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


#Datasets for further comparison
cured.data1<-Filter.sp
cured.data2<-Filter.doy2
save(cured.data1, cured.data2, file="data/curations.RData")

