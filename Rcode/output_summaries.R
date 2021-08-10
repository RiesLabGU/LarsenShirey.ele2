#Model outputs
library(tidyverse)

load("output/term_model_results.RData")
load("output/onset_model_results.RData")

sp.responses<-bind_rows(ev.on.spresp,we.on.spresp,ev.term.spresp,we.term.spresp)

sp.resp.2<-sp.responses %>%
  group_by(spname, agg, cur, metric) %>%
  mutate(group=paste(spname, agg, cur, metric, sep=".")) %>%
  pivot_wider(names_from=effort, values_from=resp) %>%
  mutate(same=ifelse(`yes`==`no`,1,0), diff=yes-no)

sp.resp.agg<-sp.responses %>%
  group_by(spname, cur, metric, effort) %>%
  mutate(group=paste(spname, cur, metric, effort, sep=".")) %>%
  pivot_wider(names_from=agg, values_from=resp) %>%
  mutate(same=ifelse(`Lat*yr`==`Lat*yr*elev`,1,0), diff=`Lat*yr*elev`-`Lat*yr`)

sp.resp.est<-sp.responses %>%
  group_by(spname, cur, agg, effort) %>%
  mutate(group=paste(spname, cur, agg, effort, sep=".")) %>%
  pivot_wider(names_from=metric, values_from=resp) %>%
  mutate(on.same=ifelse(`ev.on`==`we.on`,1,0), off.same=ifelse(`ev.off`==`we.off`,1,0))

(on.est1<-sp.resp.est %>%
  group_by(agg, on.same) %>%
  tally()  %>%
    filter(!is.na(on.same)) %>%
  group_by(agg) %>%
  mutate(nagg=sum(n), prop=n/nagg))

(off.est1<-sp.resp.est %>%
    group_by(agg, off.same) %>%
    tally()  %>%
    filter(!is.na(off.same)) %>%
    group_by(agg) %>%
    mutate(nagg=sum(n), prop=n/nagg))




sp.r.1<-sp.resp.2 %>%
  filter(cur==2) %>%
  group_by(metric, agg) %>%
  summarize(p.changed=(1-sum(same, na.rm=T)/length(same)))

sp.r.2<-sp.resp.2 %>%
  filter(same==0) %>%
  group_by(metric, agg, same,yes,no) %>%
  tally()


table(sp.resp.2$agg,sp.resp.2$`same`, sp.resp.2$metric)
table(sp.resp.2$agg,sp.resp.2$yes, sp.resp.2$metric)
table(sp.resp.2$agg,sp.resp.2$no, sp.resp.2$metric)


