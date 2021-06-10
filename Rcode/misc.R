#Unclear if we will need any of this


```
<h2> Sample effort vs. latitude model </h2>
  <p> This function uses single regression to model the correlation between latitudinal band and sampling effort metric for a given dataset, and returns the phenometric, effort metric, and model results: coefficient sign, coefficient p-value, and model rsquared.</p>
  ```{r function for lat-effort correlation, echo=T}

lat.effort = function(df1,       #phenometric dataset     
                      emetric) { #effort metric to use
  eff.model<-summary(lm(as.formula(paste(emetric,"~1+rndLat",sep="")), data=df1))
  pmetric<-df1$metric[1]
  pval<-eff.model$coefficients["rndLat",'Pr(>|t|)']
  sign<-sign(eff.model$coefficients["rndLat",'Estimate'])
  rsq<-round(as.numeric(eff.model$adj.r.squared),3)
  return(c(pmetric,emetric,sign,pval,rsq))
} 


```

```{r analytical functions}


