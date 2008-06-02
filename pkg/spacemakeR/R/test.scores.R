"test.scores" <-
function(scores,listw,nsim=999) {
alter<-ifelse(scores$values>0,"greater","less")
vec2test<-1:dim(scores$vectors)[2]
moran.w <- apply(cbind(vec2test),1, function(x) moran.mc(scores$vectors[,x],listw=listw,nsim=nsim,alternative=alter[x]))
stat<-as.numeric(unlist(lapply(moran.w,function(x) x$statistic)))
pval<-as.numeric(unlist(lapply(moran.w,function(x) x$p.value)))
res<-cbind(stat,pval)
names(res)<-c("Moran","pval")
row.names(res)<-vec2test
return(data.frame(res))
}