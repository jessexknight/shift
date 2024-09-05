source('utils.r')
source('sim/plot.r')
library('ggplot2')
library('reshape2')
# -----------------------------------------------------------------------------
# config
x = list(
  dtz=c(1,3,7,15),
  tsc=c(7,30,90))
iRR  = 3
tmax = 360
tRR.exp = function(t,iRR,tsc,dtz){
  RR = 1+(iRR-1)*exp(-t/tsc)
}
tRR.ramp = function(t,iRR,tsc,dtz){
  n  = tsc/dtz
  RR = seq(iRR,1,(1-iRR)/(n+1))[-1]
  RR = c(RR,rep(1,len(t)-len(RR)))
}
tRR.rect = function(t,iRR,tsc,dtz){
  n  = tsc/dtz
  RR = 1+(iRR-1)*c(rep(1,floor(n)),n-floor(n))
  RR = c(RR,rep(1,len(t)-len(RR)))
}
# -----------------------------------------------------------------------------
# run & plot
for (type in c('rect','ramp','exp')){
  tRR.fun = get(str('tRR.',type))
  # generate tRR vectors
  tRR = grid.apply(x,function(dtz,tsc){
    t = seq(dtz/2,tmax,dtz)
    RR = tRR.fun(t,iRR,tsc,dtz)
    cbind(dtz=dtz,tsc=tsc,t=t,RR=RR,cRR=1+cumsum(RR-1)*dtz)
  },.par=FALSE)
  tRR = melt(as.data.frame(do.call(rbind,tRR)),m=c('RR','cRR'))
  # compute cum RR ratios vs dtz = 1
  lab = aggregate(value~dtz+tsc,subset(tRR,variable='cRR'),last)
  lab = cbind(lab[order(lab$dtz),],variable='cRR')
  lab$y = max(lab$value)*as.numeric(factor(lab$dtz))/(1+len(x$dtz))
  lab$ratio = signif(lab$value / lab$value[lab$dtz==1],3)
  # -----------------------------------------------------------------------------
  # plot
  g = ggplot(tRR,aes(x=t,y=value,color=factor(dtz))) +
    facet_grid(variable~tsc,scales='free',labeller='label_both') +
    geom_text(data=lab,aes(x=tmax,y=y,label=ratio),
      size=3,hjust=1,show.legend=FALSE) +
    scale_color_viridis_d() +
    geom_step(alpha=.5) +
    geom_point(size=.1) +
    labs(x='time (days)',color='timestep\n(days)')
  g = plot.clean(g)
  plot.save('toy','tRR',str('tRR.',type),w=8,h=4)
}
