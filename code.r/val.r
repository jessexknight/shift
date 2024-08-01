source('meta.r')

# =============================================================================
# config

uid = Sys.Date()

key.vars = c('age','vio.n','dep.now','dep.past','haz.now','haz.past','ptr.n','ptr.tot')

vals = list(
  # base rates
  'Ri.all'=list(save=NULL,vars=key.vars),
  # RR age
  'aRR.vio'=list(save=c('aRR.vio'),  vars=c('vio.n1y'),  strat='age.10'),
  'aRR.dep'=list(save=c('aRR.dep_o'),vars=c('dep_o.a1y'),strat='age.10'),
  'aRR.haz'=list(save=c('aRR.haz_o'),vars=c('haz_o.a1y'),strat='age.10'),
  'aRR.ptr'=list(save=c('aRR.ptr_o'),vars=c('ptr_o.n1y'),strat='age.10'),
  # basic RR
  'RR.dep_o.dep_p'=list(save='RR.dep_o.dep_p',vars='dep_o.a1y',strat='dep.past'),
  'RR.haz_o.haz_p'=list(save='RR.haz_o.haz_p',vars='haz_o.a1y',strat='haz.past'),
  'RR.haz_o.dep_w'=list(save='RR.haz_o.dep_w',vars='haz_o.a1y',strat='dep.now'),
  'RR.haz_x.dep_w'=list(save='RR.haz_x.dep_w',vars='haz_x.a1y',strat='dep.now'),
  'RR.ptr_o.dep_w'=list(save='RR.ptr_o.dep_w',vars='ptr_o.n1y',strat='dep.now'),
  'RR.ptr_o.haz_w'=list(save='RR.ptr_o.haz_w',vars='ptr_o.n1y',strat='haz.now'),
  'RR.ptr_x.dep_w'=list(save='RR.ptr_x.dep_w',vars='ptr_x.n1y',strat='dep.now'),
  'RR.ptr_x.haz_w'=list(save='RR.ptr_x.haz_w',vars='ptr_x.n1y',strat='haz.now'),
  # transient RR
  'tRR.dep_o.vio_z'=list(save=c('iRR.dep_o.vio_z','tsc.dep_o.vio_z'),vars='dep_o.a1y',strat='vio.a1y'),
  'tRR.dep_x.vio_z'=list(save=c('iRR.dep_x.vio_z','tsc.dep_x.vio_z'),vars='dep_x.a1y',strat='vio.a1y'),
  'tRR.haz_o.vio_z'=list(save=c('iRR.haz_o.vio_z','tsc.haz_o.vio_z'),vars='haz_o.a1y',strat='vio.a1y'),
  'tRR.haz_x.vio_z'=list(save=c('iRR.haz_x.vio_z','tsc.haz_x.vio_z'),vars='haz_x.a1y',strat='vio.a1y'),
  'tRR.ptr_o.vio_z'=list(save=c('iRR.ptr_o.vio_z','tsc.ptr_o.vio_z'),vars='ptr_o.n1y',strat='vio.a1y'),
  # cumulative RR
  'nRR.dep_o.vio_n'=list(save=c('mRR.dep_o.vio_n','nsc.dep_o.vio_n'),vars='dep_o.a1y',strat='vio.n'),
  'nRR.haz_o.vio_n'=list(save=c('mRR.haz_o.vio_n','nsc.haz_o.vio_n'),vars='haz_o.a1y',strat='vio.n'),
  'nRR.ptr_o.vio_n'=list(save=c('mRR.ptr_o.vio_n','nsc.ptr_o.vio_n'),vars='ptr_o.n1y',strat='vio.n'),
  # duration RR
  'dRR.dep_x.dep_u'=list(save=c('dsc.dep_x.dep_u'),vars='dep_x.a1y',strat='dep.u'),
  'dRR.haz_x.haz_u'=list(save=c('dsc.haz_x.haz_u'),vars='haz_x.a1y',strat='haz.u'),
  # full null
  'null'=list(save=NULL,vars=key.vars,'Ri\\.m$'=0)
)
for (v in names(vals)){ vals[[v]]$name = v }

# =============================================================================
# run & plot

val.run = function(name,vars,strat='.',...){
  Ps = lapply(1:7,get.pars,n=333,
    null=ulist('Ri\\.m$'=NULL,...))
  Is = sim.runs(Ps)
  Is = Is[Is$age<amax,]
  g = val.plot(Is,vars,strat)
  plot.save('val',uid,name,h=3,w=1+3*len(vars))
}

val.plot = function(Is,vars,strat='.'){
  # plot the densities for multiple (7) seeds using:
  # - boxplot if var is binary
  # - line+ribbon otherwise
  # pre-compute group-wise densities b/c no ggplot support
  g  = c('seed',strat) # grouping variables
  Is = cbind(Is,.='')[c(g,vars)]
  if (ulen(Is[[strat]]) > 5){
    Is[[strat]] = q.cut(Is[[strat]],0:5/5) }
  Im = rbind.lapply(vars,function(var){
    x  = as.numeric(Is[[var]]) # extract data
    br = hist(x,br=min(31,ulen(x)))$br # compute breaks
    Imx = aggregate(x,Is[g],function(xi){ # for each group
      x = sum1(hist(xi,br=br,plot=FALSE)$count) }) # compute density
    if (len(br) > 3){ # continuous
      Imv = cbind(Imx[g],var=var,d.cts=c(Imx$x),d.bin=NA,b=rep(br[-len(br)],each=nrow(Imx))) }
    else { # binary
      Imv = cbind(Imx[g],var=var,d.bin=Imx$x[,2],d.cts=NA,b=1) }
  })
  g = ggplot(Im,aes(x=b,y=as.numeric(d.cts),
      color = as.factor(.data[[strat]]),
      fill  = as.factor(.data[[strat]]))) +
    facet_wrap('~var',scales='free',ncol=len(vars)) +
    stat_summary(geom='ribbon',fun.min=min,fun.max=max,alpha=.3,color=NA) +
    stat_summary(geom='line',fun=median) +
    geom_boxplot(aes(y=as.numeric(d.bin),group=interaction(b,.data[[strat]])),
      alpha=.3,outlier.alpha=1,outlier.shape=3) +
    labs(x='value',y='density',color=strat,fill=strat) +
    scale_x_continuous(expand=c(.1,.1)) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    ylim(c(0,NA))
  g = plot.clean(g)
}

# =============================================================================
# main

for (val in vals){
  do.call(val.run,val) }