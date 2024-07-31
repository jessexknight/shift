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
  Ism = melt(cbind(Is,.=''),measure=vars,id=c('seed','i',strat))
  bins = round(min(len(unique(Ism$value)),sqrt(len(unique(Is$i)))))
  g = ggplot(Ism,aes(
      x = as.numeric(value),
      y = after_stat(density),
      color = as.factor(.data[[strat]]))) +
    facet_wrap('~variable',scales='free',ncol=len(vars)) +
    stat_bin(bins=bins,pad=FALSE,position='identity',geom='path') +
    labs(x='value',y='density',color=strat) +
    scale_color_viridis_d()
  g = plot.clean(g)
}

# =============================================================================
# main

for (val in vals){
  do.call(val.run,val) }
