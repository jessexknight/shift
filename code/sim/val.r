source('sim/meta.r')

# =============================================================================
# config

v      = cli.arg('v',     1:23)
n.pop  = cli.arg('n.pop', 333)
n.seed = cli.arg('n.seed',7)
key.vars = c('age',
  'vio.nt',
  'dep.now','dep.past',
  'haz.now','haz.past',
  'ptr.nw','ptr.nt')

# -----------------------------------------------------------------------------
# RR scenarios

RR3 = c(1,3,9)
RR2 = c(1,9)
ts2 = c(30,180)
ns2 = c(10,100)

val.RR = list(
  # RR age
  aRR.vio=list(save='aRR.vio',  vars='vio.n1y',  strat='age.10'),
  aRR.dep=list(save='aRR.dep_o',vars='dep_o.a1y',strat='age.10'),
  aRR.haz=list(save='aRR.haz_o',vars='haz_o.a1y',strat='age.10'),
  aRR.ptr=list(save='aRR.ptr_o',vars='ptr_o.n1y',strat='age.10',among=quote(p1y.sex.act)),
  # basic RR
  RR.dep_o.dep_p=list(gpar=list('RR.dep_o.dep_p'=RR3),vars='dep_o.a1y',strat='p1y.dep.past',among=quote(!p1y.dep.now)),
  RR.haz_o.haz_p=list(gpar=list('RR.haz_o.haz_p'=RR3),vars='haz_o.a1y',strat='p1y.haz.past',among=quote(!p1y.haz.now)),
  RR.haz_o.dep_w=list(gpar=list('RR.haz_o.dep_w'=RR3),vars='haz_o.a1y',strat='p1y.dep.now', among=quote(!p1y.haz.now)),
  RR.haz_x.dep_w=list(gpar=list('RR.haz_x.dep_w'=RR3),vars='haz_x.a1y',strat='p1y.dep.now', among=quote( p1y.haz.now)),
  RR.ptr_o.dep_w=list(gpar=list('RR.ptr_o.dep_w'=RR3),vars='ptr_o.n1y',strat='p1y.dep.now', among=quote( p1y.sex.act)),
  RR.ptr_o.haz_w=list(gpar=list('RR.ptr_o.haz_w'=RR3),vars='ptr_o.n1y',strat='p1y.haz.now', among=quote( p1y.sex.act)),
  RR.ptr_x.dep_w=list(gpar=list('RR.ptr_x.dep_w'=RR3),vars='ptr_x.n1y',strat='p1y.dep.now', among=quote( p1y.sex.act)),
  RR.ptr_x.haz_w=list(gpar=list('RR.ptr_x.haz_w'=RR3),vars='ptr_x.n1y',strat='p1y.haz.now', among=quote( p1y.sex.act)),
  # transient RR
  tRR.dep_o.vio_zf=list(gpar=list('iRR.dep_o.vio_zf'=RR2,'tsc.dep_o.vio_zf'=ts2),vars='dep_o.a3m',strat='vio.a3m',among=quote(!p3m.dep.now)),
  tRR.dep_x.vio_zf=list(gpar=list('iRR.dep_x.vio_zf'=RR2,'tsc.dep_x.vio_zf'=ts2),vars='dep_x.a3m',strat='vio.a3m',among=quote( p3m.dep.now)),
  tRR.haz_o.vio_zf=list(gpar=list('iRR.haz_o.vio_zf'=RR2,'tsc.haz_o.vio_zf'=ts2),vars='haz_o.a3m',strat='vio.a3m',among=quote(!p3m.haz.now)),
  tRR.haz_x.vio_zf=list(gpar=list('iRR.haz_x.vio_zf'=RR2,'tsc.haz_x.vio_zf'=ts2),vars='haz_x.a3m',strat='vio.a3m',among=quote( p3m.haz.now)),
  tRR.ptr_o.vio_zf=list(gpar=list('iRR.ptr_o.vio_zf'=RR2,'tsc.ptr_o.vio_zf'=ts2),vars='ptr_o.n3m',strat='vio.a3m',among=quote( p3m.sex.act)),
  # cumulative RR
  nRR.dep_o.vio_nt=list(gpar=list('mRR.dep_o.vio_nt'=RR2,'nsc.dep_o.vio_nt'=ns2),vars='dep_o.a1y',strat='p1y.vio.nt.c',among=quote(!p1y.dep.now)),
  nRR.haz_o.vio_nt=list(gpar=list('mRR.haz_o.vio_nt'=RR2,'nsc.haz_o.vio_nt'=ns2),vars='haz_o.a1y',strat='p1y.vio.nt.c',among=quote(!p1y.haz.now)),
  nRR.ptr_o.vio_nt=list(gpar=list('mRR.ptr_o.vio_nt'=RR2,'nsc.ptr_o.vio_nt'=ns2),vars='ptr_o.n1y',strat='p1y.vio.nt.c',among=quote( p1y.sex.act)),
  # duration RR
  dRR.dep_x.dep_u=list(gpar=list('dsc.dep_x.dep_u'=ts2),vars='dep_x.a1y',strat='p1y.dep.dur.c',among=quote(p1y.dep.now)),
  dRR.haz_x.haz_u=list(gpar=list('dsc.haz_x.haz_u'=ts2),vars='haz_x.a1y',strat='p1y.haz.dur.c',among=quote(p1y.haz.now))
)
for (name in names(val.RR)){
  val.RR[[name]]$null = ulist('Ri\\.m$'=NULL,save=val.RR[[name]]$save)
  val.RR[[name]]$srvs = c(srv.val.RR)
  val.RR[[name]]$name = name
}

# -----------------------------------------------------------------------------
# base scenarios

val.base = list(
  'base'=list(vars=key.vars)
)
for (name in names(val.base)){
  val.base[[name]]$srvs = c(srv.base,srv.ptr)
  val.base[[name]]$name = name
}

# =============================================================================
# run & plot

val.run = function(name,vars,strat='.',among=quote(TRUE),srvs=NULL,gpar=list(case='base'),...){
  status(2,'val.run: ',name,' @ ',n.seed*prod(lens(gpar)))
  Ps = grid.apply(c(list(seed=1:n.seed),gpar),get.pars,n.pop=n.pop,...)
  Q = srv.apply(sim.runs(Ps),srvs=srvs,p.vars=names(gpar))
  Q = subset(Q,age < amax & eval(among))
  for (var in vars){
    g = val.plot(Q,var,strat,names(gpar))
    plot.save('val',uid,str(name,'--',var),h=3,w=1+3*prod(lens(gpar)))
  }
}

val.plot = function(Q,var,strat,gpar){
  g = c('seed','gpar',strat) # grouping variables
  Q = cbind(Q,gpar=apply(Q[gpar],1,list.str),.='')[c(g,var)]
  x = as.numeric(Q[[var]]) # extract values
  b = breaks(x) # compute common breaks
  cts = ulen(x) > 2 # is var continuous or binary
  Qx = aggregate(x,Q[g],function(xg){ # for each group
    x = sum1(hist(xg,breaks=b,right=FALSE,plot=FALSE)$count) }) # compute density
  if (cts){ Qp = cbind(p=c(Qx$x),b=rep(b[-len(b)],each=nrow(Qx))) } # continuous
  else    { Qp = cbind(p=Qx$x[,2],b=1) } # binary
  Qp = cbind(Qx[g],Qp)
  g = ggplot(Qp,aes(x=b,y=100*as.numeric(p),
      color = as.factor(.data[[strat]]),
      fill  = as.factor(.data[[strat]]))) +
    facet_wrap('~gpar',scales='free',ncol=ulen(Q$gpar)) +
    labs(y='proportion (%)',color=strat,fill=strat) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    ylim(c(0,NA)) +
    ggtitle(var)
  if (cts){ g = plot.clean(g) + xlab('value') +
    stat_summary(geom='ribbon',fun.min=min,fun.max=max,alpha=.3,color=NA) +
    stat_summary(geom='line',fun=median) }
  else { g = plot.clean(g,axis.text.x=element_blank()) + xlab('') +
    geom_boxplot(aes(group=interaction(b,gpar,.data[[strat]])),
      alpha=.3,outlier.alpha=1,outlier.shape=3) }
}

# =============================================================================
# main

vals = c(val.RR,val.base)[v] # 1:22,23:23
for (val in vals){ do.call(val.run,val,quote=TRUE) }
