source('sim/meta.r')

# =============================================================================
# config

v      = cli.arg('v',     1:24)
n.pop  = cli.arg('n.pop', 333)
n.seed = cli.arg('n.seed',7)
key.vars = c('age',
  'vio.nt',
  'dep.now','dep.past',
  'haz.now','haz.past',
  'ptr.nw','ptr.nt')

# -----------------------------------------------------------------------------
# RR scenarios

RR1 = c(8)
RR2 = c(1,8)
RR3 = c(1,4,8)
ts2 = c(30,90)
ds2 = c(90,720)
ns2 = c(10,100)
gvio = list(vio.Ri.m   = .002)
gdep = list(dep_o.Ri.m = .002, dep_x.Ri.m = .001)
ghaz = list(haz_o.Ri.m = .002, haz_x.Ri.m = .001)
gptr = list(ptr_o.Ri.m = .01,  ptr_x.Ri.m = .02,  ptr.max.m = 3)

val.RR = list(
  # RR age [1:4]
  aRR.vio=list(save='aRR.vio',  vars='vio.n1y',  strat='age.10'),
  aRR.dep=list(save='aRR.dep_o',vars='dep_o.a1y',strat='age.10'),
  aRR.haz=list(save='aRR.haz_o',vars='haz_o.a1y',strat='age.10'),
  aRR.ptr=list(save='aRR.ptr_o',vars='ptr_o.n1y',strat='age.10',among=quote(p1y.sex.act)),
  # basic RR [5:12]
  RR.dep_o.dep_p=list(gpar=ulist(gdep,     RR.dep_o.dep_p=  RR3),vars='dep_o.a1y',strat='p1y.dep.past',among=quote(!p1y.dep.now)),
  RR.haz_o.haz_p=list(gpar=ulist(ghaz,     RR.haz_o.haz_p=  RR3),vars='haz_o.a1y',strat='p1y.haz.past',among=quote(!p1y.haz.now)),
  RR.haz_o.dep_w=list(gpar=ulist(ghaz,gdep,RR.haz_o.dep_w=  RR3),vars='haz_o.a1y',strat='p1y.dep.now', among=quote(!p1y.haz.now)),
  RR.haz_x.dep_w=list(gpar=ulist(ghaz,gdep,RR.haz_x.dep_w=1/RR3),vars='haz_x.a1y',strat='p1y.dep.now', among=quote( p1y.haz.now)),
  RR.ptr_o.dep_w=list(gpar=ulist(gptr,gdep,RR.ptr_o.dep_w=1/RR3),vars='ptr_o.n1y',strat='p1y.dep.now', among=quote( p1y.sex.act)),
  RR.ptr_o.haz_w=list(gpar=ulist(gptr,ghaz,RR.ptr_o.haz_w=  RR3),vars='ptr_o.n1y',strat='p1y.haz.now', among=quote( p1y.sex.act)),
  RR.ptr_x.dep_w=list(gpar=ulist(gptr,gdep,RR.ptr_x.dep_w=  RR3),vars='ptr_x.n1y',strat='p1y.dep.now', among=quote( p1y.sex.act)),
  RR.ptr_x.haz_w=list(gpar=ulist(gptr,ghaz,RR.ptr_x.haz_w=  RR3),vars='ptr_x.n1y',strat='p1y.haz.now', among=quote( p1y.sex.act)),
  # transient RR [13:17]
  tRR.dep_o.vio_zr=list(gpar=ulist(gdep,gvio,iRR.dep_o.vio_zr=  RR2,tsc.dep_o.vio_zr=ts2),vars='dep_o.a3m',strat='vio.a3m',among=quote(!p3m.dep.now)),
  tRR.dep_x.vio_zr=list(gpar=ulist(gdep,gvio,iRR.dep_x.vio_zr=1/RR2,tsc.dep_x.vio_zr=ts2),vars='dep_x.a3m',strat='vio.a3m',among=quote( p3m.dep.now)),
  tRR.haz_o.vio_zr=list(gpar=ulist(ghaz,gvio,iRR.haz_o.vio_zr=  RR2,tsc.haz_o.vio_zr=ts2),vars='haz_o.a3m',strat='vio.a3m',among=quote(!p3m.haz.now)),
  tRR.haz_x.vio_zr=list(gpar=ulist(ghaz,gvio,iRR.haz_x.vio_zr=1/RR2,tsc.haz_x.vio_zr=ts2),vars='haz_x.a3m',strat='vio.a3m',among=quote( p3m.haz.now)),
  tRR.ptr_o.vio_zr=list(gpar=ulist(gptr,gvio,iRR.ptr_o.vio_zr=  RR2,tsc.ptr_o.vio_zr=ts2),vars='ptr_o.n3m',strat='vio.a3m',among=quote( p3m.sex.act)),
  # cumulative RR [18:20]
  nRR.dep_o.vio_nt=list(gpar=ulist(gdep,gvio,mRR.dep_o.vio_nt=RR2,nsc.dep_o.vio_nt=ns2),vars='dep_o.a1y',strat='p1y.vio.nt.c',among=quote(!p1y.dep.now)),
  nRR.haz_o.vio_nt=list(gpar=ulist(ghaz,gvio,mRR.haz_o.vio_nt=RR2,nsc.haz_o.vio_nt=ns2),vars='haz_o.a1y',strat='p1y.vio.nt.c',among=quote(!p1y.haz.now)),
  nRR.ptr_o.vio_nt=list(gpar=ulist(gptr,gvio,mRR.ptr_o.vio_nt=RR2,nsc.ptr_o.vio_nt=ns2),vars='ptr_o.n1y',strat='p1y.vio.nt.c',among=quote( p1y.sex.act)),
  # duration RR [21:22]
  dRR.dep_x.dep_u=list(gpar=ulist(gdep,dsc.dep_x.dep_u=ds2),vars='dep_x.a1y',strat='p1y.dep.dur.c',among=quote(p1y.dep.now)),
  dRR.haz_x.haz_u=list(gpar=ulist(ghaz,dsc.haz_x.haz_u=ds2),vars='haz_x.a1y',strat='p1y.haz.dur.c',among=quote(p1y.haz.now))
)
for (name in names(val.RR)){
  val.RR[[name]]$all.Ri.shape = 1e6 # reduce heterogeneity
  val.RR[[name]]$null = ulist('Ri\\.m$'=NULL,save=val.RR[[name]]$save)
  val.RR[[name]]$srvs = c(srv.val.RR)
  val.RR[[name]]$name = name
}

# -----------------------------------------------------------------------------
# base scenarios

val.base = list(
  base=list(vars=key.vars),
  dtz=list(vars=key.vars,gpar=ulist(dtz=c(3,7,14)),strat='dtz')
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
  facet = setdiff(names(gpar)[lens(gpar)>1],strat)
  fixed = list.str(Ps[[1]][setdiff(filter.names(Ps[[1]],'Ri\\.m$'),c(facet,strat))],sig=3)
  for (var in vars){
    g = val.plot(Q,var,strat,facet,fixed,name)
    plot.save('val',uid,'.debug',str(name,'--',var),h=3,w=2+3*prod(lens(gpar[facet])))
  }
}

val.plot = function(Q,var,strat,facet,fixed,name){
  g = c('seed','facet',strat) # grouping variables
  Q = cbind(Q,facet=apply(Q[facet],1,list.str,sig=3),.='')[c(g,var)]
  x = as.numeric(Q[[var]]) # extract values
  b = breaks(x) # compute common breaks
  cts = ulen(x) > 2 # is var continuous or binary
  Qx = aggregate(x,Q[g],function(xg){ # for each group
    x = sum1(hist(xg,breaks=b,right=FALSE,plot=FALSE)$count) }) # compute density
  if (cts){ Qp = cbind(p=c(Qx$x),b=rep(b[-len(b)],each=nrow(Qx))) } # continuous
  else    { Qp = cbind(p=Qx$x[,2],b=1) } # binary
  Qp = cbind(Qx[g],Qp)
  slab = str(fixed,'\n\n',strat)
  g = ggplot(Qp,aes(x=b,y=100*as.numeric(p),
      color = as.factor(.data[[strat]]),
      fill  = as.factor(.data[[strat]]))) +
    facet_wrap('~facet',scales='fixed',ncol=ulen(Q$facet)) +
    labs(y='proportion (%)',x=var,color=slab,fill=slab) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    ylim(c(0,NA)) +
    ggtitle(name)
  if (cts){ g = plot.clean(g) +
    stat_summary(geom='ribbon',fun.min=min,fun.max=max,alpha=.3,color=NA) +
    stat_summary(geom='line',fun=median) }
  else { g = plot.clean(g,axis.text.x=element_blank()) +
    geom_violin(aes(group=interaction(b,facet,.data[[strat]])),
      alpha=.3,scale='width',bw=2,draw_quantiles=1:3/4) }
}

# =============================================================================
# main

vals = c(val.RR,val.base)[v] # 1:22,23:24
for (val in vals){ do.call(val.run,val,quote=TRUE) }
