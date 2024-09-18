source('sim/meta.r')

# =============================================================================
# config

v      = cli.arg('v',     1:23)
n.pop  = cli.arg('n.pop', 333)
n.seed = cli.arg('n.seed',7)

# -----------------------------------------------------------------------------

P0 = list(
  n.pop = n.pop,
  vio.Ri.m   = .001,
  dep_o.Ri.m = .001, dep_x.Ri.m = .002,
  haz_o.Ri.m = .001, haz_x.Ri.m = .002,
  ptr_o.Ri.m = .01,  ptr_x.Ri.m = .02,  ptr.max.m = 3,
  all.Ri.shape = 1e6,
  ptr.Ri.shape = 1e6,
  null = 'xRR')

# =============================================================================
# scenarios

val.Ri = list(
  Ri.m=list(pars=ulist(P0),evts=evts[1:7])
)

# -----------------------------------------------------------------------------

RR1 = c(3)
RR2 = c(1,3)
RR3 = c(1,2,3)
ts2 = c(35,70)
ds2 = c(70,700)
ns2 = c(3,10)
age5 = c(10,20,30,40,50)
aRR5 = c( 1, 5, 2, 4, 3)/3

age.10.fun   = function(Q){ floor(Q$age.1/10)*10 }
vio.nt.c.fun = function(Q){ int.cut(Q$vio.nt,c(0,ns2)) }

val.RR = list(
  # RR age [1:4]
  aRR.vio=list(pars=ulist(P0,aRR.vio.RRs  =aRR5,aRR.vio.ages  =age5,aRR.shape='const'),evts='vio',  strat='age.10',x.cols=list(age.10=age.10.fun)),
  aRR.dep=list(pars=ulist(P0,aRR.dep_o.RRs=aRR5,aRR.dep_o.ages=age5,aRR.shape='const'),evts='dep_o',strat='age.10',x.cols=list(age.10=age.10.fun)),
  aRR.haz=list(pars=ulist(P0,aRR.haz_o.RRs=aRR5,aRR.haz_o.ages=age5,aRR.shape='const'),evts='haz_o',strat='age.10',x.cols=list(age.10=age.10.fun)),
  aRR.ptr=list(pars=ulist(P0,aRR.ptr_o.RRs=aRR5,aRR.ptr_o.ages=age5,aRR.shape='const'),evts='ptr_o',strat='age.10',x.cols=list(age.10=age.10.fun)),
  # basic RR [5:12]
  RR.dep_o.dep_p=list(pars=ulist(P0,RR.dep_o.dep_p=  RR3),evts='dep_o',strat='dep.past'),
  RR.haz_o.haz_p=list(pars=ulist(P0,RR.haz_o.haz_p=  RR3),evts='haz_o',strat='haz.past'),
  RR.haz_o.dep_w=list(pars=ulist(P0,RR.haz_o.dep_w=  RR3),evts='haz_o',strat='dep.now'),
  RR.haz_x.dep_w=list(pars=ulist(P0,RR.haz_x.dep_w=1/RR3),evts='haz_x',strat='dep.now'),
  RR.ptr_o.dep_w=list(pars=ulist(P0,RR.ptr_o.dep_w=1/RR3),evts='ptr_o',strat='dep.now'),
  RR.ptr_o.haz_w=list(pars=ulist(P0,RR.ptr_o.haz_w=  RR3),evts='ptr_o',strat='haz.now'),
  RR.ptr_x.dep_w=list(pars=ulist(P0,RR.ptr_x.dep_w=  RR3),evts='ptr_x',strat='dep.now'),
  RR.ptr_x.haz_w=list(pars=ulist(P0,RR.ptr_x.haz_w=  RR3),evts='ptr_x',strat='haz.now'),
  # transient RR [13:17]
  tRR.dep_o.vio_zr=list(pars=ulist(P0,iRR.dep_o.vio_zr=  RR2,tsc.dep_o.vio_zr=ts2,tRR.shape='step'),evts='dep_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.dep_x.vio_zr=list(pars=ulist(P0,iRR.dep_x.vio_zr=1/RR2,tsc.dep_x.vio_zr=ts2,tRR.shape='step'),evts='dep_x',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.haz_o.vio_zr=list(pars=ulist(P0,iRR.haz_o.vio_zr=  RR2,tsc.haz_o.vio_zr=ts2,tRR.shape='step'),evts='haz_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.haz_x.vio_zr=list(pars=ulist(P0,iRR.haz_x.vio_zr=1/RR2,tsc.haz_x.vio_zr=ts2,tRR.shape='step'),evts='haz_x',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.ptr_o.vio_zr=list(pars=ulist(P0,iRR.ptr_o.vio_zr=  RR2,tsc.ptr_o.vio_zr=ts2,tRR.shape='step'),evts='ptr_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  # cumulative RR [18:20]
  nRR.dep_o.vio_nt=list(pars=ulist(P0,mRR.dep_o.vio_nt=RR2,nsc.dep_o.vio_nt=ns2,nRR.shape='step'),evts='dep_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.fun)),
  nRR.haz_o.vio_nt=list(pars=ulist(P0,mRR.haz_o.vio_nt=RR2,nsc.haz_o.vio_nt=ns2,nRR.shape='step'),evts='haz_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.fun)),
  nRR.ptr_o.vio_nt=list(pars=ulist(P0,mRR.ptr_o.vio_nt=RR2,nsc.ptr_o.vio_nt=ns2,nRR.shape='step'),evts='ptr_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.fun)),
  # duration RR [21:22]
  dRR.dep_x.dep_u=list(pars=ulist(P0,dsc.dep_x.dep_u=ds2,dRR.shape='step'),evts='dep_x',strat='dep_o.dt.c',e.dts=list(dep_o=ds2)),
  dRR.haz_x.haz_u=list(pars=ulist(P0,dsc.haz_x.haz_u=ds2,dRR.shape='step'),evts='haz_x',strat='haz_o.dt.c',e.dts=list(haz_o=ds2))
)

# =============================================================================
# run & plot

val.run = function(name,pars,evts,strat='case',e.dts=NULL,x.cols=NULL){
  pars = val.par.split(pars)
  status(2,'val.run: ',name,' @ ',n.seed*pars$n.var)
  Ps = grid.apply(ulist(pars$var,seed=1:n.seed),get.pars,pars$fix,.par=FALSE)
  Ms = sim.runs(Ps)
  Y = rate.datas(Ms,p.vars=names(pars$var),e.dts=e.dts,x.cols=x.cols)
  R = rbind.lapply(evts,rate.est,Y=Y,strat=c('seed',strat,names(pars$var)))
  Q = srv.apply(Ms,srvs=c(srv.base,def.args(srv.e.dts,e.dts=e.dts)),p.vars=names(pars$var),x.cols=x.cols)
  for (evt in evts){
    val.plot.rate(name,R,evt,pars,strat)
    val.plot.prev(name,Q,evt,pars,strat)
  }
}

val.par.split = function(par){
  b = lens(par) > 1 & !grepl(names(null.sets$aRR),names(par))
  p = list(var=par[b],fix=par[!b],n.var=prod(lens(par[b])))
}

val.plot.rate = function(name,R,evt,pars,strat,t1y=364){
  R = subset(R,variable==evt)
  R$facet = apply(R[names(pars$var)],1,list.str,sig=3)
  R.ref = pars$fix[[str(evt,'.Ri.m')]]
  g = ggplot(R,aes(x='',y=t1y*value,color=as.factor(.data[[strat]]))) +
    geom_point(data=data.frame(value=R.ref),shape=9,color='red') +
    ylab('Rate (per year)')
  g = val.plot.label(g,R,c('facet',strat),value=max(R$value)+R.ref,
    function(Ri){ signif(median(Ri$value/R.ref),3) })
  g = val.plot.label(g,R,c('facet',strat),value=0,
    function(Ri){ str(rmed(Ri$ne),'\n',rmed(Ri$dt/t1y),'\n') })
  g = val.plot.finish(g,c(name,evt,'rate'),pars,strat)
}

val.plot.prev = function(name,Q,evt,pars,strat){
  vars = switch(substr(evt,1,3),
    vio = c('vio.nt', 'vio.dt'),
    dep = c('dep.now','dep.past'),
    haz = c('haz.now','haz.past'),
    ptr = c('ptr.nt', 'ptr.nw'))
  Q$facet = apply(Q[names(pars$var)],1,list.str,sig=3)
  Q = melt(Q,measure=setdiff(vars,strat))
  Q = maggregate(formula(str('value~seed+variable+facet+',strat)),Q,
    function(x){ c(p=mean(x),s=sum(x),n=len(x)) })
  g = ggplot(Q,aes(x='',y=value.p,color=as.factor(.data[[strat]]))) +
    ylab('Value (population mean)')
  g = val.plot.label(g,Q,c('variable','facet',strat),value.p=0,
    function(Qi){ str(rmed(Qi$value.s),'\n',rmed(Qi$value.n),'\n') })
  g = val.plot.finish(g,c(name,evt,'prev'),pars,strat,nr=ulen(Q$variable))
}

val.plot.label = function(g,X,strat,lab.fun,...){
  X.lab = rbind.lapply(split(X,X[strat]),function(Xi){
    cbind(Xi[1,strat],label=lab.fun(Xi),...) })
  g = g + geom_text(data=X.lab,aes(label=label),size=2.5,
    position=position_dodge(width=.75),show.legend=FALSE)
}

val.plot.finish = function(g,name,pars,strat,nrow=1){
  g = g + geom_boxplot(outlier.shape=1,outlier.alpha=1) +
    facet_grid('variable~facet',scales='free') +
    guides(color=guide_legend(ncol=2)) +
    labs(x='',color=str(list.str(pars$fix,sig=3),'\n\n',strat)) +
    ggtitle(str(name,collapse=': ')) +
    scale_color_viridis_d() +
    ylim(c(0,NA))
  g = plot.clean(g,legend.title=element_text(size=9))
  plot.save('val',uid,str(name,collapse='--'),w=2.5+2*pars$n.var,h=1+2*nrow)
}

rmed = function(x){ round(median(x)) }

# =============================================================================
# main

vals = c(val.Ri,val.RR)[v]
for (name in names(vals)){ vals[[name]]$name = name }
for (val in vals){ do.call(val.run,val,quote=TRUE) }
