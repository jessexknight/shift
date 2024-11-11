source('sim/meta.r')

# =============================================================================
# config

v = cli.arg('v',0)
n.pop = cli.arg('n.pop',333)
n.seed = cli.arg('n.seed',7)
fig.dir = 'val'

# -----------------------------------------------------------------------------

P0 = list(
  n.pop = n.pop,
  vio.Ri.my   = 1.0,
  dep_o.Ri.my = .10, dep_x.Ri.my = .50,
  haz_o.Ri.my = .20, haz_x.Ri.my = .50,
  ptr_o.Ri.my = 2.0, ptr_x.Ri.my = 1.0, ptr.max.m = 3,
  vio.Ri.cv   = .001,
  dep.Ri.cv   = .001,
  haz.Ri.cv   = .001,
  ptr.Ri.cv   = .001,
  null = 'xRR')

# =============================================================================
# scenarios

dtz = c(1,3,10)
RR1 = c(3)
RR2 = c(1,3)
RR3 = c(1,2,3)
ts2 = c(35,70)
ds2 = c(70,700)
ns2 = c(3,10)
age5 = c(10,20,30,40,50)
aRR5 = c( 1, 5, 2, 4, 3)/3

age.10.fun   = function(Q){ int.cut(Q$age.1,seq(10,50,10)) }
vio.nt.c.fun = function(Q){ int.cut(Q$vio.nt,c(0,ns2)) }

vals = list(
  # base rates
  Ri.m     =list(pars=ulist(P0),evts=evts[1:7]),
  ptr.n.pop=list(pars=ulist(P0,n.pop=c(100,300,900)),evts=c('ptr_o','ptr_x')),
  ptr.max  =list(pars=ulist(P0,ptr.max.m=c(1,3,30)), evts=c('ptr_o','ptr_x')),
  # age RR
  aRR.vio=list(pars=ulist(P0,aRR.vio.RRs  =aRR5,aRR.vio.ages  =age5,aRR.shape='const'),evts='vio',  strat='age.10',x.cols=list(age.10=age.10.fun)),
  aRR.dep=list(pars=ulist(P0,aRR.dep_o.RRs=aRR5,aRR.dep_o.ages=age5,aRR.shape='const'),evts='dep_o',strat='age.10',x.cols=list(age.10=age.10.fun)),
  aRR.haz=list(pars=ulist(P0,aRR.haz_o.RRs=aRR5,aRR.haz_o.ages=age5,aRR.shape='const'),evts='haz_o',strat='age.10',x.cols=list(age.10=age.10.fun)),
  aRR.ptr=list(pars=ulist(P0,aRR.ptr_o.RRs=aRR5,aRR.ptr_o.ages=age5,aRR.shape='const'),evts='ptr_o',strat='age.10',x.cols=list(age.10=age.10.fun)),
  # basic RR
  RR.dep_o.dep_p=list(pars=ulist(P0,RR.dep_o.dep_p=  RR3),evts='dep_o',strat='dep.past'),
  RR.haz_o.haz_p=list(pars=ulist(P0,RR.haz_o.haz_p=  RR3),evts='haz_o',strat='haz.past'),
  RR.haz_o.dep_w=list(pars=ulist(P0,RR.haz_o.dep_w=  RR3),evts='haz_o',strat='dep.now'),
  RR.haz_x.dep_w=list(pars=ulist(P0,RR.haz_x.dep_w=1/RR3),evts='haz_x',strat='dep.now'),
  RR.ptr_o.dep_w=list(pars=ulist(P0,RR.ptr_o.dep_w=1/RR3),evts='ptr_o',strat='dep.now'),
  RR.ptr_o.haz_w=list(pars=ulist(P0,RR.ptr_o.haz_w=  RR3),evts='ptr_o',strat='haz.now'),
  RR.ptr_x.dep_w=list(pars=ulist(P0,RR.ptr_x.dep_w=  RR3),evts='ptr_x',strat='dep.now'),
  RR.ptr_x.haz_w=list(pars=ulist(P0,RR.ptr_x.haz_w=  RR3),evts='ptr_x',strat='haz.now'),
  # transient RR
  tRR.dep_o.vio_zr=list(pars=ulist(P0,iRR.dep_o.vio_zr=  RR2,tsc.dep_o.vio_zr=ts2,tRR.shape='step'),evts='dep_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.dep_x.vio_zr=list(pars=ulist(P0,iRR.dep_x.vio_zr=1/RR2,tsc.dep_x.vio_zr=ts2,tRR.shape='step'),evts='dep_x',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.haz_o.vio_zr=list(pars=ulist(P0,iRR.haz_o.vio_zr=  RR2,tsc.haz_o.vio_zr=ts2,tRR.shape='step'),evts='haz_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.haz_x.vio_zr=list(pars=ulist(P0,iRR.haz_x.vio_zr=1/RR2,tsc.haz_x.vio_zr=ts2,tRR.shape='step'),evts='haz_x',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.ptr_o.vio_zr=list(pars=ulist(P0,iRR.ptr_o.vio_zr=  RR2,tsc.ptr_o.vio_zr=ts2,tRR.shape='step'),evts='ptr_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  # cumulative RR
  nRR.dep_o.vio_nt=list(pars=ulist(P0,mRR.dep_o.vio_nt=RR2,nsc.dep_o.vio_nt=ns2,nRR.shape='step'),evts='dep_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.fun)),
  nRR.haz_o.vio_nt=list(pars=ulist(P0,mRR.haz_o.vio_nt=RR2,nsc.haz_o.vio_nt=ns2,nRR.shape='step'),evts='haz_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.fun)),
  nRR.ptr_o.vio_nt=list(pars=ulist(P0,mRR.ptr_o.vio_nt=RR2,nsc.ptr_o.vio_nt=ns2,nRR.shape='step'),evts='ptr_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.fun)),
  # duration RR
  dRR.dep_x.dep_u=list(pars=ulist(P0,dsc.dep_x.dep_u=ds2,dRR.shape='step'),evts='dep_x',strat='dep_o.dt.c',e.dts=list(dep_o=ds2)),
  dRR.haz_x.haz_u=list(pars=ulist(P0,dsc.haz_x.haz_u=ds2,dRR.shape='step'),evts='haz_x',strat='haz_o.dt.c',e.dts=list(haz_o=ds2)),
  # timestep stuff
  Ri.m.dtz     =list(pars=ulist(P0,dtz=dtz),evts=evts[1:7]),
  tRR.dep_o.dtz=list(pars=ulist(P0,dtz=dtz,iRR.dep_o.vio_zr=  RR1,tsc.dep_o.vio_zr=ts2,tRR.shape='step'),evts='dep_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.dep_x.dtz=list(pars=ulist(P0,dtz=dtz,iRR.dep_x.vio_zr=1/RR1,tsc.dep_x.vio_zr=ts2,tRR.shape='step'),evts='dep_x',strat='vio.dt.c',e.dts=list(vio=ts2)),
  # combo RR
  RR2.haz_o=list(pars=ulist(P0,RR.haz_o.haz_p=  RR1,RR.haz_o.dep_w=RR1,aggr.rate=c('add','mult')),evts='haz_o',strat=c('haz.past','dep.now')),
  RR2.ptr_o=list(pars=ulist(P0,RR.ptr_o.dep_w=1/RR1,RR.ptr_o.haz_w=RR1,aggr.rate=c('add','mult')),evts='ptr_o',strat=c('dep.now','haz.now')))

# =============================================================================
# run & plot

val.run = function(name,pars,evts,strat='.',e.dts=NULL,x.cols=NULL){
  status(2,'val.run: ',name)
  Ps = get.pars.grid(pars)
  pa = attributes(Ps)
  Ms = sim.runs(Ps)
  Y = rate.datas(Ms,p.vars=names(pa$var),e.dts=e.dts,x.cols=x.cols)
  R = rbind.lapply(evts,rate.est,Y=Y,strat=c('seed',strat,names(pa$var)))
  Q = srv.apply(Ms,srvs=c(srv.base,def.args(srv.e.dts,e.dts=e.dts)),p.vars=names(pa$var),x.cols=x.cols)
  for (evt in evts){
    val.plot.rate(name,R,pa,strat,evt)
    val.plot.prev(name,Q,pa,strat,evt)
  }
}

val.plot.rate = function(name,R,pa,strat,evt){
  ref = pa$fix[[str(evt,'.Ri.my')]]/365
  g = plot.rate(R,evt=evt,strat=strat,facet=names(pa$var),ref=ref)
  g = add.info(g,list.str(pa$fix,sig=3,rnd=9))
  plot.save(g,fig.dir,uid,str(c(name,evt,'rate'),collapse='--'))
}

val.plot.prev = function(name,Q,pa,strat,evt){
  g = plot.mean(Q,evt=evt,strat=strat,facet=names(pa$var))
  g = add.info(g,list.str(pa$fix,sig=3,rnd=9))
  plot.save(g,fig.dir,uid,str(c(name,evt,'prev'),collapse='--'))
}

# =============================================================================
# main

for (name in names(vals)){ vals[[name]]$name = name }
for (val in vals[v]){ do.call(val.run,val,quote=TRUE) }
