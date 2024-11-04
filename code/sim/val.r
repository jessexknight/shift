source('sim/meta.r')

# =============================================================================
# config

v      = cli.arg('v',     0)
n.pop  = cli.arg('n.pop', 333)
n.seed = cli.arg('n.seed',7)
fig.dir = 'val'

# -----------------------------------------------------------------------------

P0 = list(
  n.pop = n.pop,
  vio.Ri.my   = .50,
  dep_o.Ri.my = .10, dep_x.Ri.my = 1.0,
  haz_o.Ri.my = .10, haz_x.Ri.my = 1.0,
  ptr_o.Ri.my = 3.0, ptr_x.Ri.my = 3.0, ptr.max.m = 3,
  all.Ri.shape = 1e6,
  ptr.Ri.shape = 1e6,
  null = 'xRR')

# =============================================================================
# scenarios

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
  # base rates [1]
  Ri.m=list(pars=ulist(P0),evts=evts[1:7]),
  # RR age [2:5]
  aRR.vio=list(pars=ulist(P0,aRR.vio.RRs  =aRR5,aRR.vio.ages  =age5,aRR.shape='const'),evts='vio',  strat='age.10',x.cols=list(age.10=age.10.fun)),
  aRR.dep=list(pars=ulist(P0,aRR.dep_o.RRs=aRR5,aRR.dep_o.ages=age5,aRR.shape='const'),evts='dep_o',strat='age.10',x.cols=list(age.10=age.10.fun)),
  aRR.haz=list(pars=ulist(P0,aRR.haz_o.RRs=aRR5,aRR.haz_o.ages=age5,aRR.shape='const'),evts='haz_o',strat='age.10',x.cols=list(age.10=age.10.fun)),
  aRR.ptr=list(pars=ulist(P0,aRR.ptr_o.RRs=aRR5,aRR.ptr_o.ages=age5,aRR.shape='const'),evts='ptr_o',strat='age.10',x.cols=list(age.10=age.10.fun)),
  # basic RR [6:13]
  RR.dep_o.dep_p=list(pars=ulist(P0,RR.dep_o.dep_p=  RR3),evts='dep_o',strat='dep.past'),
  RR.haz_o.haz_p=list(pars=ulist(P0,RR.haz_o.haz_p=  RR3),evts='haz_o',strat='haz.past'),
  RR.haz_o.dep_w=list(pars=ulist(P0,RR.haz_o.dep_w=  RR3),evts='haz_o',strat='dep.now'),
  RR.haz_x.dep_w=list(pars=ulist(P0,RR.haz_x.dep_w=1/RR3),evts='haz_x',strat='dep.now'),
  RR.ptr_o.dep_w=list(pars=ulist(P0,RR.ptr_o.dep_w=1/RR3),evts='ptr_o',strat='dep.now'),
  RR.ptr_o.haz_w=list(pars=ulist(P0,RR.ptr_o.haz_w=  RR3),evts='ptr_o',strat='haz.now'),
  RR.ptr_x.dep_w=list(pars=ulist(P0,RR.ptr_x.dep_w=  RR3),evts='ptr_x',strat='dep.now'),
  RR.ptr_x.haz_w=list(pars=ulist(P0,RR.ptr_x.haz_w=  RR3),evts='ptr_x',strat='haz.now'),
  # transient RR [14:18]
  tRR.dep_o.vio_zr=list(pars=ulist(P0,iRR.dep_o.vio_zr=  RR2,tsc.dep_o.vio_zr=ts2,tRR.shape='step'),evts='dep_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.dep_x.vio_zr=list(pars=ulist(P0,iRR.dep_x.vio_zr=1/RR2,tsc.dep_x.vio_zr=ts2,tRR.shape='step'),evts='dep_x',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.haz_o.vio_zr=list(pars=ulist(P0,iRR.haz_o.vio_zr=  RR2,tsc.haz_o.vio_zr=ts2,tRR.shape='step'),evts='haz_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.haz_x.vio_zr=list(pars=ulist(P0,iRR.haz_x.vio_zr=1/RR2,tsc.haz_x.vio_zr=ts2,tRR.shape='step'),evts='haz_x',strat='vio.dt.c',e.dts=list(vio=ts2)),
  tRR.ptr_o.vio_zr=list(pars=ulist(P0,iRR.ptr_o.vio_zr=  RR2,tsc.ptr_o.vio_zr=ts2,tRR.shape='step'),evts='ptr_o',strat='vio.dt.c',e.dts=list(vio=ts2)),
  # cumulative RR [19:21]
  nRR.dep_o.vio_nt=list(pars=ulist(P0,mRR.dep_o.vio_nt=RR2,nsc.dep_o.vio_nt=ns2,nRR.shape='step'),evts='dep_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.fun)),
  nRR.haz_o.vio_nt=list(pars=ulist(P0,mRR.haz_o.vio_nt=RR2,nsc.haz_o.vio_nt=ns2,nRR.shape='step'),evts='haz_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.fun)),
  nRR.ptr_o.vio_nt=list(pars=ulist(P0,mRR.ptr_o.vio_nt=RR2,nsc.ptr_o.vio_nt=ns2,nRR.shape='step'),evts='ptr_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.fun)),
  # duration RR [22:23]
  dRR.dep_x.dep_u=list(pars=ulist(P0,dsc.dep_x.dep_u=ds2,dRR.shape='step'),evts='dep_x',strat='dep_o.dt.c',e.dts=list(dep_o=ds2)),
  dRR.haz_x.haz_u=list(pars=ulist(P0,dsc.haz_x.haz_u=ds2,dRR.shape='step'),evts='haz_x',strat='haz_o.dt.c',e.dts=list(haz_o=ds2)),
  # partnership dynamics [24:26]
  n.pop  =list(pars=ulist(P0,n.pop=c(100,300,900)),evts=c('ptr_o','ptr_x')),
  ptr.max=list(pars=ulist(P0,ptr.max.m=c(1,3,30)), evts=c('ptr_o','ptr_x')),
  dtz.ptr=list(pars=ulist(P0,dtz=1:7),             evts=c('ptr_o','ptr_x'))
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
    val.plot.rate(name,R,pars,strat,evt)
    val.plot.prev(name,Q,pars,strat,evt)
  }
}

val.par.split = function(par){
  # split up variable vs fixed par, & count total var combo
  b = lens(par) > 1 & !grepl(names(null.sets$aRR),names(par))
  p = list(var=par[b],fix=par[!b],n.var=prod(lens(par[b])))
}

val.plot.rate = function(name,R,pars,strat,evt){
  ref = pars$fix[[str(evt,'.Ri.my')]]/365
  g = plot.rate(R,evt=evt,strat=strat,facet=names(pars$var),ref=ref)
  g = add.info(g,pars$fix)
  plot.save(g,fig.dir,uid,str(c(name,evt,'rate'),collapse='--'))
}

val.plot.prev = function(name,Q,pars,strat,evt){
  g = plot.prev(Q,evt=evt,strat=strat,facet=names(pars$var))
  g = add.info(g,pars$fix)
  plot.save(g,fig.dir,uid,str(c(name,evt,'prev'),collapse='--'))
}

# =============================================================================
# main

for (name in names(vals)){ vals[[name]]$name = name }
for (val in vals[v]){ do.call(val.run,val,quote=TRUE) }
