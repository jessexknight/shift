source('sim/meta.r')
source('sim/fit.r')
uid = '2025-04-25'

# -----------------------------------------------------------------------------
# targets / outcomes

T0 = list(
  list(id='dep.now',  type='prop',mu=NA,se=NA,w=1,vo='dep.now'),
  list(id='dep.past', type='prop',mu=NA,se=NA,w=1,vo='dep.past'),
  list(id='dep.pt.30',type='prop',mu=NA,se=NA,w=1,vo='dep.pt>.30',vsub=TRUE),
  list(id='dep.pt.10',type='prop',mu=NA,se=NA,w=1,vo='dep.pt>.10',vsub=TRUE),
  list(id='dep.pt.03',type='prop',mu=NA,se=NA,w=1,vo='dep.pt>.03',vsub=TRUE),
  list(id='dep.pt.01',type='prop',mu=NA,se=NA,w=1,vo='dep.pt>.01',vsub=TRUE),
  list(id='dep.ne.0', type='prop',mu=NA,se=NA,w=1,vo='dep.ne==0',vsub=TRUE),
  list(id='dep.ne.1', type='prop',mu=NA,se=NA,w=1,vo='dep.ne==1',vsub=TRUE),
  list(id='dep.ne.2+',type='prop',mu=NA,se=NA,w=1,vo='dep.ne>1', vsub=TRUE),
  list(id='dep.Ro',   type='pois',mu=NA,se=NA,w=1,vo='dep.past',vt='dep.tto'))
T = list()
ags = 10
for (Ti in T0){ id = Ti$id
  T[[Ti$id]] = do.call(gen.targ,Ti)
  for (a in seq(amin,amax-ags,ags)){
    Ti$id  = str(id,':',a)
    Ti$sub = str('age >= ',a,' & age < ',a+ags)
    T[[Ti$id]] = do.call(gen.targ,Ti)
  }
}

# -----------------------------------------------------------------------------
# default params

P0 = list(
  dtz = cli.arg('dtz',45), # final: 7
  n.pop = cli.arg('n.pop',10000), # final: 10000
  seed = 1:cli.arg('n.seed',7), # final: 100
  n.dur = 1,
  het.distr = 'lnorm',
  dRR.shape = 'exp',
  dep_o.Ri.my = .01,
  dep_x.Ri.my = .50,
  dep.Ri.het  = 0,
  dep.cov     = 0,
  RR.dep_o.dep_p = 1,
  dsc.dep_x.dep_u = Inf,
  run = get.run.par('dep',u=FALSE))
t1y = add.pars.time(P0,P0$dtz)$t1y

# -----------------------------------------------------------------------------
# param grid & run sims

PG = list(
  dep_o.Ri.my     = c(.001,.002,.003,.005,.01,.02,.03,.05,.10),
  dep_x.Ri.my     = c(.100,.200,.300,.500, 1 , 2 , 3 , 5 ,10 ),
  dep.Ri.het      = c(0,.1,.2,.3,.5,1,2,3,5),
  dep.cov         = c(-.9,-.6,-.3, 0,+.3,+.6,+.9),
  RR.dep_o.dep_p  = 1 + c(0,.1,.2,.3,.5,1,2,3,5),
  dsc.dep_x.dep_u = c(Inf,50,30,20,10,5,3,2,1) * t1y / log(2))

grid.path = function(p){
  hash.path(ulist(P0,PG[p]),'data','sim','depr',uid)
}

run.grid = function(p=NULL){
  Y.age = fit.run.grid(PG[p],T,P0,srvs=srv.extra)
  Y.all = subset(Y.age,!grepl(':',id))
  save.rda(Y.age,grid.path(p),'Y.age')
  save.rda(Y.all,grid.path(p),'Y.all')
}

load.grid = function(p=NULL,f=NULL){
  Y = load.rda(grid.path(p),'Y')
  Y$value = pmax(Y$value * 100, .1)
  Y$Ro  = Y$dep_o.Ri.my * 100
  Y$Rx  = Y$dep_x.Ri.my * 100
  Y$het = Y$dep.Ri.het
  Y$cor = Y$dep.cov
  Y$RRp = Y$RR.dep_o.dep_p
  Y$RRu = if.null(Y$dsc.dep_x.dep_u,Inf) * log(2) / 360
  Y[f] = lapply(Y[f],as.factor)
  return(Y)
}

# -----------------------------------------------------------------------------
# plot

# TODO

# -----------------------------------------------------------------------------
# main

run.grid(c(1,2))
run.grid(c(1,2,3,4))
run.grid(c(1,2,5,6))
