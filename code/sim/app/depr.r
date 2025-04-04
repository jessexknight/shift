source('sim/meta.r')
source('sim/fit.r')
uid = '2025-04-04'

# -----------------------------------------------------------------------------
# targets / outcomes

T = name.list(key='id',
  gen.targ(id='dep.now', type='prop',mu=NA,se=NA,w=1,vo='dep.now'),
  gen.targ(id='dep.past',type='prop',mu=NA,se=NA,w=1,vo='dep.past'))
ags = 5
for (v in names(T)){
  for (a in seq(amin,amax-ags,ags)){
    id = str(v,':',a)
    T[[id]] = gen.targ(id=id,type='prop',mu=NA,se=NA,w=1,vo=v,
      among=str('age >= ',a,' & age < ',a+ags))
}}

# -----------------------------------------------------------------------------
# default params

P0 = list(
  dtz    = cli.arg('dtz',     45),
  n.pop  = cli.arg('n.pop', 1000),
  n.seed = cli.arg('n.seed', 100),
  n.dur  = 1,
  het.distr = 'lnorm',
  dRR.shape = 'exp',
  dep_o.Ri.my = .01,
  dep_x.Ri.my = .20,
  dep.Ri.het  = 0,
  dep.cov     = 0,
  RR.dep_o.dep_p = 1,
  dsc.dep_x.dep_u = Inf,
  run = get.run.par('dep',u=FALSE))

# -----------------------------------------------------------------------------
# param grid & run sims

PG = list(
  dep_o.Ri.my     = c(.001,.002,.005,.01,.02,.05,.10),
  dep_x.Ri.my     = c(.100,.200,.500, 1 , 2 , 5 ,10 ),
  dep.Ri.het      = c(0,.1,.3,1,3),
  dep.cov         = c(-.5, 0,+.5),
  RR.dep_o.dep_p  = 1 + c(0,.1,.3,1,3),
  dsc.dep_x.dep_u = c(Inf,30,10,3,1)*360/log(2)
)

grid.path = function(p){
  hash.path(ulist(P0,PG[p]),'data','sim','depr',uid)
}

run.grid = function(p=NULL){
  Y = fit.run.grid(PG[p],T,P0)
  save.rda(Y,grid.path(p),'Y')
}

load.grid = function(p=NULL,f=NULL){
  Y = load.rda(grid.path(p),'Y')
  Y$value = Y$value * 100
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
