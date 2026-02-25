source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')
uid = '2026-02-26'
.k     = cli.arg('.k','RRo.rev.base')
.b     = cli.arg('.b', 1)
.nb    = cli.arg('.nb',1)
.debug = cli.arg('.debug',0)

# e = exposure "dep" (abuse)
# o = outcome  "haz" (depression)

# -----------------------------------------------------------------------------
# params & grid

P0 = list(
  pop.type = 'open',
  n.pop = xdf(1000,10000),
  n.dur = 1, dtz = 7,
  het.distr = 'gamma',
  run = get.run.par(c('dep','haz'),u=0))

G = name.list(key='i',
  list(i='RRo', id='RR.haz_o.dep_w',v=1,   vg=  c(1,2,3,4,8)),
  list(i='RRx', id='RR.haz_x.dep_w',v=1,   vg=1/c(1,2,3,4,8)),
  list(i='eRo', id='dep_o.Ri.my',   v=.03, vg=c(.003,.01,.03,.1)),
  list(i='eRx', id='dep_x.Ri.my',   v=.3,  vg=c(.03,.1,.3,1)),
  list(i='eHo', id='dep_o.Ri.het',  v=0,   vg=c(0,.3,1,3)),
  list(i='eHx', id='dep_x.Ri.het',  v=0,   vg=c(0,.3,1,3)),
  list(i='ecv', id='dep.cov',       v=0,   vg=c(-.5,0,+.5)),
  list(i='ep',  id='dep.prev',      v=.1,  vg=c(.03,.10,.30)),
  list(i='oRo', id='haz_o.Ri.my',   v=.03, vg=c(.003,.01,.03,.1)),
  list(i='oRx', id='haz_x.Ri.my',   v=.3,  vg=c(.03,.1,.3,1)),
  list(i='oHo', id='haz_o.Ri.het',  v=0,   vg=c(0,.3,1,3)),
  list(i='oHx', id='haz_x.Ri.het',  v=0,   vg=c(0,.3,1,3)),
  list(i='ocv', id='haz.cov',       v=0,   vg=c(-.5,0,+.5)),
  list(i='seed',id='seed',          v=NA,  vg=xdf(1:7,1:21)),
  list(i='ek',  id='e.case', v='rev',vg=NA),
  list(i='ok',  id='o.case', v='rev',vg=NA))

ids = lapply(G,`[[`,'id')
G0 = lapply(G,`[[`,'v')
Gi = function(i,...){ ulist(G0,lapply(G[i],`[[`,'vg'),...) }
PG = function(Gk,...){ ulist(P0,set.names(Gk,ids[names(Gk)]),...) }

Gk = list()
Gk$RRo.rev.base = Gi(ek='rev',c('seed','RRo'))
Gk$RRx.rev.base = Gi(ek='rev',c('seed','RRx'))
Gk$RR2.rev.base = Gi(ek='rev',c('seed','RRo','RRx'))
Gk$RRo.rev.eRo  = Gi(ek='rev',c('seed','RRo','eRo','eHo'))
Gk$RRo.rev.eRx  = Gi(ek='rev',c('seed','RRo','eRx','eHx'))
Gk$RRo.rev.eR2  = Gi(ek='rev',c('seed','RRo','eRo','eRx','ecv'),eHo=1,eHx=1)
Gk$RRo.rev.oRo  = Gi(ek='rev',c('seed','RRo','oRo','oHo'))
Gk$RRo.rev.oRx  = Gi(ek='rev',c('seed','RRo','oRx','oHx'))
Gk$RRo.rev.oR2  = Gi(ek='rev',c('seed','RRo','oRo','oRx','ocv'),oHo=1,oHx=1)
Gk$RRo.rev.2Rx  = Gi(ek='rev',c('seed','RRo','eRx','eHx','oRx','oHx'))
Gk$RRo.fix.base = Gi(ek='fix',c('seed','RRo'))
Gk$RRo.fix.eRo  = Gi(ek='fix',c('seed','RRo','ep'))
Gk$RRo.irr.base = Gi(ek='irr',c('seed','RRo'))
Gk$RRo.irr.eRo  = Gi(ek='irr',c('seed','RRo','eRo','eHo'))
Gk$RRo.irr.oRo  = Gi(ek='irr',c('seed','RRo','oRo','oHo'))
Gk$RRo.irr.2Ro  = Gi(ek='irr',c('seed','RRo','eRo','eHo','oRo','oHo'))
Gk$RRo.rev.lhs  = Gi(ek='rev',seed=c(1,1e9),lhs=xdf(1e1,1e4),
  c('RRo','eRo','eRx','eHo','eHx','ecv','oRo','oRx','oHo','oHx','ocv'))
# for (k in names(Gk)){ status(3,k,': ',prod(lens(Gk[[k]]))) } # for hpc gen

apply.case = function(P,eps=1e-12){
  if (P$e.case=='fix'){
    P$init.inds = function(I,P){
      I$dep_o.Ri = ifelse(runif(P$n.tot) < P$dep.prev,Inf,0)
      return(I) }}
  if (P$e.case!='rev'){ P$dep_x.Ri.m = eps }
  if (P$o.case!='rev'){ P$haz_x.Ri.m = eps }
  return(P)
}

get.lhs = function(Gi,seed=666){
  set.seed(seed)
  is = lens(Gi) > 1
  Gi[is] = as.data.frame(qunif(
    lhs::randomLHS(Gi$lhs,sum(is)),
    rep(unlist(lapply(Gi[is],min)),each=Gi$lhs),
    rep(unlist(lapply(Gi[is],max)),each=Gi$lhs) ))
  return(Gi)
}

# -----------------------------------------------------------------------------
# targets / outcomes

T = name.list(key='id',
  gen.targ(id='e.w',   type='prop',vo='dep.now' ),
  gen.targ(id='e.p',   type='prop',vo='dep.past'),
  gen.targ(id='o.w',   type='prop',vo='haz.now' ),
  gen.targ(id='o.p',   type='prop',vo='haz.past'),
  gen.targ(id='or.ww', type='OR',  ve='dep.now', vo='haz.now' ),
  gen.targ(id='or.wp', type='OR',  ve='dep.now', vo='haz.past'),
  gen.targ(id='or.pw', type='OR',  ve='dep.past',vo='haz.now' ),
  gen.targ(id='or.pp', type='OR',  ve='dep.past',vo='haz.past'),
  gen.targ(id='pr.ww', type='PR',  ve='dep.now', vo='haz.now' ),
  gen.targ(id='pr.pp', type='PR',  ve='dep.now', vo='haz.past'),
  gen.targ(id='pr.pw', type='PR',  ve='dep.past',vo='haz.now' ),
  gen.targ(id='pr.wp', type='PR',  ve='dep.past',vo='haz.past'))

# -----------------------------------------------------------------------------
# run sims & save/load

grid.path = function(k,.save=FALSE){
  path = hash.path(PG(Gk[[k]]),'data','sim','mass',uid,k,.save=.save)
}

run.one = function(...,.par=0){
  P1 = PG(list(...),fun=apply.case)
  Ps = get.pars.grid(P1,.par=.par)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  Q  = srv.apply(Ms,.par=.par)
  Y  = srv.targs(Q,T)
  Y[c('seed','targ.mu','targ.se','ll')] = NULL
  row.names(Y) = NULL
  return(Y)
}

run.grid = function(k){
  lhs = len(Gk[[k]]$lhs)
  Gi = { if (lhs) get.lhs(Gk[[k]]) else Gk[[k]] }
  Y = grid.apply(Gi,run.one,.grid=!lhs,.batch=.b,.nbatch=.nb,
    .rbind=1,.cbind=1,.log=3)
  save.rds(Y,grid.path(k,.save=TRUE),str('b',.nb),str('Y.',.b))
}

merge.batch = function(k){
  Y = rbind.lapply(1:.nb,function(b){
    Yb = load.rds(grid.path(k),str('b',.nb),str('Y.',b)) })
  save.rds(Y,grid.path(k),'Y')
}

# -----------------------------------------------------------------------------
# plots

# TODO

# -----------------------------------------------------------------------------
# main

# run.grid(.k)
# merge.batch(.k)

