source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')
uid = '2025-07-30'
.k  = cli.arg('.k','o6')
.b  = cli.arg('.b', 1)
.nb = cli.arg('.nb',100)

# -----------------------------------------------------------------------------
# targets / outcomes

T = name.list(key='id',
  gen.targ(id='dep.now',    type='prop',vo='dep.now'),
  gen.targ(id='dep.past',   type='prop',vo='dep.past'),
  gen.targ(id='haz.now',    type='prop',vo='haz.now'),
  gen.targ(id='haz.past',   type='prop',vo='haz.past'),
  gen.targ(id='dep.haz.or', type='OR',  ve='dep.now',vo='haz.now',ao1=FALSE),
  gen.targ(id='dep.haz.pr', type='PR',  ve='dep.now',vo='haz.now',ao1=FALSE),
  gen.targ(id='dep.haz.aor',type='OR',  ve='dep.now',vo='haz.now',ao1=FALSE,va1='age'),
  gen.targ(id='dep.haz.apr',type='PR',  ve='dep.now',vo='haz.now',ao1=FALSE,va1='age'))

# -----------------------------------------------------------------------------
# params & grid

P0 = list(
  dtz   =   cli.arg('dtz',     45), # final: 7
  n.pop =   cli.arg('n.pop',10000), # final: 10000
  seed  = 1:cli.arg('n.seed',  11), # final: 41
  n.dur = 1,
  null  = 'xRR',
  het.distr = 'lnorm',
  dep_o.Ri.my = 0, haz_o.Ri.my = 0,
  dep_x.Ri.my = 0, haz_x.Ri.my = 0,
  dep.Ri.het  = 0, haz.Ri.het  = 0,
  dep.cov     = 0, haz.cov     = 0,
  RR.haz_o.dep_w = 1,
  RR.haz_x.dep_w = 1,
  run = get.run.par(c('dep','haz'),u=FALSE))

PG0 = list(
  RR.haz_o.dep_w = 1*c(1,2,4,8),
  RR.haz_x.dep_w = 1/c(1,2,4,8),
  dep_o.Ri.my = c(.01,.03,.1,.3), haz_o.Ri.my = c(.01,.03,.1,.3),
  dep_x.Ri.my = c(0,.1,.3,1,3),   haz_x.Ri.my = c(0,.1,.3,1,3),
  dep.Ri.het  = c(0,.3,1,3),      haz.Ri.het  = c(0,.3,1,3),
  dep.cov     = c(-.5,0,+.5),     haz.cov     = c(-.5,0,+.5))

PGk = list(
  o6 = PG0[1:6],
  eo6 = PG0[c(1:4,7:8)],
  eu6 = ulist(PG0[1:6],dep.Ri.het=1,haz.Ri.het=1),
  en6 = ulist(PG0[1:6],dep.Ri.het=1,haz.Ri.het=1,dep.cov=-.5,haz.cov=-.5),
  ep6 = ulist(PG0[1:6],dep.Ri.het=1,haz.Ri.het=1,dep.cov=+.5,haz.cov=+.5))

grid.path = function(k,.save=FALSE){
  hash.path(ulist(P0,PGk[[k]],set=k),'data','sim','aor',uid,k,.save=.save)
}

# -----------------------------------------------------------------------------
# run & save/load

run.grid = function(k){
  Y = fit.run.grid(PGk[[k]],T,P0,.batch=.b,.nbatch=.nb)
  save.rda(Y,grid.path(k,.save=TRUE),str('b',.nb),str('Y.',.b))
}

recut.rda = function(k){
  # join batches + split by target id for speed
  Y = rbind.lapply(1:.nb,function(b){ load.rda(grid.path(k),str('b',.nb),str('Y.',b)) })
  for (i in unique(Y$id)){ save.rda(subset(Y,id==i),grid.path(k),str('Y.',i)) }
}

# -----------------------------------------------------------------------------
# main

# run.grid(.k)
# recut.rda(.k)
