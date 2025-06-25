source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')

# -----------------------------------------------------------------------------
# config

uid   = '2025-06-25'
seed  = cli.arg('seed',   666)
n.sam = cli.arg('n.sam',10000) # final: 100000
.b    = cli.arg('.b', 1)
.nb   = cli.arg('.nb',1)

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
T = sub.targs(T,sub.targ.age,ags=10)

# -----------------------------------------------------------------------------
# params & priors

P0 = list(
  dtz   = cli.arg('dtz',        45), # final: 7
  n.pop = cli.arg('n.pop',   10000), # final: 10000
  seed  = 1:cli.arg('n.seed',    3), # final: 7
  n.dur = 1,
  null  = 'xRR',
  run = get.run.par(c('dep','haz'),u=FALSE))

F = fpar.set(
  gen.fpar(id='dep_o.Ri.my',   lo=.01,up=.1),
  gen.fpar(id='haz_o.Ri.my',   lo=.01,up=.1),
  gen.fpar(id='dep_x.Ri.my',   lo=.3, up= 3),
  gen.fpar(id='haz_x.Ri.my',   lo=.3, up= 3),
  gen.fpar(id='dep.Ri.het',    lo= 0, up= 5),
  gen.fpar(id='haz.Ri.het',    lo= 0, up= 5),
  gen.fpar(id='dep.cov',       lo=-1, up=+1),
  gen.fpar(id='haz.cov',       lo=-1, up=+1),
  gen.fpar(id='RR.haz_o.dep_w',lo= 1, up=10),
  gen.fpar(id='RR.haz_x.dep_w',lo=.1, up= 1))

data.path = function(.save=FALSE){
  info = list(n.sam=n.sam,seed=seed,P0=P0,F=lapply(F$pars,unclass),T=T)
  hash.path(info,'data','sim','lhs',uid,n.sam,.save=.save)
}

# -----------------------------------------------------------------------------
# main

S = fpar.sam(F,n=n.sam,seed=seed)
Y = fit.run.batch(S,T,P0,.batch=.b,.nbatch=.nb)
save.rda(Y,data.path(.save=TRUE),str('b',.nb),str('Y.',.b))
