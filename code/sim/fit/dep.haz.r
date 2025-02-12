source('sim/meta.r')
source('sim/fit.r')
source('sim/mass.r')

# -----------------------------------------------------------------------------
# config

n.sam  = cli.arg('n.sam', 10000)
n.pop  = cli.arg('n.pop',   400)
n.seed = cli.arg('n.seed',    3)

# base params
P0 = list(
  n.pop = n.pop,
  n.dur = 1,
  null = 'xRR',
  dep_o.Ri.my = .05, haz_o.Ri.my = .05,
  dep_x.Ri.my = .05, haz_x.Ri.my = .05,
  dep.Ri.cv   =   0, haz.Ri.cv   =   0,
  dep.cov     =   0, haz.cov     =   0,
  run = get.run.par(c('dep','haz'),u=FALSE))

# params to fit
F = list(
  'dep_o.Ri.my'    = list(m= -1.5, lo= -3, up=  0, tx=e10, itx=log10),
  'haz_o.Ri.my'    = list(m= -1.5, lo= -3, up=  0, tx=e10, itx=log10),
  'dep_x.Ri.my'    = list(m= -1.5, lo= -3, up=  0, tx=e10, itx=log10),
  'haz_x.Ri.my'    = list(m= -1.5, lo= -3, up=  0, tx=e10, itx=log10),
  'dep.Ri.cv'      = list(m= 0,    lo=  0, up=  9),
  'haz.Ri.cv'      = list(m= 0,    lo=  0, up=  9),
  'dep.cov'        = list(m= 0,    lo= -1, up= +1),
  'haz.cov'        = list(m= 0,    lo= -1, up= +1),
  'RR.haz_o.dep_w' = list(m= 0,    lo= -1, up= +1, tx=e10, itx=log10),
  'RR.haz_x.dep_w' = list(m= 0,    lo= -1, up= +1, tx=e10, itx=log10))

# targets
T = list(
  'dep.now'    = list(type='prop',t.arg=list(p = .15,n=100),w=1,v='dep.now'),
  'dep.past'   = list(type='prop',t.arg=list(p = .45,n=100),w=1,v='dep.past'),
  'haz.now'    = list(type='prop',t.arg=list(p = .25,n=100),w=1,v='haz.now'),
  'haz.past'   = list(type='prop',t.arg=list(p = .60,n=100),w=1,v='haz.past'),
  'dep.haz.or' = list(type='OR',  t.arg=list(OR=3.00,n=100),w=1,ve='dep.now',vo='haz.now',va1='age',ao1=FALSE))

# -----------------------------------------------------------------------------
# main

fid = list.str(list(F=len(F),T=len(T),H=n.sam,n=n.pop,s=n.seed),def='',join='.')
status(1,'fit: dep.haz @ ',fid)
S = lhs.sample(F,n.sam)
L = fit.runs(S,T,P0=P0,seed=1:n.seed,aggr=FALSE)
save.csv(cbind(S,L),'data','sim','fit','dep.haz',uid,fid)
