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
F = name.list(
  gen.fpar('dep_o.Ri.my',   min=-3,max= 0,tx=e10,itx=log10),
  gen.fpar('haz_o.Ri.my',   min=-3,max= 0,tx=e10,itx=log10),
  gen.fpar('dep_x.Ri.my',   min=-3,max= 0,tx=e10,itx=log10),
  gen.fpar('haz_x.Ri.my',   min=-3,max= 0,tx=e10,itx=log10),
  gen.fpar('dep.Ri.cv',     min= 0,max= 9),
  gen.fpar('haz.Ri.cv',     min= 0,max= 9),
  gen.fpar('dep.cov',       min=-1,max=+1),
  gen.fpar('haz.cov',       min=-1,max=+1),
  gen.fpar('RR.haz_o.dep_w',min=-1,max=+1,tx=e10,itx=log10),
  gen.fpar('RR.haz_x.dep_w',min=-1,max=+1,tx=e10,itx=log10))

# targets
T = name.list(
  gen.targ('dep.now',   type='prop',mu= .15,se=.05,w=1,vo='dep.now'),
  gen.targ('dep.past',  type='prop',mu= .45,se=.05,w=1,vo='dep.past'),
  gen.targ('haz.now',   type='prop',mu= .25,se=.05,w=1,vo='haz.now'),
  gen.targ('haz.past',  type='prop',mu= .60,se=.05,w=1,vo='haz.past'),
  gen.targ('dep.haz.or',type='OR',mu=log(3),se=.50,w=1,ve='dep.now',vo='haz.now',va1='age',ao1=FALSE))


# -----------------------------------------------------------------------------
# main

fid = list.str(list(F=len(F),T=len(T),H=n.sam,n=n.pop,s=n.seed),def='',join='.')
status(1,'fit: dep.haz @ ',fid)
S = fpar.lhs(F,n.sam)
Y = fit.runs(S,T,P0=P0,seed=1:n.seed)
save.csv(Y,'data','sim','fit','dep.haz',uid,fid)
