source('sim/meta.r')
source('sim/fit.r')
source('sim/mass.r')

# -----------------------------------------------------------------------------
# config

h.init = cli.arg('h.init',   10)
n.iter = cli.arg('n.iter',  100)
n.pop  = cli.arg('n.pop',  1000)
n.seed = cli.arg('n.seed',    7)

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
F = fpar.set(
  gen.fpar(id='dep_o.Ri.my',   lo=-3,up= 0,tr=e10),
  gen.fpar(id='haz_o.Ri.my',   lo=-3,up= 0,tr=e10),
  gen.fpar(id='dep_x.Ri.my',   lo=-3,up= 0,tr=e10),
  gen.fpar(id='haz_x.Ri.my',   lo=-3,up= 0,tr=e10),
  gen.fpar(id='dep.Ri.cv',     lo= 0,up= 9),
  gen.fpar(id='haz.Ri.cv',     lo= 0,up= 9),
  gen.fpar(id='dep.cov',       lo=-1,up=+1),
  gen.fpar(id='haz.cov',       lo=-1,up=+1),
  gen.fpar(id='RR.haz_o.dep_w',lo=-1,up=+1,tr=e10),
  gen.fpar(id='RR.haz_x.dep_w',lo=-1,up=+1,tr=e10))

# targets
T = name.list(key='id',
  gen.targ(id='dep.now',   type='prop',mu= .15,se=.05,w=1,vo='dep.now'),
  gen.targ(id='dep.past',  type='prop',mu= .45,se=.05,w=1,vo='dep.past'),
  gen.targ(id='haz.now',   type='prop',mu= .25,se=.05,w=1,vo='haz.now'),
  gen.targ(id='haz.past',  type='prop',mu= .60,se=.05,w=1,vo='haz.past'),
  gen.targ(id='dep.haz.or',type='OR',mu=log(3),se=.50,w=1,ve='dep.now',vo='haz.now',va1='age',ao1=FALSE))

# -----------------------------------------------------------------------------
# main

# Y = fit.run(fpar.sam(n=1,F,tr=TRUE),T,P0,seed=1:n.seed); print(Y) # DEBUG
fid = list.str(list(F=len(F),T=len(T),n=n.pop,s=n.seed,h=h.init,i=n.iter),def='',join='.')
status(1,'fit: dep.haz @ ',fid)
O = opt.run(F,T=T,P0=P0,n.seed=n.seed,h.init=h.init,n.iter=n.iter)
save.rda(O,'data','sim','fit','dep.haz',uid,str('opt.',fid))
print(O)
