source('meta.r')

Ps = lapply(1:7,get.pars)
Ss = sim.runs(Ps)
Qs = srv.map(Ss)

print(summary(Qs[Qs$age<amax,c(
  'vio.n',
  'dep.now','dep.past',
  'haz.now','haz.past',
  'ptr.n','ptr.tot',
  'age'
)]))
