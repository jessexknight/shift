source('meta.r')

Ps = lapply(1:7,get.pars)
Is = sim.runs(Ps)

print(summary(Is[Is$age<amax,c(
  'vio.n','vio.tot',
  'dep.now','dep.past',
  'haz.now','haz.past',
  'ptr.n','ptr.tot',
  'age'
)]))
