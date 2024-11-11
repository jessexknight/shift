source('sim/meta.r')

Ps = get.pars.grid(seed=1:7)
Ms = sim.runs(Ps)
Q  = srv.apply(Ms)

print(summary(Q[,c(
  'vio.nt',
  'dep.now','dep.past',
  'haz.now','haz.past',
  'ptr.nw','ptr.nt',
  'age'
)]))
