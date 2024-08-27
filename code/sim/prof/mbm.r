source('sim/meta.r')
options(width=256)

.verb = 0
n = list(rep=3,seed=7)
f = file.path('.log','prof',sprintf('mbm_s%dr%d_%s.out',
  n$seed,n$rep,uid))

mbm = microbenchmark::microbenchmark(
  {sim.runs(lapply(1:n$seed,get.pars,n.pop=30))},
  {sim.runs(lapply(1:n$seed,get.pars,n.pop=100))},
  {sim.runs(lapply(1:n$seed,get.pars,n.pop=300))},
  {sim.runs(lapply(1:n$seed,get.pars,n.pop=1000))},
times=n$rep)

print(mbm)
sink(f)
print(mbm)
sink()
