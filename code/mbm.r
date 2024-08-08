source('meta.r')
options(width=256)

N = list(rep=3,seed=7)
f = file.path('.tmp','mbm',sprintf('mbm_s%dr%d_%s.out',
  N$seed,N$rep,Sys.Date()))

mbm = microbenchmark::microbenchmark(
  {sim.runs(lapply(1:N$seed,get.pars,n=30))},
  {sim.runs(lapply(1:N$seed,get.pars,n=100))},
  {sim.runs(lapply(1:N$seed,get.pars,n=300))},
  {sim.runs(lapply(1:N$seed,get.pars,n=1000))},
times=N$rep)

print(mbm)
sink(f)
print(mbm)
sink()
