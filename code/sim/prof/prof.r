source('sim/meta.r')

.verb = 0
n  = list(seed=3)
f  = file.path('.log','prof',sprintf('prof_s%d_%s.html',
  n$seed,uid))

out = profvis::profvis({
  Ps = lapply(1:n$seed,get.pars)
  Ms = sim.runs(Ps,.par=FALSE)
})
htmlwidgets::saveWidget(out,f)
