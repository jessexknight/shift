source('sim.r')
source('plot.r')

Ps = lapply(1:7,get.pars)
Is = sim.runs(Ps,.par=TRUE)
