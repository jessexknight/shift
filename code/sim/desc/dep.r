source('sim/meta.r')

# -----------------------------------------------------------------------------
# config

n.pop  = cli.arg('n.pop', 1000)
n.seed = cli.arg('n.seed',  21)

# base params
P0 = list(
  n.pop = n.pop,
  n.dur = 1,
  null = 'aRR',
  het.distr = 'gamma',
  dRR.shape = 'exp',
  dep_o.Ri.my = .02,
  dep_x.Ri.my = .10,
  dep.Ri.het  = 0.0,
  dep.cov     = 0.0,
  RR.dep_o.dep_p = 1,
  dsc.dep_x.dep_u = Inf,
  run = get.run.par('dep',u=FALSE))

# -----------------------------------------------------------------------------
# aggr

srv.aggr = function(Ms,grid,out='dep.now',fun=mean,age.10=FALSE){
  Q = srv.apply(Ms,p.vars=names(grid))
  f = formula(str(out,'~',str(names(grid),collapse='+'),ifelse(age.10,'+age.10','')))
  Qa = aggregate(f,Q,fun)
}

# -----------------------------------------------------------------------------
# plot

# TODO

# -----------------------------------------------------------------------------
# main

grid = list(
  # dep_o.Ri.my = c(.001,.002,.005,.01,.02,.05,.1),
  # dep_x.Ri.my = c(.01 ,.02 ,.05 ,.1 ,.2 ,.5 ,1 ),
  dep_o.Ri.my = c(.001,.003,.01,.03,.1),
  dep_x.Ri.my = c(.01 ,.03 ,.1 ,.3 ,1 ),
  dep.Ri.het  = c( 0,.1,.3,1 ,3 ),
  dep.cov     = c(-.9,  0,+.9))
Ps = get.pars.grid(ulist(P0,grid),seed=1:n.seed)
Ms = sim.runs(Ps)
save(Ms,file='Ms.o5.x5.h5.c3.s21.n1000.rda')
