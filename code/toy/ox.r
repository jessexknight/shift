# simple 2-state model with open pop to examine influence of:
# - R.o.m: mean rate of onset
# - R.x.m: mean rate of recovery
# - R.cv:  CoV of individual rates (both o & x)
# - R.cor: correlation among individual rates

source('utils.r')
source('sim/plot.r')

amin = 10
amax = 60
adur = amax-amin

def.pars = function(seed=NULL,...){
  P = list(seed=seed)
  P$dt = 7
  P$n.dur = 1
  P$n.pop = 1e3
  P$R.o.m = .1
  P$R.x.m = .1
  P$R.cv  = .001
  P$R.cor = +.5
  P = ulist(P,...)
}

init.inds = function(P){
  n = P$n.pop * (1+P$n.dur)
  age = runif(n,min=amin-P$n.dur*adur,max=amax)
  Ri.ox = copula(n,
    covs = P$R.cor,
    qfuns = list(Ri.o=qgamma,Ri.x=qgamma),
    Ri.o = list(shape=1/P$R.cv^2,scale=P$R.o.m/365*P$R.cv^2),
    Ri.x = list(shape=1/P$R.cv^2,scale=P$R.x.m/365*P$R.cv^2) )
  I = data.frame(
    i = seq(n),
    age = age,
    state = FALSE,
    Ri.ox)
}

run.sim = function(P){
  set.seed(P$seed)
  I = init.inds(P)
  ts = seq(0,P$n.dur*adur*365,P$dt)
  dtp <<- function(R){ p = 1-exp(-R*P$dt) }
  for (t in ts){
    # age & select active
    I$age = I$age + P$dt/365
    i = which(I$age > amin & I$age <= amax)
    # state change
    s = I$state[i]
    u = runif(s)
    I$state[i][!s & u < dtp(I$Ri.o[i])] = TRUE
    I$state[i][ s & u < dtp(I$Ri.x[i])] = FALSE
  }
  return(I)
}

run.sims = function(n.seed=7,...){
  print(list.str(list(...),join='   '))
  I = rbind.lapply(1:n.seed,function(seed){
    P = def.pars(seed=seed,...)
    I = run.sim(P)
    I = subset(I,age > amin & age <= amax)
    I$age.10 = int.cut(I$age,seq(10,50,10))
    I = cbind(I,seed=seed,...)
  })
}

plot.mean = function(I,y='state',...){
  f = str(y,' ~ ',str(unname(c(1,'seed',...)),collapse=' + '))
  I = aggregate(formula(f),I,mean)
  g = ggplot(I,aes.string(y=y,...)) +
    geom_boxplot() +
    scale_color_viridis_d() +
    ylim(c(0,NA))
  g = plot.clean(g)
}

# MWE
# g = plot.mean(run.sims(),color='age.10')

# run grid
X = list(
  R.o.m = c(.01,.03,.1,.3,1),
  R.x.m = c(.01,.03,.1,.3,1),
  R.cv  = c(.1,.3,1,3,5),
  R.cor = c(-.9,-.3,0,+.3,+.9))
I = do.call(rbind,grid.apply(X,run.sims))
I[names(X)] = lapply(I[names(X)],as.factor)
# plot grid
I = aggregate(state ~ R.o.m + R.x.m + R.cv + R.cor, I, mean)
g = ggplot(I,aes(y=state,color=R.cor,x=R.cv,
    group=interaction(R.o.m,R.x.m,R.cor))) +
  facet_grid('R.o.m ~ R.x.m',scales='free',labeller=label_both) +
  scale_color_viridis_d() +
  geom_line()
g = plot.clean(g)


