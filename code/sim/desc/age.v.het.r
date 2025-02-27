source('sim/meta.r')
source('sim/fit.r')

# -----------------------------------------------------------------------------
# params & config

var    = cli.arg('var',  'dep')
n.pop  = cli.arg('n.pop', 1000)
n.seed = cli.arg('n.seed',  21)

def.aRR = function(P){
  age = seq(1,amax)-amin
  P$aRR.dep_o = exp(-if.null(P$aRR.dep_o.decay,0)*age)
  P$aRR.haz_o = exp(-if.null(P$aRR.haz_o.decay,0)*age)
  P$dep_x.Ri.my = 0
  P$haz_x.Ri.my = 0
  return(P) }

# base params
P0 = list(
  n.pop = n.pop,
  n.dur = 1,
  null = 'aRR',
  het.distr = 'gamma',
  dRR.shape = 'exp',
  dep_o.Ri.my = .02, haz_o.Ri.my = .02,
  dep_x.Ri.my = .10, haz_x.Ri.my = .10,
  dep.Ri.het  = 0.0, haz.Ri.het  = 0.0,
  dep.cov     = 0.0, haz.cov     = 0.0,
  fun = def.aRR,
  run = get.run.par(var,u=FALSE))

# -----------------------------------------------------------------------------
# targets & fit

def.targs.dep = function(){
  T = name.list( # McGrath2023
    # gen.targ('dep.past.10',type='prop',mu=.003,se=.00100,vo='dep.past',among='age<=10'),
    gen.targ('dep.past.15',type='prop',mu=.018,se=.00135,vo='dep.past',among='age<=15'),
    gen.targ('dep.past.20',type='prop',mu=.052,se=.00241,vo='dep.past',among='age<=20'),
    gen.targ('dep.past.25',type='prop',mu=.076,se=.00276,vo='dep.past',among='age<=25'),
    gen.targ('dep.past.30',type='prop',mu=.095,se=.00300,vo='dep.past',among='age<=30'),
    gen.targ('dep.past.35',type='prop',mu=.112,se=.00329,vo='dep.past',among='age<=35'),
    gen.targ('dep.past.40',type='prop',mu=.130,se=.00351,vo='dep.past',among='age<=40'),
    gen.targ('dep.past.45',type='prop',mu=.146,se=.00367,vo='dep.past',among='age<=45'),
    gen.targ('dep.past.50',type='prop',mu=.161,se=.00383,vo='dep.past',among='age<=50'),
    gen.targ('dep.past.55',type='prop',mu=.175,se=.00418,vo='dep.past',among='age<=55'),
    gen.targ('dep.past.60',type='prop',mu=.183,se=.00428,vo='dep.past',among='age<=60'))
}

def.targs.haz = function(){
  T = name.list( # McGrath2023
    # gen.targ('haz.past.15',type='prop',mu=.001,se=.00100,vo='haz.past',among='age<=10'),
    gen.targ('haz.past.15',type='prop',mu=.017,se=.00153,vo='haz.past',among='age<=15'),
    gen.targ('haz.past.20',type='prop',mu=.108,se=.00396,vo='haz.past',among='age<=20'),
    gen.targ('haz.past.25',type='prop',mu=.146,se=.00455,vo='haz.past',among='age<=25'),
    gen.targ('haz.past.30',type='prop',mu=.172,se=.00491,vo='haz.past',among='age<=30'),
    gen.targ('haz.past.35',type='prop',mu=.185,se=.00493,vo='haz.past',among='age<=35'),
    gen.targ('haz.past.40',type='prop',mu=.195,se=.00502,vo='haz.past',among='age<=40'),
    gen.targ('haz.past.45',type='prop',mu=.204,se=.00517,vo='haz.past',among='age<=45'),
    gen.targ('haz.past.50',type='prop',mu=.209,se=.00522,vo='haz.past',among='age<=50'),
    gen.targ('haz.past.55',type='prop',mu=.212,se=.00527,vo='haz.past',among='age<=55'),
    gen.targ('haz.past.60',type='prop',mu=.214,se=.00525,vo='haz.past',among='age<=60'))
}

err.fun = function(x){
  Y = verb.wrap(fit.run(Si=x,T=T,P0=P0,seed=1:n.seed,.par=TRUE),2)
  # print(round(aggregate(targ.mu~name,Y,mean)$targ.mu,3)) # DEBUG
  # print(round(aggregate( est.mu~name,Y,mean)$ est.mu,3)) # DEBUG
  ll = -mean(aggregate(ll~seed,Y,sum)$ll)
}

# -----------------------------------------------------------------------------
# main

cases = list( # MAN fits
  dep = list(
    ref  = c(dep_o.Ri.my=.01,aRR.dep_o.decay=.000,dep.Ri.het=0.0),
    age  = c(dep_o.Ri.my=.01,aRR.dep_o.decay=.012,dep.Ri.het=0.0),
    het  = c(dep_o.Ri.my=.01,aRR.dep_o.decay=.000,dep.Ri.het=1.2)),
  haz = list(
    ref  = c(haz_o.Ri.my=.013,aRR.haz_o.decay=.000,haz.Ri.het=0.0),
    age  = c(haz_o.Ri.my=.035,aRR.haz_o.decay=.130,haz.Ri.het=0.0),
    het  = c(haz_o.Ri.my=.060,aRR.haz_o.decay=.000,haz.Ri.het=3.0))
)[[var]]

T = list(dep=def.targs.dep(),haz=def.targs.haz())[[var]]

Y = rbind.lapply(names(cases),function(i){
  Yi = fit.run(cases[[i]],T=T,P0=P0,seed=1:n.seed,.par=TRUE)
  Yi = cbind(Yi,case=i,age=substr(Yi$name,10,11))
},.par=FALSE)
g = plot.targs(Y,x=age,color=case,geom=geom_boxplot) +
  labs(y=str(var,'_past'),x='Age') + clr.map.d()
ggsave(plot=g,file=str(var,'.pdf'),w=8,h=5)
