source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')

# -----------------------------------------------------------------------------
# config

P0 = list(
  dtz   =   cli.arg('dtz',     60),
  n.pop =   cli.arg('n.pop',  400),
  seed  = 1:cli.arg('n.seed', 100),
  n.dur = 1,
  null  = 'xRR',
  dep_o.Ri.my = .02, haz_o.Ri.my = .02,
  dep_x.Ri.my = .20, haz_x.Ri.my = .20,
  dep.Ri.cv   =   0, haz.Ri.cv   =   0,
  dep.cov     =   0, haz.cov     =   0,
  run = get.run.par(c('dep','haz'),u=FALSE))

T = name.list(key='id',
  gen.targ(id='dep.now',   type='prop',mu=NA,se=NA,w=1,vo='dep.now'),
  gen.targ(id='dep.past',  type='prop',mu=NA,se=NA,w=1,vo='dep.past'),
  gen.targ(id='haz.now',   type='prop',mu=NA,se=NA,w=1,vo='haz.now'),
  gen.targ(id='haz.past',  type='prop',mu=NA,se=NA,w=1,vo='haz.past'),
  gen.targ(id='dep.haz.or',type='OR',  mu=NA,se=NA,w=1,ve='dep.now',vo='haz.now',va1='age',ao1=FALSE),
  gen.targ(id='dep.haz.pr',type='PR',  mu=NA,se=NA,w=1,ve='dep.now',vo='haz.now',va1='age',ao1=FALSE))

PG = list(
  RR.haz_o.dep_w=signif(2^seq( 0,+3,.1),3),
  RR.haz_x.dep_w=signif(2^seq(-3, 0,.1),3))

path = hash.path(ulist(P0,PG),'data','sim','mass','dep.haz',uid)

# -----------------------------------------------------------------------------
# plotting

add.interval = function(Y,v=3,ci=.95){
  Y$int = exp(Y$est.mu+qnorm(.5-ci/2)*Y$est.se) < v &
          exp(Y$est.mu+qnorm(.5+ci/2)*Y$est.se) > v
  return(Y)
}

plot.tile = function(Y,id,...,out='value',aggr=mean){
  Y = Y[Y$id==id,]
  Y$RRo = as.factor(Y$RR.haz_o.dep_w)
  Y$RRx = as.factor(Y$RR.haz_x.dep_w)
  f = formula(aggr.form(out,c('RRo','RRx',...)))
  g = ggplot(aggregate(f,Y,aggr),aes(x=RRo,y=RRx)) +
    geom_tile(aes.string(fill=out),color=NA) + coord_fixed() +
    labs(x='\nRR of Drinking Onset while Depressed',
         y='\nRR of Drinking Recovery while Depressed')
  g = plot.clean(g,
    axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
    legend.key.height=unit(1,'null'))
}

# -----------------------------------------------------------------------------
# main

Y = fit.run.grid(PG,T,P0); save.rda(Y,path,'Y')
Y = add.interval(load.rda(path,'Y'),v=3)

g = plot.tile(Y,'haz.now')              + clr.map.c(option='viridis',limits=c(0,.2)); print(g)
g = plot.tile(Y,'dep.haz.or')           + clr.map.c(option='plasma', limits=c(1,9));  print(g)
g = plot.tile(Y,'dep.haz.or',out='int') + clr.map.c(option='inferno',limits=c(0,1));  print(g)
