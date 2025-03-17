source('sim/mass/core.r')

# -----------------------------------------------------------------------------
# config

Ro      = cli.arg('Ro',.03)
Rx      = cli.arg('Rx',.03)
RR.step = cli.arg('RR.step',1)
n.seed  = cli.arg('n.seed',7)
n.pop   = cli.arg('n.pop',400)

P0 = list(
  n.pop = n.pop,
  n.dur = 1+1,
  null = 'xRR',
  dep_o.Ri.my =  Ro, haz_o.Ri.my =  Ro,
  dep_x.Ri.my =  Rx, haz_x.Ri.my =  Rx,
  dep.Ri.cv   =   0, haz.Ri.cv   =   0,
  dep.cov     =   0, haz.cov     =   0,
  run = get.run.par(c('dep','haz'),u=FALSE))

grid = list(
  RR.haz_o.dep_w=signif(2^seq( 0,+3,RR.step),3),
  RR.haz_x.dep_w=signif(2^seq(-3, 0,RR.step),3))
pv = names(grid)

mfuns = list(OR=def.args(glm.elu,family=binomial(link='logit'),ctx=exp))

fid = root.path(create=TRUE,'data','sim','mass','dep.haz',uid,
  sprintf('s%d.g%d.Ro%d.Rx%d',n.seed,prod(lens(grid)),Ro*100,Rx*100))

# -----------------------------------------------------------------------------
# run sim & save outputs

status(3,'running: ',fid)

Ms = sim.runs(get.pars.grid(ulist(P0,grid),seed=1:n.seed))
save(Ms,file=str(fid,'.Ms.rda'))

Qs = par.lapply(-365*0:6,srv.apply,Ms=Ms,p.vars=pv)
A = rbind.lapply(Qs,mass.calc,mfuns=mfuns,
  ve='dep.now',vo='haz.now',va1='age',by=pv,vs=pv,ao1=FALSE)
save(A,file=str(fid,'.A.rda'))

# -----------------------------------------------------------------------------
# plotting

A$RRo = as.factor(A$RR.haz_o.dep_w)
A$RRx = as.factor(A$RR.haz_x.dep_w)
A$I3  = A$lower < 3 & A$upper > 3
b = seq(0,100,10)
g = ggplot(aggregate(I3~RRo+RRx,A,mean),aes(x=RRo,y=RRx,fill=100*I3)) +
  geom_tile(color=NA) + coord_fixed() +
  labs(x    = 'Base Rates of Onset (per year)',
       y    = 'Base Rates of Recovery (per year)',
       fill = 'Proportion\nof 95% CI\ncontaining\nOR = 3 (%)\n') +
  scale_fill_viridis_b(option='inferno',
    limits=range(b),breaks=b,labels=ifelse(b%%20,'',b))
ggsave(str(fid,'.tile.pdf'),w=5,h=4)
