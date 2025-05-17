source('sim/meta.r')
source('sim/fit.r')
uid = '2025-05-17'

# -----------------------------------------------------------------------------
# targets / outcomes

T0 = list(
  list(id='dep.now',  type='prop',mu=NA,se=NA,w=1,vo='dep.now'),
  list(id='dep.past', type='prop',mu=NA,se=NA,w=1,vo='dep.past'),
  list(id='dep.pt.30',type='prop',mu=NA,se=NA,w=1,vo='dep.pt>.30',vsub=TRUE),
  list(id='dep.pt.10',type='prop',mu=NA,se=NA,w=1,vo='dep.pt>.10',vsub=TRUE),
  list(id='dep.pt.03',type='prop',mu=NA,se=NA,w=1,vo='dep.pt>.03',vsub=TRUE),
  list(id='dep.pt.01',type='prop',mu=NA,se=NA,w=1,vo='dep.pt>.01',vsub=TRUE),
  list(id='dep.ne.0', type='prop',mu=NA,se=NA,w=1,vo='dep.ne==0',vsub=TRUE),
  list(id='dep.ne.1', type='prop',mu=NA,se=NA,w=1,vo='dep.ne==1',vsub=TRUE),
  list(id='dep.ne.2+',type='prop',mu=NA,se=NA,w=1,vo='dep.ne>1', vsub=TRUE),
  list(id='dep.um.3m',type='prop',mu=NA,se=NA,w=1,vo='dep.um>t1y/4',vsub=TRUE,sub='dep.past'),
  list(id='dep.um.6m',type='prop',mu=NA,se=NA,w=1,vo='dep.um>t1y/2',vsub=TRUE,sub='dep.past'),
  list(id='dep.um.1y',type='prop',mu=NA,se=NA,w=1,vo='dep.um>t1y',  vsub=TRUE,sub='dep.past'),
  list(id='dep.um.2y',type='prop',mu=NA,se=NA,w=1,vo='dep.um>t1y*2',vsub=TRUE,sub='dep.past'),
  list(id='dep.um.5y',type='prop',mu=NA,se=NA,w=1,vo='dep.um>t1y*5',vsub=TRUE,sub='dep.past'),
  list(id='dep.eRo',  type='pois',mu=NA,se=NA,w=1,vo='dep.past',vt='dep.tto'))
T = list()
ags = 10
alos = seq(amin,amax-ags,ags)
aall = str(amin,'-',amax-1)
for (Ti in T0){ id = Ti$id
  sub = ifelse(is.null(Ti$sub),'',str(Ti$sub,' & '))
  T[[Ti$id]] = do.call(gen.targ,Ti)
  for (a in alos){
    Ti$id  = str(id,':',a)
    Ti$sub = str(sub,'age >= ',a,' & age < ',a+ags)
    T[[Ti$id]] = do.call(gen.targ,Ti)
  }
}

# -----------------------------------------------------------------------------
# default params

P0 = list(
  dtz = cli.arg('dtz',45), # final: 7
  n.pop = cli.arg('n.pop',10000), # final: 10000
  seed = 1:cli.arg('n.seed',7), # final: 100
  n.dur = 1,
  dep_o.Ri.my = .01,
  dep_x.Ri.my = 1,
  het.distr = 'lnorm',
  dep.Ri.het = 0,
  dep.cov = 0,
  dRR.shape = 'pow',
  dsc.dep_x.dep_u = Inf,
  RR.dep_o.dep_p = 1,
  run = get.run.par('dep',u=FALSE))
t1y = add.pars.time(P0,P0$dtz)$t1y

# -----------------------------------------------------------------------------
# param grid & run sims

PG = list(
  dep_o.Ri.my     = seq(.01,.1,.01),
  dep_x.Ri.my     = seq(.5,3,.5),
  dep.Ri.het      = c(0,.1,.2,.3,.5,1,2,3,4,5),
  dep.cov         = c(-.9,-.6,-.3, 0,+.3,+.6,+.9),
  RR.dep_o.dep_p  = 1+c(0,.1,.2,.3,.5,1,2,3,4,5),
  dsc.dep_x.dep_u = round(t1y*c(Inf,1,.8,.6,.4,.2)))
PGk = list(
  hom = PG[c(1,2)],
  het = PG[c(1,2,3,4)],
  int = PG[c(1,2,5,6)])

grid.path = function(k,.save=FALSE){
  hash.path(ulist(P0,PGk[[k]]),'data','sim','depr',uid,.save=.save)
}

run.grid = function(k){
  Y = fit.run.grid(PGk[[k]],T,P0,srvs=srv.extra)
  Y = cbind(Y,col.split(Y$id,':',c('out','age')))
  Y$age = add.na(int.cut(Y$age,alos,up=amax),aall)
  Y[c('targ.mu','targ.se','ll')] = NULL
  save.rda(Y,grid.path(k,.save=TRUE),'Y')
  for (o in unique(Y$out)){
    save.rda(subset(Y,out==o),grid.path(k),str('Y.',o))
  }
}

# -----------------------------------------------------------------------------
# plot setup

ext = '.png'; font = 'Alegreya Sans'
plot.1o = list(w1=2.0,h1=2,wo=2,ho=1) # plot size
ymm = list(dep.now=c(01,12),dep.past=c(05,50)) # gray rects
cmap = lapply(c(het='cividis',Ro='plasma',Rx='viridis'), # colormaps
  function(o){ clr.map.d(option=o) })

load.grid = function(k,out='dep.now',f=NULL,age=FALSE){
  Y = load.rda(grid.path(k),str('Y.',out))
  if (!age){ Y = subset(Y,age==aall) }
  iR = Y$type == 'pois' # rate out rows
  c3 = c('value','lower','upper') # out value cols
  Y[ iR,c3] = Y[ iR,c3]*100*t1y # rates per 100 PY
  Y[!iR,c3] = Y[!iR,c3]*100     # props as %
  Y$Ro  = round(Y$dep_o.Ri.my*100,1) # per 100 PY
  Y$Rx  = round(Y$dep_x.Ri.my*100,1) # per 100 PY
  Y$het = Y$dep.Ri.het        # shorthand
  Y$cor = Y$dep.cov           # shorthand
  Y$RRp = Y$RR.dep_o.dep_p    # shorthand
  Y$RRu = if.null(Y$dsc.dep_x.dep_u,Inf)/t1y
  Y[f] = lapply(Y[f],as.factor) # Y[f] -> factors
  return(Y)
}

# grid subsets for plotting
PGi = list(
  Ro  = seq(1,10,2),
  Rx  = seq(50,300,50),
  het = c(0,.1,.2,.5,1,2,5),
  cor = c(-.6,-.3, 0,+.3,+.6),
  RRp = 1+c(0,.1,.2,.5,1,2,5),
  RRu = c(Inf,1,.8,.6,.4,.2))
PGii = list(Ro=c(2,5,8),Rx=c(100,200,300),cor=c(-.6,0,+.6))
PGiii = list(Ro=c(3,7),Rx=c(100,200))
cor.lab = c('–0.9'=-.9,'–0.6'=-.6,'–0.3'=-.3,'0'=0,
            '+0.3'=+.3,'+0.6'=+.6,'+0.9'=+.9)

# aes labels
l = list(
  Ro  = 'Mean~onset rate~(per 100 PY)',
  Rx  = 'Mean~recov rate~(per 100 PY)',
  het = 'Rate CV',
  cor = 'Rate~correlation',
  RRp = 'RR relapse~vs onset',
  RRu = 'Recovery~waning scale~(years)',
  age = 'Age (years)',
  dep.now  = 'Current~depression~prevalence (%)',
  dep.past = 'Lifetime~depression~prevalence (%)')

# aes label utils
hom = function(s){ gsub('Mean~(.)','\\U\\1',s,perl=TRUE) }
grp = function(s){ gsub('~','\n',s) }
axi = function(s){ gsub('~',' ',s) }
fct = function(s){ ss = strsplit(axi(s),' \\(|\\)')[[1]]; ss[len(ss)+1] = '';
  str.lab(str(' ',ss[1],': '),str(' ',ss[2])) }

# -----------------------------------------------------------------------------
# plot core

plot.prev = function(g,out,ylim,tx){
  xmin = list(identity=-Inf,log10=0)[[tx]]
  ymmo = if.null(ymm[[out]],c(-Inf,+Inf))
  mask = def.args(annotate,'rect',alpha=1/2,fill='#ccc',xmin=xmin,xmax=+Inf)
  g = g + mask(ymin=-Inf,ymax=ymmo[1]) + mask(ymax=+Inf,ymin=ymmo[2])
}

plot.core = function(g,out='dep.now',tx='identity',ribbon=1/5,ylim=c(0,15),ci=.95){
  g = plot.prev(g,out,ylim,tx)
  g = plot.clean(g,font=font,legend.spacing=unit(0,'mm')) +
    stat_summary(geom='ribbon',color=NA,alpha=ribbon,
      fun.min=qfun((1-ci)/2),fun.max=qfun(1-(1-ci)/2)) +
    stat_summary(fun=median,geom='line') +
    scale_x_continuous(trans=tx) +
    coord_cartesian(ylim=ylim) +
    ylab(axi(l[[out]]))
}

# -----------------------------------------------------------------------------
# main

run.grid('hom')
run.grid('het')
run.grid('int')
