source('sim/meta.r')
source('sim/fit.r')
uid = '2025-05-15'

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
  list(id='dep.Ro',   type='pois',mu=NA,se=NA,w=1,vo='dep_o.Ri',vt='.',sub='!dep.now'),
  list(id='dep.Ro.n', type='pois',mu=NA,se=NA,w=1,vo='dep_o.Ri',vt='.',sub='!dep.now & !dep.past'),
  list(id='dep.Ro.p', type='pois',mu=NA,se=NA,w=1,vo='dep_o.Ri',vt='.',sub='!dep.now &  dep.past'),
  list(id='dep.Rx',   type='pois',mu=NA,se=NA,w=1,vo='dep_x.Ri',vt='.',sub=' dep.now'),
  list(id='dep.Rx.n', type='pois',mu=NA,se=NA,w=1,vo='dep_x.Ri',vt='.',sub=' dep.now & !dep.past'),
  list(id='dep.Rx.p', type='pois',mu=NA,se=NA,w=1,vo='dep_x.Ri',vt='.',sub=' dep.now &  dep.past'),
  list(id='dep.eRo.n',type='pois',mu=NA,se=NA,w=1,vo='dep.past',vt='dep.tto'))
T = list()
ags = 10
for (Ti in T0){ id = Ti$id
  sub = ifelse(is.null(Ti$sub),'',str(Ti$sub,' & '))
  T[[Ti$id]] = do.call(gen.targ,Ti)
  for (a in seq(amin,amax-ags,ags)){
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
  dep_o.Ri.my     = .005 + seq(0,.1,.01),
  dep_x.Ri.my     = seq(.5,3,.5),
  dep.Ri.het      = c(0,.1,.2,.3,.5,1,2,3,5),
  dep.cov         = c(-.9,-.6,-.3, 0,+.3,+.6,+.9),
  RR.dep_o.dep_p  = 1 + c(0,.1,.2,.3,.5,1,2,3,5),
  dsc.dep_x.dep_u = c(Inf,1,.8,.6,.4,.2))

grid.path = function(p){
  hash.path(ulist(P0,PG[p]),'data','sim','depr',uid)
}

run.grid = function(p=NULL){
  Y.age = fit.run.grid(PG[p],T,P0,srvs=srv.extra,i.vars=c('dep_o.Ri','dep_x.Ri'))
  Y.all = subset(Y.age,!grepl(':',id))
  save.rda(Y.age,grid.path(p),'Y.age')
  save.rda(Y.all,grid.path(p),'Y.all')
}

# -----------------------------------------------------------------------------
# plot setup

ext = '.png'; font = 'Alegreya Sans'
plot.1o = list(w1=2.0,h1=2,wo=2,ho=1.2) # plot size
ymm = list(dep.now=c(01,20),dep.past=c(03,60)) # gray rects
cmap = lapply(c(het='cividis',Ro='plasma',Rx='viridis'), # colormaps
  function(o){ clr.map.d(option=o) })

load.grid = function(p=NULL,f=NULL,age=FALSE){
  Y = load.rda(grid.path(p),ifelse(age,'Y.age','Y.all'))
  Y[c('targ.mu','targ.se','ll','t')] = NULL # rm cols
  Y = cbind(Y,col.split(Y$id,':',c('out','age.10'))) # id -> out,age.10
  Y$age.10 = add.na(int.cut(as.numeric(Y$age.10),seq(10,50,10),up=60),'10-59')
  iR = which(Y$type == 'pois') # rate out rows
  c3 = c('value','lower','upper') # out value cols
  Y[ iR,c3] = Y[ iR,c3] * 100 * t1y # rates per 100 PY
  Y[-iR,c3] = Y[-iR,c3] * 100       # props as %
  Y$Ro  = Y$dep_o.Ri.my * 100 # per 100 PY
  Y$Rx  = Y$dep_x.Ri.my * 100 # per 100 PY
  Y$het = Y$dep.Ri.het        # shorthand
  Y$cor = Y$dep.cov           # shorthand
  Y$RRp = Y$RR.dep_o.dep_p    # shorthand
  Y$RRu = if.null(Y$dsc.dep_x.dep_u,Inf)
  Y[f] = lapply(Y[f],as.factor) # Y[f] -> factors
  return(Y)
}

# grid subsets for plotting
PGi = list(
  Ro  = .5 + seq(0,10,2),
  Rx  = seq(50,300,50),
  het = c(0,.1,.2,.5,1,2,5),
  cor = c(-.6,-.3, 0,+.3,+.6),
  RRp = 1 + c(0,.1,.3,1,3),
  RRu = c(Inf,1,.8,.6,.4,.2))
# PGii = TODO
# PGiii = TODO
cor.lab = c('–0.9'=-.9,'–0.6'=-.6,'–0.3'=-.3,'0'=0,
            '+0.3'=+.3,'+0.6'=+.6,'+0.9'=+.9)

# aes labels
l = list(
  Ro  = 'Mean~onset rate~(per 100 PY)',
  Rx  = 'Mean~recov rate~(per PY)',
  het = 'Rate CV',
  cor = 'Rate~correlation',
  RRp = 'HR–1 relapse vs onset',
  RRu = 'Recovery rate half-life (years)',
  age = 'Age (years)',
  dep.now  = 'Current~depression~prevalence (%)',
  dep.past = 'Lifetime~depression~prevalence (%)',
  dep.Ro   = 'Observed~onset rate~(per 100 PY)')

# aes label utils
hom = function(s){ gsub('Mean~(.)','\\U\\1',s,perl=TRUE) }
grp = function(s){ gsub('~','\n',s) }
axi = function(s){ gsub('~',' ',s) }
fct = function(s){ ss = strsplit(axi(s),' \\(|\\)')[[1]]; ss[len(ss)+1] = '';
  str.lab(str(' ',ss[1],': '),str(' ',ss[2])) }

# -----------------------------------------------------------------------------
# plot core

plot.prev = function(g,out,ylim,by=NULL,ty='log10'){
  by = if.null(by,seq(0,15,1))
  xmm = layer_scales(g)$x$range$range
  mask = def.args(annotate,'rect',alpha=1/2,fill='gray',xmin=xmm[1],xmax=xmm[2])
  g = g + scale_y_continuous(trans=ty,breaks=by,labels=by,minor=NULL) +
    mask(ymin=ylim[1],ymax=if.null(ymm[[out]][1],NA)) +
    mask(ymax=ylim[2],ymin=if.null(ymm[[out]][2],NA))
}

plot.core = function(g,out='dep.now',tx=NULL,ribbon=1/5,ylim=c(0,15),by=NULL,ci=.95){
  g = plot.prev(g,out,ylim,by)
  g = plot.clean(g,font=font) +
    stat_summary(geom='ribbon',color=NA,alpha=ribbon,
      fun.min=qfun((1-ci)/2),fun.max=qfun(1-(1-ci)/2)) +
    stat_summary(fun=median,geom='line') +
    scale_x_continuous(trans=tx) +
    coord_cartesian(ylim=ylim) +
    ylab(axi(l[[out]]))
}

# -----------------------------------------------------------------------------
# main

run.grid(c(1,2))
run.grid(c(1,2,3,4))
run.grid(c(1,2,5,6))
