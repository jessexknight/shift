source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')
uid = '2025-07-25'
.k  = cli.arg('.k','o6')
.b  = cli.arg('.b', 1)
.nb = cli.arg('.nb',1)

# -----------------------------------------------------------------------------
# targets / outcomes

T = name.list(key='id',
  gen.targ(id='dep.now',    type='prop',vo='dep.now'),
  gen.targ(id='dep.past',   type='prop',vo='dep.past'),
  gen.targ(id='haz.now',    type='prop',vo='haz.now'),
  gen.targ(id='haz.past',   type='prop',vo='haz.past'),
  gen.targ(id='dep.haz.or', type='OR',  ve='dep.now',vo='haz.now',ao1=FALSE),
  gen.targ(id='dep.haz.pr', type='PR',  ve='dep.now',vo='haz.now',ao1=FALSE),
  gen.targ(id='dep.haz.aor',type='OR',  ve='dep.now',vo='haz.now',ao1=FALSE,va1='age'),
  gen.targ(id='dep.haz.apr',type='PR',  ve='dep.now',vo='haz.now',ao1=FALSE,va1='age'))

# -----------------------------------------------------------------------------
# params & grid

P0 = list(
  dtz   =   cli.arg('dtz',     45), # final: 7
  n.pop =   cli.arg('n.pop',10000), # final: 10000
  seed  = 1:cli.arg('n.seed',  41), # final: 41
  n.dur = 1,
  null  = 'xRR',
  het.distr = 'lnorm',
  dep_o.Ri.my = .04, haz_o.Ri.my = .02,
  dep_x.Ri.my = 1,   haz_x.Ri.my = .333,
  dep.Ri.het  = 0,   haz.Ri.het  = 0,
  dep.cov     = 0,   haz.cov     = 0,
  RR.haz_o.dep_w = 1,
  RR.haz_x.dep_w = 1,
  run = get.run.par(c('dep','haz'),u=FALSE))

PG = list(
  RR.haz_o.dep_w = signif(2^seq( 0,+3,.5),3),
  RR.haz_x.dep_w = signif(2^seq(-3, 0,.5),3),
  dep_o.Ri.my = c(.02,.04,.06), haz_o.Ri.my = c(.01,.02,.03),
  dep_x.Ri.my = c(.5,1,1.5),    haz_x.Ri.my = c(.167,.333,.500),
  dep.Ri.het  = c(0,1,3),       haz.Ri.het  = c(0,1,3),
  dep.cov     = c(-.5,0,+.5),   haz.cov     = c(-.5,0,+.5),
  dep.Ri.het  = 1,              haz.Ri.het  = 1) # HACK

PGk = list(
  o2  = PG[c(1:2)],       # fix all (hom)
  e2  = PG[c(1:2,11:12)], # fix all (het)
  o6  = PG[c(1:6)],       # fix hom
  e6  = PG[c(1:6,11:12)], # fix het, no cor
  e8  = PG[c(1:8)],       # vary het, no cor
  e10 = PG[c(1:10)])      # vary het, vary cor

grid.path = function(k,.save=FALSE){
  hash.path(ulist(P0,PGk[[k]],set=k),'data','sim','cseb',uid,.save=.save)
}

# -----------------------------------------------------------------------------
# run & save/load

run.grid = function(k){
  Y = fit.run.grid(PGk[[k]],T,P0,.batch=.b,.nbatch=.nb)
  save.rda(Y,grid.path(k,.save=TRUE),str('b',.nb),str('Y.',.b))
}

recut.rda = function(k){
  # join batches + split by target id for speed
  Y = rbind.lapply(1:.nb,function(b){ load.rda(grid.path(k),str('b',.nb),str('Y.',b)) })
  for (i in unique(Y$id)){ save.rda(subset(Y,id==i),grid.path(k),str('Y.',i)) }
}

load.grid = function(k,id='dep.haz.aor',f=NULL){
  Y = load.rda(grid.path(k),str('Y.',id))
  c3 = c('value','lower','upper') # out value cols
  iP = Y$type == 'prop' # prop out rows
  Y[iP,c3] = Y[iP,c3]*100 # props as %
  iZ = Y$lower==0 & Y$upper==Inf # degenerate glm
  Y[iZ,c3] = NA # TODO: better fix?
  y = function(v){ if.null(Y[[v]],P0[[v]]) }
  Y$dRo  = round(y('dep_o.Ri.my')*100,1); Y$hRo  = round(y('haz_o.Ri.my')*100,1) # per 100 PY
  Y$dRx  = round(y('dep_x.Ri.my')*100,1); Y$hRx  = round(y('haz_x.Ri.my')*100,1) # per 100 PY
  Y$dhet = y('dep.Ri.het');               Y$hhet = y('haz.Ri.het') # shorthand
  Y$dcor = y('dep.cov');                  Y$hcor = y('haz.cov')    # shorthand
  Y$RRo  = y('RR.haz_o.dep_w'); # shorthand
  Y$RRx  = y('RR.haz_x.dep_w'); # shorthand
  Y[f] = lapply(Y[f],as.factor) # Y[f] -> factors
  return(Y)
}

# PG for subsets after scaling
PG4 = list(RRo=c(1,2,4,8),RRx=1/c(1,2,4,8))
PG0 = list(dRo=4,dRx=100,hRo=2,hRx=33.3,dhet=0,hhet=0,dcor=0,hcor=0,RRo=1,RRx=1)

# -----------------------------------------------------------------------------
# plot config

ext = '.png'; font = 'Alegreya Sans'
plot.1o = list(w1=2,h1=2,wo=2,ho=1) # plot size
clrs = list( # HACK
  RRo=c('#600','#903','#c06','#f6c'), # red
  RRx=c('#063','#096','#0c9','#3fc'), # green
  dRo=c('#06c','#39f','#6cf'),hRo=c('#639','#96c','#c9f'), # blue, orange
  dRx=c('#930','#c60','#f90'),hRx=c('#960','#c90','#fc0')) # purple, yellow
cmap = lapply(clrs,function(v){ clr.map.m(values=v) })

l = list( # aes labels
  dRo='Mean~depression~onset rate~(per 100 PY)', hRo='Mean~drinking~onset rate~(per 100 PY)',
  dRx='Mean~depression~recov rate~(per 100 PY)', hRx='Mean~drinking~recov rate~(per 100 PY)',
  dhet='depression~rates~heterogeneity~(CV)',    hhet='drinking~rates~heterogeneity~(CV)',
  bhet='depression~& drinking~rates~heterogeneity~(CV)',
  dcor='depression~rates~correlation',           hcor='drinking~rates~correlation',
  dep.now='Current~depression~prevalence (%)',   haz.now='Current~drinking~prevalence (%)',
  dep.past='Lifetime~depression~prevalence (%)', haz.past='Lifetime~drinking~prevalence (%)',
  dep.haz.or='OR', dep.haz.aor='OR',
  dep.haz.pr='PR', dep.haz.apr='PR',
  RRo='IRR of~drinking~onset while~depressed',
  RRx='IRR of~drinking~recov while~depressed',
  type='Measure',adj='Age adjust')
# aes label tweaks
lg = function(v){ gsub('~','\n',l[[v]]) } # group: ~ -> wrap
la = function(v){ gsub('~',' ',l[[v]]) }  # axis: ~ -> space
lf = function(v,j=' '){ # facet: ~ -> space; insert value before (units)
  ss = strsplit(la(v),'\\(|\\)')[[1]]; ss[len(ss)+1] = '';
  str.lab(str('  ',ss[1],':',j),str('  ',ss[2])) }

# -----------------------------------------------------------------------------
# plot core

plot.line = function(g,id='dep.haz.aor',ty='log2',tx='log2',
  ribbon=.1,ci=.95,yi=0,dy=1,ymm=c(.75,8)){
  g = plot.clean(g,font=font) +
    geom_abline(intercept=yi,slope=dy,color='gray',lty='22') +
    stat_summary(geom='ribbon',color=NA,alpha=ribbon,
      fun.min=qfun((1-ci)/2),fun.max=qfun(1-(1-ci)/2)) +
    stat_summary(geom='line',fun=mean) +
    scale_x_continuous(trans=tx) +
    scale_y_continuous(trans=ty) +
    coord_cartesian(ylim=ymm) +
    ylab(la(id))
}

plot.tile = function(g,id='dep.haz.aor',ymm=c(NA,NA)){
  g = plot.clean(g,font=font) +
    stat_summary_2d(fun=mean) +
    labs(fill=lg(id)) +
    coord_fixed()
}

# -----------------------------------------------------------------------------
# plots

plot.base.meas = function(k){
  Y = subset(rbind.lapply(names(T)[5:8],load.grid,k=k),RRx==1)
  Y$adj = factor(1+grepl('haz.a',Y$id),1:2,c('No','Yes'))
  g = ggplot(Y,aes(x=RRo,y=value,color=type,fill=type,lty=adj,alpha=adj))
  g = plot.line(g) + clr.map.m(values=c(clrs$RRo[3],'#666')) +
    labs(y='Measure',x=la('RRo'),color=lg('type'),fill=lg('type'),alpha=lg('adj'),lty=lg('adj')) +
    scale_linetype_manual(values=c('solid','44')) + scale_alpha_manual(values=c(1/3,3/3))
  plot.save(g,'cseb',uid,k,'base.meas',ext=ext)
}

plot.base.line = function(k,x='o'){
  RRa = list(o='RRo',x='RRx')[[x]] # axis
  RRg = list(o='RRx',x='RRo')[[x]] # group (color)
  sub = function(Y){ df.sub(Y,str(RRg,'%in% PG4$',RRg)) }
  Y = sub(load.grid(k,f=RRg)) # default output: dep.haz.aor
  g = ggplot(df.sub(Y,str(RRg,'==1')),aes.string(y='value',x=RRa,color=RRg,fill=RRg))
  g = plot.line(g,dy=1-2*(x=='x')) + labs(x=la(RRa)) + clr.map.m(values=clrs[[RRa]][3],guide='none')
  g = add.info(g,info=str(rep(' ',17),collapse='')) # HACK
  plot.save(g,'cseb',uid,k,str('base.aor.',x),ext=ext) # RRa only (RRg==1)
  g = ggplot(Y,aes.string(y='value',x=RRa,color=RRg,fill=RRg))
  g = plot.line(g,dy=1-2*(x=='x')) + cmap[[RRg]] + labs(x=la(RRa),color=lg(RRg),fill=lg(RRg))
  plot.save(g,'cseb',uid,k,str('base.aor.',x,2),ext=ext) # RRa (axis) & RRg (color)
  for (id in c('dep.now','haz.now')){ # alternate outputs
    Y = sub(load.grid(k,id=id,f=RRg))
    g = ggplot(Y,aes.string(y='value',x=RRa,color=RRg,fill=RRg))
    g = plot.line(g,id=id,dy=NA,ymm=c(0,8),tx='identity',ty='identity') + cmap[[RRg]] + labs(x=la(RRa),color=lg(RRg),fill=lg(RRg))
    plot.save(g,'cseb',uid,k,str('base.',id,'.',x),ext=ext) # RRa (axis) & RRg (color)
  }
}

plot.base.tile = function(k){
  Y = load.grid(k,f=c('RRo','RRx'))
  g = ggplot(Y,aes(z=value,x=RRo,y=RRx))
  g = plot.tile(g) + clr.map.c(option='inferno',trans='log2',limits=c(1,8)) + labs(x=la('RRo'),y=la('RRx'))
  plot.save(g,'cseb',uid,k,'base.aor.tile',ext=ext)
}

plot.vs.rates = function(k){
  Y = load.grid(k)
  Rs = c('dRo','dRx','hRo','hRx')
  for (Ri in Rs){
    Y$f = factor(Y[[Ri]])
    sub = list.str(PG0[c('RRx',Rs[Rs!=Ri])],def='==',join=' & ')
    g = ggplot(df.sub(Y,sub),aes(y=value,x=RRo,color=f,fill=f))
    g = plot.line(g) + cmap[[Ri]] + labs(x=la('RRo'),color=lg(Ri),fill=lg(Ri))
    plot.save(g,'cseb',uid,k,str('aor.',Ri),ext=ext)
  }
}

# -----------------------------------------------------------------------------
# main

# run.grid(.k)
# recut.rda(.k)

# k2='o2'; k6='o6' # hom
# k2='e2'; k6='e6' # het
# plot.base.meas(k=k2)
# plot.base.line(k=k2,x='o')
# plot.base.line(k=k2,x='x')
# plot.base.tile(k=k2)
# plot.vs.rates (k=k6)
