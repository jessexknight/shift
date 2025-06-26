source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')
uid = '2025-05-30'
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
  seed  = 1:cli.arg('n.seed',  10), # final: 10
  n.dur = 1,
  null  = 'xRR',
  het.distr = 'lnorm',
  dep_o.Ri.my = .04, haz_o.Ri.my = .04,
  dep_x.Ri.my = 1,   haz_x.Ri.my = 1,
  dep.Ri.het  = 0,   haz.Ri.het  = 0,
  dep.cov     = 0,   haz.cov     = 0,
  RR.haz_o.dep_w = 1,
  RR.haz_x.dep_w = 1,
  run = get.run.par(c('dep','haz'),u=FALSE))

PG = list(
  dep_o.Ri.my = c(.02,.04,.06), haz_o.Ri.my = c(.02,.04,.06),
  dep_x.Ri.my = c(1,2,3),       haz_x.Ri.my = c(1,2,3),
  dep.Ri.het  = c(0,1,3),       haz.Ri.het  = c(0,1,3),
  dep.cov     = c(-.5,0,+.5),   haz.cov     = c(-.5,0,+.5),
  RR.haz_o.dep_w = signif(2^seq( 0,+3,.5),3),
  RR.haz_x.dep_w = signif(2^seq(-3, 0,.5),3))

PGk = list(
  ref = P0[names(PG[9:10])], # P0
  fix = PG[9:10],        # Ri* ~ fixed (IRR only)
  hom = PG[c(1:4,9:10)], # Ri* ~ homog (no heter + cor)
  het = PG[c(1:6,9:10)], # Ri* ~ heter + uncor
  cor = PG)              # Ri* ~ heter + cor

grid.path = function(k,.save=FALSE){
  hash.path(ulist(P0,PGk[[k]],k=k),'data','sim','cseb',uid,.save=.save)
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
  y = function(k){ if.null(Y[[k]],P0[[k]]) }
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
PG3 = list(RRo=c(1,2,  8),RRx=1/c(1,2,  8))
PG2 = list(RRo=c(  2,  8),RRx=1/c(  2,  8))
PG1 = list(RRo=c(    4  ),RRx=1/c(    4  ))
PG0 = list(dRo=4,dRx=100,hRo=4,hRx=100,dhet=0,hhet=0,dcor=0,hcor=0,RRo=1,RRx=1)
PGi = list(PG1,PG2,PG3,PG4)

# -----------------------------------------------------------------------------
# plot config

ext = '.png'; font = 'Alegreya Sans'
plot.1o = list(w1=2,h1=2,wo=2,ho=1) # plot size
cmap = lapply(c(RRo='cividis',RRx='viridis',Ri='plasma',het='inferno'),
  function(o){ clr.map.d(option=o) })

l = list( # aes labels
  dRo='Mean MD~onset rate~(per 100 PY)', hRo='Mean HD~onset rate~(per 100 PY)',
  dRx='Mean MD~recov rate~(per 100 PY)', hRx='Mean HD~recov rate~(per 100 PY)',
  dhet='MD rate~CV',                     hhet='HD rate~CV',
  dcor='MD rate~correlation',            hcor='HD rate~correlation',
  dep.now='Current MD~prevalence (%)',   haz.now='Current HD~prevalence (%)',
  dep.past='Lifetime MD~prevalence (%)', haz.past='Lifetime HD~prevalence (%)',
  dep.haz.or='OR', dep.haz.aor='OR',
  dep.haz.pr='PR', dep.haz.apr='PR',
  RRo='IRR of HD~onset while~MD',
  RRx='IRR of HD~recov while~MD')
# aes label tweaks
lg = function(v){ gsub('~','\n',l[[v]]) } # group: ~ -> wrap
la = function(v){ gsub('~',' ',l[[v]]) }  # axis: ~ -> space
lf = function(v,j=' '){ # facet: ~ -> space; insert value before (units)
  ss = strsplit(la(v),'\\(|\\)')[[1]]; ss[len(ss)+1] = '';
  str.lab(str('  ',ss[1],':',j),str('  ',ss[2])) }

# -----------------------------------------------------------------------------
# plot core

plot.line = function(g,id='dep.haz.aor',ty='log2',tx='identity',
  ribbon=.1,ci=.95,yi=0,dy=NA,ymm=c(NA,NA)){
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

plot.basic.line = function(x='o'){
  RRa = list(o='RRo',x='RRx')[[x]] # axis
  RRg = list(o='RRx',x='RRo')[[x]] # group (color)
  sub = function(Y){ df.sub(Y,str(RRg,'%in% PG4$',RRg)) }
  Y = sub(load.grid('fix',f=RRg)) # default output: dep.haz.aor
  g = ggplot(df.sub(Y,str(RRg,'==1')),aes.string(y='value',x=RRa))
  g = plot.line(g,tx='log2',dy=1-2*(x=='x'),ymm=c(.75,8)) + labs(x=la(RRa))
  plot.save(g,'cseb',uid,str('basic.aor.',x),ext=ext,size=c(3,3)) # RRa only (RRg==1)
  g = ggplot(Y,aes.string(y='value',x=RRa,color=RRg,fill=RRg))
  g = plot.line(g,tx='log2') + cmap[[RRg]] + labs(x=la(RRa),color=lg(RRg),fill=lg(RRg))
  plot.save(g,'cseb',uid,str('basic.aor.',x,2),ext=ext) # RRa (axis) & RRg (color)
  for (id in c('dep.now','haz.now')){ # alternate outputs
    Y = sub(load.grid('fix',id=id,f=RRg))
    g = ggplot(Y,aes.string(y='value',x=RRa,color=RRg,fill=RRg))
    g = plot.line(g,id=id,ty='identity',ymm=c(0,7)) + cmap[[RRg]] + labs(x=la(RRa),color=lg(RRg),fill=lg(RRg))
    plot.save(g,'cseb',uid,str('basic.',id,'.',x),ext=ext) # RRa (axis) & RRg (color)
  }
}

plot.basic.tile = function(){
  Y = load.grid('fix',f=c('RRo','RRx'))
  g = ggplot(Y,aes(z=value,x=RRo,y=RRx))
  g = plot.tile(g) + clr.map.c(option='inferno',trans='log2',limits=c(1,8)) + labs(x=la('RRo'),y=la('RRx'))
  plot.save(g,'cseb',uid,'basic.aor.tile',ext=ext)
}

plot.hom.f = function(f=1){
  # plot aor vs *Rx (axis) & *Ro (color), with f ~ {1,2,3,4} facets for RRx & RRo
  Y = subset(load.grid('hom'), RRo %in% PGi[[f]]$RRo & RRx %in% PGi[[f]]$RRx)
  g = ggplot(subset(Y,hRo==PG0$hRo & hRx==PG0$hRx),aes(y=value,x=dRx,color=factor(dRo),fill=factor(dRo)))
  g = plot.line(g,ymm=c(.5,16)) + cmap$Ri + labs(x=la('dRx'),color=lg('dRo'),fill=lg('dRo'))
  if (f > 1){ g = g + facet_grid('RRo~RRx',labeller=labeller(.rows=lf('RRo'),.cols=lf('RRx'))) }
  plot.save(g,'cseb',uid,str('hom.d.f',f),ext=ext)
  g = ggplot(subset(Y,dRo==PG0$dRo & dRx==PG0$dRx),aes(y=value,x=hRx,color=factor(hRo),fill=factor(hRo)))
  g = plot.line(g,ymm=c(.5,16)) + cmap$Ri + labs(x=la('hRx'),color=lg('hRo'),fill=lg('hRo'))
  if (f > 1){ g = g + facet_grid('RRo~RRx',labeller=labeller(.rows=lf('RRo'),.cols=lf('RRx'))) }
  plot.save(g,'cseb',uid,str('hom.h.f',f),ext=ext)
}

plot.het.1 = function(){
  # plot aor vs het, with het in dep only / haz only / both equally (fix base rates)
  Y = subset(load.grid('het'), RRo %in% PG1$RRo & RRx %in% PG1$RRx &
             dRo==PG0$dRo & dRx==PG0$dRx & hRo==PG0$hRo & hRx==PG0$hRx)
  Y$het = pmax(Y$dhet,Y$hhet)
  Y = rbind(
    cbind(subset(Y,hhet==0),   hg=factor('MD rates')),
    cbind(subset(Y,dhet==0),   hg=factor('HD rates')),
    cbind(subset(Y,dhet==hhet),hg=factor('both')))
  g = ggplot(Y,aes(y=value,x=het,lty=hg))
  g = plot.line(g,ymm=c(1,16)) +
      scale_linetype_manual(values=c('11','61','4111')) +
      labs(x='Base rate heterogeneity (CV)',lty='Heterog. in:')
  plot.save(g,'cseb',uid,'het.1',ext=ext)
}

plot.het.2 = function(){
  # plot aor vs het for all combos of dep/haz het (fix base rates)
  Y = subset(load.grid('het'), RRo %in% PG2$RRo & RRx %in% PG2$RRx &
    dRo==PG0$dRo & dRx==PG0$dRx & hRo==PG0$hRo & hRx==PG0$hRx)
  g = ggplot(Y,aes(y=value,x=dhet,color=factor(hhet),fill=factor(hhet)))
  g = plot.line(g,ymm=c(.5,32)) + cmap$het + labs(x=la('dhet'),color=lg('hhet'),fill=lg('hhet')) +
      facet_grid('RRo~RRx',labeller=labeller(.rows=lf('RRo'),.cols=lf('RRx')))
  plot.save(g,'cseb',uid,'het.d.f2',ext=ext)
  g = ggplot(Y,aes(y=value,x=hhet,color=factor(dhet),fill=factor(dhet)))
  g = plot.line(g,ymm=c(.5,32)) + cmap$het + labs(x=la('hhet'),color=lg('dhet'),fill=lg('dhet')) +
      facet_grid('RRo~RRx',labeller=labeller(.rows=lf('RRo'),.cols=lf('RRx')))
  plot.save(g,'cseb',uid,'het.h.f2',ext=ext)
}

plot.het.tile = function(){
  # plot (mean) aor vs *RR (axes) as tile with *het as facets
  Y = subset(load.grid('het',f=c('RRo','RRx')),
    dRx==PG0$dRx & hRx==PG0$hRx & dRo==PG0$dRo & hRo==PG0$hRo)
  Y = aggregate(value~RRo+RRx+dhet+hhet,Y,mean)
  Y$value = Y$value / Y$value[Y$dhet==0 & Y$hhet==0]
  g = ggplot(Y,aes(z=value,x=RRo,y=RRx))
  g = plot.tile(g) + clr.map.v(trans='log10',limits=c(1/4,4/1)) +
    facet_grid('dhet~hhet',labeller=labeller(.cols=lf('hhet'),.rows=lf('dhet'))) +
    labs(x=la('RRo'),y=la('RRx'),fill='OR (CV>0)\nvs\nOR (CV=0)')
  plot.save(g,'cseb',uid,'het.tile',ext=ext,size=c(5,4))
}

# -----------------------------------------------------------------------------
# main

# run.grid('ref')
# run.grid('fix')
# run.grid('hom')
# run.grid('het')
# run.grid('cor')
# recut.rda(TODO)

plot.basic.line('o')
plot.basic.line('x')
plot.basic.tile()
plot.hom.f(1)
plot.hom.f(2)
plot.het.1()
plot.het.2()
plot.het.tile()
