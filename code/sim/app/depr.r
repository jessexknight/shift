source('sim/meta.r')
source('sim/fit.r')
uid = '2025-05-17'

# -----------------------------------------------------------------------------
# targets / outcomes

T = name.list(key='id',
  gen.targ(id='dep.now',  type='prop',vo='dep.now'),
  gen.targ(id='dep.past', type='prop',vo='dep.past'),
  gen.targ(id='dep.pt.30',type='prop',vo='dep.pt>.30',vsub=TRUE),
  gen.targ(id='dep.pt.10',type='prop',vo='dep.pt>.10',vsub=TRUE),
  gen.targ(id='dep.pt.03',type='prop',vo='dep.pt>.03',vsub=TRUE),
  gen.targ(id='dep.pt.01',type='prop',vo='dep.pt>.01',vsub=TRUE),
  gen.targ(id='dep.ne.0', type='prop',vo='dep.ne==0',vsub=TRUE),
  gen.targ(id='dep.ne.1', type='prop',vo='dep.ne==1',vsub=TRUE),
  gen.targ(id='dep.ne.2+',type='prop',vo='dep.ne>1', vsub=TRUE),
  gen.targ(id='dep.um.3m',type='prop',vo='dep.um>t1y/4',vsub=TRUE,sub='dep.past'),
  gen.targ(id='dep.um.6m',type='prop',vo='dep.um>t1y/2',vsub=TRUE,sub='dep.past'),
  gen.targ(id='dep.um.1y',type='prop',vo='dep.um>t1y',  vsub=TRUE,sub='dep.past'),
  gen.targ(id='dep.um.2y',type='prop',vo='dep.um>t1y*2',vsub=TRUE,sub='dep.past'),
  gen.targ(id='dep.um.5y',type='prop',vo='dep.um>t1y*5',vsub=TRUE,sub='dep.past'),
  gen.targ(id='dep.eRo',  type='pois',vo='dep.past',vt='dep.tto'))
T = sub.targs(T,sub.targ.age,ags=10)

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
  ext = PG[c(1,2,5,6)])

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
ymm = list(raw=c(01,12),dep.now=c(01,12),dep.past=c(05,50)) # gray rects
cmap = lapply(c(het='plasma',Ro='cividis',Rx='viridis'), # colormaps
  function(o){ clr.map.d(option=o) })
cmap$dsc = clr.map.d(option='plasma',direction=-1)

load.grid = function(k,out='dep.now',f=NULL,age=FALSE){
  Y = load.rda(grid.path(k),str('Y.',out))
  if (!age){ Y = subset(Y,age==aall) }
  iR = Y$type == 'pois' # rate out rows
  c3 = c('value','lower','upper') # out value cols
  Y[ iR,c3] = Y[ iR,c3]*100*t1y # rates per 100 PY
  Y[!iR,c3] = Y[!iR,c3]*100     # props as %
  Y$Ro  = round(Y$dep_o.Ri.my*100,1) # per 100 PY
  Y$Rx  = round(Y$dep_x.Ri.my*100,1) # per 100 PY
  Y$het = if.null(Y$dep.Ri.het,0)            # shorthand
  Y$cor = if.null(Y$dep.cov,0)               # shorthand
  Y$RRp = if.null(Y$RR.dep_o.dep_p,1)        # shorthand
  Y$dsc = if.null(Y$dsc.dep_x.dep_u,Inf)/t1y # shorthand
  Y[f] = lapply(Y[f],as.factor) # Y[f] -> factors
  return(Y)
}

# grid subsets for plotting
PGi = list(
  Ro  = seq(1,10,2),
  Rx  = seq(50,300,50),
  het = c(0,1,2,3,4,5),
  cor = c(-.6,-.3, 0,+.3,+.6),
  RRp = 1+c(0,1,2,3,4,5),
  dsc = c(Inf,1,.8,.6,.4,.2))
PGii = list(Ro=c(2,5,8),Rx=c(100,200,300),cor=c(-.6,0,+.6))
PG2 = list(Ro=c(3,7),Rx=c(100,200))
PG1 = list(Ro=5,Rx=150,het=0,cor=0,RRp=1,dsc=Inf)
cor.lab = c('–0.9'=-.9,'–0.6'=-.6,'–0.3'=-.3,'0'=0,
            '+0.3'=+.3,'+0.6'=+.6,'+0.9'=+.9)
dsc.max = 1.2

# aes labels
l = list(
  Ro  = 'Mean~onset rate~(per 100 PY)',
  Rx  = 'Mean~recov rate~(per 100 PY)',
  het = 'Rate CV',
  cor = 'Rate~correlation',
  RRp = 'RR relapse~vs onset',
  dsc = 'Recovery~waning~half-life~(years)',
  age = 'Age (years)',
  dep.now  = 'Current~MDD~prevalence (%)',
  dep.past = 'Lifetime~MDD~prevalence (%)')
l$rel = 'Relative~MDD~prevalence'
l$raw = l$dep.now # HACK

# aes label utils
hom = function(s){ gsub('Mean~(.)','\\U\\1',s,perl=TRUE) }
grp = function(s){ gsub('~','\n',s) }
axi = function(s){ gsub('~',' ',s) }
fct = function(s){ ss = strsplit(axi(s),' \\(|\\)')[[1]]; ss[len(ss)+1] = '';
  str.lab(str(' ',ss[1],': '),str(' ',ss[2])) }
exact.fun = list(
  dep.now  = function(o,x){ k=o+x; P = 100*(adur-(1-exp(-adur*k))/k)*o/adur/k },
  dep.past = function(o,x){ k=o;   P = 100*(adur-(1-exp(-adur*k))/k)*o/adur/k })

run.exact = function(Y,n=1e5){
  Y = subset(Y,seed==1)
  qf = het.funs[[P0$het.distr]]$q
  Ye = rbind.lapply(1:nrow(Y),function(i){ Yi = Y[i,]
    fi = exact.fun[[Yi$out]]
    R = copula(n,covs=if.null(Yi$dep.cov,0),qfuns=list(o=qf,x=qf),
      o=list(m=Yi$dep_o.Ri.my,het=if.null(Yi$dep.Ri.het,0)),
      x=list(m=Yi$dep_x.Ri.my,het=if.null(Yi$dep.Ri.het,0)))
    Yi$value = mean(fi(R[,1],R[,2]))
    return(Yi)
  })
}

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

add.exact = function(g,Y){
  g = g + geom_point(data=run.exact(Y),shape=21,fill='red',size=1)
}

# -----------------------------------------------------------------------------
# plot apps

plot.eRo = function(){
  Y = subset(load.grid('het',out='dep.eRo',f='het',age=TRUE),Rx==100 & cor==0)
  g = ggplot(subset(Y,het %in% PGi$het),aes(x=Ro,y=value,color=het,fill=het))
  g = plot.core(g,'dep.eRo',ylim=c(0,10)) + cmap$het +
    facet_wrap('age',labeller=fct(l$age)) +
    labs(color=grp(l$het),fill=grp(l$het),
      x='Model input mean onset rate (per 100 PY)',
      y='Observed mean onset rate (per 100 PY)')
  plot.save(g,'depr',uid,'eRo',ext=ext)
}

plot.ever = function(){
  Y = subset(load.grid('het',out='dep.past',f='het'),Rx==100 & cor==0)
  g = ggplot(subset(Y,het %in% PGi$het),aes(x=Ro,y=value,color=het,fill=het))
  g = plot.core(g,'dep.past',ylim=c(0,80)) + cmap$het +
    labs(x=axi(l$Ro),color=grp(l$het),fill=grp(l$het)) +
    geom_abline(intercept=0,slope=25,color='gray',lty='11')
  plot.save(g,'depr',uid,'ever',ext=ext)
  plot.save(add.exact(g,Y),'depr',uid,'ever.v',ext=ext)
}

plot.edur = function(k,clr){
  umap = list('3m'=.25,'6m'=.5,'1y'=1,'2y'=2,'5y'=5)
  Y = rbind.lapply(str('dep.um.',names(umap)),load.grid,k=k,f=clr) # load
  Y = subset(Y,Ro==1 & cor==0 & RRp==1) # subset
  Y$edur = as.numeric(umap[gsub('dep.um.','',Y$out)]) # Y$out -> Y$edur
  Y = rbind(Y,df.ow(subset(Y,edur==1),edur=0,value=100)) # append dummy (edur=0)
  g = ggplot(df.sub(Y,str('Rx %in% PGii$Rx & ',clr,' %in% PGi$',clr)),
    aes.string(x='edur',y='value',color=clr,fill=clr))
  g = plot.core(g,'edur',ylim=c(0,100)) + cmap[[clr]] +
    facet_wrap('Rx',labeller=fct(l$Rx)) +
    labs(color=grp(l[[clr]]),fill=grp(l[[clr]]),
      x='Time sime onset (years)',
      y='Proportion still depressed (%)')
  plot.save(g,'depr',uid,str('edur.',clr),ext=ext)
}

plot.now.hom = function(){
  Y = subset(load.grid('hom'))
  g = ggplot(subset(Y,Rx %in% PGi$Rx),aes(x=Ro,y=value,color=factor(Rx),fill=factor(Rx)))
  g = plot.core(g) + cmap$Rx + labs(x=axi(l$Ro),color=grp(l$Rx),fill=grp(l$Rx))
  plot.save(g,'depr',uid,'now.hom.o',ext=ext)
  g = ggplot(subset(Y,Ro %in% PGi$Ro),aes(x=Rx,y=value,color=factor(Ro),fill=factor(Ro)))
  g = plot.core(g) + cmap$Ro + labs(x=axi(l$Rx),color=grp(l$Ro),fill=grp(l$Ro))
  plot.save(g,'depr',uid,'now.hom.x',ext=ext)
}

plot.now.het = function(){
  Y = subset(load.grid('het',f='Rx'),cor==0)
  g = ggplot(subset(Y,Ro %in% PGii$Ro),aes(x=het,y=value,color=Rx,fill=Rx))
  g = plot.core(g) + facet_wrap('Ro',labeller=fct(l$Ro)) + cmap$Rx +
    labs(x=axi(l$het),color=grp(l$Rx),fill=grp(l$Rx))
  plot.save(g,'depr',uid,'now.het',ext=ext)
}

plot.now.cor = function(){
  Y = subset(load.grid('het'))
  Y$cor.f = factor(Y$cor,cor.lab,names(cor.lab))
  Y1 = subset(Y,cor %in% PGii$cor & Ro %in% PGii$Ro & Rx %in% PGii$Rx)
  g = ggplot(Y1,aes(x=het,y=value,color=factor(Rx),fill=factor(Rx),lty=cor.f))
  g = plot.core(g,ribbon=0) + facet_wrap('Ro',labeller=fct(l$Ro)) + cmap$Rx +
    stat_summary(fun=median,geom='line',lty='solid',size=2/3) + # HACK
    scale_linetype_manual(values=c('11','solid','31')) +
    labs(x=axi(l$het),color=grp(l$Rx),fill=grp(l$Rx),lty=grp(l$cor))
  plot.save(g,'depr',uid,'now.cor.x',ext=ext)
  Y2 = subset(Y,het %in% PGi$het & Ro %in% PG2$Ro & Rx %in% PG2$Rx)
  g = ggplot(Y2,aes(x=cor,y=value,color=factor(het),fill=factor(het)))
  g = plot.core(g) + cmap$het +
    facet_grid('Rx~Ro',labeller=labeller(.rows=fct(l$Rx),.cols=fct(l$Ro))) +
    scale_x_continuous(breaks=cor.lab,labels=names(cor.lab)) +
    labs(x=axi(l$cor),color=grp(l$het),fill=grp(l$het))
  plot.save(g,'depr',uid,'now.cor.cor',ext=ext)
  plot.save(add.exact(g,Y2),'depr',uid,'now.cor.cor.v',ext=ext)
}

plot.now.RRp = function(alt='raw'){
  Y = subset(load.grid('ext',f='Rx'),dsc==Inf)
  if (alt=='rel'){ Y$value = Y$value / Y$value[Y$RRp==PG1$RRp] }
  g = ggplot(subset(Y,Rx %in% PGi$Rx & Ro %in% PGii$Ro),aes(x=RRp,y=value,color=Rx,fill=Rx))
  g = plot.core(g,alt) + facet_wrap('Ro',labeller=fct(l$Ro)) + cmap$Rx +
    labs(x=axi(l$RRp),color=grp(l$Rx),fill=grp(l$Rx))
  plot.save(g,'depr',uid,str('now.RRp.',alt),ext=ext)
}

plot.now.dsc = function(alt='raw'){
  f = function(dsc){ pmin(dsc.max,dsc) }
  Y = subset(load.grid('ext',f='Rx'),RRp==1)
  if (alt=='rel'){ Y$value = Y$value / Y$value[Y$dsc==PG1$dsc] }
  g = ggplot(subset(Y,Rx %in% PGi$Rx & Ro %in% PGii$Ro),aes(x=f(dsc),y=value,color=Rx,fill=Rx))
  g = plot.core(g,alt) + facet_wrap('Ro',labeller=fct(l$Ro)) + cmap$Rx +
    scale_x_reverse(breaks=f(PGi$dsc),labels=PGi$dsc) + # HACK
    labs(x=axi(l$dsc),color=grp(l$Rx),fill=grp(l$Rx))
  plot.save(g,'depr',uid,str('now.dsc.',alt),ext=ext)
}

# -----------------------------------------------------------------------------
# main

run.grid('hom')
run.grid('het')
run.grid('ext')

plot.eRo()
plot.ever()
plot.edur('het','het')
plot.edur('ext','dsc')
plot.now.hom()
plot.now.het()
plot.now.cor()
plot.now.RRp('raw')
plot.now.RRp('rel')
plot.now.dsc('raw')
plot.now.dsc('rel')
