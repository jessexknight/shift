source('utils.r')
source('sim/plot.r')

# -----------------------------------------------------------------------------
# config + plot stuff

ext = '.png'
uid = '2025-03-07'
plot.1o = list(w1=1.5,h1=2,wo=2,ho=1.2)

l = list( # labels
  Rom  = 'Mean onset rate (per 100 PY)',
  Rxm  = 'Mean recov rate (per 100 PY)',
  cRom = 'Mean onset rate\n(per 100 PY)',
  cRxm = 'Mean recov rate\n(per 100 PY)',
  Rcv  = 'Rate CV',
  cor  = 'Rate\ncorrelation',
  now  = 'Proportion depressed (now)',
  ever = 'Proportion depressed (ever)')

plot.common = function(g,cmap,yref=Inf){
  g = g + clr.map.d(option=cmap) +
    geom_hline(yintercept=yref,color='#999',lty='22') +
    scale_y_continuous(trans='log10',limits=c(.001,1)) +
    scale_x_continuous(trans='log10')
  g = plot.clean(g)
}

plot.ever = function(){
  grid = list(Rom=leq(-3,-1,.1,sig=2),Rcv=c(0,.1,.3,1,3),Rxm=.1,cor=0)
  X = grid.apply(grid,model,.rbind=TRUE)
  g = ggplot(X,aes(x=Rom,y=ever,color=as.factor(Rcv))) +
    labs(y=l$ever,x=l$Rom,color=l$Rcv) +
    geom_line()
  g = plot.common(g,cmap='inferno',yref=.10)
  plot.save(g,'cseb',uid,'ever',ext=ext)
}

plot.now.hom = function(){
  grid = list(Rcv=0,cor=0,
    Rom=leq(-3,-1,.1,sig=2),
    Rxm=c(.01,.03,.1,.3,1))
  X = grid.apply(grid,model,.rbind=TRUE)
  g = ggplot(X,aes(x=Rom,y=now,color=as.factor(Rxm))) +
    labs(y=l$now,x=l$Rom,color=l$cRxm) +
    geom_line()
  g = plot.common(g,cmap='viridis',yref=.02)
  plot.save(g,'cseb',uid,'now.hom',ext=ext)
}

plot.now.het = function(cor=FALSE){
  grid = list(Rcv=leq(-2,.5,.2),cor=c(-.5,0,+.5),
    Rom=c(.001,.003,.01,.03,.1),
    Rxm=c(.01, .03, .1, .3, 1))
  if (!cor){ grid$cor = 0 }
  X = grid.apply(grid,model,.rbind=TRUE)
  X$Onset = X$Rom
  g = ggplot(X,aes(x=Rcv,y=now,color=as.factor(Rxm),lty=as.factor(cor))) +
    facet_grid('~Onset',labeller=label_both) +
    scale_linetype_manual(values=c('11','solid','31')[-cor:+cor+2]) +
    labs(y=l$now,x=l$Rcv,color=l$cRxm,lty=l$cor) +
    geom_line()
  g = plot.common(g,cmap='viridis',yref=.02)
  plot.save(g,'cseb',uid,ifelse(cor,'now.het.cor','now.het'),ext=ext)
}

# -----------------------------------------------------------------------------
# analysis

leq = function(...,sig=1,z=NULL){ signif(c(z,10^seq(...)),sig) }

prop = function(Ro,Ru,dt=50){
  # analytic expr: % depressed now/ever, assume p(0) = 0
  p = mean((dt+(exp(-dt*Ru)-1)/Ru)*Ro/Ru/dt)
}

model = function(Rom,Rxm,Rcv,cor,dist='lnorm',n=1e5){
  cv2 = max(Rcv^2,1e-9) # HACK
  Rox = switch(dist, # sampling rates Ro,Rx ~ copula
  'gamma' = copula(n,covs=cor,qfuns=list(o=qgamma,x=qgamma),
    o=list(shape=1/cv2,scale=Rom*cv2),
    x=list(shape=1/cv2,scale=Rxm*cv2)),
  'lnorm' = copula(n,covs=cor,qfuns=list(o=qlnorm,x=qlnorm),
    o=list(meanlog=log(Rom/sqrt(1+cv2)),sdlog=sqrt(log(1+cv2))),
    x=list(meanlog=log(Rxm/sqrt(1+cv2)),sdlog=sqrt(log(1+cv2)))),
  'weibull' = copula(n,covs=cor,qfuns=list(o=qweibull,x=qweibull),
    o=fit.weibull(Rom,cv2),
    x=fit.weibull(Rxm,cv2)))
  data.frame(dist=dist,Rom=Rom,Rxm=Rxm,Rcv=Rcv,cor=cor,
    now  = prop(Rox[,1],Rox[,1]+Rox[,2]),
    ever = prop(Rox[,1],Rox[,1]))
}

# -----------------------------------------------------------------------------
# main

plot.ever()
plot.now.hom()
plot.now.het(cor=FALSE)
plot.now.het(cor=TRUE)
