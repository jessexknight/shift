source('utils.r')
library('ggplot2')
library('reshape2')

ext = '.png'
distrs = c(
  'Gamma'    = 'gamma',
  'Log-Norm' = 'lnorm',
  'Weibull'  = 'weibull',
  'Normal'   = 'norm')

sample.copula = function(n,m,cv,distr,cor,seed=666){
  set.seed(seed)
  qfun = het.funs[[distr]]$q
  pars = list(m=m,het=max(cv,1e-6))
  X = as.data.frame(copula(n,cor,list(u=qfun,v=qfun),u=pars,v=pars))
}

plot.1d = function(){
  q = seq(0,4,.001) + 1e-16 # HACK
  G = list(cv=c(0,.1,.3,1,3),distr=distrs[1:3])
  X = grid.apply(G,function(cv,distr){
    f = het.funs[[distr]]
    data.frame(check.names=0,q=q,
      'Probability Density'   =f$d(q,1,cv),
      'Cumulative Probability'=f$p(q,1,cv))
  },.rbind=1,.cbind=1)
  X$distr = factor(X$distr,distrs,names(distrs))
  X = melt(X,id=c('q',names(G)))
  g = ggplot(X,aes(x=q,y=pmin(6,value),color=factor(cv))) +
    facet_grid('variable~distr',scales='free_y') +
    scale_x_continuous(breaks=0:4) +
    clr.map.d(option='plasma') +
    geom_line() +
    labs(y='Probability',x='Hazard Rate Ratio',color='σ')
  g = plot.clean(g)
  plot.save(g,'toy','het.1d',ext=ext,size=c(7,4))
}

plot.2d = function(distr,tx='raw'){
  bs = list(raw=seq(0,20,5),log=10^seq(-5,1,3))[[tx]]
  tf = list(raw='identity',log='log10')[[tx]]
  G = list(n=1e4,m=1,cv=c(.5,1,2),distr=distr,cor=c(-.5,0,+.5))
  X = grid.apply(G,sample.copula,.rbind=1,.cbind=1)
  X$cv  = str('σ: ',X$cv)
  X$cor = str('ρ: ',X$cor)
  g = ggplot(X,aes(x=u,y=v)) +
    facet_grid('cor~cv') +
    geom_point(size=.2,alpha=.1) +
    scale_x_continuous(trans=tf,breaks=bs,labels=bs) +
    scale_y_continuous(trans=tf,breaks=bs,labels=bs) +
    coord_fixed(xlim=c(1e-6,20),ylim=c(1e-6,20))
  g = plot.clean(g)
  plot.save(g,'toy',str('het.2d.',distr,'.',tx),ext=ext,size=c(6,6))
}

plot.cor = function(){
  meths = c('Pearson (Raw)','Pearson (Log)','Spearman')
  G = list(n=1e4,m=1,cv=0:5,cor=seq(-.9,+.9,.1),distr=distrs)
  Y = grid.apply(G,function(...){
    X = sample.copula(...)
    L = log(X)
    Yi = data.frame(meth=meths,ecor=c(
      cor(X$u,X$v,method='pearson'),
      cor(L$u,L$v,method='pearson'),
      cor(X$u,X$v,method='spearman')))
  },.rbind=1,.cbind=1)
  Y$distr = factor(Y$distr,G$distr,names(G$distr))
  Y$cv = str('σ: ',Y$cv)
  g = ggplot(Y,aes(x=cor,y=ecor,color=meth,lty=meth)) +
    facet_grid('distr~cv') +
    geom_line() +
    coord_fixed() +
    scale_x_continuous(lim=c(-1,+1),breaks=c(-1,0,+1)) +
    scale_y_continuous(lim=c(-1,+1),breaks=c(-1,0,+1)) +
    scale_color_manual(values=c('#0cc','#066','#c06')) +
    scale_linetype_manual(values=c('solid','31','11')) +
    labs(x='Input Correlation (ρ)',y='Sample Correlation',
      lty='Measure',color='Measure')
  g = plot.clean(g)
  plot.save(g,'toy','het.2d.cor',ext=ext,size=c(8,5))
}

plot.1d()
plot.cor()
plot.2d('gamma',  tx='log')
plot.2d('gamma',  tx='raw')
plot.2d('lnorm',  tx='log')
plot.2d('lnorm',  tx='raw')
plot.2d('weibull',tx='log')
plot.2d('weibull',tx='raw')
