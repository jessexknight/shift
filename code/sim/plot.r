library('ggplot2')
library('reshape2')

plot.save = function(g,...,size=NULL,ext='.pdf'){
  if (missing(size)){ size = plot.size(g) }
  fname = str(root.path('out','fig',...,create=TRUE),ext)
  status(3,'saving: ',fname)
  ggsave(fname,w=size[1],h=size[2],device=cairo_pdf)
}

plot.size = function(g,w1=2,h1=2,wo=2.5,ho=1){
  # define size from facet grid
  layout = ggplot_build(g)$layout$layout
  size = c(w=wo+w1*max(layout$COL),h=ho+h1*max(layout$ROW))
}

plot.clean = function(g,...){
  g = g + theme_light() + theme(...,
    legend.title=element_text(size=9),
    strip.background=element_rect(fill='gray85'),
    strip.text.x=element_text(color='black'),
    strip.text.y=element_text(color='black'))
}

evt.vars = function(evt){
  if (is.null(evt)){ return(NULL) }
  vars = switch(substr(evt,1,3),
    vio = c('vio.nt'),
    dep = c('dep.now','dep.past'),
    haz = c('haz.now','haz.past'),
    ptr = c('ptr.nt', 'ptr.nw'))
}

plot.common = function(X,...,geom='box'){
  g = ggplot(X,aes(...)) + switch(geom,
      box    = geom_boxplot(outlier.shape=1,outlier.alpha=1,alpha=0),
      violin = geom_violin(draw_quantiles=1:3/4,alpha=0)) +
    facet_grid(var~facet,scales='free') +
    scale_color_viridis_d(guide='none') +
    ylim(c(0,NA))
  g = plot.clean(g)
}

plot.prev = function(Q,vars=NULL,evt=NULL,facet=NULL,strat='.',wts='.'){
  vars = c(vars,evt.vars(evt)) # select vars
  Q$facet = facet.label(Q[facet]) # clean-up facet label
  Q$strat = do.call(interaction,Q[strat]) # strat
  Q = melt(Q,measure=vars,var='var') # wide -> long
  Q = aggr.prev(Q,g=c('seed','var','facet','strat'),w=wts) # agg
  g = plot.common(Q,y=p,x=strat,color=strat) +
    labs(y='Value (population mean)',x=strat,color=strat)
  g = add.label(g,Q,vs='strat',p=0,
    function(Qi){ str(rmed(Qi$k),'\n',rmed(Qi$n),'\n') })
}

plot.rate = function(R,evt=NULL,facet=NULL,strat='.',ref=NA){
  R = subset(R,var==evt) # select evt
  R$facet = facet.label(R[facet]) # clean-up facet label
  R$strat = as.factor(R[[strat]]) # strat
  R.ref = data.frame(value=ref,strat=len(levels(R$strat))/2+.5)
  g = plot.common(R,y=value*365,x=strat,color=strat) +
    labs(y='Rate (per year)',x=strat,color=strat) +
    geom_point(data=R.ref,shape=9,color='red')
  g = add.label(g,R,vs='strat',value=0,
    function(Ri){ str(rmed(Ri$ne),'\n',rmed(Ri$dt/365),'\n') })
  g = add.label(g,R,vs='strat',value=max(R$value)+ref,
    function(Ri){ signif(median(Ri$value/ref),3) })
}

add.label = function(g,X,label.fun,vs=NULL,...){
  strat = c('var','facet',vs)
  X.label = rbind.lapply(split(X,X[strat]),
    function(Xi){ cbind(Xi[1,strat],label=label.fun(Xi),...) })
  g = g + geom_text(data=X.label,aes(label=label),size=2.5,
    position=pos,show.legend=FALSE)
}

add.info = function(g,x){
  info = list.str(x,sig=3,rnd=9)
  g = g + geom_point(aes(alpha=NA),x=0,y=0) +
    scale_alpha_manual(name=info,na.value=0,labels='')
}

facet.label = function(X){
  f = apply(X,1,list.str,sig=3,rnd=9)
  f = factor(f,unique(f))
}

aggr.prev = function(X,y='value',g='.',w='.'){
  Xa = rbind.lapply(split(X,X[g]),function(Xi){
    Xia = Xi[1,g]
    Xia$k = sum(Xi[[w]]*Xi[[y]])
    Xia$n = sum(Xi[[w]])
    Xia$p = Xia$k / Xia$n
    return(Xia)
  })
}

pos  = position_dodge(width=.75)
rmed = function(x){ round(median(x)) }
