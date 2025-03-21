library('ggplot2')
library('reshape2')

evt.vars = function(evt){
  if (is.null(evt)){ return(NULL) }
  vars = switch(substr(evt,1,3),
    vio = c('vio.nt','vio.past'),
    dep = c('dep.now','dep.past'),
    haz = c('haz.now','haz.past'),
    ptr = c('ptr.nt', 'ptr.nw'))
}

plot.common = function(X,...,geom='box'){
  g = ggplot(X,aes(...)) + switch(geom,
      box    = geom_boxplot(outlier.shape=1,outlier.alpha=1,alpha=0),
      violin = geom_violin(draw_quantiles=1:3/4,alpha=0)) +
    facet_grid('var~facet',scales='free') +
    scale_color_viridis_d(guide='none') +
    ylim(c(0,NA))
  g = plot.clean(g)
}

plot.mean = function(Q,vars=NULL,evt=NULL,facet=NULL,strat='.',wts='.'){
  vars = c(vars,evt.vars(evt)) # select vars
  Q$facet = facet.label(Q[facet]) # clean-up facet label
  Q$strat = interac(Q[strat]) # strat
  Q = melt(Q,measure=vars,var='var') # wide -> long
  Q = aggr.mean(Q,g=c('seed','var','facet','strat'),w=wts) # agg
  g = plot.common(Q,y=value,x=strat,color=strat) +
    labs(y='Value (population mean)',x=ilab(strat),color=ilab(strat))
  g = add.label(g,Q,vs='strat',value=0,
    function(Qi){ str(rmed(Qi$num),'\n',rmed(Qi$den),'\n') })
}

plot.rate = function(R,evt=NULL,facet=NULL,strat='.',ref=NA){
  R = subset(R,var==evt) # select evt
  R$facet = facet.label(R[facet]) # clean-up facet label
  R$strat = interac(R[strat]) # strat
  R.ref = data.frame(value=ref,strat=mean(as.numeric(R$strat)))
  g = plot.common(R,y=value*365,x=strat,color=strat) +
    labs(y='Rate (per year)',x=ilab(strat),color=ilab(strat)) +
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

add.info = function(g,info,size=8){
  g = g + guides(custom=guide_custom(grid::textGrob(
    label=info,gp=grid::gpar(fontsize=size,lineheight=1))))
}

facet.label = function(X){
  f = apply(X,1,list.str,sig=3,rnd=9,join='\n')
  f = factor(f,unique(f))
}

aggr.mean = function(X,y='value',g='.',w='.'){
  Xa = rbind.lapply(split(X,X[g]),function(Xi){
    Xia = Xi[1,g]
    Xia$num = sum(Xi[[w]]*Xi[[y]])
    Xia$den = sum(Xi[[w]])
    Xia$value = Xia$num / Xia$den
    return(Xia)
  })
}

pos  = position_dodge(width=.75)
rmed = function(x){ round(median(x)) }
ilab = function(x){ str(x,collapse=.isep) }
aes.string = function(...){ suppressWarnings(aes_string(...)) }
