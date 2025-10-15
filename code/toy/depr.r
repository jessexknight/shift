source('utils.r')
library('ggplot2')
library('reshape2')
# config ======================================================================
uid = Sys.Date()
age = list(max=60,u=.1);  age$vec=seq(age$u,age$max,age$u); age$n=len(age$vec)
dur = list(max= 5,u=.01); dur$vec=seq(dur$u,dur$max,dur$u); dur$n=len(dur$vec)
types = list(
  arr='Applied age HRR',
  drr='Applied duration HRR',
  haz='Hazard rate',
  inc='Incidence rate',
  srv='Cumulative survival',
  cum='Cumulative incidence',
  hrr='Apparent frailty HRR')
plot.1o = list(w1=2,h1=1.5,wo=1.5,ho=1)
.sub = 'type!="inc"'
.ext = '.png'
# helpers ---------------------------------------------------------------------
diffe = function(x,d0=NULL,dn=NULL){ c(d0,diff(x),dn) }
clean.df = function(X,sub=.sub){
  if.null(X,return(X))
  X = df.sub(X,sub)
  X$type = factor(X$type,names(types),types)
  return(X) }
# aao =========================================================================
load.aao = function(plot=0){
  E = load.csv('data','pub','McGrath2023'); m = age$max # HACK
  E = subset(E, age<=m & var=='dep')
  E = rbind.lapply(.par=0,split(E,E$sex),function(Ei){
    b = Ei$type=='haz'
    Ei$sv = ifelse(b,'HR','CI')
    Ei$sk = 'Extracted'
    h = Ei$value[ b]
    F = Ei$value[!b]
    Ei = rbind(Ei,
      df.ow(Ei[b,],sv='HR',sk='Converted',type='inc',value=h*exp(-cumsum(h))),
      df.ow(Ei[b,],sv='HR',sk='Converted',type='cum',value=1-exp(-cumsum(h))),
      df.ow(Ei[b,],sv='CI',sk='Converted',type='inc',value=diffe(F,dn=NA)),
      df.ow(Ei[b,],sv='CI',sk='Converted',type='haz',value=diffe(-log(1-F),dn=NA)))
  })
  if (plot){
    g = plot.aao(E,f='type~sex',lty=sk,color=sv,sub='TRUE') +
      scale_linetype_manual(values=c('21','solid')) +
        labs(lty='',color='Source')
    plot.save(g,'depr',uid,'aao.McGrath2023',ext=.ext) }
  E = df.ow(subset(E,sv=='HR'),src='data',ref='McGrath2023')
}
# model -----------------------------------------------------------------------
model.aao = function(m,cv,am=0,as=Inf,wa=Inf,s=1,...,n=1e4,dist='gamma'){
  if (as==Inf){ aHR = s * exp(-wa*(age$vec-am)) * (age$vec > am) }
  else        { aHR = s * exp(-wa*(age$vec-am)) / (1 + exp(-as*(age$vec-am))) }
  H = het.funs[[dist]]$r(n,m=m,het=cv)
  S = matrix(1,n,1+age$n)
  for (a in 1:age$n){ S[,a+1] = S[,a] * (1-H*aHR[a]*age$u) }
  F = 1-colMeans(S[,-1])
  f = diffe(F,dn=NA)/age$u
  h = f/(1-F)
  M = data.frame(src='model',age=age$vec,type=NA,value=NA,
    m=m,cv=cv,am=am,as=as,wa=wa,n=n,dist=dist,...)
  M = rbind(
    df.ow(M,type='arr',value=aHR),
    df.ow(M,type='haz',value=h),
    df.ow(M,type='inc',value=f),
    df.ow(M,type='cum',value=F),
    df.ow(M,type='hrr',value=pmin(1,h/m/aHR)))
}
plot.aao = function(E,M=NULL,f='type',...,sub=.sub){
  g = ggplot(clean.df(E,sub),aes(x=age,y=value,...)) +
    facet_grid(f,scales='free_y') +
    labs(x='Age (years)',y='Value') +
    geom_line(data=clean.df(M,sub)) +
    geom_line() +
    lims(x=c(0,age$max),y=c(0,NA))
  g = plot.clean(g)
}
# main ------------------------------------------------------------------------
main.aao = function(){
  E0 = cbind(load.aao(),cv=NA)
  sex = list(
    'female' = list(m=.014,cv=0:3,wa=c(.03,.02,.00,-.01),s=c(1.00,1.00,1.20,1.20)),
    'male'   = list(m=.008,cv=0:3,wa=c(.025,.02,.01,.00),s=c(0.90,0.90,1.00,1.20)))
  cfg = list(
    'Frailty'                       = list(am=10,as=Inf,s=1,wa=0),
    'Frailty + ramp'                = list(am=12,as=0.4,s=1,wa=0),
    'Frailty + ramp + fixed waning' = list(am=12,as=0.4,s=1,wa=.02),
    'Frailty + ramp + fit waning'   = list(am=12,as=0.4))
  for (s in names(sex)){
    E = subset(E0,sex==s)
    M = rbind.lapply(names(cfg),function(k){
      par.k = ulist(sex[[s]],cfg[[k]],k=factor(k))
      M.k = grid.apply(par.k,model.aao,.rbind=1,.grid=0,.par=0)
    })
    g = plot.aao(E,M,lty=src,color=factor(cv),f='type~k') +
      clr.map.d(option='plasma',na.value='#0cc') +
      scale_linetype_manual(values=c('21','solid')) +
      labs(lty='Source',color='Frailty\nSD (σ)') + lims(x=c(0,age$max))
    g = add.sublabs(g,clean.df(M[c('k','type')]),'tl')
    plot.save(g,'depr',uid,str('aao.',s),ext=.ext)
  }
}
# dur =========================================================================
prep.dur = function(plot=1){
  E = load.csv('data','pub','dep.edur.refs')
  E = subset(E,inc==1)
  E$g = cumsum(E$month==0)
  E$ref = factor(E$ref,unique(E$ref))
  E$epi = ifelse(E$epi=='1+','1+',ifelse(E$epi=='1','1','2+'))
  E = rbind.lapply(split(E,E$g),function(Ei){
    m = seq(0,max(Ei$month))
    Ei = df.ow(Ei[1,],month=m,p=approx(Ei$month,Ei$p,m)$y)
  })
  lapply(split(E,E$ref),function(Ei){
    Ai = aggr.dur(Ei)
    qqq = approx(Ai$value,Ai$month,3:1/4,ties='mean')$y
    cat(round(qqq,1),'@',str(Ei$ref[1]),'\n')
  })
  A = cbind(aggr.dur(E,'g'),ref=factor('*'),type='srv')
  if (plot){ plot.dur.refs(rbind.fill(E,A)) }
  save.csv(A,'data','pub','dep.edur.aggr')
}
aggr.dur = function(E,g='g',v='value'){
  E$s = E$p * E$n
  E$e = -unlist(lapply(split(E$s,E[g]),diffe,dn=0))
  A = aggregate(cbind(s=s,e=e)~month,E,sum)
  A[v] = c(1,cumprod(1-A$e/A$s)[-nrow(A)])
  A[c('e','s')] = NULL
  return(A)
}
plot.dur.refs = function(E){
  E = rbind(cbind(subset(E,month <= 24),f=factor('First 24 months')),
    cbind(E,f=factor('All available follow-up')))
  g = ggplot(E,aes(x=month/12)) +
    facet_wrap('f',ncol=1,scales='free_x') +
    geom_point(data=subset(E,ref!='*'),aes(y=100*p,shape=epi,color=ref)) +
    geom_step(data=subset(E,ref=='*'),aes(y=100*value),color='black') +
    scale_shape_manual(values=c(1,0,2)) +
    labs(x='Time since onset (years)',
         y='Proportion still depressed (%)',
         color='Reference',shape='Episode')
  g = plot.clean(g)
  plot.save(g,'depr',uid,'dur.refs',ext=.ext,size=c(7,6))
}
load.dur = function(){
  E = load.csv('data','pub','dep.edur.aggr')
  F = 1-E$value
  E = rbind(
    df.ow(E,type='cum',value=F),
    df.ow(E,type='inc',value=12*diffe(F,d0=0)),
    df.ow(E,type='haz',value=12*diffe(F,d0=0)/(1-F)))
  E = cbind(E,src='data')
}
# model -----------------------------------------------------------------------
model.dur = function(m,cv,wd=Inf,s=1,...,n=1e4,dist='gamma'){
  dHR = s * na.to.num(wd/(wd+dur$vec),1)
  H = het.funs[[dist]]$r(n,m=m,het=cv)
  S = matrix(1,n,1+dur$n)
  for (d in 1:dur$n){ S[,d+1] = S[,d] * (1-H*dHR[d]*dur$u) }
  F = 1-colMeans(S[,-1])
  f = diffe(F,dn=NA)/dur$u
  h = f/(1-F)
  M = data.frame(src='model',month=dur$vec*12,type=NA,value=NA,
    m=m,cv=cv,wd=wd,n=n,dist=dist,...)
  M = rbind(
    df.ow(M,type='drr',value=dHR),
    df.ow(M,type='haz',value=h),
    df.ow(M,type='inc',value=f),
    df.ow(M,type='cum',value=F),
    df.ow(M,type='hrr',value=pmin(1,h/m/dHR)))
}
plot.dur = function(E,M,f='type',...,sub=.sub){
  g = ggplot(clean.df(E,sub),aes(x=month/12,y=value,...)) +
    facet_grid(f,scales='free_y') +
    labs(x='Time since onset (years)',y='Value') +
    geom_line(data=clean.df(M,sub)) +
    geom_line()
  g = plot.clean(g)
}
# main ------------------------------------------------------------------------
main.dur = function(prep=0){
  if (prep){ prep.dur() }
  E = cbind(load.dur(),cv=NA)
  cfg = list(
    'Frailty only'           = list(m=2.8,cv=0:3/2,wd=Inf),
    'Frailty + fixed waning' = list(m=2.8,cv=0:3/2,wd=.5),
    'Frailty + fit waning'   = list(m=2.8,cv=0:3/2,wd=c(.3,.45,4,Inf)))
  M = rbind.lapply(names(cfg),function(k){
    par.k = ulist(cfg[[k]],k=factor(k))
    M.k = grid.apply(par.k,model.dur,.rbind=1,.grid=0,.par=0)
  })
  g = plot.dur(E,M,lty=src,color=factor(cv),f='type~k') +
    clr.map.d(option='plasma',na.value='#0cc') +
    scale_linetype_manual(values=c('11','solid')) +
    labs(lty='Source',color='Frailty\nSD (σ)') + lims(x=c(0,dur$max))
  g = add.sublabs(g,clean.df(M[c('k','type')]),'bl')
  plot.save(g,'depr',uid,'dur.main',ext=.ext)
}
# main ========================================================================
main.aao()
main.dur()
