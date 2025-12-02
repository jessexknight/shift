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
plot.1o = list(w1=2.1,h1=1.6,wo=1.3,ho=1)
.ext = '.png'
# helpers ---------------------------------------------------------------------
diffe = function(x,d0=NULL,dn=NULL){ c(d0,diff(x),dn) }
clean.df = function(X,ys){
  if.null(X,return(X))
  X = subset(X,type %in% ys)
  y = add.enum(types[ys],'i')
  X$type = factor(X$type,names(y),y)
  return(X) }
# aao =========================================================================
load.aao = function(plot=0){
  E = load.csv('data','pub','McGrath2023'); m = age$max # HACK
  E = subset(E, age<=m & var=='dep')
  E = rbind.fill(E,cbind(aggregate(value~type+age,E,mean),sex='pooled'))
  # print(aggregate(value~type+sex,subset(E,age>=10),mean))
  E = rbind.lapply(.par=0,split(E,E$sex),function(Ei){
    b = Ei$type=='haz'
    Ei$sv = ifelse(b,'HR','CI')
    Ei$sk = as.factor('Extracted')
    h = Ei$value[ b]
    F = Ei$value[!b]
    Ei = rbind(Ei,
      df.ow(Ei[b,],sv='HR',sk='Converted',type='inc',value=h*exp(-cumsum(h))),
      df.ow(Ei[b,],sv='HR',sk='Converted',type='cum',value=1-exp(-cumsum(h))),
      df.ow(Ei[b,],sv='CI',sk='Converted',type='inc',value=diffe(F,dn=NA)),
      df.ow(Ei[b,],sv='CI',sk='Converted',type='haz',value=diffe(-log(1-F),dn=NA)))
  })
  E$sex = factor(E$sex,c('female','male','pooled'))
  E$sv  = factor(E$sv,c('HR','CI'))
  if (plot){
    g = plot.aao(E,ys=c('haz','inc','cum'),f='type~sex',alpha=sk,color=sv) +
      scale_alpha_manual(values=c(1,.2)) +
      labs(alpha='',color='Source') +
      theme(legend.position='top')
    plot.save(g,'depr',uid,'aao.McGrath2023',ext=.ext,size=c(7,6)) }
  E = df.ow(subset(E,sv=='HR'),src='data',ref='McGrath2023')
}
# model -----------------------------------------------------------------------
model.aao = function(m,cv,am=0,as=Inf,wa=Inf,s=1,...,n=1e4,dist='gamma'){
  if (as==Inf){ aHR = s * exp(-wa*(age$vec-am)) * (age$vec > am) }
  else        { aHR = s * exp(-wa*(age$vec-am)) / (1 + exp(-as*(age$vec-am))) }
  H = het.funs[[dist]]$r(n,m=m,het=cv)
  S = matrix(1,n,1+age$n)
  for (a in 1:age$n){ S[,a+1] = S[,a] * exp(-H*aHR[a]*age$u) }
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
plot.aao = function(E,M=NULL,ys=NULL,f='type',...){
  ys = if.null(ys,c('arr','haz','cum','hrr'))
  g = ggplot(clean.df(E,ys),aes(x=age,y=value,...)) +
    facet_grid(f,scales='free_y') +
    labs(x='Age (years)',y='Value') +
    geom_line(data=clean.df(M,ys)) +
    geom_line() +
    lims(x=c(0,age$max),y=c(0,NA))
  g = plot.clean(g)
}
# main ------------------------------------------------------------------------
main.aao = function(){
  E0 = cbind(load.aao(),cv=NA)
  sex = list(
    'pooled' = list(m=.011,cv=0:3,wa=c(.027,.02,.005,-.02),s=c(0.95,0.95,1.10,1.40)),
    'female' = list(m=.014,cv=0:3,wa=c(.030,.02,.000,-.04),s=c(1.00,1.00,1.20,1.50)),
    'male'   = list(m=.008,cv=0:3,wa=c(.025,.02,.010,-.01),s=c(0.90,0.90,1.00,1.20)))
  cfg = list(
    '(a) Frailty'                      = list(am=10,as=Inf,s=1,wa=0),
    '(b) Frailty + ramp'               = list(am=12,as=0.4,s=1,wa=0),
    '(c) Frailty + ramp + decay'       = list(am=12,as=0.4,s=1,wa=.02),
    '(d) Frailty + ramp + decay (fit)' = list(am=12,as=0.4))
  for (s in names(sex)){
    E = subset(E0,sex==s)
    M = rbind.lapply(names(cfg),function(k){
      par.k = ulist(sex[[s]],cfg[[k]],k=factor(k))
      M.k = grid.apply(par.k,model.aao,.rbind=1,.grid=0,.par=0)
    })
    M$value[M$value>1.5] = NA # HACK
    g = plot.aao(E,M,lty=src,color=factor(cv),f='type~k') +
      clr.map.d(option='rocket',na.value='#0cf',end=.7) +
      scale_linetype_manual(values=c('21','solid')) +
      labs(lty='Source',color='Onset\nfrailty\nSD (σ)') + lims(x=c(0,age$max))
    plot.save(g,'depr',uid,str('aao.',s),ext=.ext)
  }
}
# dur =========================================================================
prep.dur = function(plot=0){
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
    Ai = aggr.srv(Ei)
    qqq = approx(Ai$value,Ai$month,3:1/4,ties='mean')$y
    cat(round(qqq,1),'@',str(Ei$ref[1]),'\n')
  })
  A = cbind(aggr.srv(E),ref='*')
  if (plot){
    g = plot.srv.refs(E,A,dt=1/12,t0='onset',shape=epi) + labs(shape='Episode')
    plot.save(g,'depr',uid,'dur.refs',ext=.ext,size=c(7,5)) }
  save.csv(A,'data','pub','dep.edur.aggr')
}
aggr.srv = function(E,g='g',v='value'){
  E$s = E$p * E$n
  E$e = -unlist(lapply(split(E$s,E[g]),diffe,dn=0))
  A = aggregate(cbind(s=s,e=e)~month+type,E,sum)
  A[v] = c(1,cumprod(1-A$e/A$s)[-nrow(A)])
  A[c('e','s')] = NULL
  return(A)
}
plot.srv.refs = function(E,A,dt=1,t0,...){
  E = clean.df(E,c('srv','haz'))
  A = rbind(A,df.ow(A,type='haz',value=diffe(-A$value,dn=NA)/dt))
  A = clean.df(A,c('srv','haz'))
  g = ggplot(E,aes(x=month/12)) +
    facet_grid('type',scales='free_y') +
    geom_point(aes(y=p,color=ref,...)) +
    geom_step(data=A,aes(y=value),color='black') +
    scale_shape_manual(values=c(1,0,2)) +
    scale_x_continuous(breaks=2*0:10) +
    labs(x=str('Time since ',t0,' (years)'),
         y='Value',color='Reference')
  g = plot.clean(g) + lims(y=c(0,NA))
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
  for (d in 1:dur$n){ S[,d+1] = S[,d] * exp(-H*dHR[d]*dur$u) }
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
plot.dur = function(E,M,ys=NULL,f='type',...){
  ys = if.null(ys,c('drr','haz','cum','hrr'))
  g = ggplot(clean.df(E,ys),aes(x=month/12,y=value,...)) +
    facet_grid(f,scales='free_y') +
    labs(x='Time since onset (years)',y='Value') +
    geom_line(data=clean.df(M,ys)) +
    geom_line()
  g = plot.clean(g)
}
# main ------------------------------------------------------------------------
main.dur = function(prep=0){
  if (prep){ prep.dur() }
  E = cbind(load.dur(),cv=NA)
  cfg = list(
    '(a) Frailty only'          = list(m=2.8,cv=0:3/2,wd=Inf),
    '(b) Frailty + decay'       = list(m=2.8,cv=0:3/2,wd=.5),
    '(c) Frailty + decay (fit)' = list(m=2.8,cv=0:3/2,wd=c(.3,.45,4,Inf)))
  M = rbind.lapply(names(cfg),function(k){
    par.k = ulist(cfg[[k]],k=factor(k))
    M.k = grid.apply(par.k,model.dur,.rbind=1,.grid=0,.par=0)
  })
  g = plot.dur(E,M,lty=src,color=factor(cv),f='type~k') +
    clr.map.d(option='mako',na.value='#f60',end=.7) +
    scale_linetype_manual(values=c('11','solid')) +
    labs(lty='Source',color='Recovery\nfrailty\nSD (σ)') + lims(x=c(0,dur$max))
  plot.save(g,'depr',uid,'dur.main',ext=.ext)
}
# int =========================================================================
prep.int = function(plot=0){
  E = load.csv('data','pub','dep.idur.refs')
  E = subset(E,inc==1)
  E$g = cumsum(E$month==0)
  E$ref = factor(E$ref,unique(E$ref))
  E = rbind.lapply(split(E,E$g),function(Ei){
    m = seq(0,max(Ei$month),6)
    mo = na.to.num(Ei$offset[1])
    f  = approxfun(Ei$month,Ei$p)
    Ei = df.ow(Ei[1,],month=m,p=f(m),po=f(m)/f(mo))
  })
  lapply(split(E,E$ref),function(Ei){
    Ai = aggr.srv(Ei)
    qqq = approx(Ai$value,Ai$month,3:1/4,ties='mean')$y
    cat(round(qqq,1),'@',str(Ei$ref[1]),'\n')
  })
  A = cbind(aggr.srv(E),ref='*')
  A$inc = -log(A$value)*12/A$month
  print(subset(A,month %in% c(6,12,24,60,120,180)))
  if (plot){
    g = plot.srv.refs(E,A,dt=1/2,t0='recovery')
    plot.save(g,'depr',uid,'int.refs',ext=.ext,size=c(7,5)) }
  save.csv(A,'data','pub','dep.idur.aggr')
}

# main ========================================================================

load.aao(plot=1)
prep.dur(plot=1)
prep.int(plot=1)
main.aao()
main.dur()
