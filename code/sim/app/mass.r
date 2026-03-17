source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')
uid = '2026-03-13'
.k     = cli.arg('.k','RR2.ad.base')
.b     = cli.arg('.b', 1)
.nb    = cli.arg('.nb',1)
.debug = cli.arg('.debug',0)

# e = exposure "dep" (abuse)
# o = outcome  "haz" (depression)

# -----------------------------------------------------------------------------
# params & grid

P0 = list(
  pop.type = 'open',
  n.pop = xdf(1000,10000),
  n.dur = 1, dtz = 7,
  het.distr = 'gamma',
  run = get.run.par(c('dep','haz'),u=0))

z = 1e-12
G = name.list(key='i',
  list(i='RRo', id='RR.haz_o.dep_w',v=1,   vg=  c(1,2,3,4,8)),
  list(i='RRx', id='RR.haz_x.dep_w',v=1,   vg=1/c(1,2,3,4,8)),
  list(i='eRo', id='dep_o.Ri.my',   v=.01, vg=c(.003,.01,.03,.1)),
  list(i='eRx', id='dep_x.Ri.my',   v=1,   vg=c(z,.03,.1,.3,1,3)),
  list(i='eHo', id='dep_o.Ri.het',  v=2,   vg=c(0,1,2,3)),
  list(i='eHx', id='dep_x.Ri.het',  v=1,   vg=c(0,1,2,3)),
  list(i='ecv', id='dep.cov',       v=-.5, vg=c(-.5,0,+.5)),
  list(i='ep',  id='dep.prev',      v=.1,  vg=c(.03,.1,.2,.3)),
  list(i='oRo', id='haz_o.Ri.my',   v=.01, vg=c(.003,.01,.03,.1)),
  list(i='oRx', id='haz_x.Ri.my',   v=3,   vg=c(z,.03,.1,.3,1,3)),
  list(i='oHo', id='haz_o.Ri.het',  v=2,   vg=c(0,1,2,3)),
  list(i='oHx', id='haz_x.Ri.het',  v=1,   vg=c(0,1,2,3)),
  list(i='ocv', id='haz.cov',       v=-.5, vg=c(-.5,0,+.5)),
  list(i='seed',id='seed',          v=NA,  vg=xdf(1:7,1:21)),
  list(i='ek',  id='exp.case',v='adult',vg=NA))

Gid = lapply(G,`[[`,'id')
G0 = lapply(G,`[[`,'v')
Gi = function(i,...){ ulist(G0,lapply(G[c('seed',i)],`[[`,'vg'),...) }
PG = function(Gk,...){ ulist(P0,set.names(Gk,Gid[names(Gk)]),...) }

Gk = list()
# childhood exposure
Gk$RRo.ch.base = Gi(ek='child',c('RRo'))
Gk$RRx.ch.base = Gi(ek='child',c('RRx'))
Gk$RR2.ch.base = Gi(ek='child',c('RRo','RRx'))
Gk$RRo.ch.ep   = Gi(ek='child',c('RRo','ep'))
Gk$RRo.ch.oRo  = Gi(ek='child',c('RRo','oRo','oHo'))
Gk$RRo.ch.oRx  = Gi(ek='child',c('RRo','oRx','oHx'))
Gk$RRo.ch.oR2  = Gi(ek='child',c('RRo','oRo','oRx','ocv'),oHo=2,oHx=1)
# adulthood exposure
Gk$RRo.ad.base = Gi(ek='adult',c('RRo'))
Gk$RRx.ad.base = Gi(ek='adult',c('RRx'))
Gk$RR2.ad.base = Gi(ek='adult',c('RRo','RRx'))
Gk$RRo.ad.eRo  = Gi(ek='adult',c('RRo','eRo','eHo'))
Gk$RRo.ad.eRx  = Gi(ek='adult',c('RRo','eRx','eHx'))
Gk$RRo.ad.eR2  = Gi(ek='adult',c('RRo','eRo','eRx','ecv'),eHo=2,eHx=1)
Gk$RRo.ad.oRo  = Gi(ek='adult',c('RRo','oRo','oHo'))
Gk$RRo.ad.oRx  = Gi(ek='adult',c('RRo','oRx','oHx'))
Gk$RRo.ad.oR2  = Gi(ek='adult',c('RRo','oRo','oRx','ocv'),oHo=2,oHx=1)
Gk$RRo.ad.lhs  = Gi(ek='adult',seed=c(1,1e9),lhs=xdf(1e1,1e5),
  c('RRo','eRo','eRx','eHo','eHx','ecv','oRo','oRx','oHo','oHx','ocv'))
# for (k in names(Gk)){ status(3,k,': ',prod(lens(Gk[[k]]))) } # for hpc gen

apply.case = function(P){
  if (P$exp.case=='child'){
    P$init.inds = function(I,P){
      I$dep_o.Ri = ifelse(runif(P$n.tot) < P$dep.prev,Inf,0)
      I$dep_x.Ri = 0
      return(I) }}
  return(P)
}

get.lhs = function(Gi,seed=666){
  set.seed(seed)
  is = lens(Gi) > 1
  Gi[is] = as.data.frame(qunif(
    lhs::randomLHS(Gi$lhs,sum(is)),
    rep(unlist(lapply(Gi[is],min)),each=Gi$lhs),
    rep(unlist(lapply(Gi[is],max)),each=Gi$lhs) ))
  return(Gi)
}

# -----------------------------------------------------------------------------
# targets / outcomes

T = name.list(key='id',
  gen.targ(id='e.w',   type='prop',vo='dep.now' ),
  gen.targ(id='e.p',   type='prop',vo='dep.past'),
  gen.targ(id='o.w',   type='prop',vo='haz.now' ),
  gen.targ(id='o.p',   type='prop',vo='haz.past'),
  gen.targ(id='or.ww', type='OR',  ve='dep.now', vo='haz.now' ),
  gen.targ(id='or.wp', type='OR',  ve='dep.now', vo='haz.past'),
  gen.targ(id='or.pw', type='OR',  ve='dep.past',vo='haz.now' ),
  gen.targ(id='or.pp', type='OR',  ve='dep.past',vo='haz.past'),
  gen.targ(id='pr.ww', type='PR',  ve='dep.now', vo='haz.now' ),
  gen.targ(id='pr.pp', type='PR',  ve='dep.now', vo='haz.past'),
  gen.targ(id='pr.pw', type='PR',  ve='dep.past',vo='haz.now' ),
  gen.targ(id='pr.wp', type='PR',  ve='dep.past',vo='haz.past'),
  gen.targ(id='or.wwa',type='OR',  ve='dep.now', vo='haz.now' ,va1='age'),
  gen.targ(id='or.wpa',type='OR',  ve='dep.now', vo='haz.past',va1='age'),
  gen.targ(id='or.pwa',type='OR',  ve='dep.past',vo='haz.now' ,va1='age'),
  gen.targ(id='or.ppa',type='OR',  ve='dep.past',vo='haz.past',va1='age'),
  gen.targ(id='pr.wwa',type='PR',  ve='dep.now', vo='haz.now' ,va1='age'),
  gen.targ(id='pr.ppa',type='PR',  ve='dep.now', vo='haz.past',va1='age'),
  gen.targ(id='pr.pwa',type='PR',  ve='dep.past',vo='haz.now' ,va1='age'),
  gen.targ(id='pr.wpa',type='PR',  ve='dep.past',vo='haz.past',va1='age'))

Tid = list(
  ORx = filter.names(T,'or...a$'),
  PRx = filter.names(T,'pr...a$'),
  XRx = filter.names(T,'.r...a$'),
  XRw = filter.names(T,'.r.wwa'))

# -----------------------------------------------------------------------------
# run sims & save/load

grid.path = function(k,.save=FALSE){
  path = hash.path(PG(Gk[[k]]),'data','sim','mass',uid,k,.save=.save)
}

run.one = function(...,.par=0){
  P1 = PG(list(...),fun=apply.case)
  Ps = get.pars.grid(P1,.par=.par)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  Q  = srv.apply(Ms,.par=.par)
  Y  = srv.targs(Q,T)
  Y[c('seed','targ.mu','targ.se','ll')] = NULL
  row.names(Y) = NULL
  return(Y)
}

run.grid = function(k){
  lhs = len(Gk[[k]]$lhs)
  GR = { if (lhs) get.lhs(Gk[[k]]) else Gk[[k]] }
  Y = grid.apply(GR,run.one,.grid=!lhs,.batch=.b,.nbatch=.nb,
    .rbind=1,.cbind=1,.log=3)
  save.rds(Y,grid.path(k,.save=TRUE),str('b',.nb),str('Y.',.b))
}

merge.batch = function(k){
  Y = rbind.lapply(1:.nb,function(b){
    Yb = load.rds(grid.path(k),str('b',.nb),str('Y.',b)) })
  save.rds(Y,grid.path(k),'Y')
}

load.grid = function(k,i='or.wwa',f=NULL){
  Y = load.rds(grid.path(k),'Y')
  Y = subset(Y,id %in% i)
  Y[c('ve','vo','te','to','dt')] = NULL
  v = c('eRo','eHo','eRx','eHx','oRo','oHo','oRx','oHx')
  Y[v][Y[v]<=z] = 0 # HACK
  Y$bias     = ifelse(Y$type=='prop',NA,Y$value/(Y$RRo/Y$RRx))
  Y$bias.adj = ifelse(Y$type=='prop',NA,(Y$value-1)/(Y$RRo/Y$RRx-1))
  Y$mass = factor(substr(Y$id,1,2),names(fl$mass),fl$mass)
  Y$erep = factor(substr(Y$id,4,4),names(fl$rep),fl$rep)
  Y$orep = factor(substr(Y$id,5,5),names(fl$rep),fl$rep)
  Y$ek   = factor(Y$ek,names(fl$ek),fl$ek)
  Y$RRx = round(Y$RRx,3)
  Y[f] = lapply(Y[f],as.factor)
  return(Y)
}

fl = list( # factor levels
  mass = c(or='OR',pr='PR'),
  rep = c(w='current',p='lifetime'),
  ek = c(adult='adulthood',child='childhood'))

reps = c('erep','orep')

# -----------------------------------------------------------------------------
# exact math @ ek=ch

avec = seq(.1,adur,.1)

efun = list(
  prev.age = function(o,x){ k=o+x; p = colMeans((1-exp(-outer(k,avec)))*o/k) },
  prev = function(o,x){ k=o+x; p = mean((adur-(1-exp(-k*adur))/k)*o/adur/k) },
  OR = function(p0,p1){ p1*(1-p0)/p0/(1-p1) },
  PR = function(p0,p1){ p1/p0 })

run.ch.exact = function(Y,n=1e4,age=FALSE){
  Y = subset(Y,seed==1)
  qf = het.funs[[P0$het.distr]]$q
  Ye = rbind.lapply(1:nrow(Y),function(i){ Yi=Y[i,]
    R = copula(n,covs=Yi$ocv,qfuns=list(o=qf,x=qf),
      o=list(m=Yi$oRo+z,het=Yi$oHo+z),
      x=list(m=Yi$oRx+z,het=Yi$oHx+z))
    R = round(R,12)+z # HACK
    prev.fun = ifelse(age,efun$prev.age,efun$prev)
    fx = switch(str(Yi$orep),current=1,lifetime=0)
    p0 = prev.fun(R[,1],       fx*R[,2])
    p1 = prev.fun(R[,1]*Yi$RRo,fx*R[,2]*Yi$RRx)
    if (age){ Yi = df.ow(Yi,a=avec,p0=p0,p1=p1,OR=efun$OR(p0,p1),PR=efun$PR(p0,p1)) }
    else {    Yi = df.ow(Yi,value=efun[[Yi$type]](p0,p1)) }
  })
}

# -----------------------------------------------------------------------------
# plot utils

labels = list(
  mass = 'Measure of~association',
  bias = 'Bias~vs~HR',
  ek   = 'Exposure',
  OR   = 'OR:~abuse and~depression',
  PR   = 'PR:~abuse and~depression',
  RRo  = 'HR: depression~onset~while abused',
  RRx  = 'HR: depression~recovery~while abused',
  iRRx = '1/HR: depression recovery~while abused',
  ep   = 'Childhood~abuse~prevalence~(%)',
  op   = 'Depression~prevalence~(%)',
  eRo  = 'Abuse~onset rate~(per 100 PY)',
  eRx  = 'Abuse~recovery rate~(per 100 PY)',
  oRo  = 'Depression~onset rate~(per 100 PY)',
  oRx  = 'Depression~recovery rate~(per 100 PY)',
  eHo  = 'Abuse~onset~frailty SD',
  eHx  = 'Abuse~recovery~frailty SD',
  oHo  = 'Depression~onset~frailty SD',
  oHx  = 'Depression~recovery~frailty SD',
  erep = 'Abuse~reporting',
  orep = 'Depression~reporting',
  age  = 'Age~(years)')

ll = function(i,grp=0){
  if (is.null(i)) return(i)
  gsub('~',ifelse(grp,'\n',' '),if.null(labels[[i]],'')) }

fct = function(s,enum=NULL){
  ss = strsplit(gsub('~',' ',s),' \\(|\\)')[[1]]; ss[len(ss)+1] = '';
  str.lab(str(' ',ss[1],': '),str(' ',ss[2]),enum=enum) }

fct_grid = function(x='.',y='.',ex=NULL,ey=NULL){
  facet_grid(str(y,'~',x),labeller=labeller(
    .cols=fct(labels[[x]],enum=ex),
    .rows=fct(labels[[y]],enum=ey))) }

sublabs = def.args(add.sublabs,fmt='i',dx=.5,size=3,family='Alegreya Sans')

cmap = lapply(list(RRo='viridis',RRx='inferno',ep='mako',
  eRo='mako',  eHo='mako',  eRx='mako',  eHx='mako',
  oRo='rocket',oHo='rocket',oRx='rocket',oHx='rocket'),
  function(o){ clr.map.d(option=o,end=.7) })
cmap$mass = clr.map.m(c('#c06','#0cc'))
cmap$null = clr.map.m('#000')

ltys = lapply(list(
    v2=c('solid','22'),
    v3=c('solid','31','11'),
    v4=c('solid','41','21','11')),
  function(v){ scale_linetype_manual(values=v) })

scales = list(
  mass = scale.y.cts(breaks=seq(0,10, 2),limits=c(0,10)),
  RRo  = scale.x.cts(breaks=seq(0, 8, 2),limits=c(0, 8)),
  bias = scale.y.cts(breaks=c(.03,.1,.3,1,3),limits=c(.03,3),trans='log10'))
scales$iRRx = scales$RRo
scales$OR = scales$PR = scales$mass

plot.core = function(x,y,clr=NULL,lty=NULL,da=1,ra=1/5,ci=.95){ list(
  scales[[x]],scales[[y]],cmap[[if.null(clr,'null')]],
  geom_hline(lty='11',color='#999',yintercept=1),
  geom_abline(lty='11',color='#999',alpha=da),
  labs(x=ll(x),y=ll(y),lty=ll(lty,1),color=ll(clr,1),fill=ll(clr,1)),
  stat_summary(geom='ribbon',color=NA,alpha=ra,
    fun.min=qfun((1-ci)/2),fun.max=qfun(1-(1-ci)/2)),
  stat_summary(geom='line',fun=mean),
  plot.clean(font='Alegreya Sans',legend.spacing.y=unit(-1,'mm'))
)}

add.ch.exact = function(Y){
  geom_point(data=run.ch.exact(Y),shape=21,fill='#fc0',size=1)
}

plot.1o = list(w1=2,h1=1.6,wo=1.5,ho=1)

plot.save.i = function(g,...,size=NULL,ext='.png'){
  plot.save(g,'mass',uid,...,ext=ext,size=size)
}

# -----------------------------------------------------------------------------
# objective plots

plot.obj.1 = function(){
  Y = rbind(load.grid('RR2.ad.base',i=Tid$XRw),
            load.grid('RR2.ch.base',i=Tid$XRw))
  g = ggplot(subset(Y,RRx==1),aes(x=RRo,y=value,lty=ek,color=mass,fill=mass)) +
    plot.core('RRo','mass','mass','ek')
  plot.save.i(g,'RRo.base')
  g = ggplot(subset(Y,RRo==1),aes(x=1/RRx,y=value,lty=ek,color=mass,fill=mass)) +
    plot.core('iRRx','mass','mass','ek')
  plot.save.i(g,'RRx.base')
  Y$RRx = as.factor(Y$RRx)
  g = ggplot(subset(Y,RRx!=.333),aes(x=RRo,y=value,lty=ek,color=RRx,fill=RRx)) +
    plot.core('RRo','OR','RRx','ek')
  plot.save.i(g,'RR2.base')
}

plot.obj.2 = function(){
  Y = rbind(load.grid('RRo.ad.base',i=Tid$XRx,f=reps),
            load.grid('RRo.ch.base',i=Tid$XRx,f=reps))
  g = ggplot(Y,aes(x=RRo,y=value,lty=ek,color=mass,fill=mass)) +
    fct_grid('erep','orep') + sublabs(Y[reps]) +
    plot.core('RRo','mass','mass','ek')
  plot.save.i(g,'RRo.reps')
}

plot.obj.3 = function(){
  for (k in c('ch.ep','ch.oRo','ch.oRx',
    'ad.eRo','ad.eRx','ad.oRo','ad.oRx')){
    R = substr(k,4,6);   iR = str('as.factor(100*',R,')')
    H = gsub('R','H',R); iH = str('interaction(mass,',H,')')
    if (R=='ep'){ H = NULL; iH = 'mass' }
    Y = subset(load.grid(str('RRo.',k),i=Tid$XRx,f=c(reps,H)),RRo==8)
    g = ggplot(Y,aes.string(x=iR,y='bias.adj',lty='mass',color=H,fill=H,group=iH)) +
      fct_grid('erep','orep') + sublabs(Y[reps]) + ylab('Bias vs onset HR') +
      plot.core(R,'bias',H,'mass',da=0)
    plot.save.i(g,str('bias.',k))
  }
}

plot.ch.valid = function(){
  # TODO: something funky @ or.wpa & oHo=0
  Y = load.grid('RRo.ch.oRo',i=Tid$XRw)
  g = ggplot(Y,aes(x=RRo,y=value,color=as.factor(100*oRo),fill=as.factor(100*oRo))) +
    fct_grid('oHo','mass') + plot.core('RRo','mass','oRo')
  plot.save.i(g + add.ch.exact(Y),'valid.ch.oRo')
  Y = load.grid('RRo.ch.oRx',i=Tid$XRw)
  g = ggplot(Y,aes(x=RRo,y=value,color=as.factor(100*oRx),fill=as.factor(100*oRx))) +
    fct_grid('oHx','mass') + plot.core('RRo','mass','oRx')
  plot.save.i(g + add.ch.exact(Y),'valid.ch.oRx')
}

plot.ch.age.i = function(slug,orep='current',fac='RRo',clr=NULL,...,mm=c(1,16)){
  ids = c(p1='exposed',p0='unexposed',OR='OR',PR='PR')
  GR = Gi(c(clr,fac),ek='child',seed=1,orep=orep,...)
  Y = run.ch.exact(expand.grid(GR),age=1)
  Y = melt(Y,m=names(ids),var='id')
  Y$type = ifelse(Y$id %in% fl$mass,ll('mass'),ll('op'))
  Y$id   = factor(Y$id,names(ids),ids)
  Y$oRx[Y$oRx<=z] = 0 # HACK
  Y[c('oRo','oRx')] = 100*Y[c('oRo','oRx')]
  Y$value = Y$value * ifelse(Y$type==ll('op'),100,1)
  Y[[clr]] = as.factor(Y[[clr]])
  g = ggplot(Y,aes.string(x='a+amin',y='value',lty='id',color=clr)) +
    facet_grid(str('type~',fac),lab=labeller(.cols=fct(labels[[fac]])),scales='free_y') +
    scale_linetype_manual(values=c(exposed='31',unexposed='13',OR='solid',PR='22')) +
    ggh4x::scale_y_facet(type==ll('op'),lim=c(0,100)) +
    ggh4x::scale_y_facet(type==ll('mass'),lim=mm,trans='log2') +
    geom_hline(data=subset(Y,id=='OR'),aes(yintercept=RRo/RRx),lty='11',color='#999') +
    labs(x=ll('age'),y='Value',lty='Variable',color=ll(clr,1)) +
    cmap[[if.null(clr,'null')]] + plot.clean(font='Alegreya Sans') +
    geom_line()
  plot.save.i(g,str('age.',slug))
}

plot.ch.age = function(){
  v = list(RR=4,oRo=c(.003,.01,.03),oRx=c(z,.1,1),oHo=c(0,.5,1,1.5))
  labels$oRx <<- gsub('Depression','Depr',labels$oRx) # HACK
  plot.ch.age.i('RRo.R2.hom',RRo=v$RR,  clr='oRo',fac='oRx',oRo=v$oRo,oRx=v$oRx,oHo=0,    oHx=0)
  plot.ch.age.i('RRo.R2.het',RRo=v$RR,  clr='oRo',fac='oRx',oRo=v$oRo,oRx=v$oRx,oHo=1,    oHx=1)
  plot.ch.age.i('RRx.R2.hom',RRx=1/v$RR,clr='oRo',fac='oRx',oRo=v$oRo,oRx=v$oRx,oHo=0,    oHx=0)
  plot.ch.age.i('RRx.R2.het',RRx=1/v$RR,clr='oRo',fac='oRx',oRo=v$oRo,oRx=v$oRx,oHo=1,    oHx=1)
  plot.ch.age.i('RRo.oHo',   RRo=v$RR,  clr='oHo',fac='oRx',oRo=.03,  oRx=v$oRx,oHo=v$oHo,oHx=1)
}

# -----------------------------------------------------------------------------
# main

# run.grid(.k)
# merge.batch(.k)

# plot.obj.1()
# plot.obj.2()
# plot.obj.3()
# plot.ch.valid()
# plot.ch.age()
