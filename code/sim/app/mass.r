source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')
uid = '2026-04-04'
.k     = cli.arg('.k','RR2.pet.base')
.b     = cli.arg('.b', 1)
.nb    = cli.arg('.nb',1)
.debug = cli.arg('.debug',0)

# e = exposure "dep" {fixed,period} or "vio" {event}
# o = outcome  "haz"

# -----------------------------------------------------------------------------
# params & grid

P0 = list(
  pop.type = 'open',
  n.pop = xdf(1000,10000),
  n.dur = 1, dtz = 7,
  het.distr = 'gamma',
  tRR.shape = 'step')

z = 1e-12
G = name.list(key='i',
  list(i='exp', id='exp.type',  def='period',grid=NA), # {fixed,period,event}
  list(i='eff', id='eff.type',  def='trans', grid=NA), # {first,trans}
  list(i='RRo', id='RR.haz_o.dep_w',def=1,   grid=  c(1,2,3,4,8)),
  list(i='RRx', id='RR.haz_x.dep_w',def=1,   grid=1/c(1,2,3,4,8)),
  list(i='eRo', id='dep_o.Ri.my',   def=.01, grid=c(.003,.01,.03,.1)),
  list(i='eRx', id='dep_x.Ri.my',   def=1,   grid=c(z,.03,.1,.3,1,3)),
  list(i='eHo', id='dep_o.Ri.het',  def=2,   grid=c(0,1,2,3)),
  list(i='eHx', id='dep_x.Ri.het',  def=1,   grid=c(0,1,2,3)),
  list(i='ecv', id='dep.cov',       def=-.5, grid=c(-.5,0,+.5)),
  list(i='ep',  id='dep.prev',      def=.1,  grid=c(.03,.1,.2,.3)),
  list(i='oRo', id='haz_o.Ri.my',   def=.01, grid=c(.003,.01,.03,.1)),
  list(i='oRx', id='haz_x.Ri.my',   def=3,   grid=c(z,.03,.1,.3,1,3)),
  list(i='oHo', id='haz_o.Ri.het',  def=2,   grid=c(0,1,2,3)),
  list(i='oHx', id='haz_x.Ri.het',  def=1,   grid=c(0,1,2,3)),
  list(i='ocv', id='haz.cov',       def=-.5, grid=c(-.5,0,+.5)),
  list(i='seed',id='seed',          def=NA,  grid=xdf(1:7,1:21)))

Gid = lapply(G,`[[`,'id')
G0 = lapply(G,`[[`,'def')
Gi = function(i,...){ ulist(G0,lapply(G[c('seed',i)],`[[`,'grid'),...) }
PG = function(G1,...){ ulist(P0,set.names(G1,Gid[names(G1)]),...) }

Gk = list()
# exposure: fixed, effect: {first / transient}
Gk$RR2.fix.base = Gi(exp='fixed', eff='trans',c('RRo','RRx'))
Gk$RRo.fix.base = Gi(exp='fixed', eff='trans',c('RRo'))
Gk$RRx.fix.base = Gi(exp='fixed', eff='trans',c('RRx'))
Gk$RRo.fix.ep   = Gi(exp='fixed', eff='trans',c('RRo','ep'))
Gk$RRo.fix.oRo  = Gi(exp='fixed', eff='trans',c('RRo','oRo','oHo'))
Gk$RRo.fix.oRx  = Gi(exp='fixed', eff='trans',c('RRo','oRx','oHx'))
# exposure: period, effect: transient
Gk$RR2.pet.base = Gi(exp='period',eff='trans',c('RRo','RRx'))
Gk$RRo.pet.base = Gi(exp='period',eff='trans',c('RRo'))
Gk$RRx.pet.base = Gi(exp='period',eff='trans',c('RRx'))
Gk$RRo.pet.eRo  = Gi(exp='period',eff='trans',c('RRo','eRo','eHo'))
Gk$RRo.pet.eRx  = Gi(exp='period',eff='trans',c('RRo','eRx','eHx'))
Gk$RRo.pet.oRo  = Gi(exp='period',eff='trans',c('RRo','oRo','oHo'))
Gk$RRo.pet.oRx  = Gi(exp='period',eff='trans',c('RRo','oRx','oHx'))
# exposure: period, effect: first
Gk$RR2.pef.base = Gi(exp='period',eff='first',c('RRo','RRx'))
Gk$RRo.pef.base = Gi(exp='period',eff='first',c('RRo'))
Gk$RRx.pef.base = Gi(exp='period',eff='first',c('RRx'))
Gk$RRo.pef.eRo  = Gi(exp='period',eff='first',c('RRo','eRo','eHo'))
Gk$RRo.pef.eRx  = Gi(exp='period',eff='first',c('RRo','eRx','eHx'))
Gk$RRo.pef.oRo  = Gi(exp='period',eff='first',c('RRo','oRo','oHo'))
Gk$RRo.pef.oRx  = Gi(exp='period',eff='first',c('RRo','oRx','oHx'))
# exposure: event, effect: transient
Gk$RR2.evt.base = Gi(exp='event', eff='trans',c('RRo','RRx'))
Gk$RRo.evt.base = Gi(exp='event', eff='trans',c('RRo'))
Gk$RRx.evt.base = Gi(exp='event', eff='trans',c('RRx'))
Gk$RRo.evt.eRo  = Gi(exp='event', eff='trans',c('RRo','eRo','eHo'))
Gk$RRo.evt.eRx  = Gi(exp='event', eff='trans',c('RRo','eRx','eHx'))
Gk$RRo.evt.oRo  = Gi(exp='event', eff='trans',c('RRo','oRo','oHo'))
Gk$RRo.evt.oRx  = Gi(exp='event', eff='trans',c('RRo','oRx','oHx'))
# for (k in names(Gk)){ status(3,k,': ',prod(lens(Gk[[k]]))) } # for hpc gen

apply.types = function(P){
  P$run = get.run.par(c('dep','haz'),u=0)
  P$mass.var = switch(P$eff.type,first='dep.past',trans='dep.now')
  if (P$exp.type=='period'){}
  if (P$exp.type=='fixed'){
    P$init.inds = function(I,P){
      I$dep_o.Ri = ifelse(runif(P$n.tot) < P$dep.prev,Inf,0)
      I$dep_x.Ri = 0
      return(I) }}
  if (P$exp.type=='event'){
    P$run = get.run.par(c('vio','haz'),u=0)
    tsc = P$t1y*switch(P$eff.type,first=adur,trans=min(adur,1/P$dep_x.Ri.my))
    P$vio.Ri.my        = P$dep_o.Ri.my
    P$vio.Ri.het       = P$dep_o.Ri.het
    P$iRR.haz_o.vio_zr = P$RR.haz_o.dep_w
    P$tsc.haz_o.vio_zr = tsc
    P$iRR.haz_x.vio_zr = P$RR.haz_x.dep_w
    P$tsc.haz_x.vio_zr = tsc }
  return(P)
}

# -----------------------------------------------------------------------------
# customize rate funs

rate.vio = function(P,J,aj,z){
  R = aggr.rate( # among all
      J$vio.Ri # base rate
    # skip: RR age, tRR vio, nRR vio
); return(R) }

rate.dep_o = function(P,J,R,aj,z){
  j = which(!J$dep.now)
  R[j] = aggr.rate( # among not dep
      J$dep_o.Ri[j] # base rate
    , map.tRR(P$tRRu.dep_o.vio_zr,J$vio.zr[j],z) # tRR vio
    # skip: RR age, dep past, nRR vio
); return(R) }

rate.dep_x = function(P,J,R,aj,z){
  j = which(J$dep.now)
  R[j] = aggr.rate( # among dep
      J$dep_x.Ri[j] # base rate
    , map.tRR(P$tRRu.dep_x.vio_zr,J$vio.zr[j],z) # tRR vio
    # skip: RR dep dur
); return(R) }

rate.haz_o = function(P,J,R,aj,z){
  j = which(!J$haz.now)
  R[j] = aggr.rate( # among not haz
      J$haz_o.Ri[j] # base rate
    , (1 + P$RRu.haz_o.dep_w * J[[P$mass.var]][j]) # RR dep now
    , map.tRR(P$tRRu.haz_o.vio_zr,J$vio.zr[j],z) # tRR vio
    # skip: RR age, RR haz past, nRR vio
); return(R) }

rate.haz_x = function(P,J,R,aj,z){
  j = which(J$haz.now)
  R[j] = aggr.rate( # among haz
      J$haz_x.Ri[j] # base rate
    , (1 + P$RRu.haz_x.dep_w * J[[P$mass.var]][j]) # RR dep now
    , map.tRR(P$tRRu.haz_x.vio_zr,J$vio.zr[j],z) # tRR vio
    # skip: RR haz dur
); return(R) }

# -----------------------------------------------------------------------------
# targets / outcomes

T = name.list(key='id',
  gen.targ(id='px:ew',    type='prop',vo='dep.now' ),
  gen.targ(id='px:ep',    type='prop',vo='dep.past'),
  gen.targ(id='px:ow',    type='prop',vo='haz.now' ),
  gen.targ(id='px:op',    type='prop',vo='haz.past'),
  gen.targ(id='or:ew.ow', type='OR',  ve='dep.now', vo='haz.now' ),
  gen.targ(id='or:ew.op', type='OR',  ve='dep.now', vo='haz.past'),
  gen.targ(id='or:ep.ow', type='OR',  ve='dep.past',vo='haz.now' ),
  gen.targ(id='or:ep.op', type='OR',  ve='dep.past',vo='haz.past'),
  gen.targ(id='pr:ew.ow', type='PR',  ve='dep.now', vo='haz.now' ),
  gen.targ(id='pr:ew.op', type='PR',  ve='dep.now', vo='haz.past'),
  gen.targ(id='pr:ep.ow', type='PR',  ve='dep.past',vo='haz.now' ),
  gen.targ(id='pr:ep.op', type='PR',  ve='dep.past',vo='haz.past'),
  gen.targ(id='aor:ew.ow',type='OR',  ve='dep.now', vo='haz.now' ,va1='age'),
  gen.targ(id='aor:ew.op',type='OR',  ve='dep.now', vo='haz.past',va1='age'),
  gen.targ(id='aor:ep.ow',type='OR',  ve='dep.past',vo='haz.now' ,va1='age'),
  gen.targ(id='aor:ep.op',type='OR',  ve='dep.past',vo='haz.past',va1='age'),
  gen.targ(id='apr:ew.ow',type='PR',  ve='dep.now', vo='haz.now' ,va1='age'),
  gen.targ(id='apr:ew.op',type='PR',  ve='dep.now', vo='haz.past',va1='age'),
  gen.targ(id='apr:ep.ow',type='PR',  ve='dep.past',vo='haz.now' ,va1='age'),
  gen.targ(id='apr:ep.op',type='PR',  ve='dep.past',vo='haz.past',va1='age'))

Tid = list()
Tid$aOR.xx = filter.names(T,'aor:e..o.')
Tid$aPR.xx = filter.names(T,'apr:e..o.')
Tid$aXR.xx = filter.names(T,'a.r:e..o.')
Tid$aXR.ww = filter.names(T,'a.r:ew.ow')

srv.event = function(P,Q,E,t){
  Q$vio.tf = sapply(E$vio,any.dt,t=t,dt=P$tsc.haz_o.vio_zr)
  Q$vio.1y = sapply(E$vio,any.dt,t=t,dt=P$t1y)
  Q$vio.6m = sapply(E$vio,any.dt,t=t,dt=P$t1y/2)
  Q$vio.3m = sapply(E$vio,any.dt,t=t,dt=P$t1y/4)
  Q$vio.1m = sapply(E$vio,any.dt,t=t,dt=P$t1y/12)
  return(Q)
}

targs.event = function(T){
  recs = c('tf','1y','6m','3m','1m')
  T = unlist(rec=0,lapply(T,function(Ti){
    Ti$arg = lapply(Ti$arg,sub,pat='dep',rep='vio')
    if (grepl('ew',Ti$id)){
      Ts = lapply(recs,function(rec){
        Ti$arg = lapply(Ti$arg,sub,pat='vio.now',rep=str('vio.',rec))
        Ti$id  = sub('ew',str('e',rec),Ti$id)
        return(Ti)
      })}
    else { Ts = list(Ti) }
  }))
  T = set.names(T,lapply(T,`[[`,'id'))
}

# -----------------------------------------------------------------------------
# run sims & save/load

grid.path = function(k,.save=FALSE){
  path = hash.path(PG(Gk[[k]]),'data','sim','mass',uid,k,.save=.save)
}

run.one = function(...,.par=0){
  P1 = PG(list(...),fun.pre=apply.types)
  Ps = get.pars.grid(P1,.par=.par)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  Q  = srv.apply(Ms,srvs=srv.event,.par=.par)
  if (P1$exp.type=='event') { T = targs.event(T) }
  Y  = srv.targs(Q,T)
  Y[c('seed','targ.mu','targ.se','ll')] = NULL
  row.names(Y) = NULL
  return(Y)
}

run.grid = function(k){
  Y = grid.apply(Gk[[k]],run.one,.grid=1,.batch=.b,.nbatch=.nb,
    .rbind=1,.cbind=1,.log=3)
  save.rds(Y,grid.path(k,.save=TRUE),str('b',.nb),str('Y.',.b))
}

merge.batch = function(k){
  Y = rbind.lapply(1:.nb,function(b){
    Yb = load.rds(grid.path(k),str('b',.nb),str('Y.',b)) })
  save.rds(Y,grid.path(k),'Y')
}

load.grid = function(k,i='aor:ew.ow',f=NULL){
  Y = load.rds(grid.path(k),'Y')
  Y = subset(Y,id %in% i)
  Y[c('ve','vo','te','to','dt')] = NULL
  v = c('eRo','eHo','eRx','eHx','oRo','oHo','oRx','oHx')
  Y[v][Y[v]<=z] = 0 # HACK
  Y$bias     = ifelse(Y$type=='prop',NA,Y$value/(Y$RRo/Y$RRx))
  Y$bias.adj = ifelse(Y$type=='prop',NA,(Y$value-1)/(Y$RRo/Y$RRx-1))
  Y$mass = factor(Y$type,names(fl$mass),fl$mass)
  Y$erep = factor(gsub('.*\\:e|\\.o.','',Y$id),names(fl$rep),fl$rep)
  Y$orep = factor(gsub('.*\\:e.*\\.o','',Y$id),names(fl$rep),fl$rep)
  Y$exp = factor(Y$exp,names(fl$exp),fl$exp)
  Y$eff = factor(Y$eff,names(fl$eff),fl$eff)
  Y$RRx = round(Y$RRx,3)
  Y[f] = lapply(Y[f],as.factor)
  return(Y)
}

fl = list() # factor levels
fl$mass = c(OR='OR',PR='PR')
fl$exp  = c(fixed='fixed',period='period',event='event')
fl$eff  = c(trans='transient',first='first-ever')
fl$rep  = c(w='current',p='lifetime',tf='effect duration',
  '1y'='past year','6m'='past 6 months','3m'='past 3 months','1m'='past month')

reps = c('erep','orep')

# -----------------------------------------------------------------------------
# exact math @ exp.type=fixed

avec = seq(.1,adur,.1)

efun = list(
  prev.age = function(o,x){ k=o+x; p = colMeans((1-exp(-outer(k,avec)))*o/k) },
  prev = function(o,x){ k=o+x; p = mean((adur-(1-exp(-k*adur))/k)*o/adur/k) },
  OR = function(p0,p1){ p1*(1-p0)/p0/(1-p1) },
  PR = function(p0,p1){ p1/p0 })

run.fix.exact = function(Y,n=1e4,age=FALSE){
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
  exp  = 'Exposure~type',
  eff  = 'Effect~type',
  RRo  = 'HR θ:~outcome onset~while exposed',
  RRx  = 'HR φ:~outcome recovery~while exposed',
  iRRx = '1/HR 1/φ:~outcome recovery~while exposed',
  ep   = 'Fixed~exposure~prevalence~(%)',
  op   = 'Outcome~prevalence~(%)',
  eRo  = 'Mean~exposure~onset rate λ~(per 100 PY)',
  eRx  = 'Mean~exposure~recovery rate γ~(per 100 PY)',
  oRo  = 'Mean~outcome~onset rate μ~(per 100 PY)',
  oRx  = 'Mean~outcome~recovery rate η~(per 100 PY)',
  eHo  = 'Exposure~onset~frailty SD~σλ',
  eHx  = 'Exposure~recovery~frailty SD~σγ',
  oHo  = 'Outcome~onset~frailty SD~σμ',
  oHx  = 'Outcome~recovery~frailty SD~ση',
  erep = 'Exposure~reporting',
  orep = 'Outcome~reporting',
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

cmap = lapply(list(RRo='cividis',RRx='cividis',ep='viridis',
  eRo='viridis',eHo='viridis',eRx='mako',  eHx='mako',
  oRo='inferno',oHo='inferno',oRx='rocket',oHx='rocket'),
  function(o){ clr.map.d(option=o,end=.8) })
cmap$exp  = clr.map.m(c('#c06','#0cc','#f90'))
cmap$null = clr.map.m('#000')

ltys = lapply(list(
    v2=c('solid','22'),
    v3=c('solid','31','11'),
    v4=c('solid','41','21','11')),
  function(v){ scale_linetype_manual(values=v) })

scales = list(
  mass = scale.y.cts(breaks=seq(0,8,2),limits=c(0,8)),
  RRo  = scale.x.cts(breaks=seq(0,8,2),limits=c(0,8)),
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

add.fix.exact = function(Y){
  geom_point(data=run.fix.exact(Y),shape=21,fill='#fc0',size=1)
}

plot.1o = list(w1=2,h1=1.6,wo=1.5,ho=1)

plot.save.i = function(g,...,size=NULL,ext='.png'){
  plot.save(g,'mass',uid,...,ext=ext,size=size)
}

# -----------------------------------------------------------------------------
# objective plots

plot.obj.1 = function(){
  Y = load.grid('RR2.pet.base',i=Tid$aXR.ww)
  g = ggplot(subset(Y,RRx==1),aes(x=RRo,y=value,lty=mass)) +
    plot.core('RRo','mass',lty='mass')
  plot.save.i(g,'RRo.pet.base')
  g = ggplot(subset(Y,RRo==1),aes(x=1/RRx,y=value,lty=mass)) +
    plot.core('iRRx','mass',lty='mass')
  plot.save.i(g,'RRx.pet.base')
  Y$RRx = as.factor(Y$RRx)
  g = ggplot(subset(Y,RRx!=.333),aes(x=RRo,y=value,lty=mass,color=RRx,fill=RRx)) +
    plot.core('RRo','mass','RRx','mass') + labs(lty='Measure')
  plot.save.i(g,'RR2.pet.base')
  Y = load.grid('RR2.fix.base',i=Tid$aXR.ww,f='RRx')
  g = ggplot(subset(Y,RRx!=.333),aes(x=RRo,y=value,lty=mass,color=RRx,fill=RRx)) +
    plot.core('RRo','mass','RRx','mass') + labs(lty='Measure')
  plot.save.i(g,'RR2.fix.base')
}

plot.obj.2 = function(){
  Y.pet = load.grid('RRo.pet.base',i=Tid$aXR.xx,f=reps)
  Y.fix = load.grid('RRo.fix.base',i=Tid$aXR.xx,f=reps)
  Y.evt = load.grid('RRo.evt.base',i=gsub('ew','etf',Tid$aXR.xx),f=reps)
  Y.all = rbind(Y.pet,Y.fix,Y.evt)
  g = ggplot(Y.pet,aes(x=RRo,y=value,lty=mass)) +
    fct_grid('erep','orep') +
    plot.core('RRo','mass','exp','mass')
  plot.save.i(g + sublabs(Y.pet[reps]),'RRo.pet.reps')
  g = g + Y.all + aes(color=exp,fill=exp)
  plot.save.i(g + sublabs(Y.all[reps]),'RRo.exp.reps')
}

plot.obj.3 = function(){
  for (k in c('fix.ep','fix.oRo','fix.oRx',
    'evt.eRo','evt.eRx','evt.oRo','evt.oRx',
    'pet.eRo','pet.eRx','pet.oRo','pet.oRx')){
    exp = substr(k,1,3); ids = Tid$aXR.xx
    R = substr(k,5,7);   iR = str('as.factor(100*',R,')')
    H = gsub('R','H',R); iH = str('interaction(mass,',H,')')
    if (R=='ep'){ H = NULL; iH = 'mass' }
    if (exp=='evt'){ ids = gsub('ew','etf',ids) }
    Y = subset(load.grid(str('RRo.',k),i=ids,f=c(reps,H)),RRo==8)
    if (Y$exp[1]=='fixed'){ Y = subset(Y,erep=='lifetime') }
    g = ggplot(Y,aes.string(x=iR,y='bias.adj',lty='mass',color=H,fill=H,group=iH)) +
      fct_grid('erep','orep') + sublabs(Y[reps]) + ylab('Bias vs onset HR') +
      plot.core(R,'bias',H,'mass',da=0)
    plot.save.i(g,str('bias.',k))
  }
}

plot.fix.valid = function(){
  # TODO: something funky @ aor.wp & oHo=0
  Y = load.grid('RRo.fix.oRo',i=Tid$aXR.ww)
  g = ggplot(Y,aes(x=RRo,y=value,color=as.factor(100*oRo),fill=as.factor(100*oRo))) +
    fct_grid('oHo','mass') + plot.core('RRo','mass','oRo')
  plot.save.i(g + add.fix.exact(Y),'valid.fix.oRo')
  Y = load.grid('RRo.fix.oRx',i=Tid$aXR.ww)
  g = ggplot(Y,aes(x=RRo,y=value,color=as.factor(100*oRx),fill=as.factor(100*oRx))) +
    fct_grid('oHx','mass') + plot.core('RRo','mass','oRx')
  plot.save.i(g + add.fix.exact(Y),'valid.fix.oRx')
}

plot.fix.age.i = function(slug,orep='current',fac='RRo',clr=NULL,...,mm=c(1,16)){
  ids = c(p1='exposed',p0='unexposed',OR='OR',PR='PR')
  GR = Gi(c(clr,fac),exp='fixed',seed=1,orep=orep,...)
  Y = run.fix.exact(expand.grid(GR),age=1)
  Y = melt(Y,m=names(ids),var='id')
  Y$type = ifelse(Y$id %in% fl$mass,ll('mass'),ll('op'))
  Y$id   = factor(Y$id,names(ids),ids)
  Y$oRx[Y$oRx<=z] = 0 # HACK
  Y[c('oRo','oRx')] = 100*Y[c('oRo','oRx')]
  Y$value = Y$value * ifelse(Y$type==ll('op'),100,1)
  Y[[clr]] = as.factor(Y[[clr]])
  g = ggplot(Y,aes.string(x='a+amin',y='value',lty='id',color=clr)) +
    facet_grid(str('type~',fac),scales='free_y',labeller=labeller(
      .cols=fct(labels[[fac]]),.rows=def.args(add.enum,fmt='i'))) +
    scale_linetype_manual(values=c(exposed='31',unexposed='13',OR='solid',PR='22')) +
    ggh4x::scale_y_facet(type==ll('op'),lim=c(0,100)) +
    ggh4x::scale_y_facet(type==ll('mass'),lim=mm,trans='log2') +
    geom_hline(data=subset(Y,id=='OR'),aes(yintercept=RRo/RRx),lty='11',color='#999') +
    labs(x=ll('age'),y='Value',lty='Variable',color=ll(clr,1)) +
    cmap[[if.null(clr,'null')]] + plot.clean(font='Alegreya Sans') +
    geom_line()
  plot.save.i(g,str('age.',slug))
}

plot.fix.age = function(){
  v = list(RR=4,oRo=c(.003,.01,.03),oRx=c(z,.1,1),oHo=c(0,.5,1,1.5))
  labels$oRx <<- gsub('Mean~outcome~recovery','Outcome~recov',labels$oRx) # HACK
  plot.fix.age.i('RRo.R2.hom',RRo=v$RR,  clr='oRo',fac='oRx',oRo=v$oRo,oRx=v$oRx,oHo=0,    oHx=0)
  plot.fix.age.i('RRo.R2.het',RRo=v$RR,  clr='oRo',fac='oRx',oRo=v$oRo,oRx=v$oRx,oHo=1,    oHx=1)
  plot.fix.age.i('RRx.R2.hom',RRx=1/v$RR,clr='oRo',fac='oRx',oRo=v$oRo,oRx=v$oRx,oHo=0,    oHx=0)
  plot.fix.age.i('RRx.R2.het',RRx=1/v$RR,clr='oRo',fac='oRx',oRo=v$oRo,oRx=v$oRx,oHo=1,    oHx=1)
  plot.fix.age.i('RRo.oHo',   RRo=v$RR,  clr='oHo',fac='oRx',oRo=.03,  oRx=v$oRx,oHo=v$oHo,oHx=1)
}

# -----------------------------------------------------------------------------
# main

# run.grid(.k)
# merge.batch(.k)

# plot.obj.1()
# plot.obj.2()
# plot.obj.3()
# plot.fix.valid()
# plot.fix.age()
