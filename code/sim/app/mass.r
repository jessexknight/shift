source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')
uid = '2026-02-26'
.k     = cli.arg('.k','RRo.rev.base')
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

G = name.list(key='i',
  list(i='RRo', id='RR.haz_o.dep_w',v=1,   vg=  c(1,2,3,4,8)),
  list(i='RRx', id='RR.haz_x.dep_w',v=1,   vg=1/c(1,2,3,4,8)),
  list(i='eRo', id='dep_o.Ri.my',   v=.03, vg=c(.003,.01,.03,.1)),
  list(i='eRx', id='dep_x.Ri.my',   v=.3,  vg=c(.03,.1,.3,1)),
  list(i='eHo', id='dep_o.Ri.het',  v=0,   vg=c(0,.3,1,3)),
  list(i='eHx', id='dep_x.Ri.het',  v=0,   vg=c(0,.3,1,3)),
  list(i='ecv', id='dep.cov',       v=0,   vg=c(-.5,0,+.5)),
  list(i='ep',  id='dep.prev',      v=.1,  vg=c(.03,.1,.2,.3)),
  list(i='oRo', id='haz_o.Ri.my',   v=.03, vg=c(.003,.01,.03,.1)),
  list(i='oRx', id='haz_x.Ri.my',   v=.3,  vg=c(.03,.1,.3,1)),
  list(i='oHo', id='haz_o.Ri.het',  v=0,   vg=c(0,.3,1,3)),
  list(i='oHx', id='haz_x.Ri.het',  v=0,   vg=c(0,.3,1,3)),
  list(i='ocv', id='haz.cov',       v=0,   vg=c(-.5,0,+.5)),
  list(i='seed',id='seed',          v=NA,  vg=xdf(1:7,1:21)),
  list(i='ek',  id='e.case', v='rev',vg=NA),
  list(i='ok',  id='o.case', v='rev',vg=NA))

Gid = lapply(G,`[[`,'id')
G0 = lapply(G,`[[`,'v')
Gi = function(i,...){ ulist(G0,lapply(G[c('seed',i)],`[[`,'vg'),...) }
PG = function(Gk,...){ ulist(P0,set.names(Gk,Gid[names(Gk)]),...) }

Gk = list()
# fixed exposure
Gk$RRo.fix.base = Gi(ek='fix',c('RRo'))
Gk$RRo.fix.ep   = Gi(ek='fix',c('RRo','ep'))
Gk$RRo.fix.oRo  = Gi(ek='fix',c('RRo','oRo','oHo'))
Gk$RRo.fix.oRx  = Gi(ek='fix',c('RRo','oRx','oHx'))
Gk$RRo.fix.oR2  = Gi(ek='fix',c('RRo','oRo','oRx','ocv'),oHo=1,oHx=1)
# irreversible exposure
Gk$RRo.irr.base = Gi(ek='irr',c('RRo'))
Gk$RRo.irr.eRo  = Gi(ek='irr',c('RRo','eRo','eHo'))
Gk$RRo.irr.oRo  = Gi(ek='irr',c('RRo','oRo','oHo'))
Gk$RRo.irr.oRx  = Gi(ek='irr',c('RRo','oRx','oHx'))
Gk$RRo.irr.oR2  = Gi(ek='irr',c('RRo','oRo','oRx','ocv'),oHo=1,oHx=1)
# reversible exposure
Gk$RRo.rev.base = Gi(ek='rev',c('RRo'))
Gk$RRx.rev.base = Gi(ek='rev',c('RRx'))
Gk$RR2.rev.base = Gi(ek='rev',c('RRo','RRx'))
Gk$RRo.rev.eRo  = Gi(ek='rev',c('RRo','eRo','eHo'))
Gk$RRo.rev.eRx  = Gi(ek='rev',c('RRo','eRx','eHx'))
Gk$RRo.rev.eR2  = Gi(ek='rev',c('RRo','eRo','eRx','ecv'),eHo=1,eHx=1)
Gk$RRo.rev.oRo  = Gi(ek='rev',c('RRo','oRo','oHo'))
Gk$RRo.rev.oRx  = Gi(ek='rev',c('RRo','oRx','oHx'))
Gk$RRo.rev.oR2  = Gi(ek='rev',c('RRo','oRo','oRx','ocv'),oHo=1,oHx=1)
Gk$RRo.rev.2Rx  = Gi(ek='rev',c('RRo','eRx','eHx','oRx','oHx'))
Gk$RRo.rev.lhs  = Gi(ek='rev',seed=c(1,1e9),lhs=xdf(1e1,1e4),
  c('RRo','eRo','eRx','eHo','eHx','ecv','oRo','oRx','oHo','oHx','ocv'))
# for (k in names(Gk)){ status(3,k,': ',prod(lens(Gk[[k]]))) } # for hpc gen

apply.case = function(P,eps=1e-12){
  if (P$e.case=='fix'){
    P$init.inds = function(I,P){
      I$dep_o.Ri = ifelse(runif(P$n.tot) < P$dep.prev,Inf,0)
      return(I) }}
  if (P$e.case!='rev'){ P$dep_x.Ri.m = eps }
  if (P$o.case!='rev'){ P$haz_x.Ri.m = eps }
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
  Gi = { if (lhs) get.lhs(Gk[[k]]) else Gk[[k]] }
  Y = grid.apply(Gi,run.one,.grid=!lhs,.batch=.b,.nbatch=.nb,
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
  v = c('ep','eRo','eRx','oRo','oRx')
  Y[v] = Y[v] * 100
  Y$bias     = ifelse(Y$type=='prop',NA,Y$value/(Y$RRo/Y$RRx))
  Y$bias.adj = ifelse(Y$type=='prop',NA,(Y$value-1)/(Y$RRo/Y$RRx-1))
  Y$mass = factor(substr(Y$id,1,2),names(fl$mass),fl$mass)
  Y$erep = factor(substr(Y$id,4,4),names(fl$report),fl$report)
  Y$orep = factor(substr(Y$id,5,5),names(fl$report),fl$report)
  Y$ek  = factor(Y$ek,names(fl$case),fl$case)
  Y$RRx = round(Y$RRx,3)
  Y[f] = lapply(Y[f],as.factor)
  return(Y)
}

fl = list( # factor levels
  mass = c(or='OR',pr='PR'),
  report = c(w='current',p='lifetime'),
  case = c(fix='fixed',irr='irreversible',rev='reversible'))

reps = c('erep','orep')

# -----------------------------------------------------------------------------
# plot utils

labels = list(
  mass = 'Measure of~association',
  bias = 'Bias~vs~HR',
  OR   = 'OR:~abuse and~depression',
  PR   = 'PR:~abuse and~depression',
  RRo  = 'HR:~depression~onset~while~abused',
  RRx  = 'HR:~depression~recovery~while~abused',
  iRRx = '1/HR:~depression recovery~while abused',
  ep   = 'Abuse~prevalence',
  op   = 'Depression~prevalence',
  eRo  = 'Abuse~onset rate~(per 100 PY)',
  eRx  = 'Abuse~recovery rate~(per 100 PY)',
  oRo  = 'Depression~onset rate~(per 100 PY)',
  oRx  = 'Depression~recovery rate~(per 100 PY)',
  eHo  = 'Abuse~onset~frailty SD',
  eHx  = 'Abuse~recovery~frailty SD',
  oHo  = 'Depression~onset~frailty SD',
  oHx  = 'Depression~recovery~frailty SD',
  erep = 'Abuse~reporting',
  orep = 'Depression~reporting')

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

sublabs = def.args(add.sublabs,dx=.5,size=3,family='Alegreya Sans',
  labs=tolower(as.roman(1:99)))

cmap = lapply(list(RRo='viridis',RRx='inferno',ep='mako',
  eRo='mako',  eHo='mako',  eRx='mako',  eHx='mako',
  oRo='rocket',oHo='rocket',oRx='rocket',oHx='rocket'),
  function(o){ clr.map.d(option=o,end=.7) })
cmap$mass = clr.map.m(c('#c06','#0cc'))

ltys = lapply(list(
    v2=c('solid','22'),
    v3=c('solid','31','11'),
    v4=c('solid','41','21','11')),
  function(v){ scale_linetype_manual(values=v) })

scales = list(
  mass = scale_y_continuous(breaks=seq(0,10,2 ),limits=c(0,10)),
  bias = scale_y_continuous(breaks=seq(0, 2,.5),limits=c(0, 2)),
  RRo  = scale_x_continuous(breaks=seq(0, 8,2 ),limits=c(0, 8)))
scales$OR = scales$PR = scales$mass

plot.core = function(x,y,clr=NULL,lty=NULL,da=1,ra=1/5,ci=.95){ list(
  scales[[x]],scales[[y]],cmap[[clr]],
  geom_hline(lty='11',color='#999',yintercept=1),
  geom_abline(lty='11',color='#999',alpha=da),
  labs(x=ll(x),y=ll(y),lty=ll(lty,1),color=ll(clr,1),fill=ll(clr,1)),
  stat_summary(geom='ribbon',color=NA,alpha=ra,
    fun.min=qfun((1-ci)/2),fun.max=qfun(1-(1-ci)/2)),
  stat_summary(geom='line',fun=mean),
  plot.clean(font='Alegreya Sans')
)}

add.stats.ci = function(){ list(
  stat_summary(geom='line',aes(y=lower),lty='22',lwd=1/4,fun=mean),
  stat_summary(geom='line',aes(y=upper),lty='22',lwd=1/4,fun=mean)
)}

plot.1o = list(w1=2,h1=1.6,wo=1.5,ho=1)

plot.save.i = function(g,...,size=NULL,ext='.png'){
  plot.save(g,'mass',uid,...,ext=ext,size=size)
}

# -----------------------------------------------------------------------------
# objective plots

plot.obj.1 = function(){
  Y = load.grid('RR2.rev.base',i=Tid$XRw)
  g = ggplot(subset(Y,RRx==1),aes(x=RRo,y=value,color=mass,fill=mass)) +
    plot.core('RRo','mass','mass')
  plot.save.i(g,'RRo.base')
  g = ggplot(subset(Y,RRo==1),aes(x=1/RRx,y=value,color=mass,fill=mass)) +
    plot.core('iRRx','mass','mass')
  plot.save.i(g,'RRx.base')
  Y$RRx = as.factor(Y$RRx)
  g = ggplot(subset(Y,RRx!=.333),aes(x=RRo,y=value,color=RRx,fill=RRx)) +
    plot.core('RRo','OR','RRx')
  plot.save.i(g,'RR2.base')
}

plot.obj.2 = function(){
  Y = load.grid('RRo.rev.base',i=Tid$XRx,f=reps)
  g = ggplot(Y,aes(x=RRo,y=value,color=mass,fill=mass)) +
    fct_grid('erep','orep') + sublabs(Y[reps]) +
    plot.core('RRo','mass','mass')
  plot.save.i(g,'RRo.reps')
  plot.save.i(g + add.stats.ci(),'RRo.reps.ci')
}

plot.obj.3 = function(){
  for (R in c('eRo','eRx','oRo','oRx')){ H = gsub('R','H',R)
    iH = str('interaction(mass,',H,')')
    Y = subset(load.grid(str('RRo.rev.',R),i=Tid$XRx,f=c(reps,R,H)),RRo==8)
    g = ggplot(Y,aes.string(x=R,y='bias.adj',lty='mass',color=H,fill=H,group=iH)) +
      fct_grid('erep','orep') + sublabs(Y[reps]) + ylab('Bias vs onset HR') +
      plot.core(R,'bias',H,'mass',da=0)
    plot.save.i(g,str('RRo.bias.',R))
  }
}

# -----------------------------------------------------------------------------
# main

# run.grid(.k)
# merge.batch(.k)

# plot.obj.1()
# plot.obj.2()
# plot.obj.3()
