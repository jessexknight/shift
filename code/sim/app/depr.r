source('sim/meta.r')
source('sim/fit.r')
uid = '2025-11-14'
.b     = cli.arg('.b', 1)
.nb    = cli.arg('.nb',1)
.debug = cli.arg('.debug',1)

# -----------------------------------------------------------------------------
# params & grid

P0 = list(
  seed = 1:xdf(7,21),
  dtz = xdf(15,7),
  n.pop = xdf(3000,10000),
  n.dur = 1,
  dep_o.Ri.my = .01,
  dep_x.Ri.my = 1,
  dep_o.Ri.het = 0,
  dep_x.Ri.het = 0,
  dep.cov = 0,
  het.distr = 'gamma',
  run = get.run.par('dep',u=FALSE))

PG = list(
  seed = 1:xdf(7,21),
  dep_o.Ri.my  = xdf(c(.01,.03,.05),seq(.01,.05,.01)),
  dep_x.Ri.my  = xdf(c(1,2,3),      seq(.5,3,.5)),
  dep_o.Ri.het = xdf(c(0,2,4,6),    seq(0,6,1)),
  dep_x.Ri.het = xdf(c(0,1,2,3),    seq(0,3,.5)),
  dep.cov      = xdf(c(-.5,0,+.5),  seq(-.9,+.9,.3)))

t1y = add.pars.time(P0,P0$dtz)$t1y
for (v in names(PG)){ PG[[v]] = round(PG[[v]],2) }
grid.path = hash.path(ulist(P0,PG),'data','sim','depr',uid,.save=TRUE)
YP0 = function(Y,v){ if.null(Y[[v]],P0[[v]]) }

# subsets for plotting
PGk = list(
  base = P0[names(PG)],
  past = PG[c(1,2,4)],
  hom  = PG[1:3],
  hetu = PG[1:5],
  hetc = PG[1:6])
PGi = list(mo=3,mx=150,ho=3,hx=1.5,cov=0)
PGii = list(mo=c(2,4),mx=c(100,200),ho=c(0,4),hx=c(0,2),cov=c('–0.6','0','+0.6'))
PGiii = list(mo=c(1,3,5),mx=c(50,150,250),ho=c(0,2,4,6),hx=c(0,1,2,3))

# -----------------------------------------------------------------------------
# run sims & save/load

T = name.list(key='id',
  gen.targ(id='dep.now', type='prop',vo='dep.now'),
  gen.targ(id='dep.past',type='prop',vo='dep.past'))

est.rates = function(K,...,strat=NULL,e=c('dep_o','dep_x')){
  if (nrow(K)==0){ return(NULL) } # HACK
  R = cbind(rbind.lapply(e,rate.est,K=K,strat=strat,.par=FALSE),...)
}

run.one = function(...,.par=FALSE){
  Ps = get.pars.grid(ulist(P0,...),.par=.par)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  K  = rate.datas(Ms,e.dts=list(dep_o=t1y,dep_x=t1y),.par=.par)
  K$age.10 = int.cut(K$age.1,avec(10))
  Y = rbind.fill(
    est.rates(K=K), # onset & recov: overall
    est.rates(K=K,strat='age.10'), # onset & recov: by age
    est.rates(K=K,strat='dep.past',e='dep_o'), # onset: by dep.past
    est.rates(K=K,strat=c('age.10','dep.past'),e='dep_o'), # onset: by age x dep.past
    est.rates(K=subset(K,age.1==10|dep_x.dt.c==t1y),e='dep_o',sub='t1oa'), # 1-yr o/r
    est.rates(K=subset(K,age.1==10),      e='dep_o',sub='t1oi'), # 1-yr onset
    est.rates(K=subset(K,dep_x.dt.c==t1y),e='dep_o',sub='t1ol'), # 1-yr relap
    est.rates(K=subset(K,dep_o.dt.c==t1y),e='dep_x',sub='t1x'),  # 1-yr recov
    srv.targs(subset(K,e=='tmax'),T=T,strat='age.10'), # prev: by age
    srv.targs(subset(K,e=='tmax'),T=T)) # prev: overall
  Y[c('seed','targ.mu','targ.se','ll','est.mu','est.se')] = NULL
  row.names(Y) = NULL
  return(Y)
}

run.grid = function(){
  Y = grid.apply(PG,run.one,.rbind=1,.cbind=1,.batch=.b,.nbatch=.nb,.log=3)
  save.rds(Y,grid.path,str('b',.nb),str('Y.',.b))
}

merge.rda = function(){
  Y = rbind.lapply(1:.nb,function(b){
    Yb = load.rds(grid.path,str('b',.nb),str('Y.',b)) })
  par.lapply(unique(Y$id),function(i){
    Yi = subset(Y,id==i)
    for (k in names(PGk)){
      P0Gk = ulist(P0[names(PG)],PGk[[k]])
      Yik = merge(Yi,do.call(expand.grid,P0Gk))
      save.rds(Yik,grid.path,str('Y-',i,'-',k))
    }
  })
}

load.grid = function(k,i='dep.now',a10=FALSE,f=NULL){
  # TODO: a few id = NA?
  Y = load.rds(grid.path,str('Y-',i,'-',k))
  Y = subset(Y,is.na(age.10)!=a10)
  j = Y$type=='rate'
  v3 = c('value','lower','upper')
  Y[ j,v3] = Y[ j,v3]*100*t1y # rates per 100 PY
  Y[!j,v3] = Y[!j,v3]*100     # props as %
  Y$mo  = YP0(Y,'dep_o.Ri.my')*100 # shorthand, per 100 PY
  Y$mx  = YP0(Y,'dep_x.Ri.my')*100 # shorthand, per 100 PY
  Y$ho  = YP0(Y,'dep_o.Ri.het')    # shorthand
  Y$hx  = YP0(Y,'dep_x.Ri.het')    # shorthand
  Y$cov = YP0(Y,'dep.cov')         # shorthand
  Y$cov = factor(Y$cov,fl$cov,names(fl$cov))
  Y[f] = lapply(Y[f],as.factor) # Y[f] -> factors
  return(Y)
}

# -----------------------------------------------------------------------------
# exact model

exact.fun = list(
  dep.now  = function(o,x){ k=o+x; P = 100*(adur-(1-exp(-adur*k))/k)*o/adur/k },
  dep.past = function(o,x){ k=o;   P = 100*(adur-(1-exp(-adur*k))/k)*o/adur/k })

run.exact = function(Y,n=1e5){
  Y = subset(Y,seed==1)
  qf = het.funs[[P0$het.distr]]$q
  Ye = rbind.lapply(1:nrow(Y),function(i){ Yi = Y[i,]
    R = copula(n,covs=YP0(Yi,'dep.cov'),qfuns=list(o=qf,x=qf),
      o=list(m=YP0(Yi,'dep_o.Ri.my'),het=YP0(Yi,'dep_o.Ri.het')),
      x=list(m=YP0(Yi,'dep_x.Ri.my'),het=YP0(Yi,'dep_x.Ri.het')))
    R = round(R,12)+1e-12 # HACK: numerical stability
    Yi$value = mean(exact.fun[[Yi$id]](R[,1],R[,2]))
    return(Yi)
  })
}

# -----------------------------------------------------------------------------
# plot utils

# aes labels
l = list(
  mo  = 'Mean~onset rate~(per 100 PY)',
  mx  = 'Mean~recov rate~(per 100 PY)',
  ho  = 'Onset rate~frailty σ',
  hx  = 'Recov rate~frailty σ',
  cov = 'Rate~correlation',
  age = 'Age (years)',
  fup = 'Follow-up~time',
  dep.now  = 'Current~MDD~prevalence (%)',
  dep.past = 'Lifetime~MDD~prevalence (%)')

# aes label utils
hom = function(s){ gsub('Mean~(.)','\\U\\1',s,perl=TRUE) }
grp = function(s){ gsub('~','\n',s) }
axi = function(s){ gsub('~',' ',s) }
fct = function(s){ ss = strsplit(axi(s),' \\(|\\)')[[1]]; ss[len(ss)+1] = '';
  str.lab(str(' ',ss[1],': '),str(' ',ss[2])) }

# factor labsl
fl = list(
  cov = c('–0.9'=-.9,'–0.6'=-.6,'–0.5'=-.5,'–0.3'=-.3,
     '0'=0,'+0.3'=+.3,'+0.5'=+.5,'+0.6'=+.6,'+0.9'=+.9),
  ro = list(
    'Rate: first onset' = 0,
    'Rate: any relapse' = 1,
    'Rate: any onset or relapse' = NA,
    'RR: relapse vs onset' = 'RR'),
  fup = c('all available'=0,'first year'=1))

# colormaps
cmap = lapply(list(mo='rocket',ho='rocket',mx='mako',hx='mako'),
  function(o){ clr.map.d(option=o) })

# grey rects
yyy = list(dep.now=c(01,12,20),dep.past=c(15,30,60)) # rect,rect,lim
rect = def.args(annotate,'rect',xmin=-Inf,xmax=+Inf,
  alpha=1/2,fill='#ccc',color='#ccc',lty='11')

plot.core = function(g,id,ylim,ribbon=1/5,ci=.95){
  yi = if.null(yyy[[id]],c(-Inf,+Inf,NA))
  g = g + rect(ymin=-Inf,ymax=yi[1]) + rect(ymax=+Inf,ymin=yi[2]) +
    stat_summary(geom='ribbon',color=NA,alpha=ribbon,
      fun.min=qfun((1-ci)/2),fun.max=qfun(1-(1-ci)/2)) +
    stat_summary(geom='line',fun=median) +
    coord_cartesian(ylim=c(0,yi[3])) +
    ylab(axi(l[[id]]))
}

plot.exact = function(Y){
  geom_point(data=run.exact(Y),shape=21,fill='#fc0',size=1)
}

plot.1o = list(w1=2.0,h1=2,wo=2,ho=1) # plot size

plot.save.i = function(g,...,ext='.png',font='Alegreya Sans'){
  g = plot.clean(g,font=font,legend.spacing=unit(0,'mm'))
  plot.save(g,'depr',uid,...,ext=ext)
}

# -----------------------------------------------------------------------------
# plot apps

plot.past = function(){
  Y = load.grid(k='past',i='dep.past',f='ho')
  g = ggplot(Y,aes(y=value,x=mo,color=ho,fill=ho))
  g = plot.core(g,'dep.past') + cmap$ho + labs(x=axi(l$mo),color=grp(l$ho),fill=grp(l$ho))
  plot.save.i(g + plot.exact(Y),'v','past.v')
  plot.save.i(g,'past')
}

plot.now.hom = function(){
  Y = load.grid(k='hom')
  g = ggplot(Y,aes(y=value,x=mo,color=factor(mx),fill=factor(mx)))
  g = plot.core(g,'dep.now') + cmap$mx + labs(x=axi(l$mo),color=grp(l$mx),fill=grp(l$mx))
  plot.save.i(g + plot.exact(Y),'v','now.hom.o.v')
  plot.save.i(g,'now.hom.o')
  g = ggplot(Y,aes(y=value,x=mx,color=factor(mo),fill=factor(mo)))
  g = plot.core(g,'dep.now') + cmap$mo + labs(x=axi(l$mx),color=grp(l$mo),fill=grp(l$mo))
  plot.save.i(g + plot.exact(Y),'v','now.hom.x.v')
  plot.save.i(g,'now.hom.x')
}

plot.now.hetc = function(){
  Y = subset(load.grid(k='hetc',f='cov'), mo==PGi$mo & mx==PGi$mx & cov%in%PGii$cov)
  g = ggplot(Y,aes(y=value,x=ho,color=factor(hx))) + facet_grid('. ~ cov',labeller=fct(l$cov))
  g = plot.core(g,'dep.now',ribbon=0) + cmap$hx + labs(x=axi(l$ho),color=grp(l$hx))
  plot.save.i(g,'now.hetc.o')
  g = ggplot(Y,aes(y=value,x=hx,color=factor(ho))) + facet_grid('. ~ cov',labeller=fct(l$cov))
  g = plot.core(g,'dep.now',ribbon=0) + cmap$ho + labs(x=axi(l$hx),color=grp(l$ho))
  plot.save.i(g,'now.hetc.x')
}

# -----------------------------------------------------------------------------
# main

# run.grid()
# merge.rda()

# plot.past()
# plot.now.hom()
# plot.now.hetc()
