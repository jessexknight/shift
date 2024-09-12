source('sim/meta.r')

# =============================================================================
# config

v      = cli.arg('v',     1:23)
n.pop  = cli.arg('n.pop', 333)
n.seed = cli.arg('n.seed',7)

# -----------------------------------------------------------------------------

P0 = list(
  n.pop = n.pop,
  vio.Ri.m   = .001,
  dep_o.Ri.m = .001, dep_x.Ri.m = .002,
  haz_o.Ri.m = .001, haz_x.Ri.m = .002,
  ptr_o.Ri.m = .01,  ptr_x.Ri.m = .02,  ptr.max.m = 3,
  all.Ri.shape = 1e6,
  ptr.Ri.shape = 1e6,
  null = 'xRR')

# =============================================================================
# scenarios

val.Ri = list(
  Ri.m=list(pars=ulist(P0),evts=evts[1:7])
)

# -----------------------------------------------------------------------------

RR1 = c(3)
RR2 = c(1,3)
RR3 = c(1,2,3)
ts2 = c(35,70)
ds2 = c(70,700)
ns2 = c(3,10)
age5 = c(10,20,30,40,50)
aRR5 = c( 1, 5, 2, 4, 3)/3

age.10.qc   = quote(floor(Y$age.1/10)*10)
vio.nt.c.qc = quote(int.cut(Y$vio.nt,c(0,ns2)))

val.RR = list(
  # RR age [1:4]
  aRR.vio=list(pars=ulist(P0,aRR.vio.RRs  =aRR5,aRR.vio.ages  =age5,aRR.shape='const'),evts='vio',  strat='age.10',x.cols=list(age.10=age.10.qc)),
  aRR.dep=list(pars=ulist(P0,aRR.dep_o.RRs=aRR5,aRR.dep_o.ages=age5,aRR.shape='const'),evts='dep_o',strat='age.10',x.cols=list(age.10=age.10.qc)),
  aRR.haz=list(pars=ulist(P0,aRR.haz_o.RRs=aRR5,aRR.haz_o.ages=age5,aRR.shape='const'),evts='haz_o',strat='age.10',x.cols=list(age.10=age.10.qc)),
  aRR.ptr=list(pars=ulist(P0,aRR.ptr_o.RRs=aRR5,aRR.ptr_o.ages=age5,aRR.shape='const'),evts='ptr_o',strat='age.10',x.cols=list(age.10=age.10.qc)),
  # basic RR [5:12]
  RR.dep_o.dep_p=list(pars=ulist(P0,RR.dep_o.dep_p=  RR3),evts='dep_o',strat='dep.past'),
  RR.haz_o.haz_p=list(pars=ulist(P0,RR.haz_o.haz_p=  RR3),evts='haz_o',strat='haz.past'),
  RR.haz_o.dep_w=list(pars=ulist(P0,RR.haz_o.dep_w=  RR3),evts='haz_o',strat='dep.now'),
  RR.haz_x.dep_w=list(pars=ulist(P0,RR.haz_x.dep_w=1/RR3),evts='haz_x',strat='dep.now'),
  RR.ptr_o.dep_w=list(pars=ulist(P0,RR.ptr_o.dep_w=1/RR3),evts='ptr_o',strat='dep.now'),
  RR.ptr_o.haz_w=list(pars=ulist(P0,RR.ptr_o.haz_w=  RR3),evts='ptr_o',strat='haz.now'),
  RR.ptr_x.dep_w=list(pars=ulist(P0,RR.ptr_x.dep_w=  RR3),evts='ptr_x',strat='dep.now'),
  RR.ptr_x.haz_w=list(pars=ulist(P0,RR.ptr_x.haz_w=  RR3),evts='ptr_x',strat='haz.now'),
  # transient RR [13:17]
  tRR.dep_o.vio_zr=list(pars=ulist(P0,iRR.dep_o.vio_zr=  RR2,tsc.dep_o.vio_zr=ts2,tRR.shape='step'),evts='dep_o',strat='vio.dt',e.dts=list(vio=ts2)),
  tRR.dep_x.vio_zr=list(pars=ulist(P0,iRR.dep_x.vio_zr=1/RR2,tsc.dep_x.vio_zr=ts2,tRR.shape='step'),evts='dep_x',strat='vio.dt',e.dts=list(vio=ts2)),
  tRR.haz_o.vio_zr=list(pars=ulist(P0,iRR.haz_o.vio_zr=  RR2,tsc.haz_o.vio_zr=ts2,tRR.shape='step'),evts='haz_o',strat='vio.dt',e.dts=list(vio=ts2)),
  tRR.haz_x.vio_zr=list(pars=ulist(P0,iRR.haz_x.vio_zr=1/RR2,tsc.haz_x.vio_zr=ts2,tRR.shape='step'),evts='haz_x',strat='vio.dt',e.dts=list(vio=ts2)),
  tRR.ptr_o.vio_zr=list(pars=ulist(P0,iRR.ptr_o.vio_zr=  RR2,tsc.ptr_o.vio_zr=ts2,tRR.shape='step'),evts='ptr_o',strat='vio.dt',e.dts=list(vio=ts2)),
  # cumulative RR [18:20]
  nRR.dep_o.vio_nt=list(pars=ulist(P0,mRR.dep_o.vio_nt=RR2,nsc.dep_o.vio_nt=ns2,nRR.shape='step'),evts='dep_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.qc)),
  nRR.haz_o.vio_nt=list(pars=ulist(P0,mRR.haz_o.vio_nt=RR2,nsc.haz_o.vio_nt=ns2,nRR.shape='step'),evts='haz_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.qc)),
  nRR.ptr_o.vio_nt=list(pars=ulist(P0,mRR.ptr_o.vio_nt=RR2,nsc.ptr_o.vio_nt=ns2,nRR.shape='step'),evts='ptr_o',strat='vio.nt.c',x.cols=list(vio.nt.c=vio.nt.c.qc)),
  # duration RR [21:22]
  dRR.dep_x.dep_u=list(pars=ulist(P0,dsc.dep_x.dep_u=ds2,dRR.shape='step'),evts='dep_x',strat='dep_o.dt',e.dts=list(dep_o=ds2)),
  dRR.haz_x.haz_u=list(pars=ulist(P0,dsc.haz_x.haz_u=ds2,dRR.shape='step'),evts='haz_x',strat='haz_o.dt',e.dts=list(haz_o=ds2))
)

# =============================================================================
# run & plot

val.run = function(name,pars,evts,strat='.',e.dts=NULL,x.cols=NULL){
  pars = val.par.split(pars)
  status(2,'val.run: ',name,' @ ',n.seed*pars$n.var)
  Ps = grid.apply(ulist(pars$var,seed=1:n.seed),get.pars,pars$fix,.par=FALSE)
  Ms = sim.runs(Ps)
  Y = cbind(rate.datas(Ms,p.vars=names(pars$var),e.dts=e.dts,x.cols=x.cols),.='')
  R = rbind.lapply(evts,rate.est,Y=Y,strat=c('seed',strat,names(pars$var)))
  for (evt in evts){
    g = val.plot.rate(R,evt,pars,strat)
    plot.save('val',uid,str(name,'--',evt),h=3,w=2+2.5*pars$n.var)
  }
}

val.par.split = function(par){
  b = lens(par) > 1 & !grepl(names(null.sets$aRR),names(par))
  p = list(var=par[b],fix=par[!b],n.var=prod(lens(par[b])))
}

val.plot.rate = function(R,evt,pars,strat,t1y=364){
  R = subset(R,event==evt)
  R = cbind(R,facet=apply(R[names(pars$var)],1,list.str,sig=3))
  R.ref = pars$fix[[str(evt,'.Ri.m')]]
  R.max = max(R$R)
  fix = list.str(pars$fix,sig=3)
  g = ggplot(R,aes(x='',y=t1y*R,color=as.factor(.data[[strat]]))) +
    geom_boxplot(alpha=.3,outlier.shape=1,outlier.alpha=1) +
    geom_point(data=data.frame(R=R.ref),shape=9,color='red') +
    facet_wrap('~facet',scales='fixed',ncol=ulen(R$facet)) +
    labs(y=str('Rate (per year): ',evt),x='',color=str(fix,'\n\n',strat)) +
    guides(color=guide_legend(ncol=2)) +
    ggtitle(str(name,': ',evt)) +
    scale_color_viridis_d() +
    ylim(c(0,NA))
  g = add.rate.label(g,R,c('facet',strat),R.max+R.ref,function(Ri){
    str(signif(median(Ri$R/R.ref),3)) })
  g = add.rate.label(g,R,c('facet',strat),0,function(Ri){
    str(round(median(Ri$ne)),'\n',round(median(Ri$dt)/t1y),'\n') })
  g = plot.clean(g,legend.title=element_text(size=9))
}

add.rate.label = function(g,R,strat,Ry,lab.fun){
  R.lab = rbind.lapply(split(R,R[strat]),function(Ri){
    cbind(Ri[1,strat],R=Ry,label=lab.fun(Ri)) })
  g = g + geom_text(data=R.lab,aes(label=label),size=2.5,
    position=position_dodge(width=.75),show.legend=FALSE)
}

# =============================================================================
# main

vals = c(val.Ri,val.RR)[v]
for (name in names(vals)){ vals[[name]]$name = name }
for (val in vals){ do.call(val.run,val,quote=TRUE) }
