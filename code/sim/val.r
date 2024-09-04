source('sim/meta.r')

# =============================================================================
# config

v      = cli.arg('v',     1:24)
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
  Ri.m=list(gpar=ulist(P0),evts=evts[1:7])
)

# -----------------------------------------------------------------------------

RR1 = c(3)
RR2 = c(1,3)
RR3 = c(1,2,3)
ts2 = c(30,90)
ds2 = c(90,720)
ns2 = c(10,100)

val.RR = list(
  # RR age [1:4]
  aRR.vio=list(gpar=ulist(P0),evts='vio',  strat='TODO'),
  aRR.dep=list(gpar=ulist(P0),evts='dep_o',strat='TODO'),
  aRR.haz=list(gpar=ulist(P0),evts='haz_o',strat='TODO'),
  aRR.ptr=list(gpar=ulist(P0),evts='ptr_o',strat='TODO'),
  # basic RR [5:12]
  RR.dep_o.dep_p=list(gpar=ulist(P0,RR.dep_o.dep_p=  RR3),evts='dep_o',strat='dep.past'),
  RR.haz_o.haz_p=list(gpar=ulist(P0,RR.haz_o.haz_p=  RR3),evts='haz_o',strat='haz.past'),
  RR.haz_o.dep_w=list(gpar=ulist(P0,RR.haz_o.dep_w=  RR3),evts='haz_o',strat='dep.now'),
  RR.haz_x.dep_w=list(gpar=ulist(P0,RR.haz_x.dep_w=1/RR3),evts='haz_x',strat='dep.now'),
  RR.ptr_o.dep_w=list(gpar=ulist(P0,RR.ptr_o.dep_w=1/RR3),evts='ptr_o',strat='dep.now'),
  RR.ptr_o.haz_w=list(gpar=ulist(P0,RR.ptr_o.haz_w=  RR3),evts='ptr_o',strat='haz.now'),
  RR.ptr_x.dep_w=list(gpar=ulist(P0,RR.ptr_x.dep_w=  RR3),evts='ptr_x',strat='dep.now'),
  RR.ptr_x.haz_w=list(gpar=ulist(P0,RR.ptr_x.haz_w=  RR3),evts='ptr_x',strat='haz.now'),
  # transient RR [13:17]
  tRR.dep_o.vio_zr=list(gpar=ulist(P0,iRR.dep_o.vio_zr=  RR2,tsc.dep_o.vio_zr=ts2),evts='dep_o',strat='TODO'),
  tRR.dep_x.vio_zr=list(gpar=ulist(P0,iRR.dep_x.vio_zr=1/RR2,tsc.dep_x.vio_zr=ts2),evts='dep_x',strat='TODO'),
  tRR.haz_o.vio_zr=list(gpar=ulist(P0,iRR.haz_o.vio_zr=  RR2,tsc.haz_o.vio_zr=ts2),evts='haz_o',strat='TODO'),
  tRR.haz_x.vio_zr=list(gpar=ulist(P0,iRR.haz_x.vio_zr=1/RR2,tsc.haz_x.vio_zr=ts2),evts='haz_x',strat='TODO'),
  tRR.ptr_o.vio_zr=list(gpar=ulist(P0,iRR.ptr_o.vio_zr=  RR2,tsc.ptr_o.vio_zr=ts2),evts='ptr_o',strat='TODO'),
  # cumulative RR [18:20]
  nRR.dep_o.vio_nt=list(gpar=ulist(P0,mRR.dep_o.vio_nt=RR2,nsc.dep_o.vio_nt=ns2),evts='dep_o',strat='TODO'),
  nRR.haz_o.vio_nt=list(gpar=ulist(P0,mRR.haz_o.vio_nt=RR2,nsc.haz_o.vio_nt=ns2),evts='haz_o',strat='TODO'),
  nRR.ptr_o.vio_nt=list(gpar=ulist(P0,mRR.ptr_o.vio_nt=RR2,nsc.ptr_o.vio_nt=ns2),evts='ptr_o',strat='TODO'),
  # duration RR [21:22]
  dRR.dep_x.dep_u=list(gpar=ulist(P0,dsc.dep_x.dep_u=ds2),evts='dep_x',strat='TODO'),
  dRR.haz_x.haz_u=list(gpar=ulist(P0,dsc.haz_x.haz_u=ds2),evts='haz_x',strat='TODO')
)

# =============================================================================
# run & plot

val.run = function(name,gpar,evts,strat='.'){
  status(2,'val.run: ',name,' @ ',n.seed*prod(lens(gpar)))
  Ps = grid.apply(ulist(gpar,seed=1:n.seed),get.pars)
  Ms = sim.runs(Ps)
  facet = setdiff(names(gpar)[lens(gpar)>1],strat)
  fixed = setdiff(names(gpar),c(strat,facet))
  Y = cbind(rate.datas(Ms,p.vars=facet),.='')
  R = rbind.lapply(evts,rate.est,Y=Y,strat=c('seed',strat,facet))
  for (evt in evts){
    g = val.plot.rate(R,evt,gpar,strat,facet,fixed)
    plot.save('val',uid,str(name,'--',evt),
      h=3,w=2+2.5*prod(lens(gpar[facet])))
  }
}

val.plot.rate = function(R,evt,gpar,strat,facet,fixed,t1y=364){
  R = subset(R,event==evt)
  R = cbind(R,facet=apply(R[facet],1,list.str,sig=3))
  R.ref = gpar[[str(evt,'.Ri.m')]]
  R.max = max(R$R)
  fix = list.str(gpar[fixed],sig=3)
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
