source('utils.r')
source('sim/plot.r')
# -----------------------------------------------------------------------------
# config
tmax = 50  # maximum time
dti  = .01 # integral timestep
dtw  = 7   # model timestep
nint = 10  # num integrals to plot
# time vectors
tw = seq(0,tmax,dtw)    # end of timestep windows
tx = seq(0,tmax,dti)    # for plotting
ti = seq(dti-dtw,0,dti) # for integrating (to be shifted)
# shape funs
tRR.funs = list(
  exp  = function(t,tsc=10,iRR=2){ tRR = 1+(t>0)*(iRR-1)*exp(-t/tsc) },
  ramp = function(t,tsc=30,iRR=2){ tRR = 1+(t>0)*(iRR-1)*pmax(0,1-t/tsc) },
  step = function(t,tsc=30,iRR=2){ tRR = 1+(t>0)*(iRR-1)*(t<=tsc) })
dRR.funs = list(
  exp  = function(t,tsc=10,iRR=2){ tRR = (t<=0)+(t>0)*exp(-t/tsc) },
  ramp = function(t,tsc=30,iRR=2){ tRR = (t<=0)+(t>0)*pmax(0,1-t/tsc) },
  step = function(t,tsc=30,iRR=2){ tRR = (t<=0)+(t>0)*(t<=tsc) } )
# integral function
tRR.int = function(tRR.fun,te,tw){
  # te: event time, tw = end of window
  sum(tRR.fun(ti-te+tw)*dti/dtw)} # integrate over window & normalize
tRR.adj = function(x){ c(1,x[1]-1+x[2],x[3:len(x)]) }
# -----------------------------------------------------------------------------
# integrate windows
xRR.win = function(xRR.funs){
  X.win = do.call(rbind,lapply(names(xRR.funs),function(shape){
    tRR.fun = xRR.funs[[shape]]
    tRR.w.all  = rowMeans(sapply(ti,function(tei){
                 sapply(tw,tRR.int,tRR.fun=tRR.fun,te=tei) }))
    tRR.w.mean = sapply(tw,tRR.int,tRR.fun=tRR.fun,te=-dtw/2)
    Xi = cbind(shape=shape,rbind(
      data.frame(t=tw,case='Single*',tRR=tRR.adj(tRR.w.mean)),
      data.frame(t=tw,case='Single', tRR=tRR.w.mean),
      data.frame(t=tw,case='Double', tRR=tRR.w.all) ))
}))}
xRR.int = function(xRR.funs){
  X.int = do.call(rbind,lapply(names(xRR.funs),function(shape){
    tRR.fun = xRR.funs[[shape]]
    Xi = merge(all=TRUE,cbind(i=0:nint),
      data.frame(shape=shape,case='',t=c(0,tx),tRR=c(1,tRR.fun(tx))))
}))}
# -----------------------------------------------------------------------------
# plot
xRR.plot = function(X.win,X.int){
  g = ggplot(X.win,aes(x=t-dtw/2,y=tRR,shape=case,color=case)) +
    facet_grid('shape',scales='free_y',labeller=label_both) +
    geom_hline(yintercept=1,color='gray') +
    geom_ribbon(data=X.int,aes(x=t-dtw/2+dtw*i/nint,ymin=1,ymax=tRR,fill=i,group=i),
      alpha=2/nint,color=NA,show.legend=FALSE) +
    geom_step(aes(lty=case)) +
    scale_linetype_manual(name='Integral',values=c('solid','62','22')) +
    scale_color_manual(name='Integral',values=c('#ff00ff','#990099','#000000')) +
    scale_fill_distiller(palette='Spectral') +
    scale_x_continuous(breaks=tw) +
    scale_y_continuous(breaks=seq(0,5,.5)) +
    labs(x='Time since event timestep midpoint (days)')
  g = plot.clean(g)
}
# -----------------------------------------------------------------------------
# main
for (f in c('t','d')){
  xRR.funs = get(str(f,'RR.funs'))
  g = xRR.plot(xRR.win(xRR.funs),xRR.int(xRR.funs)) +
    ylab(str(list(t='Transient',d='Duration')[[f]],' Rate Ratio (',f,'RR)'))
  plot.save('toy','tRR',str(f,'RR.int'),w=6,h=5)
}
