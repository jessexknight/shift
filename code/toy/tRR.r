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
 step = function(t,tsc=30,iRR=2){ tRR = 1+(t>0)*(iRR-1)*(t<=tsc) } )
# integral function
tRR.int = function(tRR.fun,te,tw){
  # te: event time, tw = end of window
  sum(tRR.fun(ti-te+tw)*dti/dtw)} # integrate over window & normalize
# -----------------------------------------------------------------------------
# integrate windows
X.win = do.call(rbind,lapply(names(tRR.funs),function(shape){
  tRR.fun = tRR.funs[[shape]]
  tRR.w.all  = rowMeans(sapply(ti,function(tei){
               sapply(tw,tRR.int,tRR.fun=tRR.fun,te=tei) }))
  tRR.w.mean = sapply(tw,tRR.int,tRR.fun=tRR.fun,te=-dtw/2)
  Xi = cbind(shape=shape,rbind(
      data.frame(t=tw,tRR=tRR.w.mean,case='midpoint'),
      data.frame(t=tw,tRR=tRR.w.all, case='complete') ))
}))
X.int = do.call(rbind,lapply(names(tRR.funs),function(shape){
  tRR.fun = tRR.funs[[shape]]
  Xi = merge(all=TRUE,cbind(i=0:nint),
    data.frame(shape=shape,case='',t=c(0,tx),tRR=c(1,tRR.fun(tx))))
}))
# -----------------------------------------------------------------------------
# plot
g = ggplot(X.win,aes(x=t,y=tRR,shape=case,color=case)) +
  facet_grid('shape') +
  geom_ribbon(data=X.int,aes(x=t+dtw*i/nint,ymin=1,ymax=tRR,fill=i,group=i),
    alpha=2/nint,color=NA,show.legend=FALSE) +
  geom_step(aes(lty=case)) +
  scale_linetype_manual(name='Integral',values=c('solid','33')) +
  scale_color_manual(name='Integral',values=c('purple','black','black')) +
  scale_fill_distiller(palette='Spectral') +
  scale_x_continuous(breaks=tw) +
  labs(x='Time since event (days)',y='Transient Rate Ratio (tRR)') +
  theme_light()
g = plot.clean(g)
plot.save('toy','tRR',str('tRR.int'),w=8,h=5)

