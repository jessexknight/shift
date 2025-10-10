source('utils.r')
library('ggplot2')
library('reshape2')

# config & dummy data ---------------------------------------------------------
set.seed(666)
chains = 7
iter = 5000
burn = 2000
i = seq(iter)[-seq(burn)]
options(mc.cores=chains)
data = list(
  Py = c(dep=.03),                 # outcome prevalence
  Pe = c(sex=.50,ace=.25,sle=.35), # exposure prevalence
  OR = c(sex=2.0,ace=2.8,sle=2.5), # odds ratios
  N  = 1e4)                        # total population size
X = expand.grid(lapply(data$Pe,function(P){ c(FALSE,TRUE) }))
data$Ne = len(data$Pe) # num exposures
data$Ns = 2^data$Ne    # num total strata
data$se0 = lapply(X,function(Xs){ which(!Xs) }) # indices of non-exposed
data$se1 = lapply(X,function(Xs){ which( Xs) }) # indices of exposed
# build & run model -----------------------------------------------------------
pars  = c(
  'mean prevalence'='pym',
  'SD prevalence'='pyd',
  'CV = SD / mean'='cv',
  'cumulative population'='ocps',
  'prevalence by strata'='opys')
model = rstan::stan_model('toy/cv/infer.stan',auto_write=TRUE)
fit = rstan::sampling(model,data=data,pars=pars,chains=chains,iter=iter)
# clean & plot outputs
S0 = rbind.lapply(1:chains,function(k){
  Sk = as.data.frame(fit@sim$samples[[k]])
  Sk = cbind(chain=factor(k),i=i,Sk[i,])
})
# param posterior -------------------------------------------------------------
sapply(S0[,pars[1:3]],qfun(p5)) # numeric summary
sapply(S0[,pars[1:3]],mean) # numeric summary
S = melt(S0,id=c('chain','i'),m=pars[1:3],vari='f')
S$f = factor(S$f,pars,names(pars))
g = ggplot(S,aes(x=value,color=chain)) +
  facet_wrap('f',scales='free',ncol=3) +
  geom_density() +
  clr.map.d(aes='color') +
  labs(y='Density',x='Value',color='Chain')
g = plot.clean(g,axis.text.y=element_blank())
plot.save(g,'toy','het.infer.par',size=c(7,2.5),ext='.png')
# prev dist posterior ---------------------------------------------------------
po = seq(0,1,.01)
vars = c('Density'='p','Cumulative'='cp')
cols = list(cp=filter.names(S0,'ocps'),py=filter.names(S0,'opys'))
S = rbind.lapply(1:nrow(S0),function(j){ # interp to same x
  cp = approx(c(0,S0[j,cols$py],1),c(0,S0[j,cols$cp],1),po,method='constant')$y
  p = c(diff(cp),0)
  cbind(chain=S0$chain[j],i=S0$i[j],prev=po,cp=cp,p=p)
})
S = melt(as.data.frame(S),m=vars,vari='f')
S$f = factor(S$f,vars,names(vars))
g = ggplot(S,aes(y=100*value,x=100*prev)) +
  facet_grid('f',scales='free') +
  stat_summary(fun.min=qfun(.025),fun.max=qfun(.975),geom='ribbon',alpha=1/4) +
  geom_step(data=subset(S,i==i[1]),aes(color=factor(chain)),position=dodge(1/4)) +
  stat_summary(fun=mean,geom='line') +
  clr.map.d(guide='none') +
  lims(x=c(NA,50)) +
  labs(y='Probability (%)',x='Prevalence (%)')
g = plot.clean(g)
plot.save(g,'toy','het.infer.dist',size=c(5,5),ext='.png')
