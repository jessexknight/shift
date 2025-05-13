source('utils.r')
library('ggplot2')
library('reshape2')
# config & dummy data
chains = 7
iter = 5000
burn = 1000
i = seq(iter)[-seq(burn)]
options(mc.cores=chains)
data = list(
  Py = c(dep=.03),                 # outcome prevalence
  Pe = c(sex=.50,ace=.10,oth=.25), # exposure prevalence
  OR = c(sex=2.0,ace=5.0,oth=3.0), # odds ratios
  N  = 1e4)                        # total population size
X = expand.grid(lapply(data$Pe,function(P){ c(FALSE,TRUE) }))
data$Ne = len(data$Pe) # num exposures
data$Nk = 2^data$Ne    # num total strata
data$e0 = lapply(X,function(Xe){ which(!Xe) }) # indices of non-exposed
data$e1 = lapply(X,function(Xe){ which( Xe) }) # indices of exposed
# build & run model
pars  = c('p',Mean='pym',SD='pys',CV='cv')
model = rstan::stan_model('toy/cv/infer.stan',auto_write=TRUE)
fit = rstan::sampling(model,data=data,pars=pars,chains=chains,iter=iter)
# clean & plot outputs
S0 = rbind.lapply(1:chains,function(k){
  Sk = as.data.frame(fit@sim$samples[[k]])
  Sk = cbind(chain=factor(k),i=i,Sk[i,])
})
S = melt(S0,id=c('chain','i'),vari='f')
S$f = factor(S$f,pars,names(pars))
g = ggplot(S[!is.na(S$f),],aes(x=value,color=chain)) +
  facet_wrap('f',scales='free',ncol=3) +
  geom_density()
g = plot.clean(g,axis.text.y=element_blank())
ggsave('Rplots.png',w=8,h=3)
sapply(S0[,-1],qfun(p5)) # numeric summary
