source('sim/meta.r')

# mass = measures of association

elu.names = c('estim','lower','upper')
tt = c('.te','.to') # suffix for vars measured at time of exposure or outcome
eps = 1e-6 # for clipping small & large mass

# -----------------------------------------------------------------------------
# prep survey data & calculate mass

mass.calc = function(mfuns,ve,vo,Q1,Q2,va1=NULL,va2=NULL,vs=NULL,by=NULL,ao1=TRUE){
  # compute 1+ measures of assoc for vo ~ ve + va1 + va2,
  # stratified by vs, from Q1 & Q2 (2 timepoints)
  if (missing(Q1)){ Q1 = Q2[Q2$t==min(Q2$t),] }
  if (missing(Q2)){ Q2 = Q1[Q1$t==max(Q1$t),] }
  va = c(va1,va2)       # adjust vars
  vs = c(vs,'seed')     # strat vars
  by = c(by,'i','seed') # index vars for merging
  Q = merge(by=by,suffix=tt,
    Q1[c(by,'t',ve,vo,va1)],
    Q2[c(by,'t',ve,vo,va2)])
  if (ao1){ va = c(va,str(vo,tt[1])) }
  A = rbind.lapply(split(Q,Q[vs]),function(Qi){ # strata
    Ai = rbind.lapply(names(mfuns),function(mass){ # for each mass
      elu = mfuns[[mass]](Qi,ve,vo,va) # calculate the elu
      Aim = cbind(Qi[1,vs,drop=FALSE], # strat vars
        te = Qi$t.te[1], # t1 (exposure)
        to = Qi$t.to[1], # t2 (outcome)
        dt = Qi$t.to[1] - Qi$t.te[1], # t2 - t1
        mass = mass, # mass name
        adj = str(va,collapse=', '), # adjust vars
        as.list(elu)) # estim, lower, upper
    })
  })
  A = mass.clean(A)
}

mass.clean = function(A){
  # clip elu - TODO: is this needed?
  A$estim = pmax(eps,pmin(1/eps,A$estim))
  A$lower = pmax(eps,pmin(1/eps,A$lower))
  A$upper = pmax(eps,pmin(1/eps,A$upper))
  return(A)
}

# -----------------------------------------------------------------------------
# glm (general linear model) = default mfun (function to compute mass)

glm.elu = function(Q,ve,vo,va,family,ctx,among=quote(TRUE),...){
  # run glm for: otx(vo) ~ etx(ve) + va
  # return (estim, lower, upper) = ctx(coef(ve)) with 95% CI
  Q = subset(Q,eval(parse(text=among)))
  f = glm.formula(ve,vo,va,...)
  m = no.warn(glm(f,family,Q))
  elu.coef = c(coef(m)[2],confint.default(m,2))
  elu = set.names(ctx(elu.coef),elu.names)
}

glm.formula = function(ve,vo,va,otx='',etx=''){
  # helper to generate formula; otx & etx = transforms as strings
  f = formula(str(
    otx,'(',vo,tt[2],') ~ ',
    etx,'(',ve,tt[1],') + ',
    str(c(1,va),collapse=' + ')))
}

# -----------------------------------------------------------------------------
# basic mass plot

mass.plot = function(A,...,facet=NULL,w=1,ylim=NULL){
  # plot the 95% CI for (estim, lower, upper) as 3 lineranges
  # assuming multiple "observations" of each mass
  pos = position_dodge(width=2/3)
  qfun = function(p){ def.args(quantile,p=p) }
  geom = function(...){ stat_summary(
    fun=mean,fun.min=qfun(.05),fun.max=qfun(.95),
    geom='linerange',position=pos,...)}
  g = ggplot(A,aes(...)) +
    facet_grid(facet) +
    geom_hline(yintercept=1,color='gray') +
    # geom_pointrange(aes(y=estim,ymin=lower,ymax=upper),position=pos) +
    geom(aes(y=estim),lwd=w,  alpha=1/1) +
    geom(aes(y=lower),lwd=w*3,alpha=1/3) +
    geom(aes(y=upper),lwd=w*3,alpha=1/3) +
    scale_y_continuous(trans='log10') +
    coord_cartesian(ylim=ylim) +
    labs(y='Value')
  g = plot.clean(g)
}

mass.add.ref = function(g,value,...){
  # add a reference hline (i.e. target mass value)
  data = data.frame(ulist(value=value,...))
  g = g + geom_hline(data=data,aes(yintercept=value),lty='11')
}
