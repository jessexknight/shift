source('sim/meta.r')

# -----------------------------------------------------------------------------
# prep survey data & calculate mass = measures of association

mass.calc = function(ofun,ve,vo,Q1,Q2,
    va1=NULL,va2=NULL,vs=NULL,by=NULL,among=NULL,ao1=TRUE){
  # compute a measure of assoc for vo ~ ve + va1 + va2,
  # stratified by vs, from Q1 & Q2 (2 timepoints)
  if (missing(Q1)){ Q1 = Q2[Q2$t==min(Q2$t),] }
  if (missing(Q2)){ Q2 = Q1[Q1$t==max(Q1$t),] }
  va = c(va1,va2)       # adjust vars
  vs = c(vs,'seed')     # strat vars
  by = c(by,'i','seed') # index vars for merging
  Q = merge(by=by,suffix=c('.te','.to'),
    Q1[c(by,'t',ve,vo,va1)],
    Q2[c(by,'t',ve,vo,va2)])
  if (ao1){ va = c(va,str(vo,'.te')) }
  Q = srv.sub(Q,among)
  A = rbind.lapply(split(Q,Q[vs]),function(Qi){ # strata
    out = ofun(Qi,ve,vo,va) # calculate mass
    Ai = cbind(Qi[1,vs,drop=FALSE], # strat vars
      ve = ve,
      vo = vo,
      te = Qi$t.te[1], # t1 (exposure)
      to = Qi$t.to[1], # t2 (outcome)
      dt = Qi$t.to[1] - Qi$t.te[1], # t2 - t1
      adj = str(va,collapse=', '), # adjust vars
      as.list(out)) # est.mu, est.se, value, lower, upper
  },.par=FALSE)
}

# -----------------------------------------------------------------------------
# glm (general linear model) = default ofun (function to compute mass)

glm.run = function(Q,ve,vo,va,family,...){
  # run glm for: otx(vo) ~ etx(ve) + va
  f = glm.formula(ve,vo,va,...)
  m = no.warn(glm(f,family,Q))
}

glm.formula = function(ve,vo,va,otx='',etx=''){
  # helper to generate formula; otx & etx = transforms as strings
  f = formula(str(
    otx,'(',vo,'.to) ~ ',
    etx,'(',ve,'.te) + ',
    str(c(1,va),collapse=' + ')))
}

glm.out = function(...,ctx=identity){
  # run glm & extract outputs = (est.mu, est.se, estim, lower, upper)
  # for first non-intercept coef
  m = glm.run(...)
  est = summary(m)$coef[2,1:2] # raw coef & std err
  eci = confint.default(m,2) # transformed coef & 95% CI
  out = list(
    est.mu = est[1],
    est.se = est[2],
    value = ctx(est[1]),
    lower = ctx(eci[1]),
    lower = ctx(eci[2]))
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
