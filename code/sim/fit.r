
# -----------------------------------------------------------------------------
# params

fit.par = function(name,min,max,tx=NULL,itx=NULL){
  Fi = list(
    name = name,
    min  = min,
    max  = max,
    tx   = if.null(tx,identity),
    itx  = if.null(itx,identity))
}

fit.sam.tx = function(F,S,d='tx'){
  # forward (tx) or reverse (itx) transform S0 <-> S
  # where S0 ~ unif[0,1] are quantiles of S ~ tx(unif[lo,up])
  for (Fi in F){
    fun = if.null(Fi[[d]],identity)
    S[[Fi$name]] = switch(d,
      tx  = fun(qunif(S[[Fi$name]], Fi$min,Fi$max)),
      itx = punif(fun(S[[Fi$name]]),Fi$min,Fi$max))
  }
  return(S)
}

fit.sam.lhs = function(F,n,seed=666){
  set.seed(seed)
  S0 = lhs::randomLHS(n,len(F)); colnames(S0) = names(F)
  S  = fit.sam.tx(F,cbind(id=1:n,as.data.frame(S0)))
}

# -----------------------------------------------------------------------------
# targets

fit.targ = function(name,type,mu,se,...,w=1){
  Ti = list(
    name = name,
    type = type,
    mu = mu, # estimate
    se = se, # std err
    w = w, # weight
    fun = def.args(targ.funs[[type]],...))
}

targ.calc = function(Q,ofun,...,vs=NULL){
  # compute outputs from Q for simple targets
  vs = c(vs,'t','seed')
  Y = rbind.lapply(split(Q,Q[vs]),function(Qi){ # strata
    out = ofun(Qi,...) # calculate estimate
    Yi = cbind(Qi[1,vs,drop=FALSE],as.list(out))
  },.par=FALSE)
}

prop.out = function(Q,vo){
  x = Q[[vo]]; k = sum(x); n = len(x); p = k/n
  out = list(
    est.mu = p,
    est.se = p*(1-p)/n,
    value = p,
    lower = qbeta(.025,k+.5,n-k+.5),
    upper = qbeta(.975,k+.5,n-k+.5))
}

pois.out = function(Q,vo,vt){
  k = sum(Q[[vo]]); t = sum(Q[[vt]]); p = k/t
  u = log(p); use = 1/sqrt(t*p)
  out = list(
    est.mu = u,
    est.st = use,
    value = p,
    lower = exp(u+qnorm(.025)*use),
    upper = exp(u+qnorm(.975)*use))
}

targ.funs = list(
  prop = def.args(targ.calc,ofun=prop.out),
  pois = def.args(targ.calc,ofun=pois.out),
  OR   = def.args(mass.calc,ofun=def.args(glm.out,family=binomial,ctx=exp)),
  PR   = def.args(mass.calc,ofun=def.args(glm.out,family=poisson, ctx=exp)))

# -----------------------------------------------------------------------------
# log-likelihoods

srv.targs = function(Q,T){
  # weighted log-likelihoods for each target T given survey Q
  Ys = lapply(T,function(Ti){
    Yi = Ti$fun(Q)      # estimate
    ll = targ.ll(Ti,Yi) # likelihood
    Yi = cbind(name=Ti$name,targ.mu=Ti$mu,targ.se=Ti$se,ll=ll,Yi)
  })
  Y = rbind.lapply(Ys,`[`,Reduce(intersect,lapply(Ys,colnames)))
}

targ.ll = function(Ti,Yi){
  z = (Ti$mu - Yi$est.mu) / sqrt(Ti$se^2 + pmin(Yi$est.se,Ti$se)^2)
  ll = dnorm(z,log=TRUE)
}

# -----------------------------------------------------------------------------
# fitting

fit.run = function(Si,T,P0=NULL,...,.par=FALSE){
  # get log-likelihood given sample Si, targets T, base params P0
  # i.e. get params, run model, run survey, get log-likelihoods
  Ps = get.pars.grid(ulist(P0,Si),...)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  Q  = srv.apply(Ms)
  Y  = srv.targs(Q,T)
  status(2,'fit.run:',
    ' [L] : ',sum(aggregate(ll~name,Y,mean)$ll),
    ' [S] : ',list.str(Si,join=', ',sig=3))
  Y = cbind(Si,Y)
}

fit.runs = function(S,T,...){
  # fit.run in parallel for each sample (row) in data.frame S
  Y = verb.wrap(rbind.lapply(apply(S,1,as.list),fit.run,T=T,...),2)
}
