
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

fit.par.tx = function(F,S0,dir='tx'){
  # forward (tx) or reverse (itx) transform S0 <-> S
  # where S0 ~ unif[0,1] are quantiles of S ~ tx(unif[lo,up])
  S = lapply(seqa(F),function(i){
    Fi  = F[[i]]
    fun = if.null(Fi[[dir]],identity)
    Si  = switch(dir,
      tx  = fun(qunif(S0[,i+1],Fi$min,Fi$max)),
      itx = punif(fun(S0[,i+1]),Fi$min,Fi$max))
  })
  S = cbind(id=S0[,1],as.data.frame(S,col.names=names(F)))
}

fit.par.lhs = function(F,n,seed=666){
  set.seed(seed)
  S = fit.par.tx(F,cbind(id=1:n,lhs::randomLHS(n,len(F))))
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
  z = (Ti$mu - Yi$est.mu) / sqrt(Ti$se^2 + Yi$est.se^2)
  ll = dnorm(z,log=TRUE)
}

# -----------------------------------------------------------------------------
# fitting

# TODO: finish update

fit.run = function(Si,T,P0=NULL,...,aggr=TRUE,.par=FALSE){
  # get log-likelihood given sample Si, targets T, base params P0
  # i.e. get params, run model, run survey, get log-likelihoods
  Ps = get.pars.grid(ulist(P0,Si),...)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  Q  = srv.apply(Ms)
  Ls = ll.srv(Q,T)
  Ls = c(total=sum(Ls),Ls)
  status(2,'fit.run:',
    '\n > L: ',list.str(Ls,join=', ',sig=3),
    '\n > S: ',list.str(Si,join=', ',sig=3))
  if (aggr){ return(Ls[1]) } else { return(Ls) }
}

fit.runs = function(S,T,...){ .v = .verb; .verb <<- 2 # HACK
  # fit.run in parallel for each sample (row) in data.frame S
  L = rbind.lapply(apply(S,1,as.list),fit.run,T=T,...)
  .verb <<- .v; return(L) # HACK
}
