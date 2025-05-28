stfu(library(mlrMBO))

# -----------------------------------------------------------------------------
# fitted params

gen.fpar = makeNumericParam
fpar.set = makeParamSet
length.ParamSet = function(F){ length(F$pars) }
fpar.sam = function(...,seed=666){ set.seed(seed); generateDesign(...) }

# -----------------------------------------------------------------------------
# targets

gen.targ = function(id,type,...,sub=NULL,mu=NA,se=NA,w=1){
  Ti = list(id=id,type=type,mu=mu,se=se,sub=sub,w=w,
    fun=def.args(targ.funs[[type]],sub=sub,...))
}

targ.calc = function(Q,ofun,...,vs=NULL){
  # compute outputs from Q for simple targets
  vs = c(vs,'t','seed')
  Y = rbind.lapply(split(Q,Q[vs]),function(Qi){ # strata
    out = ofun(Qi,...) # calculate estimate
    Yi = cbind(Qi[1,vs,drop=FALSE],as.list(out))
  },.par=FALSE)
}

prop.out = function(Q,vo,sub=NULL,vsub=FALSE){
  Q = df.sub(Q,sub)
  x = if (vsub){ df.sub(Q,vo)$. } else { Q[[vo]] }
  k = sum(x); n = nrow(Q); p = k/n
  out = list(
    est.mu = p,
    est.se = p*(1-p)/n,
    value = p,
    lower = qbeta(.025,k+.5,n-k+.5),
    upper = qbeta(.975,k+.5,n-k+.5))
}

pois.out = function(Q,vo,vt,sub=NULL,vt.na='act.ut'){
  Q = df.sub(Q,sub)
  k = sum(Q[[vo]])
  t = sum(if.na(Q[[vt]],Q[[vt.na]]))
  p = k/t; u = log(p); use = 1/sqrt(t*p)
  out = list(
    est.mu = u,
    est.se = use,
    value = p,
    lower = exp(u+qnorm(.025)*use),
    upper = exp(u+qnorm(.975)*use))
}

targ.funs = list(
  prop = def.args(targ.calc,ofun=prop.out),
  pois = def.args(targ.calc,ofun=pois.out),
  OR   = def.args(mass.calc,ofun=def.args(glm.out,family=binomial,ctx=exp)),
  PR   = def.args(mass.calc,ofun=def.args(glm.out,family=poisson, ctx=exp)))
  # TODO: quasipoisson above?

# -----------------------------------------------------------------------------
# log-likelihoods

srv.targs = function(Q,T,vs=NULL,aggr.seed=FALSE){
  # log-likelihoods for each target T given survey Q
  if (aggr.seed){ Q$seed = 0 }
  Ys = lapply(T,function(Ti){
    Yi = Ti$fun(Q,vs=vs) # estimate
    ll = targ.ll(Ti,Yi)  # likelihood
    Yi = cbind(id=Ti$id,type=Ti$type,targ.mu=Ti$mu,targ.se=Ti$se,ll=ll,Yi)
  })
  Y = rbind.lapply(Ys,`[`,Reduce(intersect,lapply(Ys,colnames)),.par=FALSE)
}

targ.ll = function(Ti,Yi){
  z = (Ti$mu - Yi$est.mu) / sqrt(Ti$se^2 + pmin(Yi$est.se,Ti$se)^2)
  ll = dnorm(z,log=TRUE)
}

# -----------------------------------------------------------------------------
# fitting

fit.run = function(Si,T,P0=NULL,...,srvs=NULL,aggr=FALSE,.par=TRUE,
                   p.vars=NULL,i.vars=NULL){
  # get log-likelihoods given sample Si, targets T, base params P0
  # i.e. get params, run model, run survey, get log-likelihoods
  Ps = get.pars.grid(ulist(P0,Si),...,.par=.par)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  Q  = srv.apply(Ms,srvs=srvs,.par=.par,p.vars=p.vars,i.vars=i.vars)
  Y  = srv.targs(Q,T,aggr=aggr)
  Y  = cbind(as.list(Si),Y,row.names=NULL)
}

fit.run.grid = function(PG,T,P0=NULL,srvs=NULL,aggr=FALSE,.par=TRUE,
                        p.vars=NULL,i.vars=NULL,.batch=1,.nbatch=1){
  # fit.run over the grid of params PG & rbind the results
  gs = c(lens(PG),seed=len(if.null(P0$seed,1:7)))
  status(2,'fit.run.grid: ',prod(gs),' @ ',list.str(gs))
  Y = grid.apply(PG,function(...){
    status(3,list.str(list(...)))
    Yi = verb.wrap(fit.run(Si=list(...),T=T,P0=P0,srvs=srvs,aggr=aggr,.par=FALSE,
      p.vars=p.vars,i.vars=i.vars),0)
  },.rbind=TRUE,.par=.par,.batch=.batch,.nbatch=.nbatch)
}

opt.run = function(F,T,P0=NULL,...,h.init=4,n.iter=100){
  # multi-objective optimize fitted params F given targets T using mlrMBO
  fun = function(Si){
    Y = verb.wrap(fit.run(Si,T=T,P0=P0,...,aggr=TRUE),0)
    obj = -Y$ll }
  J  = makeMultiObjectiveFunction(fn=fun,par=F,n.obj=len(T))
  C  = stfu(makeMBOControl(n.obj=len(T),y.name=names(T)))
  C  = setMBOControlInfill(C,makeMBOInfillCritDIB())
  C  = setMBOControlTermination(C,iters=n.iter)
  S0 = fpar.sam(n=h.init*len(F$pars)+1,F)
  O  = mbo(J,control=C,design=S0,show.info=TRUE)
}
