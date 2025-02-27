
# =============================================================================
# pars

get.pars = function(seed=0,...,dtz=7,case='base',null=NULL,save=NULL,fun=identity){
  P = list(case=case,seed=seed) # meta
  P = add.pars.def(P)         # default (upstream)
  P = add.pars.time(P,dtz)    # timestep-related
  P = null.pars(P,null,save)  # null some Ri,RR
  P = ulist(P,...)            # overwrite any (upstream)
  P = add.pars.cond(P)        # conditional (downstream)
  P = fun(P)                  # anything else
}

add.pars.def = function(P=NULL){
  P$run = list(vio=TRUE,dep=TRUE,haz=TRUE,ptr=TRUE,sex=TRUE)
  # pop size & duration
  P$n.pop = 1000
  P$n.dur = 1+1
  # base rates (per year)
  P$vio.Ri.my   = 0.67     # (mean) base rate: violence
  P$dep_o.Ri.my =  .02     # (mean) base rate: depression begin
  P$dep_x.Ri.my =  .10     # (mean) base rate: depression end
  P$haz_o.Ri.my =  .05     # (mean) base rate: hazdrink begin
  P$haz_x.Ri.my =  .10     # (mean) base rate: hazdrink end
  P$ptr_o.Ri.my = 1.00     # (mean) base rate: partner begin
  P$ptr_x.Ri.my = 0.50     # (mean) base rate: partner end
  P$sex.Ri.95   = c(.1,.5) # (95% CI) base rate: sex within ptr
  P$cdm.Pi.95   = c(.2,.8) # (95% CI) prob: condom use
  # base rate covariance, CoV, etc.
  P$ptr.max.m   = 3.00     # (mean) max num partners
  P$dep.cov     = -.9      # approx covariance among dep_o,dep_x
  P$haz.cov     = -.9      # approx covariance among haz_o,haz_x
  P$ptr.cov     = +.9      # approx covariance among ptr_o,ptr_x,ptr.max
  P$vio.Ri.het  = 0.5      # heterog (gamma CV): vio
  P$dep.Ri.het  = 0.5      # heterog (gamma CV): dep_o,dep_x
  P$haz.Ri.het  = 0.5      # heterog (gamma CV): haz_o,haz_x
  P$ptr.Ri.het  = 0.1      # heterog (gamma CV): ptr_o,ptr_x
  P$het.distr   = 'gamma'  # heterog distr type
  # *RR shapes & aggr
  P$aggr.rate = 'mult'
  P$aRR.shape = 'spline'
  P$tRR.shape = 'ramp'
  P$nRR.shape = 'ramp'
  P$dRR.shape = 'ramp'
  # RR: age -> *
  P$aRR.vio.ages   = c(amin,amax) # (age points) RR: age -> vio
  P$aRR.vio.RRs    = c(1.00,1.00) # (RR  points) RR: age -> vio
  P$aRR.dep_o.ages = c(amin,amax) # (age points) RR: age -> dep begin
  P$aRR.dep_o.RRs  = c(1.00,1.00) # (RR  points) RR: age -> dep begin
  P$aRR.haz_o.ages = c(amin,amax) # (age points) RR: age -> haz begin
  P$aRR.haz_o.RRs  = c(1.00,1.00) # (RR  points) RR: age -> haz begin
  P$aRR.ptr_o.ages = c(amin,  30,amax) # (age points) RR: age -> ptr begin
  P$aRR.ptr_o.RRs  = c(3.00,1.00,0.10) # (RR  points) RR: age -> ptr begin
  # RR: * -> vio
  P$iRR.vio.vio_zr   = 1     # (initial RR) transient RR: vio -> vio
  P$tsc.vio.vio_zr   = 30    # (time scale) transient RR: vio -> vio
  P$mRR.vio.vio_nt   = 1     # (max RR)  cumulative RR: vio -> vio
  P$nsc.vio.vio_nt   = 10    # (n scale) cumulative RR: vio -> vio
  # RR: * -> dep begin
  P$ RR.dep_o.dep_p  = 1     # RR: dep past -> dep begin
  P$iRR.dep_o.vio_zr = 1     # (initial RR) transient RR: vio -> dep begin
  P$tsc.dep_o.vio_zr = 30    # (time scale) transient RR: vio -> dep begin
  P$mRR.dep_o.vio_nt = 1     # (max RR)    cumulative RR: vio -> dep begin
  P$nsc.dep_o.vio_nt = 10    # (n scale)   cumulative RR: vio -> dep begin
  # RR: * -> dep end
  P$dsc.dep_x.dep_u  = Inf   # (dur scale)  duration RR: dep dur -> dep end
  P$iRR.dep_x.vio_zr = 1     # (initial RR) transient RR: vio -> dep end
  P$tsc.dep_x.vio_zr = 30    # (time scale) transient RR: vio -> dep end
  # RR: * -> haz begin
  P$ RR.haz_o.haz_p  = 1     # RR: haz past -> haz begin
  P$ RR.haz_o.dep_w  = 1     # RR: dep now -> haz begin
  P$iRR.haz_o.vio_zr = 1     # (initial RR) transient RR: vio -> haz begin
  P$tsc.haz_o.vio_zr = 30    # (time scale) transient RR: vio -> haz begin
  P$mRR.haz_o.vio_nt = 1     # (max RR)    cumulative RR: vio -> haz begin
  P$nsc.haz_o.vio_nt = 10    # (n scale)   cumulative RR: vio -> haz begin
  # RR: * -> haz end
  P$ RR.haz_x.dep_w  = 1     # RR: dep now -> haz end
  P$dsc.haz_x.haz_u  = Inf   # (dur scale)  duration RR: haz dur -> haz end
  P$iRR.haz_x.vio_zr = 1     # (initial RR) transient RR: vio -> haz end
  P$tsc.haz_x.vio_zr = 30    # (time scale) transient RR: vio -> haz end
  # RR: * -> ptr begin
  P$ RR.ptr_o.dep_w  = 1     # RR: dep now -> ptr begin
  P$ RR.ptr_o.haz_w  = 1     # RR: haz now -> ptr begin
  P$iRR.ptr_o.vio_zr = 1     # (initial RR) transient RR: vio -> ptr begin
  P$tsc.ptr_o.vio_zr = 30    # (time scale) transient RR: vio -> ptr begin
  P$mRR.ptr_o.vio_nt = 1     # (max RR)    cumulative RR: vio -> ptr begin
  P$nsc.ptr_o.vio_nt = 10    # (n scale)   cumulative RR: vio -> ptr begin
  # RR: * -> ptr end
  P$ RR.ptr_x.dep_w  = 1     # RR: dep now -> ptr end
  P$ RR.ptr_x.haz_w  = 1     # RR: haz now -> ptr end
  return(P)
}

add.pars.time = function(P,dtz){
  P$dtz  = dtz              # days in 1 timestep
  P$z3m  = round(365/dtz/4) # timesteps in 3 months
  P$z6m  = 2 * P$z3m        # timesteps in 6 months
  P$z1y  = 4 * P$z3m        # timesteps in 1 year
  P$t1y  = dtz * P$z1y      # days in 1 year
  P$t3m  = dtz * P$z3m      # days in 3 months
  P$t6m  = dtz * P$z6m      # days in 6 months
  return(P)
}

add.pars.cond = function(P){
  P$zf    = P$n.dur*adur*P$z1y # final timestep
  P$tf    = P$zf * P$dtz       # final time (days)
  P$n.tot = P$n.pop * (1+P$n.dur) # total inds needed
  P$het   = het.funs[[P$het.distr]] # funs for sampling ind rates
  P$sex.Ri.shapes = fit.beta(P$sex.Ri.95) # (shape1,shape2): sex rate
  P$cdm.Pi.shapes = fit.beta(P$cdm.Pi.95) # (shape1,shape2): condom prob
  # RR: age
  P$aRR.vio   = def.RR.age(P$aRR.vio.ages,  P$aRR.vio.RRs,  P$aRR.shape) # RR: age -> vio
  P$aRR.dep_o = def.RR.age(P$aRR.dep_o.ages,P$aRR.dep_o.RRs,P$aRR.shape) # RR: age -> dep begin
  P$aRR.haz_o = def.RR.age(P$aRR.haz_o.ages,P$aRR.haz_o.RRs,P$aRR.shape) # RR: age -> haz begin
  P$aRR.ptr_o = def.RR.age(P$aRR.ptr_o.ages,P$aRR.ptr_o.RRs,P$aRR.shape) # RR: age -> ptr begin
  # tRR: vio
  P$tRRu.vio.vio_zr   = def.tRR(P$tRR.shape,P$iRR.vio.vio_zr,  P$tsc.vio.vio_zr,  P$dtz) - 1 # tRR-1: vio -> vio
  P$tRRu.dep_o.vio_zr = def.tRR(P$tRR.shape,P$iRR.dep_o.vio_zr,P$tsc.dep_o.vio_zr,P$dtz) - 1 # tRR-1: vio -> dep begin
  P$tRRu.dep_x.vio_zr = def.tRR(P$tRR.shape,P$iRR.dep_x.vio_zr,P$tsc.dep_x.vio_zr,P$dtz) - 1 # tRR-1: vio -> dep end
  P$tRRu.haz_o.vio_zr = def.tRR(P$tRR.shape,P$iRR.haz_o.vio_zr,P$tsc.haz_o.vio_zr,P$dtz) - 1 # tRR-1: vio -> haz begin
  P$tRRu.haz_x.vio_zr = def.tRR(P$tRR.shape,P$iRR.haz_x.vio_zr,P$tsc.haz_x.vio_zr,P$dtz) - 1 # tRR-1: vio -> haz end
  P$tRRu.ptr_o.vio_zr = def.tRR(P$tRR.shape,P$iRR.ptr_o.vio_zr,P$tsc.ptr_o.vio_zr,P$dtz) - 1 # tRR-1: vio -> ptr begin
  # nRR: vio
  P$nRR.vio.vio_nt   = def.nRR(P$nRR.shape,P$mRR.vio.vio_nt,  P$nsc.vio.vio_nt,  P$t1y) # nRR: vio -> vio
  P$nRR.dep_o.vio_nt = def.nRR(P$nRR.shape,P$mRR.dep_o.vio_nt,P$nsc.dep_o.vio_nt,P$t1y) # nRR: vio -> dep begin
  P$nRR.haz_o.vio_nt = def.nRR(P$nRR.shape,P$mRR.haz_o.vio_nt,P$nsc.haz_o.vio_nt,P$t1y) # nRR: vio -> haz begin
  P$nRR.ptr_o.vio_nt = def.nRR(P$nRR.shape,P$mRR.ptr_o.vio_nt,P$nsc.ptr_o.vio_nt,P$t1y) # nRR: vio -> ptr begin
  # dRR: durs
  P$dRRu.dep_x.dep_u = def.dRR(P$dRR.shape,P$dsc.dep_x.dep_u,P$dtz,P$t1y) - 1 # dRR-1: dep dur -> dep end
  P$dRRu.haz_x.haz_u = def.dRR(P$dRR.shape,P$dsc.haz_x.haz_u,P$dtz,P$t1y) - 1 # dRR-1: dep dur -> dep end
  # pre-compute Ri.m (per day) from Ri.my (per year)
  for (x in filter.names(P,'Ri\\.my$')){
    P[[gsub('my$','m',x)]] = P[[x]] / P$t1y
  }
  # pre-compute RR-1 for all RR.*
  for (x in filter.names(P,'^RR')){
    P[[gsub('RR','RRu',x)]] = P[[x]] - 1
  }
  # for (x in filter.names(P,'^(t|n|d)RR.(?!shape)')){ plot(grepl('RRu',x)+P[[x]]); title(x) } # DEBUG
  return(P)
}

null.pars = function(P,null,save){
  P$null = null    # store config
  P$save = save    # store config
  P.save = P[save] # save exempt
  map = flist(null.sets[null]) # merge regex list
  for (re in names(map)){ # for each regex
    for (x in filter.names(P,re)){ # for each matching par
      P[[x]] = map[[re]] }} # overwrite
  P = ulist(P,P.save) # restore saved
}

null.sets = list(
  Ri  = list('Ri\\.m'=0),
  aRR = list('^aRR\\..*\\.(ages|RRs)$'=1),
   RR = list('^RR\\.'=1),
  tRR = list('^iRR\\.'=1,'^tsc\\.'=1e-12),
  nRR = list('^mRR\\.'=1,'^nsc\\.'=Inf),
  dRR = list('^dsc\\.'=Inf))
null.sets$all = flist(null.sets)      # no events
null.sets$xRR = flist(null.sets[2:6]) # only Ri
null.sets$eRR = flist(null.sets[3:6]) # only Ri & aRR

vec.pars = names(which(lens(add.pars.def()) > 1))

get.run.par = function(v,u=FALSE){
  # get run list for vars (v) + upstream
  vars = c('vio','dep','haz','ptr','sex')
  run = set.names(vars %in% v,vars)
  if (u){ run[1:max(which(run))] = TRUE }
  run = as.list(run)
}

# -----------------------------------------------------------------------------

get.pars.grid = function(pars=list(),...,seed=1:7,.par=TRUE,.grid=TRUE){
  # get Ps (list of P) for all combos of pars & ... & seed
  pars = ulist(pars,...) # merge / overwrite
  v  = lens(pars) > 1 & ! names(pars) %in% vec.pars # pars that vary
  pa = list(seed=seed,var=pars[v],fix=pars[!v],n.var=ifelse(.grid,prod,max)(lens(pars[v])))
  status(3,'get.pars: ',pa$n.var,' x ',len(seed))
  P0s = grid.apply(ulist(pa$var,seed=0),get.pars,args=pa$fix,.par=.par,.grid=.grid) # dummy seed
  Ps  = flist(lapply(P0s,function(P0){ # replicate P0s for each seed (faster)
    lapply(seed,function(s){ P0$seed = s; return(P0) }) }))
  attributes(Ps) = pa
  return(Ps)
}

# =============================================================================
# heterogeneity funs

het.funs = list(
  gamma = list( # m = mean; het = CoV (sd / mean)
    r = function(n,m,het){ cv2 = max(het^2,1e-6); x = rgamma(n,shape=1/cv2,scale=m*cv2) },
    q = function(p,m,het){ cv2 = max(het^2,1e-6); x = qgamma(p,shape=1/cv2,scale=m*cv2) }),
  weibull = list(
    r = function(n,m,het){ f = fit.weibull(m,het^2); x = rweibull(n,shape=f$shape,scale=f$scale) },
    q = function(p,m,het){ f = fit.weibull(m,het^2); x = qweibull(p,shape=f$shape,scale=f$scale) }),
  lnorm = list(
    r = function(n,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); x = rlnorm(n,meanlog=u,sdlog=s) },
    q = function(p,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); x = qlnorm(p,meanlog=u,sdlog=s) }),
  R2 = list( # m = mean; het = xR; p0 = 0.5 (fixed)
    r = function(n,m,het){ x0 = 2*m/(1+het); x = rR2(n,x0=x0,xR=het,p0=.5) },
    q = function(p,m,het){ x0 = 2*m/(1+het); x = qR2(p,x0=x0,xR=het,p0=.5) }))

# =============================================================================
# effect funs

def.RR.age = function(age,RR,shape='spline',eps=.001){
  n = len(RR)
  RR  = c(RR[1],RR,RR[n])
  age = c(age[1]-eps,age,age[n]+eps)
  age.out = seq(1,amax)
  RR.age = switch(shape,
    spline = splinefun(age,RR,method='monoH.FC')(age.out),
    linear = approx(age,RR,age.out,method='liner',rule=2)$y,
    const  = approx(age,RR,age.out,method='const',rule=2)$y)
}

# -----------------------------------------------------------------------------

def.nRR = function(shape,mRR,nsc,t1y){
  # cumulative RR: choose shape function for count: no timestep issues
  if (nsc==Inf){ return(rep(1,t1y*adur)) }
  n = 0:(t1y*adur) # nmax = all active timesteps
  nRR = 1 + (mRR-1) * switch(shape,
    exp  = 1 - exp(-n/nsc),
    ramp = pmin(1, n/nsc),
    step = n >= nsc)
}

def.tRR = function(shape,iRR,tsc,dtz){
  # transient RR: choose shape function & tmax, then integrate & adjust
  tRR.t = switch(shape,
    exp  = function(t){ 1 + (t>=0)*(iRR-1) * exp(-t/tsc) },
    ramp = function(t){ 1 + (t>=0)*(iRR-1) * pmax(0,1-t/tsc) },
    step = function(t){ 1 + (t>=0)*(iRR-1) * (t <= tsc) })
  tmax = switch(shape,exp=8*tsc,ramp=tsc,step=tsc) + dtz
  tRR = int.tRR(tRR.t,dtz,tmax)/dtz
}

def.dRR = function(shape,dsc,dtz,t1y){
  # duration RR: choose shape function & dmax, then integrate & adjust
  dRR.d = switch(shape,
    exp  = function(d){ (d<0) + (d>=0) * exp(-d/dsc) },
    ramp = function(d){ (d<0) + (d>=0) * pmax(0,1-d/dsc) },
    step = function(d){ (d<0) + (d>=0) * (d <= dsc) })
  dmax = t1y*adur # dmax = all active timesteps
  dRR = int.tRR(dRR.d,dtz,dmax)/dtz
}

int.tRR = function(tRR.t,dtz,tmax){
  # integrates tRR.t on intervals {0,dtz,...,tmax} - dtz/2
  # note: this approximates double integral via interval midpoints
  tz = seq(0,tmax,dtz) - dtz/2 # lower bounds
  tRR = sapply(tz,function(ti){ integrate(tRR.t,ti,ti+dtz)$value })
}

# -----------------------------------------------------------------------------

map.nRR = function(nRR,n,ze,z){
  # lookup nRR kernel for the number (n) of exposures
  # note: if z = ze, we average nRR[n] & nRR[n-1]
  #       since same-timestep exposure has half effect
  RR = 0.5 * (nRR[n+1] + nRR[n+na.to.num(z>ze,1)])
}

map.tRR = function(tRRu,ze,z){
  # lookup tRR kernel for now (z) given most recent event (ze)
  # note: if [z-ze+1] is out-of-bounds: NA -> 0 so RR = 1 + 0
  # note: we add +1 since tRRu[1] reflects same-timestep
  RR = 1 + na.to.num(tRRu[z-ze+1])
}
