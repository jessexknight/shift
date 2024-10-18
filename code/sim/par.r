
# =============================================================================
# pars

get.pars = function(seed=0,...,dtz=7,case='base',null=NULL,save=NULL){
  P = list(case=case,seed=seed,id=sprintf('%6d',seed))
  P = dt.pars(P,dtz)
  P$n.pop = 1000
  P$n.dur = 1+1
  P$null = null
  # base rates
  P$vio.Ri.m    = 1/P$t1y    # (mean) base rate: violence
  P$dep_o.Ri.m  = 0.01/P$t1y # (mean) base rate: depression begin
  P$dep_x.Ri.m  = 1.00/P$t1y # (mean) base rate: depression end
  P$haz_o.Ri.m  = 0.01/P$t1y # (mean) base rate: hazdrink begin
  P$haz_x.Ri.m  = 1.00/P$t1y # (mean) base rate: hazdrink end
  P$ptr_o.Ri.m  = 1/30       # (mean) base rate: partner begin
  P$ptr_x.Ri.m  = 1/30       # (mean) base rate: partner end
  P$sex.Ri.95   = c(.1,.5)   # (95% CI) base rate: sex within ptr
  P$cdm.Pi.95   = c(.2,.8)   # (95% CI) prob: condom use
  # base rate covariance, shapes, etc.
  P$ptr.max.m   = 1.50      # (mean) max num partners
  P$dep.cov     = -.9       # approx covariance among dep_o,dep_x
  P$haz.cov     = -.9       # approx covariance among haz_o,haz_x
  P$ptr.cov     = +.9       # approx covariance among ptr_o,ptr_x,ptr.max
  P$ptr.Ri.shape = 3        # (gamma shape): ptr_o,ptr_x
  P$all.Ri.shape = 1        # (gamma shape): all other base rates
  # *RR shapes
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
  P$aRR.ptr_o.ages = c(amin,amax) # (age points) RR: age -> ptr begin
  P$aRR.ptr_o.RRs  = c(1.00,1.00) # (RR  points) RR: age -> ptr begin
  # RR: * -> dep begin
  P$ RR.dep_o.dep_p  = 3     # RR: dep past -> dep begin
  P$iRR.dep_o.vio_zr = 2     # (initial RR) transient RR: vio -> dep begin
  P$tsc.dep_o.vio_zr = 30    # (time scale) transient RR: vio -> dep begin
  P$mRR.dep_o.vio_nt = 2     # (max RR)    cumulative RR: vio -> dep begin
  P$nsc.dep_o.vio_nt = 10    # (n scale)   cumulative RR: vio -> dep begin
  # RR: * -> dep end
  P$dsc.dep_x.dep_u  = P$t1y # (dur scale)  duration RR: dep dur -> dep end
  P$iRR.dep_x.vio_zr = 1/2   # (initial RR) transient RR: vio -> dep end
  P$tsc.dep_x.vio_zr = 30    # (time scale) transient RR: vio -> dep end
  # RR: * -> haz begin
  P$ RR.haz_o.haz_p  = 3     # RR: haz past -> haz begin
  P$ RR.haz_o.dep_w  = 3     # RR: dep now -> haz begin
  P$iRR.haz_o.vio_zr = 2     # (initial RR) transient RR: vio -> haz begin
  P$tsc.haz_o.vio_zr = 30    # (time scale) transient RR: vio -> haz begin
  P$mRR.haz_o.vio_nt = 2     # (max RR)    cumulative RR: vio -> haz begin
  P$nsc.haz_o.vio_nt = 10    # (n scale)   cumulative RR: vio -> haz begin
  # RR: * -> haz end
  P$ RR.haz_x.dep_w  = 1/3   # RR: dep now -> haz end
  P$dsc.haz_x.haz_u  = P$t1y # (dur scale)  duration RR: haz dur -> haz end
  P$iRR.haz_x.vio_zr = 1/2   # (initial RR) transient RR: vio -> haz end
  P$tsc.haz_x.vio_zr = 30    # (time scale) transient RR: vio -> haz end
  # RR: * -> ptr begin
  P$ RR.ptr_o.dep_w  = 1/2   # RR: dep now -> ptr begin
  P$ RR.ptr_o.haz_w  = 2     # RR: haz now -> ptr begin
  P$iRR.ptr_o.vio_zr = 1/2   # (initial RR) transient RR: vio -> ptr begin
  P$tsc.ptr_o.vio_zr = 30    # (time scale) transient RR: vio -> ptr begin
  P$mRR.ptr_o.vio_nt = 1/2   # (max RR)    cumulative RR: vio -> ptr begin
  P$nsc.ptr_o.vio_nt = 10    # (n scale)   cumulative RR: vio -> ptr begin
  # RR: * -> ptr end
  P$ RR.ptr_x.dep_w  = 2     # RR: dep now -> ptr end
  P$ RR.ptr_x.haz_w  = 2     # RR: haz now -> ptr end
  # overwrite & add conditional
  P = null.pars(P,null,save)
  P = ulist(P,...)
  P = cond.pars(P)
}

dt.pars = function(P,dtz){
  P$dtz  = dtz              # days in 1 timestep
  P$z3m  = round(364/dtz/4) # timesteps in 3 months
  P$z6m  = 2 * P$z3m        # timesteps in 6 months
  P$z1y  = 4 * P$z3m        # timesteps in 1 year
  P$t1y  = dtz * P$z1y      # days in 1 year
  P$t3m  = dtz * P$z3m      # days in 3 months
  P$t6m  = dtz * P$z6m      # days in 6 months
  return(P)
}

cond.pars = function(P){
  P$zf    = P$n.dur*adur*P$z1y # final timestep
  P$tf    = P$zf * P$dtz       # final time (days)
  P$n.tot = P$n.pop * (1+P$n.dur) # total inds needed
  P$sex.Ri.shapes = fit.beta(P$sex.Ri.95) # (shape1,shape2): sex rate
  P$cdm.Pi.shapes = fit.beta(P$cdm.Pi.95) # (shape1,shape2): condom prob
  # RR: age
  P$aRR.vio   = def.RR.age(P$aRR.vio.ages,  P$aRR.vio.RRs,  P$aRR.shape) # RR: age -> vio
  P$aRR.dep_o = def.RR.age(P$aRR.dep_o.ages,P$aRR.dep_o.RRs,P$aRR.shape) # RR: age -> dep begin
  P$aRR.haz_o = def.RR.age(P$aRR.haz_o.ages,P$aRR.haz_o.RRs,P$aRR.shape) # RR: age -> haz begin
  P$aRR.ptr_o = def.RR.age(P$aRR.ptr_o.ages,P$aRR.ptr_o.RRs,P$aRR.shape) # RR: age -> ptr begin
  # tRR: vio
  P$tRRu.dep_o.vio_zr = def.tRR(P$tRR.shape,P$iRR.dep_o.vio_zr,P$tsc.dep_o.vio_zr,P$dtz) - 1 # tRR-1: vio -> dep begin
  P$tRRu.dep_x.vio_zr = def.tRR(P$tRR.shape,P$iRR.dep_x.vio_zr,P$tsc.dep_x.vio_zr,P$dtz) - 1 # tRR-1: vio -> dep end
  P$tRRu.haz_o.vio_zr = def.tRR(P$tRR.shape,P$iRR.haz_o.vio_zr,P$tsc.haz_o.vio_zr,P$dtz) - 1 # tRR-1: vio -> haz begin
  P$tRRu.haz_x.vio_zr = def.tRR(P$tRR.shape,P$iRR.haz_x.vio_zr,P$tsc.haz_x.vio_zr,P$dtz) - 1 # tRR-1: vio -> haz end
  P$tRRu.ptr_o.vio_zr = def.tRR(P$tRR.shape,P$iRR.ptr_o.vio_zr,P$tsc.ptr_o.vio_zr,P$dtz) - 1 # tRR-1: vio -> ptr begin
  # nRR: vio
  P$nRR.dep_o.vio_nt = def.nRR(P$nRR.shape,P$mRR.dep_o.vio_nt,P$nsc.dep_o.vio_nt,P$t1y) # nRR: vio -> dep begin
  P$nRR.haz_o.vio_nt = def.nRR(P$nRR.shape,P$mRR.haz_o.vio_nt,P$nsc.haz_o.vio_nt,P$t1y) # nRR: vio -> haz begin
  P$nRR.ptr_o.vio_nt = def.nRR(P$nRR.shape,P$mRR.ptr_o.vio_nt,P$nsc.ptr_o.vio_nt,P$t1y) # nRR: vio -> ptr begin
  # dRR: durs
  P$dRRu.dep_x.dep_u = def.dRR(P$dRR.shape,P$dsc.dep_x.dep_u,P$dtz,P$t1y) - 1 # dRR-1: dep dur -> dep end
  P$dRRu.haz_x.haz_u = def.dRR(P$dRR.shape,P$dsc.haz_x.haz_u,P$dtz,P$t1y) - 1 # dRR-1: dep dur -> dep end
  # pre-compute RR-1 for all RR.*
  for (x in filter.names(P,'^RR')){
    P[[gsub('RR','RRu',x)]] = P[[x]] - 1
  }
  # for (x in filter.names(P,'^(t|n|d)RR.(?!shape)')){ plot(grepl('RRu',x)+P[[x]]); title(x) } # DEBUG
  return(P)
}

null.pars = function(P,null,save){
  P.save = P[save] # save exempt
  map = flist(null.sets[null]) # merge regex list
  for (re in names(map)){ # for each regex
    for (x in filter.names(P,re)){ # for each matching par
      P[[x]] = map[[re]] }} # overwrite
  P = ulist(P,P.save) # restore saved
}

null.sets = list(
  Ri  = list('Ri\\.m$'=0),
  RR  = list('^RR\\.'=1),
  aRR = list('^aRR\\..*\\.(ages|RRs)$'=1),
  tRR = list('^iRR\\.'=1,'^tsc\\.'=1e-12),
  nRR = list('^mRR\\.'=1,'^nsc\\.'=Inf),
  dRR = list('^dsc\\.'=Inf))
null.sets$xRR = flist(null.sets[2:6])
null.sets$all = flist(null.sets)

# =============================================================================
# effect funs

def.RR.age = function(age,RR,shape='spline',eps=.001){
  n = len(RR)
  RR  = c(RR[1],RR,RR[n])
  age = c(age[1]-eps,age,age[n]+eps)
  age.out = seq(amin,amax)
  RR.age = switch(shape,
    spline = splinefun(age,RR,method='monoH.FC')(age.out),
    linear = approx(age,RR,age.out,method='liner',rule=2)$y,
    const  = approx(age,RR,age.out,method='const',rule=2)$y)
}

# -----------------------------------------------------------------------------

def.nRR = function(shape,mRR,nsc,t1y){
  # cumulative RR: chose shape function for count: no timestep issues
  if (nsc==Inf){ return(rep(1,t1y*adur)) }
  n = 0:(t1y*adur) # nmax = all active timesteps
  nRR = 1 + (mRR-1) * switch(shape,
    exp  = 1 - exp(-n/nsc),
    ramp = pmin(1, n/nsc),
    step = n >= nsc)
}

# -----------------------------------------------------------------------------

def.tRR = function(shape,iRR,tsc,dtz){
  # transient RR: chose shape function & tmax, then integrate & adjust
  tRR.t = switch(shape,
    exp  = function(t){ 1 + (t>=0)*(iRR-1) * exp(-t/tsc) },
    ramp = function(t){ 1 + (t>=0)*(iRR-1) * pmax(0,1-t/tsc) },
    step = function(t){ 1 + (t>=0)*(iRR-1) * (t <= tsc) })
  tmax = switch(shape,exp=8*tsc,ramp=tsc,step=tsc) + dtz
  tRR = adj.tRR(int.tRR(tRR.t,dtz,tmax)/dtz)
}

def.dRR = function(shape,dsc,dtz,t1y){
  # duration RR: chose shape function & dmax, then integrate & adjust
  dRR.d = switch(shape,
    exp  = function(d){ (d<0) + (d>=0) * exp(-d/dsc) },
    ramp = function(d){ (d<0) + (d>=0) * pmax(0,1-d/dsc) },
    step = function(d){ (d<0) + (d>=0) * (d <= dsc) })
  dmax = t1y*adur # dmax = all active timesteps
  dRR = adj.tRR(int.tRR(dRR.d,dtz,dmax)/dtz)
}

int.tRR = function(tRR.t,dtz,tmax){
  # integrates tRR.t on intervals {0,dtz,...,tmax} - dtz/2
  # note: this approximates double integral via interval midpoints
  tz = seq(0,tmax,dtz) - dtz/2 # lower bounds
  tRR = sapply(tz,function(ti){ integrate(tRR.t,ti,ti+dtz)$value })
}

adj.tRR = function(tRR){
  # remove tRR[1] and add tRR[1]-1 to tRR[2]
  # why: we don't model same-timestep tRR & dRR effects
  #      but we need to conserve RR-1 mass / AUC
  n = len(tRR)
  if (n==1) return( tRR )
  if (n==2) return( tRR[1]-1+tRR[2] )
  return( c(tRR[1]-1+tRR[2],tRR[3:n]) )
}

map.tRR = function(tRRu,ze,z){
  # lookup tRR kernel for now (z) given most recent event (ze)
  # note: if [z-ze] is out-of-bounds: NA -> 0 so RR = 1 + 0
  # note: z > ze is guaranteed in sim.r via J vs I
  RR = 1 + na.to.num(tRRu[z-ze])
}
