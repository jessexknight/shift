source('meta.r')
source('utils.r')

# =============================================================================
# effect funs

fit.rr = function(.,t,flat=c(1,1)){
  # .$t: RR time control point "x" (days); .$rr: RR control point "y"
  # .$rrt: (optional) total RR-days -> re-scales final rr to match
  # t: time (days) to evaluate the RR spline at
  if (missing(t)){ t = seq(dtz,max(.$t),dtz) }
  n = len(.$t)+flat[1]
  if (flat[1]){ .$t = c(.$t[1]-eps,.$t); .$rr = c(.$rr[1],.$rr) }
  if (flat[2]){ .$t = c(.$t,.$t[n]+eps); .$rr = c(.$rr,.$rr[n]) }
  rr = splinefun(.$t,.$rr,method='monoH.FC')(t)
  if (!is.null(.$rrt)){ rr = 1 + (rr-1) * (.$rrt-1) / (dtz*sum(rr-1)) }
  list(t=t,rr=rr,rrt=1+dtz*sum(rr-1))
}

fit.rr.age = function(.){
  # convenience function
  rr.age = fit.rr(.,seq(amin,amax),flat=c(1,1))
}

get.rr.evt = function(ze,z,urr.z){
  # lookup & sum RR kernel for today (z) given prior events (ze)
  rr = 1 + sum(urr.z[z+1-ze],na.rm=TRUE)
}

# =============================================================================
# initialization funs

init.evts = function(Is){
  # initialize event vectors for each event type & individual
  E = lapply(Is$i,function(i){ numeric() })
  Es = lapply(evts,function(e){ E })
}

init.inds = function(P){
  # initialize all individuals needed for timesteps 1:zf
  nd = P$zf/z1y/adur # num sim adurs
  n  = P$n * (1+nd)  # total inds needed
  Is = data.frame(
    i = seq(n),
    age      = runif(n,min=amin-nd*adur,max=amax),
    age.act  = runif(n,min=amin,max=20),
    # partnerships
    ptr.n    = 0,
    ptr.o.r0 = rexp(n=n,rate=1/P$ptr.o.r0.m),
    ptr.max  = pmin(3,1+rgeom(n=n,prob=1/P$ptr.max.m)),
    cdm.p0   = runif(n=n,min=0,max=1),
    # structural factors & mediators
    vio.e.r0 = rexp(n=n,rate=1/P$vio.e.r0.m),
    dep.o.r0 = rexp(n=n,rate=1/P$dep.o.r0.m),
    dep.x.r0 = rexp(n=n,rate=1/P$dep.x.r0.m),
    dep.now  = FALSE,
    dep.past = FALSE,
    dep.z0   = -Inf,
    alc.o.r0 = rexp(n=n,rate=1/P$alc.o.r0.m),
    alc.x.r0 = rexp(n=n,rate=1/P$alc.x.r0.m),
    alc.now  = FALSE,
    alc.past = FALSE,
    alc.z0   = -Inf
  )
}

init.ptrs = function(P,Is,z){
  # generate partners among Is
  # mixing is "random" by splitting into 1st/2nd half
  n = nrow(Is)/2
  if (n==0){ return(NULL) }
  i1 = 0+(1:n)
  i2 = n+(1:n)
  Ks = data.frame(
    i1 = Is$i[i1], # i of individual 1
    i2 = Is$i[i2], # i of individual 2
    z0 = z, # timestep ptr starts
    zf = z + rweibull(n,shape=.5,scale=.5*P$ptr.dt.m/dtz), # timestep ptr ends
    f.sex = runif(n,.1,.5), # freq of sex in ptr
    cdm = Is$cdm.p0[i1]/2 + Is$cdm.p0[i2]/2 # condom prob (average)
  )
}

# =============================================================================
# rate & prob funs

rate.vio.e = function(P,Js,aj){
  R = ( # among all
      Js$vio.e.r0 # base rate
    * P$rr.vio.e.age[aj] # RR age
); return(R) }

rate.dep.o = function(P,Js,R,aj,z,e.vio){
  j = which(!Js$dep.now)
  R[j] = ( # among not dep
      Js$dep.o.r0[j] # base rate
    * P$rr.dep.o.age[aj[j]] # RR age
    * (1 + P$urr.dep.o.dep.p * Js$dep.past[j]) # RR dep past
    * vapply(e.vio[j],get.rr.evt,0,z,P$urr.dep.o.vio.z) # RR vio recent
); return(R) }

rate.dep.x = function(P,Js,R,z,e.vio){
  j = which(Js$dep.now)
  R[j] = ( # among dep
      Js$dep.x.r0[j] # base rate
    * 2^(dtz * (z - Js$dep.z0[j]) / P$rr.dep.x.th) # RR dep dur
    * vapply(e.vio[j],get.rr.evt,0,z,P$urr.dep.x.vio.z) # RR vio recent
); return(R) }

rate.alc.o = function(P,Js,R,aj,z,e.vio){
  j = which(!Js$alc.now)
  R[j] = ( # among not alc
      Js$alc.o.r0[j] # base rate
    * P$rr.alc.o.age[aj[j]] # RR age
    * (1 + P$urr.alc.o.alc.p * Js$alc.past[j]) # RR alc past
    * vapply(e.vio[j],get.rr.evt,0,z,P$urr.alc.o.vio.z) # RR vio recent
    * (1 + P$urr.alc.o.dep.n * Js$dep.now[j]) # RR dep now
); return(R) }

rate.alc.x = function(P,Js,R,z,e.vio){
  j = which(Js$alc.now)
  R[j] = ( # among alc
      Js$alc.x.r0[j] # base rate
    * 2^(dtz * (z - Js$alc.z0[j]) / P$rr.alc.x.th) # RR alc dur
    * vapply(e.vio[j],get.rr.evt,0,z,P$urr.alc.x.vio.z) # RR vio recent
    * (1 + P$urr.alc.x.dep.n * Js$dep.now[j]) # RR dep now
); return(R) }

rate.ptr = function(P,Js,R,aj){
  j = which(Js$ptr.n < Js$ptr.max)
  R[j] = ( # among avail
      Js$ptr.o.r0[j] # base rate
    * P$rr.ptr.o.age[aj[j]] # RR age
    * (1 + P$urr.ptr.o.dep.n * Js$dep.now[j]) # RR dep now
    * (1 + P$urr.ptr.o.alc.n * Js$alc.now[j]) # RR alc now
); return(R) }

prob.cdm = function(P,Ks,k,Is){
  R = ( # among all
      Ks$cdm[k] # base prob
    * P$rr.cdm.b.dep.n ^ (Is$dep.now[Ks$i1[k]] + Is$dep.now[Ks$i2[k]]) # RP dep
    * P$rr.cdm.b.alc.n ^ (Is$alc.now[Ks$i1[k]] + Is$alc.now[Ks$i2[k]]) # RP alc
); return(R) }

# =============================================================================
# run simulation

sim.run = function(P){
  # initialization ------------------------------------------------------------
  set.seed(P$seed)
  Is = init.inds(P)  # individuals
  Es = init.evts(Is) # events
  Ks = NULL          # partnerships
  for (z in 1:P$zf){
    # if (z %% z1y == 0) { print(z/z1y) } # DEBUG
    # age inds ----------------------------------------------------------------
    Is$age = Is$age + 1/z1y
    # break ptrs --------------------------------------------------------------
    b = Ks$zf < z
    i = c(Ks$i1[b],Ks$i2[b])
    Is$ptr.n[i] = Is$ptr.n[i] - 1
    Ks = Ks[!b,]
    # select active inds ------------------------------------------------------
    i = which(Is$age > Is$age.act & Is$age < amax)
    Js = Is[i,]       # read only copy of active
    ij = match(Js$i,Is$i) # map j -> j
    aj = floor(Js$age-amin+1) # age vector for j
    R0 = numeric(nrow(Js)) # init rate = 0 for j
    e.vio = Es$vio.e[i]
    # update vio events -------------------------------------------------------
    i = ij[which(runif(ij) < dtz * rate.vio.e(P,Js,aj))]
    Es$vio.e[i] = lapply(Es$vio.e[i],append,z)
    # update dep onset --------------------------------------------------------
    i = ij[which(runif(ij) < dtz * rate.dep.o(P,Js,R0,aj,z,e.vio))]
    Is$dep.now[i] = TRUE
    Is$dep.past[i] = TRUE
    Is$dep.z0[i] = z
    Es$dep.o[i] = lapply(Es$dep.o[i],append,z)
    # update dep recovery -----------------------------------------------------
    i = ij[which(runif(ij) < dtz * rate.dep.x(P,Js,R0,z,e.vio))]
    Is$dep.now[i] = FALSE
    Es$dep.x[i] = lapply(Es$dep.x[i],append,z)
    # update alc onset --------------------------------------------------------
    i = ij[which(runif(ij) < dtz * rate.alc.o(P,Js,R0,aj,z,e.vio))]
    Is$alc.now[i] = TRUE
    Is$alc.past[i] = TRUE
    Is$alc.z0[i] = z
    # update alc recovery -----------------------------------------------------
    i = ij[which(runif(ij) < dtz * rate.alc.x(P,Js,R0,z,e.vio))]
    Is$alc.now[i] = FALSE
    Es$alc.x[i] = lapply(Es$alc.x[i],append,z)
    # form ptrs ---------------------------------------------------------------
    i = ij[even.len(which(runif(ij) < dtz * rate.ptr(P,Js,R0,aj)))]
    Is$ptr.n[i] = Is$ptr.n[i] + 1
    Es$ptr.o[i] = lapply(Es$ptr.o[i],append,z)
    Ks = rbind(Ks,init.ptrs(P,Is[i,],z))
    # sex in ptrs -------------------------------------------------------------
    # TODO: sex should be possible > 1 per week
    # TODO: 1-exp(-dt*R) elsewhere
    k = which(runif(nrow(Ks)) < dtz * Ks$f.sex)
    cdm.b = runif(k) < prob.cdm(P,Ks,k,Is)
    i = c(Ks$i1[k],Ks$i2[k])
    Es$sex.e[i] = lapply(Es$sex.e[i],append,z)
    Es$cdm.b[i] = mapply(append,Es$cdm.b[i],cdm.b)
  }
  # lapply(Es,function(E){ e = sapply(E,len); hist(e,0:max(e)) }) # DEBUG
  Is = sim.out(Is,Es,P)
}

sim.out = function(Is,Es,P,rm.dum=TRUE){
  # clean-up simulation output
  if (rm.dum){ # remove initial dummy population
    i = which(Is$age < (P$zf/z1y+amin))
    Is = Is[i,]
    Es = lapply(Es,`[`,i)
  }
  # compute some extra variables
  Is$age.1 = floor(Is$age)     # 1-year age bins
  Is$age.5 = floor(Is$age/5)*5 # 5-year age bins
  Is$vio.tot = sapply(Es$vio.e,len) # lifetime vio
  Is$ptr.tot = sapply(Es$ptr.o,len) # lifetime ptrs
  Is$sex.tot = sapply(Es$sex.e,len) # lifetime sex
  Is$cdm.ls  = sapply(Es$cdm.b,last) # last sex cdm
  Is$cdm.m   = sapply(Es$cdm.b,mean) # lifetime cdm mean
  Is = cbind(seed=P$seed,Is)
}

sim.runs = function(Ps,.par=TRUE){
  # run.sim in parallel for each (P)arameter set in Ps
  Is = rbind.lapply(Ps,sim.run,.par=.par)
}

# =============================================================================
# params

get.pars = function(seed=0,...){
  # all rates in 1/days, durations in days
  P = list(seed=seed)
  P$n = 100
  P$zf = z1y*adur*2     # total timesteps
  P$ptr.o.r0.m     =  .01  # base rate of partner seeking (mean)
  P$ptr.max.m      = 1.50  # max num partners (mean)
  P$vio.e.r0.m     =  .001 # base rate of violence event (mean)
  P$dep.o.r0.m     =  .0001# base rate of depression onset (mean)
  P$dep.x.r0.m     =  .001 # base rate of depression recovery (mean)
  P$alc.o.r0.m     =  .0001# base rate of alcohol onset (mean)
  P$alc.x.r0.m     =  .001 # base rate of alcohol recovery (mean)
  P$ptr.dt.m       = 364   # duration of partnerships (mean)
  P$rr.dep.x.th    = 364   # half-life of RR for depression tunnel
  P$rr.alc.x.th    = 364   # half-life of RR for alcohol tunnel
  P$rr.dep.o.dep.p = 3.0   # RR of depression onset if depressed in past
  P$rr.alc.o.alc.p = 3.0   # RR of alcohol onset if alcohol in past
  P$rr.alc.o.dep.n = 2.0   # RR of alcohol onset if depressed now
  P$rr.alc.x.dep.n = 0.5   # RR of alcohol recovery if depressed now
  P$rr.ptr.o.dep.n = 1.0   # RR of partner seeking if depressed now
  P$rr.ptr.o.alc.n = 1.5   # RR of partner seeking if alcohol now
  P$rr.cdm.b.dep.n = 1.0   # RR of condom use if depressed now
  P$rr.cdm.b.alc.n =  .50  # RR of condom use if alcohol now
  P = list.update(P,...)
  P = add.pars(P)
}

add.pars = function(P){
  # define RR splines
  rr.age. = list( # RR age
    vio.e = list(t=c(15,20,30,60),rr=c(0.0,1.0,1.0,0.5)),
    dep.o = list(t=c(15,20,30,60),rr=c(0.0,1.0,1.0,0.0)),
    alc.o = list(t=c(15,20,30,60),rr=c(0.0,1.0,1.0,1.0)),
    ptr.o = list(t=c(15,20,30,60),rr=c(1.0,1.0,1.0,0.0)))
  # plot.rr(rr.age.,lapply(rr.age.,fit.rr.age)) + xlab('Age') + ylim(c(0,1))
  #   plot.save('par','rr.age',h=2.5,w=5)
  rr.vio. = list( # RR vio
    dep.o = list(t=14*(0:4),rr=1+1.0*c(1.0,0.95,0.5,0.05,0.0)),
    dep.x = list(t=14*(0:4),rr=1-0.5*c(1.0,0.95,0.5,0.05,0.0)),
    alc.o = list(t=14*(0:4),rr=1+1.0*c(1.0,0.95,0.5,0.05,0.0)),
    alc.x = list(t=14*(0:4),rr=1-0.5*c(1.0,0.95,0.5,0.05,0.0)))
  # plot.rr(rr.vio.,lapply(rr.vio.,fit.rr)) + scale_y_continuous(trans='log2')
  #   plot.save('par','rr.vio',h=2.5,w=5)
  P$rr.vio.e.age    = fit.rr.age(rr.age.$vio.e)$rr # RR of violence event by age
  P$rr.dep.o.age    = fit.rr.age(rr.age.$dep.o)$rr # RR of depression onset by age
  P$rr.alc.o.age    = fit.rr.age(rr.age.$alc.o)$rr # RR of alcohol onset by age
  P$rr.ptr.o.age    = fit.rr.age(rr.age.$ptr.o)$rr # RR of partner seeking by age
  P$urr.dep.o.vio.z = fit.rr(rr.vio.$dep.o)$rr - 1 # RR-1 of dep onset per violence event
  P$urr.dep.x.vio.z = fit.rr(rr.vio.$dep.x)$rr - 1 # RR-1 of dep recov per violence event
  P$urr.alc.o.vio.z = fit.rr(rr.vio.$alc.o)$rr - 1 # RR-1 of alc onset per violence event
  P$urr.alc.x.vio.z = fit.rr(rr.vio.$alc.x)$rr - 1 # RR-1 of alc recov per violence event
  P$urr.dep.o.dep.p = P$rr.dep.o.dep.p - 1         # RR-1 of depression onset if depressed in past
  P$urr.alc.o.alc.p = P$rr.alc.o.alc.p - 1         # RR-1 of alcohol onset if alcohol in past
  P$urr.alc.o.dep.n = P$rr.alc.o.dep.n - 1         # RR-1 of alcohol onset if depressed now
  P$urr.alc.x.dep.n = P$rr.alc.x.dep.n - 1         # RR-1 of alcohol recovery if depressed now
  P$urr.ptr.o.dep.n = P$rr.ptr.o.dep.n - 1         # RR-1 of partner seeking if depressed now
  P$urr.ptr.o.alc.n = P$rr.ptr.o.alc.n - 1         # RR-1 of partner seeking if alcohol now
  # why pre-compute RR-1: "sum_z RR_z" should be computed as "1 + sum_z (RR_z - 1)"
  # for RR_z < 1: this could yield RR_total < 0 (!) but this is rare & may be correct
  return(P)
}
