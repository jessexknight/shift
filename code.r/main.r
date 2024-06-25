source('meta.r')
source('utils.r')

# =============================================================================
# effect funs

fit.rr.age = function(a.ref,rr.ref){
  as = seq(amin,amax)
  rr.a = splinefun(a.ref,rr.ref,method='monoH.FC')(as)
}

fit.rr.dz = function(t.ref,q.ref){
  # fit decaying logistic kernel to reach q.ref by t.ref
  err.fun = function(par){
    t.par = qlogis(1-q.ref,loc=par[1],scale=par[2])
    err = sum((t.par-t.ref)^2) }
  par = optim(c(loc=0,scale=1),err.fun,lower=c(NA,1e-3),method='L-BFGS-B')$par
}

get.rr.dz = function(loc,scale,rr.tot,eps=1e-3){
  # pre-compute scaled decaying logistic kernel
  zs = 0:ceiling(qlogis(1-eps,loc=loc,scale=scale)/dtz)
  rr.z = 1-plogis(zs*dtz,loc=loc,scale=scale)
  # plot(0:z1y,1-plogis(0:z1y,loc=loc,scale=scale)); lines(zs*dtz,rr.z) # DEBUG
  rr.z = rr.z * rr.tot / sum(rr.z)
}

get.rr.evt = function(ze,z,rr.dz){
  # lookup & sum RR kernel for today (z) given prior events (ze)
  rr = sum(rr.dz[z+1-ze],na.rm=TRUE)
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
  ny = P$zf/z1y/adur # num sim years
  n  = P$n * (1+ny)        # total inds needed
  Is = data.frame(
    i = seq(n),
    age = runif(n,min=amin-ny*adur,max=amax),
    age.act = runif(n,min=amin,max=20),
    # partnerships
    ptr.n = 0,
    ptr.r0 = rexp(n=n,rate=1/P$ptr.r0.m),
    ptr.max = pmin(3,1+rgeom(n=n,prob=1/P$ptr.max.m)),
    cdm.p0 = runif(n=n,min=0,max=1),
    # structural factors & mediators
    vio.r0 = rexp(n=n,rate=1/P$vio.r0.m),
    dep.o.r0 = rexp(n=n,rate=1/P$dep.o.r0.m),
    dep.x.r0 = rexp(n=n,rate=1/P$dep.x.r0.m),
    dep.now  = FALSE,
    dep.evr  = FALSE,
    dep.z0   = NA
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
    zf = z + rweibull(n,shape=.5,scale=.5*P$ptr.dz.m), # timestep ptr ends
    f.sex = runif(n,.1,.5), # freq of sex in ptr
    cdm = Is$cdm.p0[i1]/2 + Is$cdm.p0[i2]/2 # condom prob (average)
  )
}

# =============================================================================
# run simulation

sim.run = function(P){
  # initialization ------------------------------------------------------------
  set.seed(P$seed)
  Is = init.inds(P)  # individuals
  Es = init.evts(Is) # events
  Ks = NULL          # partnerships
  for (z in 1:P$zf){
    # age inds ----------------------------------------------------------------
    Is$age = Is$age + 1/z1y
    # break ptrs --------------------------------------------------------------
    k.brk = Ks$zf < z
    i.brk = c(Ks$i1[k.brk],Ks$i2[k.brk])
    Is$ptr.n[i.brk] = Is$ptr.n[i.brk] - 1
    Ks = Ks[!k.brk,]
    # select active inds ------------------------------------------------------
    i.act = which(Is$age > Is$age.act & Is$age < amax)
    Js = Is[i.act,]       # read only copy of active
    ij = match(Js$i,Is$i) # map j -> j
    aj = floor(Js$age-amin) # age vector for j
    # update vio --------------------------------------------------------------
    i = ij[which(runif(ij) < Js$vio.r0 * dtz * exp(0
      + P$rr.vio.age[aj]
    ))]
    Es$vio[i] = lapply(Es$vio[i],append,z)
    # update dep --------------------------------------------------------------
    rr.dep.vio = sapply(Es$vio[i.act],get.rr.evt,z,P$eff.vio.dep.dz)
    # dep.o (onset)
    i = ij[which(runif(ij) < (!Js$dep.now) * Js$dep.o.r0 * dtz * exp(0
      + P$rr.dep.age[aj]
      + rr.dep.vio
    ))]
    Is$dep.now[i] = TRUE
    Is$dep.evr[i] = TRUE
    Is$dep.z0[i] = z
    Es$dep.o[i] = lapply(Es$dep.o[i],append,z)
    # dep.x (recovery)
    i = ij[which(runif(ij) < (Js$dep.now) * Js$dep.x.r0 * dtz * exp(0
      - rr.dep.vio
      + P$rr.dep.dur * (z - Js$dep.z0) * dtz
    ))]
    Is$dep.now[i] = FALSE
    Es$dep.x[i] = lapply(Es$dep.x[i],append,z)
    # form ptrs ---------------------------------------------------------------
    i = ij[even.len(which(runif(ij) < (Js$ptr.n < Js$ptr.max) * Js$ptr.r0 * dtz * exp(0
      + P$rr.ptr.age[aj]
      + P$rr.ptr.dep * Js$dep.now
    )))]
    Is$ptr.n[i] = Is$ptr.n[i] + 1
    Es$ptr.o[i] = lapply(Es$ptr.o[i],append,z)
    Ks = rbind(Ks,init.ptrs(P,Is[i,],z))
    # sex in ptrs -------------------------------------------------------------
    Xs = Ks[runif(nrow(Ks)) < Ks$f.sex,]
    cdm = runif(nrow(Xs)) < Xs$cdm * exp(0
      + P$rr.cdm.dep * (Is$dep.now[Xs$i1] + Is$dep.now[Xs$i2])
    )
    i = c(Xs$i1,Xs$i2)
    Es$sex[i] = lapply(Es$sex[i],append,z)
    Es$cdm[i] = mapply(append,Es$cdm[i],cdm)
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
  Is$ptr.tot = sapply(Es$ptr.o,len) # lifetime ptrs
  Is$sex.tot = sapply(Es$sex,len) # lifetime sex
  Is$cdm.ls  = sapply(Es$cdm,last) # last sex cdm
  Is$cdm.m   = sapply(Es$cdm,mean) # lifetime cdm mean
  Is = cbind(seed=P$seed,Is)
}

sim.runs = function(Ps){
  # run.sim in parallel for each (P)arameter set in Ps
  # Is = do.call(rbind,lapply(Ps,sim.run)) # DEBUG
  Is = do.call(rbind,parallel::mclapply(Ps,sim.run,mc.cores=7))
}

# -----------------------------------------------------------------------------
# fit effects

# print(fit.eff.dz(c(30,60),c(.5,.01)))
a.ref = c( 15, 20, 25, 30, 40, 50)
rr.age = list(
  vio = c(1.0,1.0,1.0,1.0,0.9,0.7),
  dep = c(0.5,1.0,1.0,1.0,0.7,0.3),
  ptr = c(1.0,1.0,1.0,1.0,0.8,0.5))

# =============================================================================
# main

# parameters
P = list(seed=0)
P$n = 100
P$zf = z1y*adur*2
P$ptr.r0.m    = .05
P$ptr.max.m   = 1.25
P$vio.r0.m    = .002
P$dep.o.r0.m  = .001
P$dep.x.r0.m  = .01
P$ptr.dz.m    = z1y
P$rr.vio.age  = fit.rr.age(a.ref,rr.age$vio)
P$rr.dep.age  = fit.rr.age(a.ref,rr.age$dep)
P$rr.dep.vio.dz = get.rr.dz(loc=30,scale=6.53,rr.tot=1)
P$rr.dep.dur  = -0.01
P$rr.ptr.age  = fit.rr.age(a.ref,rr.age$ptr)
P$rr.ptr.dep  = +0.7
P$rr.cdm.dep  = -0.7
# run model
Ps = lapply(1:7,function(s){ P$seed = s; P })
Is = sim.runs(Ps)
