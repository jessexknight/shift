# =============================================================================
# tiny utils

dtz  =  7 # days in 1 timestep
z1y  = 52 # timesteps in 1 year
amin = 15 # age of cohort entry
amax = 50 # age of cohort exit
adur = amax - amin # duration in cohort
evts = c('vio','dep.o','dep.x'); names(evts) = evts # event types

even.len = function(i){
  # truncate vector i to have an even length
  length(i) = length(i) - (length(i) %% 2)
  return(i)
}

# =============================================================================
# convolution effect utils

fit.eff.dz = function(q.ref,t.ref){
  # fit decaying logistic kernel to reach q.ref by t.ref
  err.fun = function(par){
    t.par = qlogis(1-q.ref,loc=par[1],scale=par[2])
    err = sum((t.par-t.ref)^2) }
  par = optim(c(loc=0,scale=1),err.fun,lower=c(NA,1e-3),method='L-BFGS-B')$par
}

get.eff.dz = function(loc,scale,e.tot,eps=1e-3){
  # pre-compute scaled decaying logistic kernel
  zs = 0:ceiling(qlogis(1-eps,loc=loc,scale=scale)/dtz)
  eff.z = 1-plogis(zs*dtz,loc=loc,scale=scale)
  # plot(0:z1y,1-plogis(0:z1y,loc=loc,scale=scale)); lines(zs*dtz,eff.z) # DEBUG
  eff.z = eff.z * e.tot / sum(eff.z)
}

get.eff.evt = function(ze,z,eff.dz){
  # lookup & sum effect kernel for today (z) given prior events (ze)
  eff = sum(eff.dz[z+1-ze],na.rm=TRUE)
}

# =============================================================================
# initialization & generation utils

init.evts = function(Is){
  # initialize event vectors for each event type & individual
  E = lapply(Is$i,function(i){ numeric() })
  Es = lapply(evts,function(e){ E })
}

init.inds = function(P,zs){
  # initialize all individuals needed for timesteps zs
  ny = length(zs)/z1y/adur # num sim years
  n  = P$n * (1+ny)        # total inds needed
  Is = data.frame(
    i = seq(n),
    age = runif(n,min=amin-ny*adur,max=amax),
    # partnerships
    ptr.n = 0,
    ptr.r0 = rexp(n=n,rate=1/P$ptr.r0.m),
    ptr.max = pmin(3,1+rgeom(n=n,prob=1/P$ptr.max.m)),
    cdm.p0 = runif(n=n,min=0,max=1),
    # structural factors & mediators
    vio.r0 = rexp(n=n,rate=1/P$vio.r0.m),
    dep.o.r0 = rexp(n=n,rate=1/P$dep.o.r0.m),
    dep.x.r0 = rexp(n=n,rate=1/P$dep.x.r0.m),
    dep.stat = 0
  )
}

gen.ptrs = function(P,Is,z){
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
    cdm = Is$cdm.p0[i1]/2 + Is$cdm.p0[i2]/2 # condom prob (average)
  )
}

# =============================================================================
# run simulation

run.sim = function(P,zs){
  # initialization ------------------------------------------------------------
  Is = init.inds(P,zs) # individuals
  Es = init.evts(Is)   # events
  Ks = NULL            # partnerships
  for (z in zs){
    # age inds ----------------------------------------------------------------
    Is$age = Is$age + 1/z1y
    # break ptrs --------------------------------------------------------------
    k.brk = Ks$zf < z
    i.brk = c(Ks$i1[k.brk],Ks$i2[k.brk])
    Is$ptr.n[i.brk] = Is$ptr.n[i.brk] - 1
    Ks = Ks[!k.brk,]
    # select active inds ------------------------------------------------------
    i.act = which(Is$age > amin & Is$age < amax)
    Js = Is[i.act,]       # read only copy of active
    ij = match(Js$i,Is$i) # map j -> j
    # update vio --------------------------------------------------------------
    i = ij[which(runif(ij) < Js$vio.r0 * dtz * exp(0))]
    Es$vio[i] = lapply(Es$vio[i],append,z)
    # update dep --------------------------------------------------------------
    b.dep = as.logical(Js$dep.stat)
    eff.vio.dep = sapply(Es$vio[i.act],get.eff.evt,z,P$eff.vio.dep.dz)
    # dep.o (onset)
    i = ij[which(runif(ij) < (!b.dep) * Js$dep.o.r0 * dtz * exp(0
      + eff.vio.dep
    ))]
    Is$dep.stat[i] = 1
    Es$dep.o[i] = lapply(Es$dep.o[i],append,z)
    # dep.x (recovery)
    i = ij[which(runif(ij) < (b.dep) * Js$dep.x.r0 * dtz * exp(0
      - eff.vio.dep
    ))]
    Is$dep.stat[i] = 0
    Es$dep.x[i] = lapply(Es$dep.x[i],append,z)
    # form ptrs ---------------------------------------------------------------
    i = ij[even.len(which(runif(ij) < (Js$ptr.n < Js$ptr.max) * Js$ptr.r0 * dtz * exp(0
    )))]
    Is$ptr.n[i] = Is$ptr.n[i] + 1
    Ks = rbind(Ks,gen.ptrs(P,Is[i,],z))
  }
  # lapply(Es,function(E){ e = sapply(E,len); hist(e,0:max(e)) }) # DEBUG
  return(Is)
}

# -----------------------------------------------------------------------------
# fit conv effects

# print(fit.eff.dz(c(.5,.01),c(30,60)))

# =============================================================================
# main

# parameters
P = list()
P$n = 1000
P$ptr.r0.m    = .05
P$ptr.max.m   = 1.25
P$vio.r0.m    = .002
P$dep.o.r0.m  = .001
P$dep.x.r0.m  = .01
P$ptr.dz.m    = z1y
P$eff.vio.dep.dz = get.eff.dz(loc=30,scale=6.53,e.tot=1)
# run model
set.seed(0)
zs = seq(z1y*adur)
Is = run.sim(P,zs)
