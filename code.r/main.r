source('meta.r')
source('utils.r')
source('plot.r')

# =============================================================================
# effect funs

fit.rr = function(.,t,flat=c(1,1)){
  if (missing(t)){ t = seq(dtz,max(.$t),dtz) }
  n = len(.$t)+flat[1]
  if (flat[1]){ .$t = c(.$t[1]-eps,.$t); .$rr = c(.$rr[1],.$rr) }
  if (flat[2]){ .$t = c(.$t,.$t[n]+eps); .$rr = c(.$rr,.$rr[n]) }
  rr = list(t=t,rr=splinefun(.$t,.$rr,method='monoH.FC')(t))
}

fit.rr.age = function(.){
  rr.age = fit.rr(.,seq(amin,amax),flat=c(1,0))
}

get.rr.evt = function(ze,z,rr.z){
  # lookup & sum RR kernel for today (z) given prior events (ze)
  rr = sum(rr.z[z+1-ze],na.rm=TRUE)
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
    aj = floor(Js$age-amin+1) # age vector for j
    # update vio --------------------------------------------------------------
    i = ij[which(runif(ij) < dtz
      * Js$vio.r0 # base rate
      * P$rr.vio.age[aj] # RR age
    )]
    Es$vio[i] = lapply(Es$vio[i],append,z)
    # update dep onset --------------------------------------------------------
    i = ij[which(runif(ij) < dtz
      * (!Js$dep.now) # among not dep
      * Js$dep.o.r0 # base rate
      * P$rr.dep.o.age[aj] # RR age
      * sapply(Es$vio[i.act],get.rr.evt,z,P$rr.dep.o.vio.z) # RR vio
    )]
    Is$dep.now[i] = TRUE
    Is$dep.evr[i] = TRUE
    Is$dep.z0[i] = z
    Es$dep.o[i] = lapply(Es$dep.o[i],append,z)
    # update dep recovery -----------------------------------------------------
    i = ij[which(runif(ij) < dtz
      * (Js$dep.now) # among dep
      * Js$dep.x.r0 # base rate
      * sapply(Es$vio[i.act],get.rr.evt,z,P$rr.dep.x.vio.z) # RR vio
      * 2^(dtz * (z - Js$dep.z0) / P$rr.dep.x.th) # RR dep dur
    )]
    Is$dep.now[i] = FALSE
    Es$dep.x[i] = lapply(Es$dep.x[i],append,z)
    # form ptrs ---------------------------------------------------------------
    i = ij[even.len(which(runif(ij) < dtz
      * (Js$ptr.n < Js$ptr.max) # among avail
      * Js$ptr.r0 # base rate
      * P$rr.ptr.age[aj] # RR age
      * (1 + P$arr.ptr.dep * Js$dep.now) # RR dep
    ))]
    Is$ptr.n[i] = Is$ptr.n[i] + 1
    Es$ptr.o[i] = lapply(Es$ptr.o[i],append,z)
    Ks = rbind(Ks,init.ptrs(P,Is[i,],z))
    # sex in ptrs -------------------------------------------------------------
    Xs = Ks[runif(nrow(Ks)) < Ks$f.sex,]
    cdm = (runif(nrow(Xs)) < 1
      * Xs$cdm # base prob
      * P$rr.cdm.dep ^ (Is$dep.now[Xs$i1] + Is$dep.now[Xs$i2]) # RR dep
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
# define rr splines

rr.age. = list(
  vio = list(t=c(15,20,30,50),rr=c(0.0,1.0,1.0,0.7)),
  dep = list(t=c(15,20,30,50),rr=c(0.0,1.0,1.0,0.3)),
  ptr = list(t=c(15,   30,50),rr=c(1.0,    1.0,0.5)))
# plot.rr(rr.age.,lapply(rr.age.,fit.rr.age)) + xlab('Age') + ylim(c(0,1))
#   plot.save('par','rr.age',h=2.5,w=5)
rr.vio. = list(
  dep.o = list(t=14*(0:4),rr=1+1.0*c(1.0,0.95,0.5,0.05,0.0)),
  dep.x = list(t=14*(0:4),rr=1-0.5*c(1.0,0.95,0.5,0.05,0.0)))
# plot.rr(rr.vio.,lapply(rr.vio.,fit.rr)) + scale_y_continuous(trans='log2')
#   plot.save('par','rr.vio',h=2.5,w=5)

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
P$rr.vio.age = fit.rr.age(rr.age.$vio)$rr
P$rr.dep.o.age = fit.rr.age(rr.age.$dep)$rr
P$rr.dep.o.vio.z = fit.rr(rr.vio.$dep.o)$rr
P$rr.dep.x.vio.z = fit.rr(rr.vio.$dep.x)$rr
P$rr.dep.x.th = 364
P$rr.ptr.age  = fit.rr.age(rr.age.$ptr)$rr
P$arr.ptr.dep = 0.5
P$rr.cdm.dep  = 0.75
# run model
Ps = lapply(1:7,function(s){ P$seed = s; P })
Is = sim.runs(Ps)
