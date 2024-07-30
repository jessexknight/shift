
# =============================================================================
# initialization funs

init.evts = function(Is){
  # initialize event vectors for each event type & individual
  E = lapply(Is$i,function(i){ numeric() })
  Es = lapply(evts,function(e){ E })
}

init.inds = function(P){
  n = P$ntot
  Is = data.frame(
    i = seq(n),
    age = runif(n,min=amin-P$ndur*adur,max=amax),
    age.act = runif(n,min=amin,max=20),
    # violence
    vio.Ri = rexp(n=n,rate=1/P$vio.Ri.m),
    vio.zf = NA,
    vio.n  = 0,
    # depression
    dep_o.Ri = rexp(n=n,rate=1/P$dep_o.Ri.m),
    dep_x.Ri = rexp(n=n,rate=1/P$dep_x.Ri.m),
    dep.now  = FALSE,
    dep.past = FALSE,
    dep.zo   = NA,
    # depression
    haz_o.Ri = rexp(n=n,rate=1/P$haz_o.Ri.m),
    haz_x.Ri = rexp(n=n,rate=1/P$haz_x.Ri.m),
    haz.now  = FALSE,
    haz.past = FALSE,
    haz.zo   = NA,
    # partnerships
    ptr_o.Ri = rexp(n=n,rate=1/P$ptr_o.Ri.m),
    ptr_x.Ri = rexp(n=n,rate=1/P$ptr_x.Ri.m),
    ptr.max = 1+rgeom(n=n,prob=1/P$ptr.max.m),
    ptr.n = 0
    # TODO: condoms
  )
}

init.ptrs = function(P,Is,i,z){
  # generate partners among Is
  # mixing is random by splitting into 1st/2nd half
  n = len(i)/2
  if (n==0){ return(NULL) }
  i1 = i[0+(1:n)]
  i2 = i[n+(1:n)]
  Ks = data.frame(
    i1 = Is$i[i1], # i of partner 1
    i2 = Is$i[i2], # i of partner 2
    zo = z # timestep ptr starts
    # TODO: f.sex, p.cdm
  )
}

# =============================================================================
# rate & prob funs

rate.to.prob = function(R){ p = 1-exp(-R*dtz) }
rate.to.bool = function(R){ b = runif(R) < (1-exp(-R*dtz)) }
rate.to.num  = function(R){ n = rpois(len(R),R*dtz) }

rate.vio = function(P,Js,aj){
  R = ( # among all
      Js$vio.Ri # base rate
    * P$aRR.vio[aj] # RR age
); return(R) }

rate.dep_o = function(P,Js,R,aj,z){
  j = which(!Js$dep.now)
  R[j] = ( # among not dep
      Js$dep_o.Ri[j] # base rate
    * P$aRR.dep_o[aj[j]] # RR age
    * (1 + P$RRu.dep_o.dep_p * Js$dep.past[j]) # RR dep past
    * map.tRR(P$tRRu.dep_o.vio_z,Js$vio.zf[j],z) # tRR vio
    * P$nRR.dep_o.vio_n[Js$vio.n[j]+1] # nRR vio
); return(R) }

rate.dep_x = function(P,Js,R,aj,z){
  j = which(Js$dep.now)
  R[j] = ( # among dep
      Js$dep_x.Ri[j] # base rate
    * map.tRR(P$dRRu.dep_x.dep_u,Js$dep.zo[j],z) # RR dep dur
    * map.tRR(P$tRRu.dep_x.vio_z,Js$vio.zf[j],z) # tRR vio
); return(R) }

rate.haz_o = function(P,Js,R,aj,z){
  j = which(!Js$haz.now)
  R[j] = ( # among not haz
      Js$haz_o.Ri[j] # base rate
    * P$aRR.haz_o[aj[j]] # RR age
    * (1 + P$RRu.haz_o.haz_p * Js$haz.past[j]) # RR haz past
    * (1 + P$RRu.haz_o.dep_w * Js$dep.now[j]) # RR dep now
    * map.tRR(P$tRRu.haz_o.vio_z,Js$vio.zf[j],z) # tRR vio
    * P$nRR.haz_o.vio_n[Js$vio.n[j]+1] # nRR vio
); return(R) }

rate.haz_x = function(P,Js,R,aj,z){
  j = which(Js$haz.now)
  R[j] = ( # among haz
      Js$haz_x.Ri[j] # base rate
    * map.tRR(P$dRRu.haz_x.haz_u,Js$haz.zo[j],z) # RR haz dur
    * (1 + P$RRu.haz_x.dep_w * Js$dep.now[j]) # RR dep now
    * map.tRR(P$tRRu.haz_x.vio_z,Js$vio.zf[j],z) # tRR vio
); return(R) }

rate.ptr_o = function(P,Js,R,aj,z){
  j = which(Js$ptr.n < Js$ptr.max)
  R[j] = ( # among avail
      Js$ptr_o.Ri[j] # base rate
    * P$aRR.ptr_o[aj[j]] # RR age
    * (1 + P$RRu.ptr_o.dep_w * Js$dep.now[j]) # RR dep now
    * (1 + P$RRu.ptr_o.haz_w * Js$haz.now[j]) # RR haz now
    * map.tRR(P$tRRu.ptr_o.vio_z,Js$vio.zf[j],z) # tRR vio
    * P$nRR.ptr_o.vio_n[Js$vio.n[j]+1] # nRR vio
); return(R) }

rate.ptr_x = function(P,Ks,Is){
  i1 = Ks$i1
  i2 = Ks$i2
  R = ( # among all (ptrs)
      0.5 * (Is$ptr_x.Ri[i1] + Is$ptr_x.Ri[i2])
    # TODO: RR age
    * (1 + P$RRu.ptr_x.dep_w * (Is$dep.now[i1] + Is$dep.now[i2])) # RR dep now
    * (1 + P$RRu.ptr_x.haz_w * (Is$haz.now[i1] + Is$haz.now[i2])) # RR haz now
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
    # ptr dissolve ------------------------------------------------------------
    b = rate.to.bool(rate.ptr_x(P,Ks,Is))
    i = c(Ks$i1[b],Ks$i2[b])
    Is$ptr.n[i] = Is$ptr.n[i] - 1
    Es$ptr_x[i] = lapply(Es$ptr_x[i],append,z)
    Ks = Ks[!b,]
    # select active inds ------------------------------------------------------
    i = which(Is$age > Is$age.act & Is$age < amax)
    Js = Is[i,] # read only copy of active
    ij = match(Js$i,Is$i) # map j -> j
    aj = floor(Js$age-amin+1) # age vector for j
    R0 = numeric(nrow(Js)) # init rate = 0 for j
    # update vio events -------------------------------------------------------
    i = ij[which(rate.to.bool(rate.vio(P,Js,aj)))]
    Is$vio.zf[i] = z
    Is$vio.n[i] = Is$vio.n[i] + 1
    Es$vio[i] = lapply(Es$vio[i],append,z)
    # update dep onset --------------------------------------------------------
    i = ij[which(rate.to.bool(rate.dep_o(P,Js,R0,aj,z)))]
    Is$dep.now[i] = TRUE
    Is$dep.past[i] = TRUE
    Is$dep.zo[i] = z
    Es$dep_o[i] = lapply(Es$dep_o[i],append,z)
    # update dep recovery -----------------------------------------------------
    i = ij[which(rate.to.bool(rate.dep_x(P,Js,R0,aj,z)))]
    Is$dep.now[i] = FALSE
    Es$dep_x[i] = lapply(Es$dep_x[i],append,z)
    # update haz onset --------------------------------------------------------
    i = ij[which(rate.to.bool(rate.haz_o(P,Js,R0,aj,z)))]
    Is$haz.now[i] = TRUE
    Is$haz.past[i] = TRUE
    Is$haz.zo[i] = z
    Es$haz_o[i] = lapply(Es$haz_o[i],append,z)
    # update haz recovery -----------------------------------------------------
    i = ij[which(rate.to.bool(rate.haz_x(P,Js,R0,aj,z)))]
    Is$haz.now[i] = FALSE
    Es$haz_x[i] = lapply(Es$haz_x[i],append,z)
    # ptr formation -----------------------------------------------------------
    i = ij[even.len(which(rate.to.bool(rate.ptr_o(P,Js,R0,aj,z))))]
    Is$ptr.n[i] = Is$ptr.n[i] + 1
    Es$ptr_o[i] = lapply(Es$ptr_o[i],append,z)
    Ks = rbind(Ks,init.ptrs(P,Is,i,z))
    # sex in ptrs -------------------------------------------------------------
    # TODO
  }
  Is = sim.out(P,Is,Es)
}

sim.out = function(P,Is,Es,rm.dum=TRUE){
  # clean-up simulation output
  if (rm.dum){ # remove initial dummy population
    i = which(Is$age < (P$zf/z1y+amin))
    Is = Is[i,]
    Es = lapply(Es,`[`,i)
  }
  # compute some extra variables
  Is$ptr.tot = sapply(Es$ptr_o,len) # lifetime ptrs
  Is = cbind(seed=P$seed,Is)
}

sim.runs = function(Ps,.par=TRUE){
  # run.sim in parallel for each (P)arameter set in Ps
  Is = rbind.lapply(Ps,sim.run,.par=.par & len(Ps)>1)
}
