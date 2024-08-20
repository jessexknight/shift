
# =============================================================================
# initialization funs

init.evts = function(I){
  # initialize event vectors for each event type & individual
  E0 = lapply(I,function(i){ numeric() })
  E  = lapply(evts,function(e){ E0 })
}

init.inds = function(P){
  n = P$n.tot
  # sample correlated parameters ----------------------------------------------
  dep = as.data.frame(copula(n,
    covs = P$dep.cov,
    qfuns = list(o.Ri=qgamma,x.Ri=qgamma),
    o.Ri = list(shape=P$all.Ri.shape,scale=P$dep_o.Ri.m/P$all.Ri.shape),
    x.Ri = list(shape=P$all.Ri.shape,scale=P$dep_x.Ri.m/P$all.Ri.shape)))
  # plot(dep,col=rgb(0,0,0,.1)) # DEBUG
  haz = as.data.frame(copula(n,
    covs = P$haz.cov,
    qfuns = list(o.Ri=qgamma,x.Ri=qgamma),
    o.Ri = list(shape=P$all.Ri.shape,scale=P$haz_o.Ri.m/P$all.Ri.shape),
    x.Ri = list(shape=P$all.Ri.shape,scale=P$haz_x.Ri.m/P$all.Ri.shape)))
  # plot(haz,col=rgb(0,0,0,.1)) # DEBUG
  ptr = as.data.frame(copula(n,
    covs = P$ptr.cov,
    qfuns = list(o.Ri=qgamma,x.Ri=qgamma,max=qgeom),
    o.Ri = list(shape=P$ptr.Ri.shape,scale=P$ptr_o.Ri.m/P$ptr.Ri.shape),
    x.Ri = list(shape=P$ptr.Ri.shape,scale=P$ptr_x.Ri.m/P$ptr.Ri.shape),
    max  = list(prob=1/P$ptr.max.m)))
  # for (i in 1:3) plot(ptr[,-i],col=rgb(0,0,0,.1)) # DEBUG
  # create main df of individuals ---------------------------------------------
  age = runif(n,min=amin-P$n.dur*adur,max=amax)
  vio.Ri = rgamma(n=n,shape=P$all.Ri.shape,scale=P$vio.Ri.m/P$all.Ri.shape)
  sex.Ri = rbeta(n=n,shape1=P$cdm.Pi.shapes[1],shape2=P$cdm.Pi.shapes[2])
  cdm.Pi = rbeta(n=n,shape1=P$cdm.Pi.shapes[1],shape2=P$cdm.Pi.shapes[2])
  I = data.frame(
    i = seq(n),
    z.born  = -age*z1y,
    age     = +age,
    age.act = runif(n,min=amin,max=20),
    # violence
    vio.Ri = vio.Ri,
    vio.zf = NA,
    vio.nt = 0,
    # depression
    dep_o.Ri = dep$o.Ri,
    dep_x.Ri = dep$x.Ri,
    dep.now  = FALSE,
    dep.past = FALSE,
    dep.zo   = NA,
    # hazardous drinking
    haz_o.Ri = haz$o.Ri,
    haz_x.Ri = haz$x.Ri,
    haz.now  = FALSE,
    haz.past = FALSE,
    haz.zo   = NA,
    # partnerships
    ptr_o.Ri = ptr$o.Ri,
    ptr_x.Ri = ptr$x.Ri,
    ptr.max  = 1+ptr$max,
    ptr.nw   = 0,
    # sex freq & condom prob
    sex.Ri = sex.Ri,
    cdm.Pi = cdm.Pi
  )
}

init.ptrs = function(P,I,i,z){
  # generate partners among I
  # mixing is random by splitting into 1st/2nd half
  n = len(i)/2
  if (n==0){ return(NULL) }
  i1 = i[0+(1:n)]
  i2 = i[n+(1:n)]
  K = data.frame(
    i1 = I$i[i1], # i of partner 1
    i2 = I$i[i2], # i of partner 2
    f.sex = I$sex.Ri[i1]/2 + I$sex.Ri[i2]/2, # sex freq
    p.cdm = I$cdm.Pi[i1]/2 + I$cdm.Pi[i2]/2, # condom prob
    zo = z # timestep ptr begins
  )
}

# =============================================================================
# rate & prob funs

rate.to.prob = function(R){ p = 1-exp(-R*dtz) }
rate.to.bool = function(R){ b = runif(R) < (1-exp(-R*dtz)) }
rate.to.num  = function(R){ n = rpois(len(R),R*dtz) }

rate.vio = function(P,J,aj){
  R = ( # among all
      J$vio.Ri # base rate
    * P$aRR.vio[aj] # RR age
); return(R) }

rate.dep_o = function(P,J,R,aj,z){
  j = which(!J$dep.now)
  R[j] = ( # among not dep
      J$dep_o.Ri[j] # base rate
    * P$aRR.dep_o[aj[j]] # RR age
    * (1 + P$RRu.dep_o.dep_p * J$dep.past[j]) # RR dep past
    * map.tRR(P$tRRu.dep_o.vio_zf,J$vio.zf[j],z) # tRR vio
    * P$nRR.dep_o.vio_nt[J$vio.nt[j]+1] # nRR vio
); return(R) }

rate.dep_x = function(P,J,R,aj,z){
  j = which(J$dep.now)
  R[j] = ( # among dep
      J$dep_x.Ri[j] # base rate
    * map.tRR(P$dRRu.dep_x.dep_u,J$dep.zo[j],z) # RR dep dur
    * map.tRR(P$tRRu.dep_x.vio_zf,J$vio.zf[j],z) # tRR vio
); return(R) }

rate.haz_o = function(P,J,R,aj,z){
  j = which(!J$haz.now)
  R[j] = ( # among not haz
      J$haz_o.Ri[j] # base rate
    * P$aRR.haz_o[aj[j]] # RR age
    * (1 + P$RRu.haz_o.haz_p * J$haz.past[j]) # RR haz past
    * (1 + P$RRu.haz_o.dep_w * J$dep.now[j]) # RR dep now
    * map.tRR(P$tRRu.haz_o.vio_zf,J$vio.zf[j],z) # tRR vio
    * P$nRR.haz_o.vio_nt[J$vio.nt[j]+1] # nRR vio
); return(R) }

rate.haz_x = function(P,J,R,aj,z){
  j = which(J$haz.now)
  R[j] = ( # among haz
      J$haz_x.Ri[j] # base rate
    * map.tRR(P$dRRu.haz_x.haz_u,J$haz.zo[j],z) # RR haz dur
    * (1 + P$RRu.haz_x.dep_w * J$dep.now[j]) # RR dep now
    * map.tRR(P$tRRu.haz_x.vio_zf,J$vio.zf[j],z) # tRR vio
); return(R) }

rate.ptr_o = function(P,J,R,aj,z){
  j = which(J$age > J$age.act & J$ptr.nw < J$ptr.max)
  R[j] = ( # among sex active & avail
      J$ptr_o.Ri[j] # base rate
    * P$aRR.ptr_o[aj[j]] # RR age
    * (1 + P$RRu.ptr_o.dep_w * J$dep.now[j]) # RR dep now
    * (1 + P$RRu.ptr_o.haz_w * J$haz.now[j]) # RR haz now
    * map.tRR(P$tRRu.ptr_o.vio_zf,J$vio.zf[j],z) # tRR vio
    * P$nRR.ptr_o.vio_nt[J$vio.nt[j]+1] # nRR vio
); return(R) }

rate.ptr_x = function(P,K,I){
  i1 = K$i1
  i2 = K$i2
  R = ( # among all (ptrs)
      0.5 * (I$ptr_x.Ri[i1] + I$ptr_x.Ri[i2])
    # TODO: RR age
    * (1 + P$RRu.ptr_x.dep_w * (I$dep.now[i1] + I$dep.now[i2])) # RR dep now
    * (1 + P$RRu.ptr_x.haz_w * (I$haz.now[i1] + I$haz.now[i2])) # RR haz now
); return(R) }

# =============================================================================
# run simulation

sim.run = function(P,rm.dum=TRUE){
  status(4,P$id)
  # initialization ------------------------------------------------------------
  set.seed(P$seed)
  I = init.inds(P) # individuals
  E = init.evts(P) # events
  K = NULL         # partnerships
  # event loop ----------------------------------------------------------------
  for (z in 1:P$zf){
    # if (z %% z1y == 0) { print(z/z1y) } # DEBUG
    # age inds ----------------------------------------------------------------
    I$age = I$age + 1/z1y
    # end ptrs ----------------------------------------------------------------
    b = rate.to.bool(rate.ptr_x(P,K,I))
    i = c(K$i1[b],K$i2[b])
    I$ptr.nw[i] = I$ptr.nw[i] - 1 # TODO: bug if any repeated i
    E$ptr_x[i] = lapply(E$ptr_x[i],append,z)
    E$ptr_u[i] = wapply(append,E$ptr_u[i],z-K$zo[b])
    K = K[!b,]
    # select active inds ------------------------------------------------------
    i = which(I$age > amin & I$age < amax)
    J = I[i,] # read only copy of active
    ij = match(J$i,I$i) # map j -> j
    aj = floor(J$age-amin+1) # age vector for j
    R0 = numeric(nrow(J)) # init rate = 0 for j
    # vio events --------------------------------------------------------------
    i = ij[which(rate.to.bool(rate.vio(P,J,aj)))]
    I$vio.zf[i] = z
    I$vio.nt[i] = I$vio.nt[i] + 1
    E$vio[i] = lapply(E$vio[i],append,z)
    # begin dep ---------------------------------------------------------------
    i = ij[which(rate.to.bool(rate.dep_o(P,J,R0,aj,z)))]
    I$dep.now[i] = TRUE
    I$dep.past[i] = TRUE
    I$dep.zo[i] = z
    E$dep_o[i] = lapply(E$dep_o[i],append,z)
    # end dep -----------------------------------------------------------------
    i = ij[which(rate.to.bool(rate.dep_x(P,J,R0,aj,z)))]
    I$dep.now[i] = FALSE
    E$dep_x[i] = lapply(E$dep_x[i],append,z)
    # begin haz ---------------------------------------------------------------
    i = ij[which(rate.to.bool(rate.haz_o(P,J,R0,aj,z)))]
    I$haz.now[i] = TRUE
    I$haz.past[i] = TRUE
    I$haz.zo[i] = z
    E$haz_o[i] = lapply(E$haz_o[i],append,z)
    # end haz -----------------------------------------------------------------
    i = ij[which(rate.to.bool(rate.haz_x(P,J,R0,aj,z)))]
    I$haz.now[i] = FALSE
    E$haz_x[i] = lapply(E$haz_x[i],append,z)
    # begin ptrs --------------------------------------------------------------
    i = ij[even.len(which(rate.to.bool(rate.ptr_o(P,J,R0,aj,z))))]
    I$ptr.nw[i] = I$ptr.nw[i] + 1
    E$ptr_o[i] = lapply(E$ptr_o[i],append,z)
    K = rbind(K,init.ptrs(P,I,i,z))
    # sex in ptrs -------------------------------------------------------------
    # TODO
  }
  # clean-up ------------------------------------------------------------------
  if (rm.dum){ # remove initial dummy population
    i = which(I$z.born > -amin*z1y)
    I = I[i,]
    E = lapply(E,`[`,i)
  }
  M = list(P=P,I=I,E=E)
}

sim.runs = function(Ps,.par=TRUE){
  # run.sim in parallel for each (P)arameter set in Ps
  status(3,'sim.runs: ',len(Ps))
  Ms = par.lapply(Ps,sim.run,.par=.par); status(4,'\n')
  return(Ms)
}
