
# =============================================================================
# pars

get.pars = function(seed=0,...,null=NULL){
  P = list(seed=seed)
  P$n = 1000
  P$zf = z1y*adur*2
  # base rates
  P$vio.Ri.m    = 1/364     # (mean) base rate: violence event
  P$dep_o.Ri.m  = 1/3640    # (mean) base rate: depression onset
  P$dep_x.Ri.m  = 1/364     # (mean) base rate: depression recovery
  P$haz_o.Ri.m  = 1/3640    # (mean) base rate: hazdrink onset
  P$haz_x.Ri.m  = 1/364     # (mean) base rate: hazdrink recovery
  P$ptr_o.Ri.m  = 1/35      # (mean) base rate: partner formation
  P$ptr_x.Ri.m  = 1/364     # (mean) base rate: partner dissolution
  P$ptr.max.m   = 1.50      # (mean) max num partners
  # RR: age -> *
  P$aRR.vio.ages   = c(amin,amax) # (age points) RR: age -> vio
  P$aRR.vio.RRs    = c(1.00,1.00) # (RR  points) RR: age -> vio
  P$aRR.dep_o.ages = c(amin,amax) # (age points) RR: age -> dep onset
  P$aRR.dep_o.RRs  = c(1.00,1.00) # (RR  points) RR: age -> dep onset
  P$aRR.haz_o.ages = c(amin,amax) # (age points) RR: age -> haz onset
  P$aRR.haz_o.RRs  = c(1.00,1.00) # (RR  points) RR: age -> haz onset
  P$aRR.ptr_o.ages = c(amin,amax) # (age points) RR: age -> ptr form
  P$aRR.ptr_o.RRs  = c(1.00,1.00) # (RR  points) RR: age -> ptr form
  # RR: * -> dep onset
  P$ RR.dep_o.dep_p = 3     # RR: dep past -> dep onset
  P$iRR.dep_o.vio_z = 2     # (initial RR) transient RR: vio -> dep onset
  P$tsc.dep_o.vio_z = 30    # (time scale) transient RR: vio -> dep onset
  P$mRR.dep_o.vio_n = 2     # (max RR)    cumulative RR: vio -> dep onset
  P$nsc.dep_o.vio_n = 10    # (n scale)   cumulative RR: vio -> dep onset
  # RR: * -> dep recov
  P$dsc.dep_x.dep_u = 364   # (dur scale)  duration RR: dep dur -> dep recov
  P$iRR.dep_x.vio_z = 1/2   # (initial RR) transient RR: vio -> dep recov
  P$tsc.dep_x.vio_z = 30    # (time scale) transient RR: vio -> dep recov
  # RR: * -> haz onset
  P$ RR.haz_o.haz_p = 3     # RR: haz past -> haz onset
  P$ RR.haz_o.dep_w = 3     # RR: dep now -> haz onset
  P$iRR.haz_o.vio_z = 2     # (initial RR) transient RR: vio -> haz onset
  P$tsc.haz_o.vio_z = 30    # (time scale) transient RR: vio -> haz onset
  P$mRR.haz_o.vio_n = 2     # (max RR)    cumulative RR: vio -> haz onset
  P$nsc.haz_o.vio_n = 10    # (n scale)   cumulative RR: vio -> haz onset
  # RR: * -> haz recov
  P$ RR.haz_x.dep_w = 1/3   # RR: dep now -> haz recov
  P$dsc.haz_x.haz_u = 364   # (dur scale)  duration RR: haz dur -> haz recov
  P$iRR.haz_x.vio_z = 1/2   # (initial RR) transient RR: vio -> haz recov
  P$tsc.haz_x.vio_z = 30    # (time scale) transient RR: vio -> haz recov
  # RR: * -> ptr form
  P$ RR.ptr_o.dep_w = 1/2   # RR: dep now -> ptr form
  P$ RR.ptr_o.haz_w = 2     # RR: haz now -> ptr form
  P$iRR.ptr_o.vio_z = 1/2   # (initial RR) transient RR: vio -> ptr form
  P$tsc.ptr_o.vio_z = 30    # (time scale) transient RR: vio -> ptr form
  P$mRR.ptr_o.vio_n = 1/2   # (max RR)    cumulative RR: vio -> ptr form
  P$nsc.ptr_o.vio_n = 10    # (n scale)   cumulative RR: vio -> ptr form
  # RR: * -> ptr dissol
  P$ RR.ptr_x.dep_w = 2     # RR: dep now -> ptr dissol
  P$ RR.ptr_x.haz_w = 2     # RR: haz now -> ptr dissol
  # overwrite & add conditional
  if (is.list(null)){ P = null.pars(P,null,null$save) }
  P = ulist(P,...)
  P = cond.pars(P)
}

cond.pars = function(P){
  P$ndur = P$zf/z1y/adur    # num sim adurs (1 is dummy)
  P$ntot = P$n * (1+P$ndur) # total inds needed
  # RR: age
  P$aRR.vio   = def.RR.age(P$aRR.vio.ages,P$aRR.vio.RRs) # RR: age -> vio
  P$aRR.dep_o = def.RR.age(P$aRR.dep_o.ages,P$aRR.dep_o.RRs) # RR: age -> dep onset
  P$aRR.haz_o = def.RR.age(P$aRR.haz_o.ages,P$aRR.haz_o.RRs) # RR: age -> haz onset
  P$aRR.ptr_o = def.RR.age(P$aRR.ptr_o.ages,P$aRR.ptr_o.RRs) # RR: age -> ptr form
  # tRR: vio
  P$tRRu.dep_o.vio_z = def.tRR.exp(P$iRR.dep_o.vio_z,P$tsc.dep_o.vio_z) - 1 # tRR-1: vio -> dep onset
  P$tRRu.dep_x.vio_z = def.tRR.exp(P$iRR.dep_x.vio_z,P$tsc.dep_x.vio_z) - 1 # tRR-1: vio -> dep recov
  P$tRRu.haz_o.vio_z = def.tRR.exp(P$iRR.haz_o.vio_z,P$tsc.haz_o.vio_z) - 1 # tRR-1: vio -> haz onset
  P$tRRu.haz_x.vio_z = def.tRR.exp(P$iRR.haz_x.vio_z,P$tsc.haz_x.vio_z) - 1 # tRR-1: vio -> haz recov
  P$tRRu.ptr_o.vio_z = def.tRR.exp(P$iRR.ptr_o.vio_z,P$tsc.ptr_o.vio_z) - 1 # tRR-1: vio -> ptr form
  # nRR: vio
  P$nRR.dep_o.vio_n = def.nRR.exp(P$mRR.dep_o.vio_n,P$nsc.dep_o.vio_n) # nRR: vio -> dep onset
  P$nRR.haz_o.vio_n = def.nRR.exp(P$mRR.haz_o.vio_n,P$nsc.haz_o.vio_n) # nRR: vio -> haz onset
  P$nRR.ptr_o.vio_n = def.nRR.exp(P$mRR.ptr_o.vio_n,P$nsc.ptr_o.vio_n) # nRR: vio -> ptr form
  # dRR: durs
  P$dRRu.dep_x.dep_u = def.dRR.exp(P$dsc.dep_x.dep_u) - 1 # RR-1: dep dur -> dep recov
  P$dRRu.haz_x.haz_u = def.dRR.exp(P$dsc.haz_x.haz_u) - 1 # RR-1: dep dur -> dep recov
  # pre-compute RR-1 for all RR.*
  for (x in filter.names(P,'^RR')){
    P[[gsub('RR','RRu',x)]] = P[[x]] - 1
  }
  # for (x in filter.names(P,'^(t|c|d)RR')){ plot(grepl('RRu',x)+P[[x]]); title(x) } # DEBUG
  return(P)
}

null.pars = function(P,null,save){
  # overwrite most pars via regex list; but save (exempt) some by name
  # to skip any default {re}, use null=list({re}=NULL)
  null.default = list(
    'Ri\\.m$'      = 0,   # base rates
    '^.?RR\\.'     = 1,   # RR, aRR, iRR, mRR
    '^tsc\\.'      = eps, # time scales
    '^(d|n)sc\\.'  = Inf) # dur & n scales
  null = ulist(null.default,null)
  P.save = P[save] # save exempt
  for (re in names(null)){ # for each regex
    if (is.null(null[[re]])){ next } # do nothing if NULL
    for (x in filter.names(P,re)){ # for each matching par
      P[[x]] = null[[re]] }} # overwrite
  P = ulist(P,P.save) # restore saved
}

# =============================================================================
# effect funs

def.RR.age = function(age,RR,eps=.001){
  n = len(RR)
  RR  = c(RR[1],RR,RR[n])
  age = c(age[1]-eps,age,age[n]+eps)
  RR.age = splinefun(age,RR,method='monoH.FC')(seq(amin,amax))
}

def.nRR.exp = function(mRR,nsc){
  n = 0:(z1y*adur) # nmax = all active timesteps
  nRR = 1 + (mRR-1) * (1-exp(-n/nsc))
}

def.dRR.exp = function(tsc){
  z = 1:(z1y*adur) # dmax = all active timesteps
  dRR = exp(-z*dtz/tsc)
}

def.tRR.exp = function(iRR,tsc,eps=.001){
  z = 1:ceiling(qexp(p=1-eps,rate=1/tsc)/dtz) # ensure we cover most AUC
  tRR = 1 + (iRR-1) * exp(-z*dtz/tsc)
}

map.tRR = function(tRRu,ze,z){
  # lookup & sum RR kernel for now (z) given prior events (ze)
  RR = 1 + na.to.num(tRRu[z+1-ze])
}

num.dz = function(zes,z,dz){
  n = sum(zes <= z & zes >= z+1-dz)
}

any.dz = function(zes,z,dz){
  b = num.dz(zes,z,dz) > 0
}
