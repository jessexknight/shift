
# =============================================================================
# pars

get.pars = function(seed=0,...){
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
  # RR: * -> dep onset
  P$ RR.dep_o.dep_p = 3     # RR: dep past -> dep onset
  P$iRR.dep_o.vio_z = 2     # (initial RR) transient RR: vio -> dep onset
  P$tsc.dep_o.vio_z = 30    # (time scale) transient RR: vio -> dep onset
  P$mRR.dep_o.vio_n = 2     # (max RR)    cumulative RR: vio -> dep onset
  P$nsc.dep_o.vio_n = 10    # (n scale)   cumulative RR: vio -> dep onset
  # RR: * -> dep recov
  P$tsc.dep_x.dep_u = 364   # (time scale)  duration RR: dep dur -> dep recov
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
  P$tsc.haz_x.haz_u = 364   # (time scale) duration RR haz dur -> haz recov
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
  P = list.update(P,...)
  P = add.pars(P)
}

add.pars = function(P){
  P$ndur = P$zf/z1y/adur    # num sim adurs (1 is dummy)
  P$ntot = P$n * (1+P$ndur) # total inds needed
  # RR: age
  P$RR.vio.age   = def.RR.age(c(10,60),c(1,1)) # RR: age -> vio
  P$RR.dep_o.age = def.RR.age(c(10,60),c(1,1)) # RR: age -> dep onset
  P$RR.haz_o.age = def.RR.age(c(10,60),c(1,1)) # RR: age -> haz onset
  P$RR.ptr_o.age = def.RR.age(c(10,60),c(1,1)) # RR: age -> ptr form
  # tRR: vio
  P$tRRu.dep_o.vio_z = def.tRR.exp(P$iRR.dep_o.vio_z,P$tsc.dep_o.vio_z) - 1 # tRR-1: vio -> dep onset
  P$tRRu.dep_x.vio_z = def.tRR.exp(P$iRR.dep_x.vio_z,P$tsc.dep_x.vio_z) - 1 # tRR-1: vio -> dep recov
  P$tRRu.haz_o.vio_z = def.tRR.exp(P$iRR.haz_o.vio_z,P$tsc.haz_o.vio_z) - 1 # tRR-1: vio -> haz onset
  P$tRRu.haz_x.vio_z = def.tRR.exp(P$iRR.haz_x.vio_z,P$tsc.haz_x.vio_z) - 1 # tRR-1: vio -> haz recov
  P$tRRu.ptr_o.vio_z = def.tRR.exp(P$iRR.ptr_o.vio_z,P$tsc.ptr_o.vio_z) - 1 # tRR-1: vio -> ptr form
  # cRR: vio
  P$cRR.dep_o.vio_n = def.cRR.exp(P$mRR.dep_o.vio_n,P$nsc.dep_o.vio_n) # cRR: vio -> dep onset
  P$cRR.haz_o.vio_n = def.cRR.exp(P$mRR.haz_o.vio_n,P$nsc.haz_o.vio_n) # cRR: vio -> haz onset
  P$cRR.ptr_o.vio_n = def.cRR.exp(P$mRR.ptr_o.vio_n,P$nsc.ptr_o.vio_n) # cRR: vio -> ptr form
  # dRR: durs
  P$dRRu.dep_x.dep_u = def.dRR.exp(P$tsc.dep_x.dep_u) - 1 # RR-1: dep dur -> dep recov
  P$dRRu.haz_x.haz_u = def.dRR.exp(P$tsc.haz_x.haz_u) - 1 # RR-1: dep dur -> dep recov
  # pre-compute RR-1 for all RR.*
  for (x in filter.names(P,'^RR')){
    P[[gsub('RR','RRu',x)]] = P[[x]] - 1
  }
  # for (x in filter.names(P,'^(t|c|d)RR')){ plot(grepl('RRu',x)+P[[x]]); title(x) } # DEBUG
  return(P)
}

# =============================================================================
# effect funs

def.RR.age = function(age,RR,eps=.001){
  n = len(RR)
  RR  = c(RR[1],RR,RR[n])
  age = c(age[1]-eps,age,age[n]+eps)
  RR.age = splinefun(age,RR,method='monoH.FC')(seq(amin,amax))
}

def.cRR.exp = function(mRR,nsc,nmax=1e4){
  n = 0:nmax # n > nmax is very bad since we will get 0 not mRR (TODO)
  cRR = 1 + (mRR-1) * (1-exp(-n/nsc))
}

def.dRR.exp = function(tsc){
  z = 1:(z1y*adur) # longest possible duration
  dRR = exp(-z*dtz/tsc)
}

def.tRR.exp = function(iRR,tsc,eps=.001){
  z = 1:ceiling(qexp(p=1-eps,rate=1/tsc)/dtz) # ensure we cover most AUC
  tRR = 1 + (iRR-1) * exp(-z*dtz/tsc)
}

map.tRR = function(tRRu,ze,z){
  # lookup & sum RR kernel for now (z) given prior events (ze)
  if (len(tRRu) < 1){ traceback(2); q() }
  RR = 1 + na.to.num(tRRu[z+1-ze])
}
