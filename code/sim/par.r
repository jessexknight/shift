
# =============================================================================
# pars

get.pars = function(seed=0,...,case='base',null=NULL){
  P = list(case=case,seed=seed,id=sprintf('%6d',seed))
  P$n.pop = 1000
  P$n.dur = 1+1
  # base rates
  P$vio.Ri.m    = 1/t1y     # (mean) base rate: violence
  P$dep_o.Ri.m  = 0.01/t1y  # (mean) base rate: depression begin
  P$dep_x.Ri.m  = 1.00/t1y  # (mean) base rate: depression end
  P$haz_o.Ri.m  = 0.01/t1y  # (mean) base rate: hazdrink begin
  P$haz_x.Ri.m  = 1.00/t1y  # (mean) base rate: hazdrink end
  P$ptr_o.Ri.m  = 1/35      # (mean) base rate: partner begin
  P$ptr_x.Ri.m  = 1/t1y     # (mean) base rate: partner end
  P$ptr.max.m   = 1.50      # (mean) max num partners
  P$dep.cov     = -.9       # approx covariance among dep_o,dep_x
  P$haz.cov     = -.9       # approx covariance among haz_o,haz_x
  P$ptr.cov     = +.9       # approx covariance among ptr_o,ptr_x,ptr.max
  P$ptr.Ri.shape = 3        # (gamma shape): ptr_o,ptr_x
  P$all.Ri.shape = 1        # (gamma shape): all other base rates
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
  P$iRR.dep_o.vio_zf = 2     # (initial RR) transient RR: vio -> dep begin
  P$tsc.dep_o.vio_zf = 30    # (time scale) transient RR: vio -> dep begin
  P$mRR.dep_o.vio_nt = 2     # (max RR)    cumulative RR: vio -> dep begin
  P$nsc.dep_o.vio_nt = 10    # (n scale)   cumulative RR: vio -> dep begin
  # RR: * -> dep end
  P$dsc.dep_x.dep_u  = t1y   # (dur scale)  duration RR: dep dur -> dep end
  P$iRR.dep_x.vio_zf = 1/2   # (initial RR) transient RR: vio -> dep end
  P$tsc.dep_x.vio_zf = 30    # (time scale) transient RR: vio -> dep end
  # RR: * -> haz begin
  P$ RR.haz_o.haz_p  = 3     # RR: haz past -> haz begin
  P$ RR.haz_o.dep_w  = 3     # RR: dep now -> haz begin
  P$iRR.haz_o.vio_zf = 2     # (initial RR) transient RR: vio -> haz begin
  P$tsc.haz_o.vio_zf = 30    # (time scale) transient RR: vio -> haz begin
  P$mRR.haz_o.vio_nt = 2     # (max RR)    cumulative RR: vio -> haz begin
  P$nsc.haz_o.vio_nt = 10    # (n scale)   cumulative RR: vio -> haz begin
  # RR: * -> haz end
  P$ RR.haz_x.dep_w  = 1/3   # RR: dep now -> haz end
  P$dsc.haz_x.haz_u  = t1y   # (dur scale)  duration RR: haz dur -> haz end
  P$iRR.haz_x.vio_zf = 1/2   # (initial RR) transient RR: vio -> haz end
  P$tsc.haz_x.vio_zf = 30    # (time scale) transient RR: vio -> haz end
  # RR: * -> ptr begin
  P$ RR.ptr_o.dep_w  = 1/2   # RR: dep now -> ptr begin
  P$ RR.ptr_o.haz_w  = 2     # RR: haz now -> ptr begin
  P$iRR.ptr_o.vio_zf = 1/2   # (initial RR) transient RR: vio -> ptr begin
  P$tsc.ptr_o.vio_zf = 30    # (time scale) transient RR: vio -> ptr begin
  P$mRR.ptr_o.vio_nt = 1/2   # (max RR)    cumulative RR: vio -> ptr begin
  P$nsc.ptr_o.vio_nt = 10    # (n scale)   cumulative RR: vio -> ptr begin
  # RR: * -> ptr end
  P$ RR.ptr_x.dep_w  = 2     # RR: dep now -> ptr end
  P$ RR.ptr_x.haz_w  = 2     # RR: haz now -> ptr end
  # overwrite & add conditional
  if (is.list(null)){ P = null.pars(P,null,null$save) }
  P = ulist(P,...)
  P = cond.pars(P)
}

cond.pars = function(P){
  P$zf    = P$n.dur*adur*z1y  # final timestep
  P$n.tot = P$n.pop * (1+P$n.dur) # total inds needed
  # RR: age
  P$aRR.vio   = def.RR.age(P$aRR.vio.ages,P$aRR.vio.RRs) # RR: age -> vio
  P$aRR.dep_o = def.RR.age(P$aRR.dep_o.ages,P$aRR.dep_o.RRs) # RR: age -> dep begin
  P$aRR.haz_o = def.RR.age(P$aRR.haz_o.ages,P$aRR.haz_o.RRs) # RR: age -> haz begin
  P$aRR.ptr_o = def.RR.age(P$aRR.ptr_o.ages,P$aRR.ptr_o.RRs) # RR: age -> ptr begin
  # tRR: vio
  P$tRRu.dep_o.vio_zf = def.tRR.exp(P$iRR.dep_o.vio_zf,P$tsc.dep_o.vio_zf) - 1 # tRR-1: vio -> dep begin
  P$tRRu.dep_x.vio_zf = def.tRR.exp(P$iRR.dep_x.vio_zf,P$tsc.dep_x.vio_zf) - 1 # tRR-1: vio -> dep end
  P$tRRu.haz_o.vio_zf = def.tRR.exp(P$iRR.haz_o.vio_zf,P$tsc.haz_o.vio_zf) - 1 # tRR-1: vio -> haz begin
  P$tRRu.haz_x.vio_zf = def.tRR.exp(P$iRR.haz_x.vio_zf,P$tsc.haz_x.vio_zf) - 1 # tRR-1: vio -> haz end
  P$tRRu.ptr_o.vio_zf = def.tRR.exp(P$iRR.ptr_o.vio_zf,P$tsc.ptr_o.vio_zf) - 1 # tRR-1: vio -> ptr begin
  # nRR: vio
  P$nRR.dep_o.vio_nt = def.nRR.exp(P$mRR.dep_o.vio_nt,P$nsc.dep_o.vio_nt) # nRR: vio -> dep begin
  P$nRR.haz_o.vio_nt = def.nRR.exp(P$mRR.haz_o.vio_nt,P$nsc.haz_o.vio_nt) # nRR: vio -> haz begin
  P$nRR.ptr_o.vio_nt = def.nRR.exp(P$mRR.ptr_o.vio_nt,P$nsc.ptr_o.vio_nt) # nRR: vio -> ptr begin
  # dRR: durs
  P$dRRu.dep_x.dep_u = def.dRR.exp(P$dsc.dep_x.dep_u) - 1 # dRR-1: dep dur -> dep end
  P$dRRu.haz_x.haz_u = def.dRR.exp(P$dsc.haz_x.haz_u) - 1 # dRR-1: dep dur -> dep end
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
  null = ulist(null.all,null)
  P.save = P[save] # save exempt
  for (re in names(null)){ # for each regex
    if (is.null(null[[re]])){ next } # do nothing if NULL
    for (x in filter.names(P,re)){ # for each matching par
      P[[x]] = null[[re]] }} # overwrite
  P = ulist(P,P.save) # restore saved
}

null.all = list(
  'Ri\\.m$'     = 0,     # base rates
  '^.?RR\\.'    = 1,     # RR, aRR, iRR, mRR
  '^tsc\\.'     = 1e-12, # time scales
  '^(d|n)sc\\.' = Inf)   # dur & num scales

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
