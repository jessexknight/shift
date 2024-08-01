
# =============================================================================
# config

p.vars = c('case','seed')
i.cols = c(
  'i',
  'age.act',
  'vio.Ri',
  'dep_o.Ri','dep_x.Ri',
  'haz_o.Ri','haz_x.Ri',
  'ptr_o.Ri','ptr_x.Ri','ptr.max'
)

# =============================================================================
# survey funs

srv.map = function(Ss,z,srvs=NULL){
  # run survey(s) on Ss at timestep z, return joined Qs
  if (missing(z)){ z = Ss[[1]]$P$zf }
  Zs = par.lapply(Ss,srv.base,z=z) # base survey (always)
  for (srv in srvs){               # other surveys (optional)
    Zs = par.lapply(Zs,srv,z=z) }
  Qs = rbind.lapply(Zs,`[[`,'Qs')  # extract Qs
}

srv.base = function(S,z){
  # base survey: recreate Qs = Is at time z <= P$zf using Es
  attach(S) # -> P, Is, Es
  Es = lapply(Es,lapply,function(zes){ zes[zes <= z] }) # clip events
  Qs = cbind(P[p.vars],z=z,Is[i.cols]) # init Qs ~= Is
  Qs$age      = Is$age - (P$zf-z)/364
  Qs$sex.act  = Qs$age > Qs$age.act
  Qs$vio.n    = sapply(Es$vio,len)
  Qs$vio.zf   = sapply(Es$vio,last)
  n_o = sapply(Es$dep_o,len)
  Qs$dep.now  = n_o > sapply(Es$dep_x,len)
  Qs$dep.past = n_o > 0
  Qs$dep.zo   = sapply(Es$dep_o,last)
  Qs$dep.u    = ifelse(Qs$dep.now,z+1-Qs$dep.zo,NA)
  n_o = sapply(Es$haz_o,len)
  Qs$haz.now  = n_o > sapply(Es$haz_x,len)
  Qs$haz.past = n_o > 0
  Qs$haz.zo   = sapply(Es$haz_o,last)
  Qs$haz.u    = ifelse(Qs$haz.now,z+1-Qs$haz.zo,NA)
  Qs$ptr.tot  = sapply(Es$ptr_o,len) # lifetime ptrs
  Qs$ptr.n    = Qs$ptr.tot - sapply(Es$ptr_x,len) # current ptrs
  # v = intersect(names(Is),names(Qs)) # DEBUG (z = P$zf)
  # print(all(Is[v]==Qs[v] | is.na(Is[v]))) # DEBUG
  Z = list(P=P,Qs=Qs,Es=Es) # sim-like output
}
