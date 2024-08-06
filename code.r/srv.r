
# =============================================================================
# config

p.vars = c('case','seed')
i.cols = c(
  'i',
  'z.born','age.act',
  'vio.Ri',
  'dep_o.Ri','dep_x.Ri',
  'haz_o.Ri','haz_x.Ri',
  'ptr_o.Ri','ptr_x.Ri','ptr.max'
)

# =============================================================================
# survey funs

srv.apply = function(Ms,z,srvs=c(srv.base)){
  # apply 1+ surveys (srvs) to sim outputs (Ms) at time (z)
  if (missing(z)){ z = Ms[[1]]$P$zf }
  Ps = lapply(Ms,`[[`,'P')
  Es = lapply(Ms,`[[`,'E')
  Qs = lapply(Ms,srv.init,z=z)
  for (srv in srvs){
    Qs = par.mapply(srv,Ps,Qs,Es,z) }
  Q = do.call(rbind,Qs)
}

srv.init = function(M,z){
  Q = cbind(M$P[p.vars],z=z,M$I[i.cols]) # init Q ~= I
  # Q = srv.base(M$P,Q,M$E,z); v = intersect(names(M$I),names(Q)) # DEBUG
  # print(all.equal(M$I[v],Q[v])) # DEBUG
}

srv.base = function(P,Q,E,z,fmt='%s'){
  E = lapply(E,lapply,clip.zes,z=z) # clip events
  Q$age      = (z-Q$z.born)/z1y
  Q$sex.act  = Q$age > Q$age.act
  Q$vio.n    = sapply(E$vio,len)
  Q$vio.nf   = sapply(E$vio,last)
  Q$dep.now  = sapply(E$dep_o,len) > sapply(E$dep_x,len)
  Q$dep.past = sapply(E$dep_o,len) > 0
  Q$dep.zo   = sapply(E$dep_o,last)
  Q$dep.u    = ifelse(Q$dep.now,z+1-Q$dep.zo,NA)
  Q$haz.now  = sapply(E$haz_o,len) > sapply(E$haz_x,len)
  Q$haz.past = sapply(E$haz_o,len) > 0
  Q$haz.zo   = sapply(E$haz_o,last)
  Q$haz.u    = ifelse(Q$haz.now,z+1-Q$haz.zo,NA)
  Q$ptr.tot  = sapply(E$ptr_o,len)
  Q$ptr.n    = Q$ptr.tot - sapply(E$ptr_x,len)
  names(Q) = sprintf(fmt,names(Q)) # format names for conflicts
  return(Q)
}

# -----------------------------------------------------------------------------

srv.val = function(P,Q,E,z){
  E = lapply(E,lapply,clip.zes,z=z) # clip events
  Q = srv.base(P,Q,E,z)
  Q$age.10 = floor(Q$age/10)*10
  # events in past 1 year
  Q$vio.n1y   = sapply(E$vio,num.dz,z,z1y)
  Q$vio.a1y   = sapply(E$vio,any.dz,z,z1y)
  Q$dep_o.a1y = sapply(E$dep_o,any.dz,z,z1y)
  Q$dep_x.a1y = sapply(E$dep_x,any.dz,z,z1y)
  Q$haz_o.a1y = sapply(E$haz_o,any.dz,z,z1y)
  Q$haz_x.a1y = sapply(E$haz_x,any.dz,z,z1y)
  Q$ptr_o.n1y = sapply(E$ptr_o,num.dz,z,z1y)
  Q$ptr_x.n1y = sapply(E$ptr_x,num.dz,z,z1y)
  # states 1 year prior
  Q = cbind(Q,srv.base(P,Q,E,z-z1y,fmt='yp.%s'))
  Q$yp.dep.u.c = int.cut(Q$yp.dep.u,z1y*c(0,1,10))
  Q$yp.haz.u.c = int.cut(Q$yp.haz.u,z1y*c(0,1,10))
  Q$yp.vio.n.c = int.cut(Q$yp.vio.n,c(0,1,10))
  return(Q)
}

# =============================================================================
# question / event funs

clip.zes = function(zes,z){ zes[zes <= z] }

num.dz = function(zes,z,dz){ n = sum(zes <= z & zes >= z+1-dz) }

any.dz = function(zes,z,dz){ b = num.dz(zes,z,dz) > 0 }
