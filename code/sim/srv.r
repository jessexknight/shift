
# =============================================================================
# config

.p.vars = c('case','seed')
.i.vars = c(
  'i',
  'z.born','age.act',
  'vio.Ri',
  'dep_o.Ri','dep_x.Ri',
  'haz_o.Ri','haz_x.Ri',
  'ptr_o.Ri','ptr_x.Ri','ptr.max'
)

# =============================================================================
# survey funs

srv.apply = function(Ms,z,srvs=c(srv.base),p.vars=NULL,i.vars=NULL){
  # apply 1+ surveys (srvs) to sim outputs (Ms) at time (z)
  status(3,'srv.apply: ',len(Ms))
  if (missing(z)){ z = Ms[[1]]$P$zf }
  Ps = lapply(Ms,`[[`,'P')
  Es = lapply(Ms,`[[`,'E')
  Qs = lapply(Ms,srv.init,z=z,p.vars=p.vars,i.vars=i.vars)
  for (srv in srvs){
    Qs = par.mapply(srv,Ps,Qs,Es,z) }
  Q = do.call(rbind,Qs)
}

srv.init = function(M,z,p.vars,i.vars){
  p.vars = unique(c(.p.vars,p.vars))
  i.vars = unique(c(.i.vars,i.vars))
  Q = cbind(M$P[p.vars],z=z,M$I[i.vars]) # init Q ~= I
  # Q = srv.base(M$P,Q,M$E,z); print(all.equal(M$I[i.vars],Q[i.vars])) # DEBUG
}

srv.base = function(P,Q,E,z,fmt='%s'){
  E = lapply(E,lapply,clip.zes,z=z) # clip events
  Q$age      = (z-Q$z.born)/z1y
  Q$sex.act  = Q$age > Q$age.act
  Q$vio.nt   = sapply(E$vio,len)
  Q$vio.zf   = sapply(E$vio,last)
  Q$dep.now  = sapply(E$dep_o,len) > sapply(E$dep_x,len)
  Q$dep.past = sapply(E$dep_o,len) > 0
  Q$dep.zo   = sapply(E$dep_o,last)
  Q$dep.dur  = ifelse(Q$dep.now,z+1-Q$dep.zo,NA)
  Q$haz.now  = sapply(E$haz_o,len) > sapply(E$haz_x,len)
  Q$haz.past = sapply(E$haz_o,len) > 0
  Q$haz.zo   = sapply(E$haz_o,last)
  Q$haz.dur  = ifelse(Q$haz.now,z+1-Q$haz.zo,NA)
  Q$ptr.nt   = sapply(E$ptr_o,len)
  Q$ptr.nw   = Q$ptr.nt - sapply(E$ptr_x,len)
  names(Q) = sprintf(fmt,names(Q)) # format names for conflicts
  return(Q)
}

srv.ptr = function(P,Q,E,z){
  E = lapply(E,lapply,clip.zes,z=z) # clip events
  Q$ptr.nw.c = int.cut(Q$ptr.nw,c(0,1,2))
  ptr.w.dur = wapply(dur.ptr.w,E$ptr_o,E$ptr_x,E$ptr_u,z) # ongo ptrs durs
  Q$ptr.p.dur.m = dtz*sapply(E$ptr_u,mean)   # past ptrs durs
  Q$ptr.w.dur.m = dtz*sapply(ptr.w.dur,mean) # ongo ptrs durs
  return(Q)
}

# -----------------------------------------------------------------------------

srv.val.RR = function(P,Q,E,z){
  E = lapply(E,lapply,clip.zes,z=z) # clip events
  Q = srv.base(P,Q,E,z)
  Q$age.10 = floor(Q$age/10)*10
  # events in past 1 year
  Q$vio.n1y   = sapply(E$vio,  num.dz,z,z1y)
  Q$vio.a1y   = sapply(E$vio,  any.dz,z,z1y)
  Q$dep_o.a1y = sapply(E$dep_o,any.dz,z,z1y)
  Q$dep_x.a1y = sapply(E$dep_x,any.dz,z,z1y)
  Q$haz_o.a1y = sapply(E$haz_o,any.dz,z,z1y)
  Q$haz_x.a1y = sapply(E$haz_x,any.dz,z,z1y)
  Q$ptr_o.n1y = sapply(E$ptr_o,num.dz,z,z1y)
  Q$ptr_x.n1y = sapply(E$ptr_x,num.dz,z,z1y)
  # states 1 year prior
  Q = cbind(Q,srv.base(P,Q,E,z-z1y,fmt='p1y.%s'))
  Q$p1y.dep.dur.c = int.cut(Q$p1y.dep.dur,z1y*c(0,1,10))
  Q$p1y.haz.dur.c = int.cut(Q$p1y.haz.dur,z1y*c(0,1,10))
  Q$p1y.vio.nt.c  = int.cut(Q$p1y.vio.nt,c(0,1,10,100))
  # events in past 3 months
  Q$vio.n3m   = sapply(E$vio,  num.dz,z,z3m)
  Q$vio.a3m   = sapply(E$vio,  any.dz,z,z3m)
  Q$dep_o.a3m = sapply(E$dep_o,any.dz,z,z3m)
  Q$dep_x.a3m = sapply(E$dep_x,any.dz,z,z3m)
  Q$haz_o.a3m = sapply(E$haz_o,any.dz,z,z3m)
  Q$haz_x.a3m = sapply(E$haz_x,any.dz,z,z3m)
  Q$ptr_o.n3m = sapply(E$ptr_o,num.dz,z,z3m)
  Q$ptr_x.n3m = sapply(E$ptr_x,num.dz,z,z3m)
  # states 3 months prior
  Q = cbind(Q,srv.base(P,Q,E,z-z3m,fmt='p3m.%s'))
  Q$p3m.dep.dur.c = int.cut(Q$p3m.dep.dur,z3m*c(0,1,10))
  Q$p3m.haz.dur.c = int.cut(Q$p3m.haz.dur,z3m*c(0,1,10))
  Q$p3m.vio.nt.c  = int.cut(Q$p3m.vio.nt,c(0,1,10,100))
  return(Q)
}

# =============================================================================
# question / event funs

clip.zes = function(zes,z){ zes[zes <= z] }

num.dz = function(zes,z,dz){ n = sum(zes <= z & zes >= z+1-dz) }

any.dz = function(zes,z,dz){ b = num.dz(zes,z,dz) > 0 }

dur.ptr.w = function(zos,zxs,us,z){ z - zos[!(zos %in% (zxs-us))] }