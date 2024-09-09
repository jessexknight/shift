
# =============================================================================
# config

.p.vars = c('case','seed')
.i.vars = c('i','t.born','age.act','ptr.max')

# =============================================================================
# survey funs

srv.apply = function(Ms,t,srvs=c(srv.base),p.vars=NULL,i.vars=NULL){
  # apply 1+ surveys (srvs) to sim outputs (Ms) at time (t)
  status(3,'srv.apply: ',len(Ms))
  if (missing(t)){ t = Ms[[1]]$P$tf }
  Ps = lapply(Ms,`[[`,'P')
  Es = lapply(Ms,`[[`,'E')
  Qs = lapply(Ms,srv.init,t=t,p.vars=p.vars,i.vars=i.vars)
  for (srv in srvs){
    Qs = par.mapply(srv,Ps,Qs,Es,t) }
  Q = do.call(rbind,Qs)
}

srv.init = function(M,t,p.vars=NULL,i.vars=NULL){
  p.vars = unique(c(.p.vars,p.vars))
  i.vars = unique(c(.i.vars,i.vars))
  Q = cbind(M$P[p.vars],t=t,M$I[i.vars]) # init Q ~= I
  # df.compare(srv.base(M$P,Q,M$E,t),M$I) # DEBUG
}

srv.base = function(P,Q,E,t,fmt='%s'){
  E = clip.evts(E,t=t)
  Q$age      = (t-Q$t.born)/P$t1y
  Q$age.1    = floor(Q$age)
  Q$sex.act  = Q$age > Q$age.act
  Q$vio.nt   = sapply(E$vio,len)
  Q$vio.tr   = sapply(E$vio,last)
  Q$dep.now  = sapply(E$dep_o,len) > sapply(E$dep_x,len)
  Q$dep.past = sapply(E$dep_o,len) > 0
  Q$dep.to   = sapply(E$dep_o,last)
  Q$dep.dur  = ifelse(Q$dep.now,t-Q$dep.to,NA)
  Q$haz.now  = sapply(E$haz_o,len) > sapply(E$haz_x,len)
  Q$haz.past = sapply(E$haz_o,len) > 0
  Q$haz.zo   = sapply(E$haz_o,last)
  Q$haz.dur  = ifelse(Q$haz.now,t-Q$haz.to,NA)
  Q$ptr.nt   = sapply(E$ptr_o,len)
  Q$ptr.nw   = Q$ptr.nt - sapply(E$ptr_x,len)
  names(Q) = sprintf(fmt,names(Q)) # format names for conflicts
  return(Q)
}

# -----------------------------------------------------------------------------

srv.val.RR = function(P,Q,E,t){
  E = clip.evts(E,t=t)
  Q = srv.base(P,Q,E,t)
  Q$age.10 = floor(Q$age/10)*10
  # events in past 1 year
  Q$vio.n1y   = sapply(E$vio,  num.dt,t,P$t1y)
  Q$vio.a1y   = sapply(E$vio,  any.dt,t,P$t1y)
  Q$dep_o.a1y = sapply(E$dep_o,any.dt,t,P$t1y)
  Q$dep_x.a1y = sapply(E$dep_x,any.dt,t,P$t1y)
  Q$haz_o.a1y = sapply(E$haz_o,any.dt,t,P$t1y)
  Q$haz_x.a1y = sapply(E$haz_x,any.dt,t,P$t1y)
  Q$ptr_o.n1y = sapply(E$ptr_o,num.dt,t,P$t1y)
  Q$ptr_x.n1y = sapply(E$ptr_x,num.dt,t,P$t1y)
  # states 1 year prior
  Q = cbind(Q,srv.base(P,Q,E,t-P$t1y,fmt='p1y.%s'))
  Q$p1y.dep.dur.c = int.cut(Q$p1y.dep.dur/P$t1y,c(0,1,5))
  Q$p1y.haz.dur.c = int.cut(Q$p1y.haz.dur/P$t1y,c(0,1,5))
  Q$p1y.vio.nt.c  = int.cut(Q$p1y.vio.nt,c(0,3,30))
  # events in past 3 months
  Q$vio.n3m   = sapply(E$vio,  num.dt,t,P$t3m)
  Q$vio.a3m   = sapply(E$vio,  any.dt,t,P$t3m)
  Q$dep_o.a3m = sapply(E$dep_o,any.dt,t,P$t3m)
  Q$dep_x.a3m = sapply(E$dep_x,any.dt,t,P$t3m)
  Q$haz_o.a3m = sapply(E$haz_o,any.dt,t,P$t3m)
  Q$haz_x.a3m = sapply(E$haz_x,any.dt,t,P$t3m)
  Q$ptr_o.n3m = sapply(E$ptr_o,num.dt,t,P$t3m)
  Q$ptr_x.n3m = sapply(E$ptr_x,num.dt,t,P$t3m)
  # states 3 months prior
  Q = cbind(Q,srv.base(P,Q,E,t-P$t3m,fmt='p3m.%s'))
  Q$p3m.dep.dur.c = int.cut(Q$p3m.dep.dur/P$t1y,c(0,1,5))
  Q$p3m.haz.dur.c = int.cut(Q$p3m.haz.dur/P$t1y,c(0,1,5))
  Q$p3m.vio.nt.c  = int.cut(Q$p3m.vio.nt,c(0,3,30))
  return(Q)
}

# =============================================================================
# rate funs

rate.datas = function(Ms,t,dt=t,...,among=quote(TRUE)){
  status(3,'rate.datas: ',len(Ms))
  if (missing( t)){  t = Ms[[1]]$P$tf }
  Y = rbind.lapply(Ms,rate.data,t=t,...)
  Y = rate.data.sub(Y,t,dt,among=among)
}

rate.data = function(M,t,p.vars=NULL,i.vars=NULL,e.dts=NULL){
  tia = function(i,a){ Q$t.born[i] + M$P$t1y * a } # i age -> time
  Q = srv.init(M,t,p.vars,i.vars)
  Y = rbind.lapply(1:nrow(Q),function(i){
    if (tia(i,amin) > t){ return(NULL) } # (unobserved)
    Ei = ulist(lapply(M$E,`[[`,i),       # add events:
      age  = tia(i,seq(amin,amax)),      # - birthdays
      act  = tia(i,Q$age.act[i]),        # - sexual activity
      tmax = rep(t,tia(i,amax) > t))     # - clip (end obs)
    for (e in names(e.dts)){             # - lagged events
      for (dt in e.dts[[e]]){
        Ei[[str(e,dt,'dt')]] = Ei[[e]]+dt }}
    ti = clip.tes(sort(do.call(c,Ei)),t) # obs event times
    ei = gsub('\\d*$','',names(ti))      # obs event names
    ein = ei[-len(ti)]   # for speed
    Yi = cbind(Q[i,],
      to = ti[-len(ti)], # period start
      tx = ti[-1],       # period end
      e  = ei[-1],       # event at period end
      vio.nt   = cumsum(ein=='vio'),
      dep.now  = cumsum(ein=='dep_o')-cumsum(ein=='dep_x'),
      dep.past = cummax(ein=='dep_o'),
      haz.now  = cumsum(ein=='haz_o')-cumsum(ein=='haz_x'),
      haz.past = cummax(ein=='haz_o'),
      ptr.nw   = cumsum(ein=='ptr_o')-cumsum(ein=='ptr_x'),
      ptr.nt   = cumsum(ein=='ptr_o'),
      age.1    = cumsum(ein=='age')+amin-1,
      sex.act  = cummax(ein=='act'),
    row.names=NULL)
    for (e in names(e.dts)){
      for (dt in e.dts[[e]]){ # pmin because *RR do not stack
        Yi[str(e,dt,'dt')] = pmin(1,cumsum(ein==e)-cumsum(ein==str(e,dt,'dt'))) }}
    return(Yi)
  },.par=FALSE)
  Y$age.10 = floor(Y$age.1/10)*10
  for (e in names(e.dts)){ # e.g. vio.dt: periods with vio in past (30,90,...) days
    cols = str(e,e.dts[[e]],'dt')
    Y[str(e,'.dt')] = factor(rowSums(Y[cols]),len(cols):0,c(sort(e.dts[[e]]),'NR'))
    Y[cols] = NULL }
  # df.compare(subset(Y,e=='tmax'),srv.base(M$P,Q,M$E,t=t)) # DEBUG
  return(Y)
}

rate.data.sub = function(Y,t,dt=t,among=quote(TRUE)){
  tx.w = t    # obs end
  to.w = t-dt # obs start
  Y = subset(Y, to <= tx.w & tx >= to.w) # observed
  Y$e [Y$tx > tx.w] = 'tmax' # clip event
  Y$tx[Y$tx > tx.w] = tx.w   # clip end
  Y$to[Y$to < to.w] = to.w   # clip start
  Y = subset(Y,among) # any other subset
}

rate.est = function(Y,e,strat='seed'){
  Y = switch(e,
    vio = Y,
    dep_o = subset(Y,dep.now==0),
    dep_x = subset(Y,dep.now==1),
    haz_o = subset(Y,haz.now==0),
    haz_x = subset(Y,haz.now==1),
    ptr_o = subset(Y,sex.act & ptr.nw < ptr.max),
    ptr_x = subset(Y,sex.act & ptr.nw > 0))
  y.split = split(1:nrow(Y),Y[strat])
  R = rbind.lapply(y.split,function(y){
    ne = sum(Y$e[y]==e)
    dt = sum(Y$tx[y]-Y$to[y])
    cbind(Y[y[1],strat,drop=FALSE],
      event=e,ne=ne,dt=dt,R=ne/dt,
      # poisson 95% CI
      R.lo=qchisq(.025,2*ne  )/dt/2,
      R.hi=qchisq(.975,2*ne+2)/dt/2)
  })
}

# =============================================================================
# question / event funs

clip.evts = function(E,t){ E = lapply(E,lapply,clip.tes,t=t) }

clip.tes = function(tes,t){ tes[tes <= t] }

num.dt = function(tes,t,dt){ n = sum(tes <= t & tes > t-dt) }

any.dt = function(tes,t,dt){ b = num.dt(tes,t,dt) > 0 }
