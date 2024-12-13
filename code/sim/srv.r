
# =============================================================================
# config

.p.vars = c('case','seed')
.i.vars = c('i','t.born','age.act','ptr.max')

# =============================================================================
# survey funs

srv.apply = function(Ms,t,srvs=c(srv.base),p.vars=NULL,i.vars=NULL,x.cols=NULL){
  # apply 1+ surveys (srvs) to sim outputs (Ms) at time (t)
  status(3,'srv.apply: ',len(Ms))
  if (missing(t)){ t = Ms[[1]]$P$tf }
  if (t <= 0){ t = Ms[[1]]$P$tf + t }
  Ps = lapply(Ms,`[[`,'P')
  Es = lapply(Ms,`[[`,'E')
  Qs = lapply(Ms,srv.init,t=t,p.vars=p.vars,i.vars=i.vars); gc()
  for (srv in srvs){
    Qs = par.mapply(srv,Ps,Qs,Es,t) }; gc()
  Q = do.call(rbind,Qs); gc()
  for (x in names(x.cols)){
    Q[[x]] = x.cols[[x]](Q) }
  return(Q)
}

srv.init = function(M,t,p.vars=NULL,i.vars=NULL){
  p.vars = unique(c(.p.vars,p.vars))
  i.vars = unique(c(.i.vars,i.vars))
  Q = cbind(M$P[p.vars],.=1,t=t,M$I[i.vars]) # init Q ~= I
  # df.compare(srv.base(M$P,Q,M$E,t),M$I) # DEBUG
}

srv.base = function(P,Q,E,t){
  E = clip.evts(E,t=t)
  Q$age      = (t - Q$t.born)/P$t1y
  Q$age.1    = floor(Q$age)
  Q$age.10   = int.cut(Q$age,seq(10,50,10))
  Q$sex.act  = Q$age > Q$age.act
  Q$vio.nt   = sapply(E$vio,len)
  Q$vio.tr   = sapply(E$vio,last)
  Q$vio.dt   = t - Q$vio.tr
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
  return(Q)
}

srv.e.dts = function(P,Q,E,t,e.dts){
  E = clip.evts(E,t=t)
  for (e in names(e.dts)){
    dt.cols = sapply(e.dts[[e]],function(dt){ sapply(E[[e]],any.dt,t,dt) })
    Q[str(e,'.dt.c')] = factor(rowSums(dt.cols),len(e.dts[[e]]):0,c(sort(e.dts[[e]]),'NR'))
  } # TODO: ^ make common fun w/ rate.data?
  # TODO: verify w/ df.compare
  return(Q)
}

# =============================================================================
# time vector summary

vec.datas = function(Ms,...){
  status(3,'vec.datas: ',len(Ms))
  V = rbind.lapply(Ms,vec.data,...)
  # Q = srv.apply(lapply(Ms,sim.sub,sub='act'))[c('seed',names(V)[ncol(V)-1:7+1])] # DEBUG
  # df.compare(V[V$t==max(V$t),],aggregate(.~seed,Q,sum)) # DEBUG
}

vec.data = function(M,strat='.',frame='g',p.vars=NULL){
  # TODO: might have off-by-one error
  p.vars = unique(c(.p.vars,p.vars,intersect(strat,names(M$P))))
  i.vars = intersect(strat,names(M$I))
  Q = cbind(M$P[p.vars],.=1,M$I[i.vars])
  t.enter = floor(M$I$t.born) + amin * M$P$t1y # enter model
  t.exit  = floor(M$I$t.born) + amax * M$P$t1y # exit model
  M$E$t0  = pmax(     1,t.enter) # left-clip enter
  M$E$tf  = pmin(M$P$tf,t.exit ) # right-clip exit
  M$E$null = lapply(1:nrow(Q),function(i){ numeric() }) # dummy
  tff = switch(frame,g=M$P$tf,i=M$P$t1y*adur) # global vs ind time
  if (frame=='i'){ M$E = lapply(M$E,function(es){ wapply(`-`,es,t.enter-1) }) }
  # print(lapply(M$E[c('t0','tf')],function(e){ summary(unlist(e)) })) # DEBUG
  i.strat = fast.split(1:nrow(Q),Q[strat])
  V = rbind.lapply(i.strat,function(i){ # for each strat
    Ei  = lapply(M$E,`[`,i) # select inds
    tfi = unlist(Ei$tf) + 1
    Vi  = data.frame(Q[i[1],],
      t = 1:tff,
      n        = tox.vec.i(Ei$t0,   tfi,     tfi,tff),
      vio.nt   = tox.vec.i(Ei$vio,  Ei$null, tfi,tff),
      dep.now  = tox.vec.i(Ei$dep_o,Ei$dep_x,tfi,tff),
      dep.past = tox.vec.i(Ei$dep_o,Ei$null, tfi,tff,t1=TRUE),
      haz.now  = tox.vec.i(Ei$haz_o,Ei$haz_x,tfi,tff),
      haz.past = tox.vec.i(Ei$haz_o,Ei$null, tfi,tff,t1=TRUE),
      ptr.nw   = tox.vec.i(Ei$ptr_o,Ei$ptr_x,tfi,tff),
      ptr.nt   = tox.vec.i(Ei$ptr_o,Ei$null, tfi,tff),
    row.names=NULL)
  })
}

tox.vec.i = function(toi,txi,tfi,tff,t1=FALSE){
  # cumsum Vi(t) = oi(t)-xi(t) on [1,tff] ensuring Vi(tfi) = 0
  # i.e. each individual's Vi(t) returns to zero upon model exit
  if (t1){ toi = ti.1(toi); txi = ti.1(txi) } # only oi[1],xi[1]
  tox.vec(to=unlist(toi),tx=unlist(wapply(tx.pad,toi,txi,tfi)),tf=tff)
}

ti.1 = function(ti){ ifelse(lens(ti),lapply(ti,`[`,1),ti) }
tx.pad = function(to,tx,tf){ c(tx,rep(tf,len(to)-len(tx))) }
tox.vec = function(to,tx,tf){ cumsum(tabulate(to,tf)-tabulate(tx,tf)) }

# =============================================================================
# rate funs

rate.datas = function(Ms,t,dt=t,...,among=quote(TRUE)){
  status(3,'rate.datas: ',len(Ms))
  if (missing(t)){ t = Ms[[1]]$P$tf }
  Y = rbind.lapply(Ms,rate.data,t=t,...); status(4,'\n')
  Y = rate.data.sub(Y,t,dt,among=among)
}

rate.data = function(M,t,p.vars=NULL,i.vars=NULL,e.dts=NULL,x.cols=NULL){
  # TODO: think carefully about if/where dtz / 2 is needed
  status(4,id=M$P$seed)
  tia = function(i,a){ Q$t.born[i] + M$P$t1y * a } # i age -> time
  Q = srv.init(M,t,p.vars,i.vars)
  Y = rbind.lapply(1:nrow(Q),function(i){
    if (tia(i,amin) > t){ return(NULL) } # (unobserved)
    Ei = ulist(lapply(M$E,`[[`,i),       # add events:
      age  = tia(i,seq(amin,amax)),      # - birthdays
      act  = tia(i,Q$age.act[i]),        # - sexual activity
      tmax = rep(t,tia(i,amax) > t))     # - clip (end obs)
    for (e in names(e.dts)){             # - lagged events
      Ei.dts = set.names(lapply(e.dts[[e]],`+`,Ei[[e]]),str(e,e.dts[[e]],'dt'))
      Ei = append(Ei,Ei.dts,which(names(Ei)==e)) }
    ti = clip.tes(sort(do.call(c,Ei)),t) # obs event times
    ei = gsub('\\d*$','',names(ti))      # obs event names
    ein = ei[-len(ei)]   # for speed
    Yi = cbind(Q[i,],
      to = ti[-len(ti)], # period start
      tx = ti[-1],       # period end
      e  = ei[-1],       # event at period end
      age.1    = cumsum(ein=='age')+amin-1,
      sex.act  = cummax(ein=='act'),
      vio.nt   = cumsum(ein=='vio'),
      dep.now  = cumsum(ein=='dep_o')-cumsum(ein=='dep_x'),
      dep.past = cummax(ein=='dep_o'),
      haz.now  = cumsum(ein=='haz_o')-cumsum(ein=='haz_x'),
      haz.past = cummax(ein=='haz_o'),
      ptr.nw   = cumsum(ein=='ptr_o')-cumsum(ein=='ptr_x'),
      ptr.nt   = cumsum(ein=='ptr_o'),
    row.names=NULL)
    for (e in names(e.dts)){
      for (dt in e.dts[[e]]){ # pmin because *RR do not stack
        Yi[str(e,dt,'dt')] = pmin(1,cumsum(ein==e)-cumsum(ein==str(e,dt,'dt'))) }}
    return(Yi)
  },.par=FALSE)
  for (e in names(e.dts)){ # e.g. vio.dt: periods with vio in past (30,90,...) days
    cols = str(e,e.dts[[e]],'dt')
    Y[str(e,'.dt.c')] = factor(rowSums(Y[cols]),len(cols):0,c(sort(e.dts[[e]]),'NR'))
    Y[cols] = NULL }
  for (x in names(x.cols)){
    Y[[x]] = x.cols[[x]](Y) }
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
  Y$w = switch(e, # person-time weight
    vio   = 1,
    dep_o = Y$dep.now==0,
    dep_x = Y$dep.now==1,
    haz_o = Y$haz.now==0,
    haz_x = Y$haz.now==1,
    ptr_o = Y$sex.act & Y$ptr.nw < Y$ptr.max,
    ptr_x = Y$ptr.nw)
  y.strat = fast.split(1:nrow(Y),Y[strat])
  R = rbind.lapply(y.strat,function(y){
    ne = sum(Y$e[y]==e)
    dt = sum((Y$tx[y] - Y$to[y]) * Y$w[y])
    cbind(Y[y[1],strat,drop=FALSE],
      var=e,ne=ne,dt=dt,value=ne/dt,
      # poisson 95% CI
      value.lo=qchisq(.025,2*ne  )/dt/2,
      value.hi=qchisq(.975,2*ne+2)/dt/2)
  })
}

# =============================================================================
# question / event funs

clip.evts = function(E,t){ E = lapply(E,lapply,clip.tes,t=t) }

clip.tes = function(tes,t){ tes[tes <= t] }

num.dt = function(tes,t,dt){ n = sum(tes <= t & tes > t-dt) }

any.dt = function(tes,t,dt){ b = num.dt(tes,t,dt) > 0 }
