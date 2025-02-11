
# -----------------------------------------------------------------------------
# param sampling

tx.sample = function(F,S0,dir='tx'){
  # forward (tx) or reverse (itx) transform S0 <-> S
  # where S0 ~ [0,1] are quantiles of S ~ tx(unif[lo,up])
  S = as.data.frame(lapply(seqa(F),function(i){
    Fi  = F[[i]]
    fun = if.null(Fi[[dir]],identity)
    Si  = switch(dir,
      tx  = fun(qunif(S0[,i],Fi$lo,Fi$up)),
      itx = punif(fun(S0[,i]),Fi$lo,Fi$up))
  }),col.names=names(F))
}

lhs.sample = function(F,n){
  S = tx.sample(F,lhs::randomLHS(n,len(F)))
}

plot.sample = function(S,color=''){
  # ggpairs plot & option to color by a factor var e.g. rank.cut(L$total)
  stfu(library('GGally',quietly=TRUE))
  grid.draw.ggmatrix <<- print
  g = ggpairs(cbind(S,color=color),map=aes(color=color,alpha=0))
  for (i in 1:g$nrow){ for (j in 1:g$ncol){ g[i,j] = g[i,j] + clr.map.d }}
  g = plot.clean(g)
}

# -----------------------------------------------------------------------------
# log-likelihoods

ll.srv = function(Q,T){
  # weighted log-likelihoods for each target T given survey Q
  # e.g. Ti$type = 'prop' uses ll.prop(Ti$t.arg,Q=Q,...)
  Ls = unlist(lapply(T,function(Ti){
    fun  = get(str('ll.',Ti$type))
    args = ulist(within(Ti,rm(type,w)),Q=Q)
    Li   = do.call(fun,args) * Ti$w
  }))
}

ll.prop = function(t.arg,Q,v){
  x = Q[[v]] # extract data (survey responses)
  p.x = mean(x); n.x = len(x);  k.x = sum(x)  # data
  p.t = t.arg$p; n.t = t.arg$n; k.t = p.t*n.t # target
  p = (k.t+k.x) / (n.t+n.x)
  z = (p.t-p.x) / sqrt(p*(1-p)*(1/n.t+1/n.x))
  ll = dnorm(z,log=TRUE)
}

# ll.pois = function(t.arg,Q,v){} # TODO

ll.OR = function(t.arg,Q,...){
  A = mass.calc(mfuns=mass.funs.ce['OR'],Q1=Q,...)
  or.x = mean(A$coef)
  se.x = mean(A$std.err)
  or.t = t.arg$OR; s = sqrt(sqrt(or.t))
  se.t = (s+1/s)*2 / sqrt(t.arg$n) # HACK
  z  = (or.t-or.x) / sqrt(se.t^2+se.x^2)
  ll = dnorm(z,log=TRUE)
}

# -----------------------------------------------------------------------------
# fitting

fit.run = function(Si,T,P0=NULL,...,aggr=TRUE,.par=FALSE){
  # get log-likelihood given sample Si, targets T, base params P0
  # i.e. get params, run model, run survey, get log-likelihoods
  Ps = get.pars.grid(ulist(P0,Si),...)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  Q  = srv.apply(Ms)
  Ls = ll.srv(Q,T)
  Ls = c(total=sum(Ls),Ls)
  status(2,'fit.run:',
    '\n > L: ',list.str(Ls,join=', ',sig=3),
    '\n > S: ',list.str(Si,join=', ',sig=3))
  if (aggr){ return(Ls[1]) } else { return(Ls) }
}

fit.runs = function(S,T,...){ .verb <<- 2 # HACK
  # fit.run in parallel for each sample (row) in data.frame S
  L = rbind.lapply(apply(S,1,as.list),fit.run,T=T,...)
}
