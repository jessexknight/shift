source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')
uid = '2026-03-03'
.debug = cli.arg('.debug',0)
.lhs = cli.arg('.lhs',xdf(1e3,1e5))
.k  = cli.arg('.k','RRo')
.nb = cli.arg('.nb',1)
.b  = cli.arg('.b', 1)

# -----------------------------------------------------------------------------
# params

P0 = list(
  pop.type = 'open',
  n.pop = xdf(1000,10000),
  n.dur = 1, dtz = 7,
  het.distr = 'gamma',
  run = get.run.par(c('dep','haz'),u=0))

PF = name.list(key='i',
  list(i='RRo', id='RR.haz_o.dep_w',def=  0,  lo=  0,up=  +1,tr='e10'),
  list(i='RRx', id='RR.haz_x.dep_w',def=  0,  lo= -1,up=   0,tr='e10'),
  list(i='eRo', id='dep_o.Ri.my',   def= -1.5,lo= -3,up= -.5,tr='e10'),
  list(i='oRo', id='haz_o.Ri.my',   def= -1.5,lo= -3,up= -.5,tr='e10'),
  list(i='eRx', id='dep_x.Ri.my',   def= -0.5,lo= -2,up= +.5,tr='e10'),
  list(i='oRx', id='haz_x.Ri.my',   def= -0.5,lo= -2,up= +.5,tr='e10'),
  list(i='eHo', id='dep_o.Ri.het',  def= -3,  lo= -3,up= +.5,tr='e10'),
  list(i='oHo', id='haz_o.Ri.het',  def= -3,  lo= -3,up= +.5,tr='e10'),
  list(i='eHx', id='dep_x.Ri.het',  def= -3,  lo= -3,up= +.5,tr='e10'),
  list(i='oHx', id='haz_x.Ri.het',  def= -3,  lo= -3,up= +.5,tr='e10'),
  list(i='ecv', id='dep.cov',       def=  0,  lo= -1,up= +1, tr='identity'),
  list(i='ocv', id='haz.cov',       def=  0,  lo= -1,up= +1, tr='identity'),
  list(i='seed',id='seed',          def=  1,  lo=  1,up=1e9, tr='identity'))

PFk = list(
  RRo  = PF[-2],
  RRx  = PF[-1],
  full = PF)

as.fpar = function(k){
  Fk = do.call(fpar.set,lapply(PFk[[k]],function(P){
    gen.fpar(id=P$id,def=P$def,lo=P$lo,up=P$up,trafo=str.fun(P$tr))
  }))
}

# -----------------------------------------------------------------------------
# targets / outcomes

T = name.list(key='id',
  gen.targ(id='e.w',   type='prop',vo='dep.now' ),
  gen.targ(id='e.p',   type='prop',vo='dep.past'),
  gen.targ(id='o.w',   type='prop',vo='haz.now' ),
  gen.targ(id='o.p',   type='prop',vo='haz.past'),
  gen.targ(id='or.ww', type='OR',  ve='dep.now', vo='haz.now' ),
  gen.targ(id='or.wp', type='OR',  ve='dep.now', vo='haz.past'),
  gen.targ(id='or.pw', type='OR',  ve='dep.past',vo='haz.now' ),
  gen.targ(id='or.pp', type='OR',  ve='dep.past',vo='haz.past'),
  gen.targ(id='pr.ww', type='PR',  ve='dep.now', vo='haz.now' ),
  gen.targ(id='pr.pp', type='PR',  ve='dep.now', vo='haz.past'),
  gen.targ(id='pr.pw', type='PR',  ve='dep.past',vo='haz.now' ),
  gen.targ(id='pr.wp', type='PR',  ve='dep.past',vo='haz.past'),
  gen.targ(id='or.wwa',type='OR',  ve='dep.now', vo='haz.now' ,va1='age'),
  gen.targ(id='or.wpa',type='OR',  ve='dep.now', vo='haz.past',va1='age'),
  gen.targ(id='or.pwa',type='OR',  ve='dep.past',vo='haz.now' ,va1='age'),
  gen.targ(id='or.ppa',type='OR',  ve='dep.past',vo='haz.past',va1='age'),
  gen.targ(id='pr.wwa',type='PR',  ve='dep.now', vo='haz.now' ,va1='age'),
  gen.targ(id='pr.ppa',type='PR',  ve='dep.now', vo='haz.past',va1='age'),
  gen.targ(id='pr.pwa',type='PR',  ve='dep.past',vo='haz.now' ,va1='age'),
  gen.targ(id='pr.wpa',type='PR',  ve='dep.past',vo='haz.past',va1='age'))

# -----------------------------------------------------------------------------
# run sims & save/load

data.path = function(k,.save=FALSE){
  info = ulist(P0,PFk[[k]],lhs=.lhs)
  path = hash.path(info,'data','sim','mass',uid,str('lhs.',k),.save=.save)
}

run.one = function(...,.par=0){
  Ps = get.pars.grid(ulist(P0,...),.par=.par)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  Q  = srv.apply(Ms,.par=.par)
  Y  = srv.targs(Q,T)
  Y[c('seed','targ.mu','targ.se','ll')] = NULL
  row.names(Y) = NULL
  return(Y)
}

run.lhs = function(k){
  S = fpar.sam(.lhs,as.fpar(k=k),tr=TRUE)
  Y = grid.apply(S,run.one,.grid=0,.rbind=1,.cbind=1,.batch=.b,.nbatch=.nb,.log=3)
  save.rds(Y,data.path(k,.save=1),str('b',.nb),str('Y.',.b))
}

merge.batch = function(k){
  Y = rbind.lapply(1:.nb,function(b){
    Yb = load.rds(data.path(k),str('b',.nb),str('Y.',b)) })
  save.rds(Y,grid.path(k),'Y')
}

# -----------------------------------------------------------------------------
# emulator

# TODO

# -----------------------------------------------------------------------------
# main

run.lhs(.k)
