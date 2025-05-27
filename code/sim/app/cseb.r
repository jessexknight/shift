source('sim/meta.r')
source('sim/mass.r')
source('sim/fit.r')
uid = '2025-05-25'

# -----------------------------------------------------------------------------
# targets / outcomes

T = name.list(key='id',
  gen.targ(id='dep.now',    type='prop',vo='dep.now'),
  gen.targ(id='dep.past',   type='prop',vo='dep.past'),
  gen.targ(id='haz.now',    type='prop',vo='haz.now'),
  gen.targ(id='haz.past',   type='prop',vo='haz.past'),
  gen.targ(id='dep.haz.or', type='OR',  ve='dep.now',vo='haz.now',ao1=FALSE),
  gen.targ(id='dep.haz.pr', type='PR',  ve='dep.now',vo='haz.now',ao1=FALSE),
  gen.targ(id='dep.haz.aor',type='OR',  ve='dep.now',vo='haz.now',ao1=FALSE,va1='age'),
  gen.targ(id='dep.haz.apr',type='PR',  ve='dep.now',vo='haz.now',ao1=FALSE,va1='age'))

# -----------------------------------------------------------------------------
# params & grid

P0 = list(
  dtz   =   cli.arg('dtz',     45), # final: 7
  n.pop =   cli.arg('n.pop', 1000), # final: 1000
  seed  = 1:cli.arg('n.seed',  21), # final: 100
  n.dur = 1,
  null  = 'xRR',
  het.distr = 'lnorm',
  dep_o.Ri.my = .02, haz_o.Ri.my = .02,
  dep_x.Ri.my = 1,   haz_x.Ri.my = 1,
  dep.Ri.cv   = 0,   haz.Ri.cv   = 0,
  dep.cov     = 0,   haz.cov     = 0,
  RR.haz_o.dep_w = 1,
  RR.haz_x.dep_w = 1,
  run = get.run.par(c('dep','haz'),u=FALSE))

PG = list(
  dep_o.Ri.my = c(.02,.04,.06), haz_o.Ri.my = c(.02,.04,.06),
  dep_x.Ri.my = c(1,2,3),       haz_x.Ri.my = c(1,2,3),
  dep.Ri.cv   = c(0,1,3),       haz.Ri.cv   = c(0,1,3),
  dep.cov     = c(-.5,0,+.5),   haz.cov     = c(-.5,0,+.5),
  RR.haz_o.dep_w = signif(2^seq( 0,+3,.5),3),
  RR.haz_x.dep_w = signif(2^seq(-3, 0,.5),3))

PGk = list(
  ref = PG[1:8],         # RR = 1
  fix = PG[9:10],        # Ri* ~ fixed
  hom = PG[c(1:4,9:10)], # Ri* ~ homog
  unc = PG[c(1:6,9:10)], # Ri* ~ heter + uncor
  het = PG)              # Ri* ~ heter + cor

grid.path = function(k,.save=FALSE){
  hash.path(ulist(P0,PGk[[k]]),'data','sim','cseb',uid,.save=.save)
}

# -----------------------------------------------------------------------------
# run & save/load

run.grid = function(k){
  Y = fit.run.grid(PGk[[k]],T,P0)
  save.rda(Y,grid.path(k,.save=TRUE),'Y')
}

load.grid = function(k,f=NULL){
  Y = load.rda(grid.path(k),'Y')
  iP = Y$type == 'prop' # prop out rows
  c3 = c('value','lower','upper') # out value cols
  Y[iP,c3] = Y[iP,c3]*100 # props as %
  y = function(k){ if.null(Y[[k]],P0[[k]]) }
  Y$d.Ro = round(y('dep_o.Ri.my')*100,1); Y$h.Ro = round(y('haz_o.Ri.my')*100,1) # per 100 PY
  Y$d.Rx = round(y('dep_x.Ri.my')*100,1); Y$h.Rx = round(y('haz_x.Ri.my')*100,1) # per 100 PY
  Y$d.het = y('dep.Ri.cv');               Y$h.het = y('haz.Ri.cv') # shorthand
  Y$d.cor = y('dep.cov');                 Y$h.cor = y('haz.cov')   # shorthand
  Y$RRo   = y('RR.haz_o.dep_w'); # shorthand
  Y$RRx   = y('RR.haz_x.dep_w'); # shorthand
  Y[f] = lapply(Y[f],as.factor) # Y[f] -> factors
  return(Y)
}

# -----------------------------------------------------------------------------
# plotting

# TODO

# -----------------------------------------------------------------------------
# main

# run.grid('ref')
# run.grid('fix')
# run.grid('hom')
# run.grid('unc')
# run.grid('het')
