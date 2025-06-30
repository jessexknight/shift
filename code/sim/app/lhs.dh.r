source('sim/meta.r')
source('sim/fit.r')

# -----------------------------------------------------------------------------
# config
uid   = '2025-06-30'
seed  = cli.arg('seed',    666)
n.sam = cli.arg('n.sam',100000) # final: 100000
.b    = cli.arg('.b', 1)
.nb   = cli.arg('.nb',1)

# -----------------------------------------------------------------------------
# params & priors

P0 = list(
  dtz   = cli.arg('dtz',        45), # final: 7
  n.pop = cli.arg('n.pop',   10000), # final: 10000
  seed  = 1:cli.arg('n.seed',    3), # final: 7
  n.dur = 1,
  null  = 'xRR',
  run = get.run.par(c('dep','haz'),u=FALSE))

F = fpar.set(
  gen.fpar(id='dep_o.Ri.my',   lo=.01,up=.1),
  gen.fpar(id='haz_o.Ri.my',   lo=.01,up=.1),
  gen.fpar(id='dep_x.Ri.my',   lo=.3, up= 3),
  gen.fpar(id='haz_x.Ri.my',   lo=.3, up= 3),
  gen.fpar(id='dep.Ri.het',    lo= 0, up= 5),
  gen.fpar(id='haz.Ri.het',    lo= 0, up= 5),
  gen.fpar(id='dep.cov',       lo=-1, up=+1),
  gen.fpar(id='haz.cov',       lo=-1, up=+1),
  gen.fpar(id='RR.haz_o.dep_w',lo= 1, up=10),
  gen.fpar(id='RR.haz_x.dep_w',lo=.1, up= 1))

data.path = function(.save=FALSE){
  info = list(n.sam=n.sam,seed=seed,P0=P0,F=lapply(F$pars,unclass))
  hash.path(info,'data','sim','lhs',uid,n.sam,.save=.save)
}

# -----------------------------------------------------------------------------
# dep:haz state stuff

ages = amin:(amax-1)
GI = expand.grid(dep.past=0:1,dep.now=0:1,haz.past=0:1,haz.now=0:1)
bi = GI$dep.now <= GI$dep.past & GI$haz.now <= GI$haz.past # possible
li = str('s',apply(GI,1,str,collapse='')) # labels
gc = names(GI) # columns

# -----------------------------------------------------------------------------
# main

run.sam = function(P0=NULL,...,.par=TRUE){
  Ps = get.pars.grid(P0,...,.par=.par)
  Ms = sim.runs(Ps,sub='act',.par=.par)
  H  = rbind.lapply(Ms,function(M){
    # build histogram by age (rows) & state (cols)
    Hi = lapply(split(floor(M$I$age-amin+1),M$I[gc]),tabulate,nbins=adur)
    Hi = data.frame(seed=M$P$seed,age=ages,set.names(Hi,li)[bi])
  },.par=.par)
}

run.lhs = function(){
  S = fpar.sam(F,n=n.sam,seed=seed)
  status(2,'run.lhs: ',nrow(S),' [',.b,'/',.nb,'] ')
  H = grid.apply(S,function(...){
    status(3,list.str(list(...),sig=4))
    Hi = verb.wrap(run.sam(...,P0=P0,.par=FALSE),0)
  },.grid=FALSE,.rbind=TRUE,.cbind=TRUE,.par=TRUE,.batch=.b,.nbatch=.nb)
  save.rda(H,data.path(.save=TRUE),str('b',.nb),str('H.',.b))
}

post.lhs = function(){
  H = rbind.lapply(1:.nb,function(b){
    load.rda(data.path(),str('b',.nb),str('H.',b)) })
  save.rda(H,data.path(),'H')
}

# run.lhs()
# post.lhs()
