# -----------------------------------------------------------------------------
# options + aliases

cli.arg = function(name,default=NA){
  args = strsplit(commandArgs(trailingOnly=TRUE),'=')
  x = args[[match(name,sapply(args,`[`,1))]][2]
  if (length(x)){ x = type.convert(x,as.is=TRUE) }
  else          { x = default }
}

options(
  stringsAsFactors=FALSE,
  showNCalls=500,
  width=200,
  warn=1)

len = length
lens = lengths
seqn = seq_len
seqa = seq_along
str = paste0
set.names = setNames
col.split = reshape2::colsplit
no.warn = suppressWarnings
stfu = suppressPackageStartupMessages

def.args = function(f,...){
  # pre-specify some args, but allow later override of named args
  args.pre = list(...)
  f.pre = function(...){
    args = c(args.pre,list(...))
    u = names(args)==''
    do.call(f,c(args[u],ulist(args[!u]))) }
}

# -----------------------------------------------------------------------------
# files + i/o

# proj.root: full path to parent of /code.r/fio.r
proj.root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]

root.path = function(...,ext='',create=FALSE){
  # e.g. root.path('abc','123',ext='.csv') returns proj.root/abc/123.csv
  path = str(file.path(proj.root,...),ext)
  if (create & !dir.exists(dirname(path))){
    dir.create(dirname(path),recursive=TRUE) }
  return(path)
}

.verb = cli.arg('.verb',4)

verb.wrap = function(code,.verb.tmp){
  # temporarily override .verb
  .verb.mem = .verb
  .verb <<- .verb.tmp
  out = force(code)
  .verb <<- .verb.mem
  return (out)
}

status = function(lvl,...,id=NULL){
  if (lvl > .verb | lvl <= 0){ return() }
  pre = list(c(rep('-',80),'\n'),'',' > ','')[[lvl]]
  end = list('\n','\n','\n','')[[lvl]]
  cat(pre,...,sprintf('%6d',id),end,sep='')
}

load.csv = function(...,ext='.csv'){
  fname = root.path(...,ext=ext)
  status(3,'load: ',fname)
  read.csv(file=fname,fileEncoding='Latin1')
}

save.csv = function(X,...,ext='.csv'){
  fname = root.path(...,ext=ext,create=TRUE)
  status(3,'save: ',fname)
  write.csv(X,file=fname,row.names=FALSE)
}

load.rda = function(...,ext='.rda'){
  fname = root.path(...,ext=ext)
  status(3,'load: ',fname)
  load(fname)
  return(X)
}

save.rda = function(X,...,ext='.rda'){
  fname = root.path(...,ext=ext,create=TRUE)
  status(3,'save: ',fname)
  save(X,file=fname)
}

load.json = function(...,ext='.json'){
  fname = root.path(...,ext=ext)
  status(3,'load: ',fname)
  rjson::fromJSON(file=fname)
}

save.json = function(X,...,ext='.json',indent=2){
  fname = root.path(...,ext=ext,create=TRUE)
  status(3,'save: ',fname)
  write(rjson::toJSON(X,indent=indent),file=fname)
}

hash.path = function(info,...,.len=11,.save=TRUE){
  hash = substr(digest::sha1(info),1,.len)
  info = ulist(info,time=str(Sys.time()))
  if (.save){ save.json(info,...,hash,'info') }
  return(file.path(...,hash))
}

# -----------------------------------------------------------------------------
# plot tools

# note: load ggplot2 elsewhere for speed

plot.save = function(g,...,size=NULL,ext='.pdf'){
  if (missing(size)){ size = plot.size(g) }
  if (ext=='.pdf'){ dev = cairo_pdf } else { dev = NULL }
  fname = root.path('out','fig',...,ext=ext,create=TRUE)
  status(3,'save: ',fname)
  ggsave(plot=g,file=fname,w=size[1],h=size[2],device=dev)
}

plot.1o = list(w1=2,h1=2,wo=1,ho=1) # width & height of each facet & offsets

plot.size = function(g,...){
  # define size from facet grid
  s = ulist(plot.1o,...)
  layout = ggplot_build(g)$layout$layout
  size = c(w=s$wo+s$w1*max(layout$COL),h=s$ho+s$h1*max(layout$ROW))
}

clr.map.c = def.args(ggplot2::scale_color_viridis_c,aes=c('color','fill'))
clr.map.d = def.args(ggplot2::scale_color_viridis_d,aes=c('color','fill'),begin=.1,end=.9)
clr.map.b = def.args(ggplot2::scale_color_viridis_b,aes=c('color','fill'))
qfun = function(p){ def.args(quantile,p=p) } # e.g. for stat_summary(...,fun=qfun(.5))

str.lab = function(pre='',post=''){
  # e.g. facet_*(...,labeller=str.lab('Y = ',' %')) -> 'Y = ... %'
  ggplot2::as_labeller(function(x){ str(pre,x,post) })
}

plot.clean = function(g,font=NULL,...){
  g = g + theme_light() + theme(...,
    text=element_text(family=font),
    strip.background=element_rect(fill='#eee'),
    strip.text.x=element_text(color='black'),
    strip.text.y=element_text(color='black'))
}

# -----------------------------------------------------------------------------
# vector tools

e10 = function(x){ 10^x } # log10 inverse

sum1 = function(x){ x/sum(x) }

na.to.num = function(x,num=0){
  # replace NA in x with num
  x[is.na(x)] = num
  return(x)
}

if.na = function(x,alt){
  ifelse(is.na(x),alt,x)
}

if.null = function(x,alt){
  if (is.null(x)){ return(alt) } else { return(x) }
}

add.na = function(x,label='<NA>'){
  x = addNA(x)
  attributes(x)$levels[len(levels(x))] = label
  return(x)
}

int.cut = function(x,low,up=Inf){
  # cut with simplified labels (assume integers)
  # e.g. int.cut(1:6,c(1,2,3,5)) -> c('1','2','3-4','3-4','5+','5+')
  high = c(low[2:len(low)]-1,up)
  labels = ifelse(low==high,low,str(low,'-',high))
  labels = gsub('-Inf','+',gsub(str('-',up,'$'),str('-',up-1),labels))
  x.cut = cut(x,breaks=c(low,Inf),labels=labels,right=FALSE)
}

rank.cut = function(x,q){
  x.rank = (rank(x)-1)/(len(x)-1)
  x.cut = cut(x.rank,breaks=c(0,q,1),right=FALSE,incl=TRUE)
}

breaks = function(x,n=30){
  b = hist(x,breaks=min(n,ulen(x)),plot=FALSE)$breaks
}

ulen = function(x){
  # e.g. ulen(c(1,1,1,2,3)) -> 3
  len(unique(x))
}

reppend = function(x,xa,n){
  append(x,rep.int(xa,n))
}

first = function(x){
  # return the first element in x or NA if len(x) == 0
  if (len(x)){ x[1] } else { NA } # redundant
}

last = function(x){
  # return the last element in x or NA if len(x) == 0
  if (len(x)){ x[len(x)] } else { NA }
}

# -----------------------------------------------------------------------------
# list tools

ulist = function(x=list(),xu=list(),...){
  # e.g. ulist(list(a=1,b=2),xu=list(a=3),b=4) -> list(a=3,b=4)
  x = c(x,xu,list(...))
  x[!duplicated(names(x),fromLast=TRUE)]
}

flist = function(x){
  # flatten list(list(a=1),list(b=2)) -> list(a=1,b=2)
  do.call(c,unname(x))
}

name.list = function(...,key='name'){
  # e.g. name.list(list(name='a',v=1),...) -> list(a=list(name='a',v=1),...)
  x = list(...)
  x = set.names(x,lapply(x,`[[`,key))
}

filter.names = function(x,re,b=TRUE){
  # e.g. filter.names(list(a1=0,a2=0,ba=0),'^a') -> c('a1','a2')
  names(x)[grepl(re,names(x),perl=TRUE)==b]
}

list.str = function(x,def=' = ',join=', ',sig=Inf,rnd=Inf){
  # e.g. list.str(list(a=1,b=2)) -> 'a = 1\nb = 2'
  f = function(x){ ifelse(is.numeric(x),signif(round(x,rnd),sig),x) }
  paste(names(x),lapply(x,f),sep=def,collapse=join)
}

.isep = ' x '
interac = function(...){ interaction(...,sep=.isep) }

# -----------------------------------------------------------------------------
# data.frame tools

df.sub = function(X,sub=NULL){
  # sub should be a character vector
  if (is.null(sub)){ return(X) }
  subset(X,eval(parse(text=sub)))
}

df.compare = function(x,y,v=NULL,cast=as.numeric){
  # check if x[v] == y[v] (for debug)
  v = if.null(v,intersect(names(x),names(y)))
  eq = all.equal(lapply(x[v],cast),lapply(y[v],cast))
  status(3,'df.compare @ ',paste(v,collapse=','),': ',
    ifelse(eq==TRUE,'OK',paste('\n',eq)))
}

df.ow = function(X,...){
  # overwrite cols in X
  as.data.frame(ulist(X,...))
}

# -----------------------------------------------------------------------------
# *apply

.cores = cli.arg('.cores',7)

par.lapply = function(...,.par=TRUE){
  if (.par && len(list(...)[[1]]) > 1){
    parallel::mclapply(...,mc.cores=.cores)
  } else {
    lapply(...)
  }
}

par.mapply = function(...,.par=TRUE){
  if (.par){
    parallel::mcmapply(...,mc.cores=.cores,SIMPLIFY=FALSE)
  } else {
    mapply(...,SIMPLIFY=FALSE)
  }
}

rbind.lapply = function(...){
  do.call(rbind,c(par.lapply(...)))
}

wapply = function(...){
  mapply(...,SIMPLIFY=FALSE)
}

grid.apply = function(x,fun,args=list(),...,
  .par=TRUE,.rbind=FALSE,.cbind=FALSE,.grid=TRUE,.batch=1,.nbatch=1){
  # e.g. grid.lapply(list(a=1:2,b=3:4),fun,c=5) runs:
  # fun(a=1,b=3,c=5), fun(a=2,b=3,c=5), fun(a=1,b=4,c=5), fun(a=2,b=4,c=5)
  # optional: split grid args into .nbatch & run .batch only
  xg = ifelse(.grid,expand.grid,as.data.frame)(x,stringsAsFactors=FALSE)
  ng = nrow(xg); gi = seqn(ng);
  grid.args   = lapply(gi,function(i){ ulist(as.list(xg[i,,drop=FALSE]),args,...) })
  grid.args   = split(grid.args,ceiling(gi*min(ng,.nbatch)/ng))[[min(ng,.batch)]]
  grid.fun    = ifelse(.cbind,function(...){ cbind(fun(...),...) },fun)
  grid.lapply = ifelse(.rbind,rbind.lapply,par.lapply)
  grid.lapply(grid.args,do.call,what=grid.fun,.par=.par)
}

fast.split = function(...){
  collapse::rsplit(...,flatten=TRUE)
}

# -----------------------------------------------------------------------------
# stats

p2 = c(lo=.025,hi=.975)
p3 = c(lo=.025,md=.5,hi=.975)
p5 = c(.025,.25,.5,.75,.975)
q2 = qnorm(p2)

fit.beta = function(qs,ps=p2){
  efun = function(par){ e = sum(abs(ps-pbeta(qs,par[1],par[2]))) }
  optim(c(1,1),efun,method='L-BFGS-B',lower=0)$par
}

fit.weibull = function(m,cv2,...){
  efun = function(k){ s = gamma(1+1/k)^2; e = ((gamma(1+2/k)-s)/s-cv2)^2 }
  k = optimize(efun,c(1e-6,1e+6))$minimum
  par = list(shape=k,scale=m/gamma(1+1/k),...)
}

# R2 = dummy 2-group distr with { p0: x0, 1-p0: x0*xR }
qR2 = function(p,x0,xR,p0=.5){ x = x0 * (1 + (xR-1) * (p > p0)) }
rR2 = function(n,x0,xR,p0=.5){ x = qR2(runif(n),x0,xR,p0) }

het.funs = list(
  # m = mean; het = CV (sd / mean)
  gamma = list(
    r = function(n,m,het){ cv2 = max(het^2,1e-9); rgamma(n,shape=1/cv2,scale=m*cv2) },
    d = function(x,m,het){ cv2 = max(het^2,1e-9); dgamma(x,shape=1/cv2,scale=m*cv2) },
    p = function(q,m,het){ cv2 = max(het^2,1e-9); pgamma(q,shape=1/cv2,scale=m*cv2) },
    q = function(p,m,het){ cv2 = max(het^2,1e-9); qgamma(p,shape=1/cv2,scale=m*cv2) }),
  weibull = list(
    r = function(n,m,het){ f = fit.weibull(m,het^2); rweibull(n,shape=f$shape,scale=f$scale) },
    d = function(x,m,het){ f = fit.weibull(m,het^2); dweibull(x,shape=f$shape,scale=f$scale) },
    p = function(q,m,het){ f = fit.weibull(m,het^2); pweibull(q,shape=f$shape,scale=f$scale) },
    q = function(p,m,het){ f = fit.weibull(m,het^2); qweibull(p,shape=f$shape,scale=f$scale) }),
  lnorm = list(
    r = function(n,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); rlnorm(n,meanlog=u,sdlog=s) },
    d = function(x,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); dlnorm(x,meanlog=u,sdlog=s) },
    p = function(q,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); plnorm(q,meanlog=u,sdlog=s) },
    q = function(p,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); qlnorm(p,meanlog=u,sdlog=s) }),
  R2 = list( # m = mean; het = xR; p0 = 0.5 (fixed)
    r = function(n,m,het){ x0 = 2*m/(1+het); x = rR2(n,x0=x0,xR=het,p0=.5) },
    q = function(p,m,het){ x0 = 2*m/(1+het); x = qR2(p,x0=x0,xR=het,p0=.5) }))

copula = function(n,covs,qfuns,...){
  # joint sample from qfuns (args in ...) with correlation (covs)
  # e.g. copula(100,0.5,list(a=qexp,b=qunif),a=list(rate=1),b=list(min=0,max=1))
  #      draws 100 samples from rexp & runif with 50% corrleation
  # final correlations are not exact due to quantile transformations
  d = len(qfuns)
  sigma = matrix(0,d,d)
  sigma[lower.tri(sigma)] = covs
  sigma = sigma + t(sigma) + diag(d)
  ms = mvtnorm::rmvnorm(n,sigma=sigma) # sample from mvn
  ps = as.list(data.frame(pnorm(ms))) # normal cdf transform
  xs = mapply(function(qfun,p,args){ # for each qfun, p vector, args
    do.call(qfun,ulist(args,p=p)) # target distr quantile transform
  },qfuns,ps,list(...))
}

aggr.form = function(y,x,ny=NULL,.log=3){
  ny = if.null(ny,if.null(names(y),y))
  f = str('cbind( ',paste(ny,y,sep=' = ',collapse=' , '),' ) ~ ',str(x,collapse=' + '))
  status(.log,'aggr: ',f)
  return(f)
}
