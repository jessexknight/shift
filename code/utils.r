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
  nwarnings=1e4,
  width=200)

len = length
lens = lengths
str = paste0

# -----------------------------------------------------------------------------
# files + i/o

# proj.root: full path to parent of /code.r/fio.r
proj.root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]

root.path = function(...,create=FALSE){
  # e.g. root.path('abc','123') returns proj.root/abc/123
  path = file.path(proj.root,...)
  if (create & !dir.exists(dirname(path))){
    dir.create(dirname(path),recursive=TRUE) }
  return(path)
}

.verb = cli.arg('.verb',4)

status = function(lvl,...){
  if (lvl > .verb){ return() }
  pre = list(c(rep('-',80),'\n'),'',' > ','')[[lvl]]
  end = list('\n','\n','\n','')[[lvl]]
  cat(pre,...,end,sep='')
}

# -----------------------------------------------------------------------------
# vector tools

sum1 = function(x){
  x/sum(x)
}

na.to.num = function(x,num=0){
  # replace NA in x with num
  x[is.na(x)] = num
  return(x)
}

int.cut = function(x,low){
  # cut with simplified labels (assume integers)
  # e.g. int.cut(1:6,c(1,2,3,5)) -> c('1','2','3-4','3-4','5+','5+')
  high = c(low[2:len(low)]-1,Inf)
  labels = gsub('-Inf','+',ifelse(low==high,low,str(low,'-',high)))
  x.cut = cut(x,breaks=c(low,Inf),labels=labels,right=FALSE)
}

breaks = function(x,n=30){
  b = hist(x,breaks=min(n,ulen(x)),plot=FALSE)$breaks
}

ulen = function(x){
  # e.g. ulen(c(1,1,1,2,3)) -> 3
  len(unique(x))
}

even.len = function(x){
  # truncate vector x to have an even length
  length(x) = len(x) - (len(x) %% 2)
  return(x)
}

reppend = function(x,xa,n){
  append(x,rep.int(xa,n))
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

filter.names = function(x,re,b=TRUE){
  # e.g. filter.names(list(a1=0,a2=0,ba=0),'^a') -> c('a1','a2')
  names(x)[grepl(re,names(x))==b]
}

list.str = function(x,def=': ',join=', '){
  # e.g. list.str(list(a=1,b=2)) -> 'a: 1, b: 2'
  paste(names(x),x,sep=def,collapse=join)
}

# -----------------------------------------------------------------------------
# *apply

.cores = cli.arg('.cores',7)

par.lapply = function(...,.par=TRUE){
  if (.par && len(list(...)[[1]]) > 1){
    parallel::mclapply(...,mc.cores=.cores) }
  else {
    lapply(...)
  }
}

par.mapply = function(...){
  parallel::mcmapply(...,mc.cores=.cores,SIMPLIFY=FALSE)
}

rbind.lapply = function(...){
  do.call(rbind,par.lapply(...))
}

wapply = function(...){
  mapply(...,SIMPLIFY=FALSE)
}

grid.apply = function(x,fun,...,.par=TRUE){
  # e.g. grid.lapply(list(a=1:2,b=3:4),fun,c=5) runs:
  # fun(a=1,b=3,c=5), fun(a=2,b=3,c=5), fun(a=1,b=4,c=5), fun(a=2,b=4,c=5)
  xg = expand.grid(x)
  args = lapply(1:nrow(xg),function(i){ c(as.list(xg[i,]),list(...)) })
  par.lapply(args,do.call,what=fun,.par=.par)
}

# -----------------------------------------------------------------------------
# stats

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
    do.call(qfun,c(list(p=p),args)) # target distr quantile transform
  },qfuns,ps,list(...))
}
