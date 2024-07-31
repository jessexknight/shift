options(width=180)

len = length
lens = lengths

# proj.root: full path to parent of /code.r/fio.r
proj.root = strsplit(file.path(getwd(),''),file.path('','code.r',''))[[1]][1]

root.path = function(...,create=FALSE){
  # e.g. root.path('abc','123') returns proj.root/abc/123
  path = file.path(proj.root,...)
  if (create & !dir.exists(dirname(path))){
    dir.create(dirname(path),recursive=TRUE) }
  return(path)
}

na.to.num = function(x,num=0){
  # replace NA in x with num
  x[is.na(x)] = num
  return(x)
}

even.len = function(x){
  # truncate vector x to have an even length
  length(x) = len(x) - (len(x) %% 2)
  return(x)
}

reppend = function(x,xa,n){
  append(x,rep.int(xa,n))
}

ulist = function(x=list(),xu=list(),...){
  # e.g. ulist(list(a=1,b=2),xu=list(a=3),b=4) -> list(a=3,b=4)
  x = c(x,xu,list(...))
  x[!duplicated(names(x),fromLast=TRUE)]
}

par.lapply = function(...,.par=TRUE){
  if (.par){
    parallel::mclapply(...,mc.cores=7)
  } else {
    lapply(...)
  }
}

rbind.lapply = function(...){
  do.call(rbind,par.lapply(...))
}

filter.names = function(x,re,b=TRUE){
  names(x)[grepl(re,names(x))==b]
}
