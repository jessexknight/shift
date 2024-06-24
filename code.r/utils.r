len = length
lens = lengths

even.len = function(x){
  # truncate vector x to have an even length
  length(x) = len(x) - (len(x) %% 2)
  return(x)
}

last = function(x){
  # get the last element in x or NA if len(x) == 0
  ifelse(len(x),tail(x,1),NA)
}

alen = function(x){
  # check if len(x) > 0
  len(x) > 0
}
