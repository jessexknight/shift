from pathos.multiprocessing import ProcessingPool as ppool

def sum1(x):
  return x/sum(x)

def mean(x):
  return sum(x)/len(x)

def pmap(fun,Xs,para=True,**kwds):
  fmap = ppool().map if para else map
  return list(fmap(lambda X: fun(X,**kwds),Xs))
