from pathos.multiprocessing import ProcessingPool as ppool
from pandas import cut
import numpy as np

def log(lvl,msg=''):
  # print with immediate output & hierarchy format
  if lvl == 0: print('-'*80+'\n'+str(msg),flush=True)
  if lvl == 1: print(            str(msg),flush=True)
  if lvl == 2: print(      ' > '+str(msg),flush=True)
  if lvl == 3: print(            str(msg),flush=True,end=' ')
  if lvl  < 0: print('',flush=True)

def sum1(x):
  return x/sum(x)

def mean(x):
  return sum(x)/len(x)

def rint(x):
  # round & cast to int
  return np.round(x).astype(int)

def gint(x,bins):
  # group (cut) & label as ints
  labs  = ['{:d}-{:d}'.format(bl,bu-1) for bl,bu in zip(bins[:-1],bins[1:])]
  labs += ['{:d}+'.format(bins[-1])]
  return cut(x,[*bins]+[np.Inf],labels=labs,right=False,precision=0)

def dictpop(x,ks):
  return {k:x.pop(k) for k in ks if k in x}

def pmap(fun,Xs,para=True,**kwds):
  # run fun in parallel for each X in Xs
  fmap = ppool().map if para else map
  return list(fmap(lambda X: fun(X,**kwds),Xs))
