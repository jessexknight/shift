import numpy as np
import pandas as pd
from utils import fio,stats
from utils import plot
import seaborn as sb

def ci_aggr(Xg):
  Xi = Xg.iloc[0]
  Xi['n'] = int(Xg['n'].mean())
  Xi['y'] = np.average(Xg['y'],weights=Xg['n'])
  Xi['y.lo'] = np.average(Xg['y.lo'],weights=Xg['n'])
  Xi['y.up'] = np.average(Xg['y.up'],weights=Xg['n'])
  return Xi

def pois_aggr(Xg,n=10000):
  Xi = Xg.iloc[0]
  Xi['n'] = int(Xg['n'].mean())
  Xi['y'] = np.average(Xg['y'],weights=Xg['n'])
  Xi['y.lo'],Xi['y.up'] = np.quantile(
    a=stats.pois(Xi['y'],rng=R).rvs(Xi['n']*n).reshape((Xi['n'],n)).mean(axis=0),
    q=(.025,.975),axis=0)
  return Xi

def load(var):
  fname = fio.rootpath('data','.private','am',var+'.a5v.csv')
  X = pd.read_csv(fname).rename(columns=dict(avg='y',lci='y.lo',uci='y.up'))
  X['var'] = var
  aggr = pois_aggr if 'ptr' in var else ci_aggr
  return X.groupby('age.5',as_index=False).apply(aggr)

def plots(X):
  for var in X:
    plot.target(X[var],x='age.5',y='y',color=plot.c1)
    plot.save('cal','targ_'+var)

R = stats.rng()
vars = ['cdm_ls','ptr_3m_all','ptr_cur','ptr_life']
X = {var:load(var) for var in vars}
plots(X)
