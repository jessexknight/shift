import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
from sim import srv
from utils import log,fio,dictpop

# sb.set_palette('inferno') # TODO
sb.set(style='whitegrid',rc={
  'patch.linewidth': 0,
  'axes.edgecolor': '#000',
  'grid.color':'#eee',
})
c1 = '#990033'

def save(*args,w=5,h=3.5):
  fig = plt.gcf()
  fig.set_size_inches((w,h))
  fig.tight_layout()
  fname = fio.genpath(fio.rootpath('out','fig',*args))+'.pdf'
  log(2,'plot.save: '+fname)
  plt.savefig(fname)
  plt.close()

aks = [d+k for d in ('x','y') for k in ('lim','label','scale','ticks')]

def clean(fun):
  # e.g. call ax.set_xlim(kwds['xlim']) after ax = fun(...)
  def deco(*args,**kwds):
    akwd = dictpop(kwds,aks)
    ax = fun(*args,**kwds)
    for k,v in akwd.items():
      getattr(ax,'set_'+k)(v)
    return ax
  return deco

@clean
def line(X,err='pi',**kwds):
  kwds.update(errorbar=err)
  return sb.lineplot(X,**kwds)

@clean
def band(X,y1,y2,legend=False,**kwds):
  X2 = pd.concat((X,X))
  X2['.tmp'] = [*X[y1],*X[y2]]
  return line(X2,y='.tmp',err=('pi',100),dashes=(0,1),legend=legend,**kwds)

@clean
def violin(X,**kwds):
  kwds.update(fill=False,inner='quart',inner_kws=dict(dashes=(1,0)))
  return sb.violinplot(X,**kwds)

def gist(X,y,b,type='line',**kwds):
  # plot grouped histogram
  plot = dict(line=line,violin=violin)[type]
  g = ['seed']+[v for v in kwds.values() if v in X]
  H = srv.gist(X,y,b,g)
  return plot(H,y=y,x='b',xlabel=y,ylabel=' ',**kwds)

def distr(Ps,key,cum=False,log=False,eps=.01,n=100):
  # plot density of P[key] for a list of Ps
  xmin = np.min([P[key].ppf(0)     for P in Ps])
  xmax = np.max([P[key].ppf(1-eps) for P in Ps])
  x = np.arange(xmin,xmax,(xmax-xmin)/n)
  if cum: d = [P[key].cdf(x) for P in Ps]
  else:   d = [P[key].pdf(x) for P in Ps]
  X  = pd.DataFrame({'x':np.tile(x,len(Ps)),'d':np.array(d).flatten()})
  return line(X,x='x',y='d',color=c1,xlabel=key,
    ylabel=('Cumulative ' if cum else '')+('Log ' if log else '')+'Density',
    yscale=('log' if log else 'linear'))
