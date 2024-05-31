import numpy as np
import pandas as pd
import sim
from utils import sum1

class Survey():
  def __init__(self,N,sfun,qfun,z0,zr=0):
    self.N = N # network
    self.sfun = sfun # sampling fun
    self.qfun = qfun # questions fun
    self.z0 = z0 # first survey
    self.zr = zr # survey interval
    self.X = pd.DataFrame() # data
    self.N.add_evt(z0,self.run)

  def __getitem__(self,k):
    return self.X[k]

  def __setitem__(self,k,v):
    self.X[k] = v

  def __str__(self):
    return str(self.X)

  def run(self,z):
    if self.zr: self.N.add_evt(z+self.zr,self.run)
    self.X = pd.concat([self.X,
      self.qfun(self.sfun(self.N.I),z=z) ])

  def save(self,fname):
    self.X.replace({False: 0, True: 1}).to_csv(fname,
      index=False,na_rep='NA')

class Meta(Survey):
  def __init__(self,Qs):
    self.X = pd.concat((Q.X for Q in Qs))

  def run(self,z):
    pass

def sfun(Is):
  # survey sampling function - TODO
  return Is

def qfun(Is,z,pars=[],**kwds):
  # survey question function
  N = Is[0].N
  kwds.update({k:N.P[k] for k in ['seed']+pars},z=z)
  kwds.update({k:[getattr(I,k) for I in Is] for k in sim.ind_attrs})
  X = pd.DataFrame(kwds)
  X['vio_a3m'] = n_vio(Is,      **sim.k3m(z)) > 0
  X['ptr_n3m'] = n_ptr(Is,'tot',**sim.k3m(z))
  return X

def nz(z,z0=-np.Inf,zf=np.Inf):
  # count length of z with values: z0 <= z < zf
  # i.e. how many events were logged in this interval
  i = np.searchsorted(z,[z0,zf],side='right')
  return i[1]-i[0]

def nze(Is,log,z0,zf):
  # count events per individual in the interval z0 <= z < zf
  return np.array([nz(I.logs[log],z0=z0,zf=zf) for I in Is])

def nzp(Is,log,what,z0,zf):
  # count periods per individual in the interval z0 <= z < zf
  # new: n started during the interval
  # end: n ended during the interval
  # cur: n active at the end of the interval
  # tot: n having any overlap with the interval
  begin = 'begin_'+log
  end   = 'end_'  +log
  if what == 'new':
    return np.array([nz(I.logs[begin],z0=z0,zf=zf) for I in Is])
  if what == 'end':
    return np.array([nz(I.logs[end],  z0=z0,zf=zf) for I in Is])
  if what == 'cur':
    return np.array([nz(I.logs[begin],zf=zf) - nz(I.logs[end],zf=zf) for I in Is])
  if what == 'tot':
    return np.array([nz(I.logs[begin],zf=zf) - nz(I.logs[end],zf=z0) for I in Is])

def n_ptr(Is,what,z0=-np.Inf,zf=np.Inf):
  return nzp(Is,log='ptr',what=what,z0=z0,zf=zf)

def n_dep(Is,what,z0=-np.Inf,zf=np.Inf):
  return nzp(Is,log='dep',what=what,z0=z0,zf=zf)

def n_vio(Is,z0=-np.Inf,zf=np.Inf):
  return nze(Is,log='vio',z0=z0,zf=zf)

def attr(Is,what):
  return np.array([getattr(I,what) for I in Is])

def gist(X,var,b,g=None):
  # grouped histogram for repeated measures
  g = ['seed'] if g is None else g
  f = lambda x: sum1(np.histogram(x,[*b,np.Inf])[0])
  H = X.groupby(g)[var].apply(f).reset_index()
  return H.explode(var).assign(b=[*b]*len(H))
