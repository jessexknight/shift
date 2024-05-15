import numpy as np
import pandas as pd
import model

class Survey():
  def __init__(self,Is):
    self.I = Is
    self.X = pd.DataFrame(
      {k:[getattr(I,k) for I in Is] for k in model.ind_attrs})

  def __getitem__(self,k):
    return self.X[k]

  def __setitem__(self,k,v):
    self.X[k] = v

  def __str__(self):
    return str(self.X)

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
