import numpy as np
from utils import stats,ppool
import sim

class Individual():
  #@profile
  def __init__(self,i,N,age,ptr_max,ptr_r0,cdm_p0,dep_r0,dep_x0,vio_r0):
    self.i = i
    self.N = N
    self.age = age
    self.ptr_max = ptr_max
    self.ptr_r0 = ptr_r0
    self.cdm_p0 = cdm_p0
    self.dep_r0 = dep_r0
    self.dep_x0 = dep_x0
    self.vio_r0 = vio_r0
    self.P = []
    self.logs = {y:[] for y in sim.logs}
    self.rpp = N.P['rng']['ptr'].poisson
    self.rip = N.P['rng']['ind'].poisson
    self.rpi = N.P['rng']['ptr'].integers
    self.dep = False
    self.cdm = False

  def __str__(self):
    return '<I:{}>'.format(self.i)

  def __repr__(self):
    return '<I:{}>'.format(self.i)

  #@profile
  def exit(self,z):
    self.logs['exit'].append(z)
    for P in [*self.P]: # BUGFIX
      P.end(z)
    self.N.I.remove(self)
    self.N.Ix.append(self)

  def set_cdm(self,z):
    if self.P:
      self.cdm = self.P[self.rpi(len(self.P))].cdm

  def get_ptr_rate(self,z):
    dage = np.abs(self.age - self.N.P['age_ptr_rmax'])
    return self.ptr_r0 * np.exp(0
      + self.N.P['age_dm:ptr_r']  * dage
      + self.N.P['vio_a3m:ptr_r'] * sim.a3m(self.logs['vio'],z)
      + self.N.P['dep_cur:ptr_r'] * self.dep
    )

  #@profile
  def n_begin_ptr(self,z):
    n = self.rpp(self.get_ptr_rate(z) * sim.dtz)
    return min(self.ptr_max - len(self.P), n)

  #@profile
  def begin_ptr(self,z,P):
    self.logs['begin_ptr'].append(z)
    self.P.append(P)

  #@profile
  def end_ptr(self,z,P):
    self.logs['end_ptr'].append(z)
    self.P.remove(P)

  def get_dep_rate(self,z):
    return self.dep_r0 * np.exp(0
      + self.N.P['vio_a3m:dep_r'] * sim.a3m(self.logs['vio'],z)
    )

  def get_dep_reco(self,z):
    return self.dep_x0 * np.exp(0
      + self.N.P['dep_dur:dep_x'] * sim.dur(self.logs['begin_dep'],z)
    )

  #@profile
  def set_dep(self,z):
    if self.dep:
      if self.rip(self.get_dep_reco(z) * sim.dtz) > 0:
        self.logs['end_dep'].append(z)
        self.dep = False
    else:
      if self.rip(self.get_dep_rate(z) * sim.dtz) > 0:
        self.logs['begin_dep'].append(z)
        self.dep = True

  def get_vio_rate(self,z):
    return self.vio_r0 * np.exp(0
    )

  #@profile
  def set_vio(self,z):
    n = self.rip(self.get_vio_rate(z) * sim.dtz)
    self.logs['vio'] += [z]*n

class Partnership():
  #@profile
  def __init__(self,I1,I2,z0,zdur):
    if I1 == I2: return None # HACK
    self.N  = I1.N
    self.I1 = I1
    self.I2 = I2
    self.z0 = z0
    self.zdur = zdur
    self.riu = self.N.P['rng']['ptr'].random
    self.set_cdm(z0)
    self.I1.begin_ptr(z0,self)
    self.I2.begin_ptr(z0,self)
    self.N.add_evt(z0+zdur,self.end)
    self.active = True

  def __str__(self):
    return '<P:{}:{}>'.format(self.I1.i,self.I2.i)

  def __repr__(self):
    return '<P:{}:{}>'.format(self.I1.i,self.I2.i)

  #@profile
  def end(self,z):
    if not self.active: return
    self.I1.end_ptr(z,self)
    self.I2.end_ptr(z,self)
    self.active = False

  #@profile
  def set_cdm(self,z):
    I1,I2 = self.I1,self.I2
    p0      = .5 * I1.cdm_p0 + .5 * I2.cdm_p0
    ddur    = self.zdur * sim.dtz - self.N.P['ptr_dur_m']
    vio_a3m = sim.a3m(I1.logs['vio'],z) + sim.a3m(I2.logs['vio'],z)
    dep_cur = I1.dep + I2.dep
    self.cdm = self.riu() < stats.plogis(0
      + stats.qlogis(p0)
      + self.N.P['ptr_dur:cdm_p'] * ddur
      + self.N.P['vio_a3m:cdm_p'] * vio_a3m
      + self.N.P['dep_cur:cdm_p'] * dep_cur
    )

class Network():
  def __init__(self,P):
    self.P = P # params
    self.E = {} # events
    self.I = [] # active
    self.Iq = [] # queued
    self.Ix = [] # exited

  #@profile
  def add_evt(self,z,evt,**kwds):
    if z not in self.E: self.E[z] = []
    self.E[z].append((evt,kwds))

  #@profile
  def run(self,zs):
    # initialize individuals
    self.gen_inds(len(zs))
    self.add_inds(self.P['n'])
    for z in zs:
      # scheduled events
      for evt,kwds in self.E.get(z,()):
        evt(z=z,**kwds)
      # every timestep events
      self.age_inds(z)
      self.update_inds(z)
      self.begin_ptrs(z)
    return self

  #@profile
  def begin_ptrs(self,z):
    I = [J for I in self.I for J in [I]*I.n_begin_ptr(z)]
    if len(I) % 2:
      I.pop(-1)
    n = int(len(I)/2)
    # TODO: add shuffle? (or proper mixing)
    list(map(Partnership,
      I[:n],I[n:],[z]*n,
      self.P['ptr_dur'].rvs(n)//sim.dtz+1,
    ))

  #@profile
  def gen_inds(self,nz):
    # pre-generate all individuals for speed
    ni0 = self.P['n']
    niz = int(stats.pois(nz*self.P['new_ind_m']).ppf(.999))
    n = ni0 + niz
    self.Iq.extend(map(Individual,
      range(n),
      [self]*n,
      self.P['age'].rvs(ni0).tolist() + [sim.amin]*niz,
      self.P['ptr_max'].rvs(n),
      self.P['ptr_r0'].rvs(n),
      self.P['cdm_p0'].rvs(n),
      self.P['dep_r0'].rvs(n),
      self.P['dep_x0'].rvs(n),
      self.P['vio_r0'].rvs(n),
    ))

  #@profile
  def add_inds(self,n):
    self.I.extend(self.Iq[:n])
    self.Iq = self.Iq[n:]

  #@profile
  def age_inds(self,z):
    self.add_inds(self.P['new_ind'].rvs())
    for I in self.I:
      I.age += sim.dtz/365
      if I.age >= 50:
        I.exit(z)

  #@profile
  def update_inds(self,z):
    for I in self.I:
      I.set_vio(z)
      I.set_dep(z)
      I.set_cdm(z)
