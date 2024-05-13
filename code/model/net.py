from utils import stats
import model

class Individual():
  def __init__(self,i,N,age,ptr_max,ptr_rate,cdm_pref,dep_rate,vio_rate):
    self.i = i
    self.N = N
    self.age = age
    self.ptr_max = ptr_max
    self.ptr_r0 = ptr_rate
    self.vio_r0 = vio_rate
    self.dep_r0 = dep_rate
    self.cdm_p0 = cdm_pref
    self.dep_state = 0
    self.P = set()
    self.logs = []
    self.rpp = N.D['rng']['ptr'].poisson
    self.active = True

  def __str__(self):
    return '<I:{}>'.format(self.i)

  def __repr__(self):
    return '<I:{}>'.format(self.i)

  def exit(self,k):
    self.logs.append({'k':k,'e':'exit'})
    for P in [*self.P]:
      P.end(k)
    self.N.I.remove(self)
    self.N.J.append(self)
    self.active = False

  def get_ptr_rate(self,k):
    return self.ptr_r0

  def n_begin_ptr(self,k):
    n = self.rpp(self.get_ptr_rate(k)*model.dtk)
    return min(self.ptr_max-len(self.P),n)*self.active

  def begin_ptr(self,k,P):
    self.logs.append({'k':k,'e':'begin_ptr'})
    self.P.add(P)

  def end_ptr(self,k,P):
    self.logs.append({'k':k,'e':'end_ptr'})
    self.P.remove(P)

class Partnership():
  def __init__(self,I1,I2,k0,dur):
    if I1 == I2: return None # HACK
    self.N  = I1.N
    self.I1 = I1
    self.I2 = I2
    self.k0 = k0
    self.dur = dur
    self.set_cdm(k0)
    self.I1.begin_ptr(k0,self)
    self.I2.begin_ptr(k0,self)
    self.N.add_evt(k0+dur,self.end)
    self.active = True

  def __str__(self):
    return '<P:{}:{}>'.format(self.I1.i,self.I2.i)

  def __repr__(self):
    return '<P:{}:{}>'.format(self.I1.i,self.I2.i)

  def end(self,k):
    if not self.active: return
    self.I1.end_ptr(k,self)
    self.I2.end_ptr(k,self)
    self.active = False

  def set_cdm(self,k):
    self.cdm = stats.plogis(0
      +.5*stats.qlogis(self.I1.cdm_p0)
      +.5*stats.qlogis(self.I2.cdm_p0))

class Network():
  def __init__(self,D):
    self.D = D  # distrs
    self.E = {} # events
    self.I = [] # active
    self.J = [] # exited
    self.imax = 0 # next I.i
    self.rip = D['rng']['ind'].poisson
    self.add_evt(0,self.begin_ptrs)
    self.add_evt(0,self.age_inds)

  def attrs(self,attr,fun=None):
    if fun is None: fun = lambda x: x
    return [fun(getattr(I,attr)) for I in self.I]

  def add_evt(self,k,evt,**kwds):
    if k not in self.E: self.E[k] = []
    self.E[k].append((evt,kwds))

  def run(self,ks):
    for k in ks:
      for evt,kwds in self.E.get(k,()):
        evt(k=k,**kwds)

  def begin_ptrs(self,k):
    self.add_evt(k+1,self.begin_ptrs)
    I = [J for I in self.I for J in [I]*I.n_begin_ptr(k)]
    if len(I) % 2:
      I.pop(-1)
    n = int(len(I)/2)
    list(map(Partnership,
      I[:n],I[n:],[k]*n,
      self.D['ptr_dur'].rvs(n)//model.dtk+1,
    ))

  def add_inds(self,n,k=0,ages=None):
    self.imax += n
    self.I.extend(map(Individual,
      range(self.imax-n,self.imax),
      [self]*n,
      [model.amin]*n if ages is None else ages,
      self.D['ptr_max'].rvs(n),
      self.D['ptr_rate'].rvs(n),
      self.D['cdm_pref'].rvs(n),
      self.D['dep_rate'].rvs(n),
      self.D['vio_rate'].rvs(n),
    ))

  def age_inds(self,k):
    self.add_evt(k+1,self.age_inds)
    n = self.rip(len(self.I)*model.dtk/365/model.adur)
    self.add_inds(n=n,k=k)
    for I in self.I:
      I.age += model.dtk/365
      if I.age > 50:
        I.exit(k)
