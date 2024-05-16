from utils import stats,ppool
import model

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
    self.P = set()
    self.logs = {y:[] for y in model.logs}
    self.rpp = N.P['rng']['ptr'].poisson
    self.rpi = N.P['rng']['ind'].poisson
    self.depressed = False
    self.active = True

  def __str__(self):
    return '<I:{}>'.format(self.i)

  def __repr__(self):
    return '<I:{}>'.format(self.i)

  #@profile
  def exit(self,z):
    self.logs['exit'].append(z)
    for P in [*self.P]:
      P.end(z)
    self.N.I.remove(self)
    self.N.J.append(self)
    self.active = False

  def get_ptr_rate(self,z):
    return self.ptr_r0

  #@profile
  def n_begin_ptr(self,z):
    n = self.rpp(self.get_ptr_rate(z)*model.dtz)
    return min(self.ptr_max-len(self.P),n)*self.active

  #@profile
  def begin_ptr(self,z,P):
    self.logs['begin_ptr'].append(z)
    self.P.add(P)

  #@profile
  def end_ptr(self,z,P):
    self.logs['end_ptr'].append(z)
    self.P.remove(P)

  def get_dep_rate(self,z):
    return self.dep_r0

  def get_dep_reco(self,z):
    return self.dep_x0

  #@profile
  def set_dep(self,z):
    if self.depressed:
      if self.rpi(self.get_dep_reco(z)*model.dtz)>0:
        self.logs['end_dep'].append(z)
        self.depressed = False
    else:
      if self.rpi(self.get_dep_rate(z)*model.dtz)>0:
        self.logs['begin_dep'].append(z)
        self.depressed = True

  def get_vio_rate(self,z):
    return self.vio_r0

  #@profile
  def set_vio(self,z):
    n = self.rpi(self.get_vio_rate(z)*model.dtz)
    self.logs['vio'] += [z]*n

class Partnership():
  #@profile
  def __init__(self,I1,I2,z0,dur):
    if I1 == I2: return None # HACK
    self.N  = I1.N
    self.I1 = I1
    self.I2 = I2
    self.z0 = z0
    self.dur = dur
    self.set_cdm(z0)
    self.I1.begin_ptr(z0,self)
    self.I2.begin_ptr(z0,self)
    self.N.add_evt(z0+dur,self.end)
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
    self.cdm = stats.plogis(0
      +.5*stats.qlogis(self.I1.cdm_p0)
      +.5*stats.qlogis(self.I2.cdm_p0))

class Network():
  def __init__(self,P):
    self.P = P  # params
    self.E = {} # events
    self.I = [] # active
    self.J = [] # exited
    self.imax = 0 # next I.i
    self.rpi = P['rng']['ind'].poisson
    self.add_inds(P['n'],
      ages=model.amin+model.adur*P['rng']['ind'].random(P['n']))

  #@profile
  def add_evt(self,z,evt,**kwds):
    if z not in self.E: self.E[z] = []
    self.E[z].append((evt,kwds))

  #@profile
  def run(self,zs):
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
    list(map(Partnership,
      I[:n],I[n:],[z]*n,
      self.P['ptr_dur'].rvs(n)//model.dtz+1,
    ))

  #@profile
  def add_inds(self,n,z=0,ages=None):
    self.imax += n
    self.I.extend(map(Individual,
      range(self.imax-n,self.imax),
      [self]*n,
      [model.amin]*n if ages is None else ages,
      self.P['ptr_max'].rvs(n),
      self.P['ptr_r0'].rvs(n),
      self.P['cdm_p0'].rvs(n),
      self.P['dep_r0'].rvs(n),
      self.P['dep_x0'].rvs(n),
      self.P['vio_r0'].rvs(n),
    ))

  #@profile
  def age_inds(self,z):
    n = self.rpi(len(self.I)*model.dtz/365/model.adur)
    self.add_inds(n=n,z=z)
    for I in self.I:
      I.age += model.dtz/365
      if I.age > 50:
        I.exit(z)

  #@profile
  def update_inds(self,z):
    for I in self.I:
      I.set_dep(z)
      I.set_vio(z)

def run_n(Ps,zs,para=True):
  frun = lambda P: Network(P).run(zs)
  fmap = ppool().map if para else map
  return list(fmap(frun,Ps))
