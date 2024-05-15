import numpy as np
import scipy.stats as ss

q2 = (.025,.975)
q3 = (.025,.5,.975)
q4 = (.025,.25,.75,.975)
q5 = (.025,.25,.5,.75,.975)

# generators

def rng(seed=None,state=None):
  R = np.random.default_rng(seed)
  if state:
    R.bit_generator.state = state
  return R

def rngs(keys,seed=None,ds=1e9):
  seedmap = lambda i: seed if seed is None else seed+int(ds)*i
  return {k:rng(seedmap(i)) for i,k in enumerate(keys)}

def state(R):
  if isinstance(R,np.random.Generator):
    return R.bit_generator.state
  else:
    return rng(state=R)

def states(Rs):
  return {k:state(R) for k,R in Rs.items()}

RNG = rng()

def drng(dfun):
  def deco(*args,rng=RNG,**kwds):
    d = dfun(*args,**kwds)
    d.random_state = rng
    return d
  return deco

def lhs(D,n,seed=None):
  qs = ss.qmc.LatinHypercube(len(D),seed=seed).random(n)
  return [{k:D[k].ppf(qk) for k,qk in zip(D,q)} for q in qs]

# discrete distrs

@drng
def bern(p):
  return ss.bernoulli(p=p)

@drng
def geom(m,z=0,rng=RNG):
  return ss.geom(p=1/m,loc=z)

@drng
def pois(m,z=0,rng=RNG):
  return ss.poisson(mu=m,loc=z)

@drng
def nbinom(m,sd,z=0,eps=1e-6,rng=RNG):
  p = min(1-eps,m/sd/sd)
  return ss.nbinom(n=m*p/(1-p),p=p,loc=z)

# continuous distrs

@drng
def unif(l=0,u=1,rng=RNG):
  return ss.uniform(loc=l,scale=u-l)

@drng
def exp(m,z=0,rng=RNG):
  return ss.expon(scale=m,loc=z)

@drng
def gamma(m,sd,z=0,rng=RNG):
  scale=sd*sd/m
  return ss.gamma(a=m/scale,scale=scale,loc=z)

@drng
def norm(m,sd,rng=RNG):
  return ss.norm(loc=m,scale=sd)

# transformations

def plogis(q):
  return 1/(1+np.exp(-q))

def qlogis(p):
  return -np.log(1/p-1)
