import itertools as it
import sys
sys.setrecursionlimit(10000)
dtz  =  7 # days per timestep
z3m  = 13 # timesteps per 3 months
z1y  = 52 # timesteps per 1 year
amin = 15 # age of model entry
amax = 50 # age of model exit
adur = 35 # years in model
a5g  = [18,20,25,30,35,40,45,50]
logs = (  # types of logs
  'begin_ptr',
  'end_ptr',
  'begin_dep',
  'end_dep',
  'vio',
  'exit',
)
ind_attrs = ( # individual attributes
  'i',
  'age',
  'ptr_max',
  'ptr_r0',
  'cdm_p0',
  'dep_r0',
  'dep_x0',
  'vio_r0',
  'dep',
  'cdm',
)
ptr_attrs = ( # partner attributes
  'z0',
  'dur',
  'cdm',
)
def k3m(z): # kwds for past 3 months from z
  return {'z0':z-z3m,'zf':z}
def k1y(z): # kwds for past 1 year from z
  return {'z0':z-z1y,'zf':z}
def a3m(log,z): # if any in past 3 months from z
  return len(log) and log[-1] >= (z - z3m)
def a1y(log,z): # if any in past 1 year from z
  return len(log) and log[-1] >= (z - z1y)
def n3m(log,z): # how many in past 3 months from z
  return len(list(it.takewhile(lambda l: l >= (z - z3m), log[::-1])))
def n1y(log,z): # how many in past 1 year from z
  return len(list(it.takewhile(lambda l: l >= (z - z1y), log[::-1])))
def dur(log,z):
  return z - log[-1] if log else 0
