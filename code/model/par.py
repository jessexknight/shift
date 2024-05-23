import numpy as np
from copy import deepcopy
import model
from utils import stats

def def_cal_distrs(seed=None):
  # sampling distributions for calibrated params
  rng = stats.rng(seed=seed)
  return {
  'ptr_max_m':  stats.unif(rng=rng,l= 1,u= 3),
  'ptr_r0_lm':  stats.unif(rng=rng,l=-5,u=-2),
  'dep_r0_lm':  stats.unif(rng=rng,l=-7,u=-5),
  'dep_x0_lm':  stats.unif(rng=rng,l=-7,u=-5),
  'vio_r0_lm':  stats.unif(rng=rng,l=-7,u=-4),
  'ptr_dur_lm': stats.unif(rng=rng,l= 6,u= 2),
  'dep_dur:dep_x': stats.unif(rng=rng,l= 0,u=1),
  'vio_a3m:dep_r': stats.unif(rng=rng,l= 0,u=6),
  'vio_a3m:ptr_r': stats.unif(rng=rng,l=-2,u=2),
  'vio_a3m:cdm_p': stats.unif(rng=rng,l=-1,u=1),
  'dep_cur:ptr_r': stats.unif(rng=rng,l= 0,u=3),
  'dep_cur:cdm_p': stats.unif(rng=rng,l=-4,u=0),
  'ptr_dur:cdm_p': stats.unif(rng=rng,l=-0.1,u=0),
  }

def get_cal_sample(seed=None):
  # get a LHS sample of calibrated params
  return stats.lhs(def_cal_distrs(seed=seed),seed=seed)

def get_cal_fixed():
  # get a fixed sample of calibrated params
  return {
  'ptr_max_m':      2.0, # 2.0
  'ptr_r0_lm':     -3.0, # .05
  'dep_r0_lm':     -7.0, # .001
  'dep_x0_lm':     -4.0, # .02
  'vio_r0_lm':     -6.0, # .002
  'ptr_dur_lm':     4.0, # 55
  'dep_dur:dep_x':  0.5, # 0.5
  'vio_a3m:dep_r':  3.0, # 3.0
  'vio_a3m:ptr_r':  0.0, # 0.0
  'vio_a3m:cdm_p':  0.0, # 0.0
  'dep_cur:ptr_r':  2.0, # 2.0
  'dep_cur:cdm_p': -2.0, # -2.0
  'ptr_dur:cdm_p': -0.01, # -0.01
  }

def get_n_all(seeds,Ps=None,**kwds):
  # repeat get_all for every seed
  if Ps is None: Ps = [get_cal_fixed()]
  return [get_all(P,seed=seed,**kwds) for P in Ps for seed in seeds]

def get_all(P,seed=None,**kwds):
  # add fixed, dependent, distr params to P for a given seed
  P = deepcopy(P)
  P['n'] = 1000
  P.update(seed=seed,**kwds)
  P['new_ind_m'] = P['n']*model.dtz/365/model.adur
  P.update(get_net_distrs(P,seed=seed))
  return P

def get_net_distrs(P,seed=None,states=None):
  rng = stats.states(states) if states else \
        stats.rngs(['ind','ptr'],seed=seed)
  return {
  'rng': rng,
  # individual-level
  'new_ind': stats.pois(rng=rng['ind'],m=P['new_ind_m']),
  'age':     stats.unif(rng=rng['ind'],l=model.amin,u=model.amax),
  'ptr_max': stats.geom(rng=rng['ind'],m=P['ptr_max_m']),
  'ptr_r0':  stats.exp (rng=rng['ind'],m=np.exp(P['ptr_r0_lm'])),
  'cdm_p0':  stats.unif(rng=rng['ind'],l=0,u=1),
  'dep_r0':  stats.exp (rng=rng['ind'],m=np.exp(P['dep_r0_lm'])),
  'dep_x0':  stats.exp (rng=rng['ind'],m=np.exp(P['dep_x0_lm'])),
  'vio_r0':  stats.exp (rng=rng['ind'],m=np.exp(P['vio_r0_lm'])),
  # partnership-level
  'ptr_dur': stats.exp(rng=rng['ptr'],m=np.exp(P['ptr_dur_lm'])),
  }
