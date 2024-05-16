import numpy as np
import model
from utils import stats

def def_cal_distrs(seed=None):
  rng = stats.rng(seed=seed)
  return {
  'ptr_max_m':  stats.unif(rng=rng,l= 1,u= 3),
  'ptr_r0_lm':  stats.unif(rng=rng,l=-5,u=-2),
  'dep_r0_lm':  stats.unif(rng=rng,l=-7,u=-5),
  'dep_x0_lm':  stats.unif(rng=rng,l=-7,u=-5),
  'vio_r0_lm':  stats.unif(rng=rng,l=-7,u=-4),
  'ptr_dur_lm': stats.unif(rng=rng,l=-6,u=-3),
  }

def get_n_sample(seeds,**kwds):
  D  = def_cal_distrs()
  Ps = stats.lhs(D,len(seeds),seed=seeds[0])
  return [get_depend(P,seed=seed,**kwds) for P,seed in zip(Ps,seeds)]

def get_depend(P,seed=None,**kwds):
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
  'age':     stats.unif(rng=rng['ind'],l=model.amin,u=model.amax),
  'new_ind': stats.pois(rng=rng['ind'],m=P['new_ind_m']),
  'ptr_max': stats.geom(rng=rng['ind'],m=P['ptr_max_m']),
  'ptr_r0':  stats.exp (rng=rng['ind'],m=np.exp(P['ptr_r0_lm'])),
  'cdm_p0':  stats.unif(rng=rng['ind'],l=0,u=1),
  'dep_r0':  stats.exp (rng=rng['ind'],m=np.exp(P['dep_r0_lm'])),
  'dep_x0':  stats.exp (rng=rng['ind'],m=np.exp(P['dep_x0_lm'])),
  'vio_r0':  stats.exp (rng=rng['ind'],m=np.exp(P['vio_r0_lm'])),
  # partnership-level
  'ptr_dur': stats.exp(rng=rng['ptr'],m=np.exp(P['dep_r0_lm'])),
  }
