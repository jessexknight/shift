import model
from utils import stats

def distrs(seed=None,states=None):
  # unit = days
  rng = stats.states(states) if states else \
        stats.rngs(['ind','ptr'],seed=seed)
  D = {'rng':rng}
  # individual-level
  D['new_ind']  = stats.pois(rng=rng['ind'],m=1/365/model.adur)
  D['ptr_max']  = stats.geom(rng=rng['ind'],m=2)
  D['cdm_pref'] = stats.unif(rng=rng['ind'],l=0,u=1)
  D['ptr_rate'] = stats.exp (rng=rng['ind'],m=1/21)
  D['dep_rate'] = stats.exp (rng=rng['ind'],m=1/70)
  D['vio_rate'] = stats.exp (rng=rng['ind'],m=1/70)
  # partnership-level
  D['ptr_dur'] = stats.geom(rng=rng['ptr'],m=21)
  # D['age_pref'] = stats.norm(rng=rng['ptr'],m=0,sd=5.1)
  return D

def init_ages(D,n):
  return model.amin + model.adur * D['rng']['ind'].random(n)
