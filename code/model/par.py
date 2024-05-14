import model
from utils import stats

def distrs(seed=None,states=None):
  # unit = days
  rng = stats.states(states) if states else \
        stats.rngs(['ind','ptr'],seed=seed)
  D = {'rng':rng}
  # individual-level
  D['new_ind'] = stats.pois(rng=rng['ind'],m=1/365/model.adur)
  D['ptr_max'] = stats.geom(rng=rng['ind'],m=2)
  D['ptr_r0']  = stats.exp (rng=rng['ind'],m=1/90)
  D['cdm_p0']  = stats.unif(rng=rng['ind'],l=0,u=1)
  D['dep_r0']  = stats.exp (rng=rng['ind'],m=1/365)
  D['dep_x0']  = stats.exp (rng=rng['ind'],m=1/365)
  D['vio_r0']  = stats.exp (rng=rng['ind'],m=1/365)
  # partnership-level
  D['ptr_dur'] = stats.exp(rng=rng['ptr'],m=90)
  # D['age_pref'] = stats.norm(rng=rng['ptr'],m=0,sd=5.1)
  return D

def init_ages(D,n):
  return model.amin + model.adur * D['rng']['ind'].random(n)
