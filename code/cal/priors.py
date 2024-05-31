from sim import par
from utils import plot
# generate 1000 parameter sets
Ps = par.get_cal_sample(1000,seed=0)
Ps = par.get_n_all([0],Ps=Ps)
# plots
plot.distr(Ps,'cdm_p0'),           plot.save('par','cdm_p0')
plot.distr(Ps,'ptr_max',cum=True), plot.save('par','ptr_max')
plot.distr(Ps,'ptr_dur',cum=True), plot.save('par','ptr_dur')
plot.distr(Ps,'ptr_r0', cum=True), plot.save('par','ptr_r0')
plot.distr(Ps,'dep_r0', cum=True), plot.save('par','dep_r0')
plot.distr(Ps,'dep_x0', cum=True), plot.save('par','dep_x0')
plot.distr(Ps,'vio_r0', cum=True), plot.save('par','vio_r0')
