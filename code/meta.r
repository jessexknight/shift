# source model code
source('utils.r')
source('par.r')
source('sim.r')
source('srv.r')
source('plot.r')

# unique ID
uid = Sys.Date()

# time stuff
dtz  =  7 # days in 1 timestep
z1y  = 52 # timesteps in 1 year
z3m  = 13 # timesteps in 3 months
z6m  = 26 # timesteps in 6 months
t1y  = dtz * z1y # days in 1 year
amin = 10 # age of cohort entry
amax = 60 # age of cohort exit
adur = amax - amin # duration in cohort
eps  = 1e-12 # a small number

# event types
evts = c(
  'vio',
  'dep_o','dep_x',
  'haz_o','haz_x',
  'ptr_o','ptr_x',
  'sex',
  'cdm')
names(evts) = evts
