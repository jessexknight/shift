
# constants
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
  'sex',  'cdm')
names(evts) = evts

# source model code
source('utils.r')
source('sim/par.r')
source('sim/sim.r')
source('sim/srv.r')
source('sim/plot.r')

# unique ID
uid = Sys.Date()
