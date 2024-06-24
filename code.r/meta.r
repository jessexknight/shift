
# time stuff
dtz  =  7 # days in 1 timestep
z1y  = 52 # timesteps in 1 year
z3m  = 13 # timesteps in 3 months
z6m  = 26 # timesteps in 6 months
amin = 15 # age of cohort entry
amax = 50 # age of cohort exit
adur = amax - amin # duration in cohort

# event types
evts = c('vio','dep.o','dep.x','ptr.o','sex','cdm')
names(evts) = evts
