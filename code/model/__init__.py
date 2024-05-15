dtz  =  7 # days per timestep
z3m  = 13 # timesteps per 3 months
z1y  = 52 # timesteps per 1 year
z3mf = lambda zf: {'z0':zf-z3m,'zf':zf}
z1yf = lambda zf: {'z0':zf-z1y,'zf':zf}
amin = 15 # age of model entry
amax = 50 # age of model exit
adur = 35 # years in model
logs = (  # types of logs
  'begin_ptr','end_ptr',
  'begin_dep','end_dep',
  'vio',
  'exit',
)
