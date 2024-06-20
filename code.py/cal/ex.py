from sim import par,net,srv,z1y,a5g
from utils import plot,pmap,graph,rint,gint
import seaborn as sb
# years
yb = 35 # burn-in
ys =  3 # survey
# function for parallel
def run(P,zs):
  N = net.Network(P)
  Q = srv.Survey(N,srv.sfun,srv.qfun,z0=(yb+1)*z1y,zr=z1y)
  N.run(zs)
  return N,Q
# setup & run
zs = range(0,1+(yb+ys)*z1y)
Ps = par.get_n_all(range(7),n=1000)
Ns,Qs = zip(*pmap(run,Ps,zs=zs))
# create Meta & add columns
M = srv.Meta(Qs)
M['year']  = rint(M['z']/z1y)
M['age_1'] = rint(M['age'])
M['age_5'] = gint(M['age'],a5g)
M.reindex()
# plots: attribs vs age & time
plot.line(M.X,y='ptr_n3m',x='age_1',hue='year',err='ci',ylim=(0,3)),  plot.save('ex','ptr_m.a1.t')
plot.line(M.X,y='vio_a3m',x='age_1',hue='year',err='ci',ylim=(0,1)),  plot.save('ex','vio_m.a1.t')
plot.line(M.X,y='dep',    x='age_1',hue='year',err='ci',ylim=(0,.5)), plot.save('ex','dep_m.a1.t')
plot.line(M.X,y='cdm',    x='age_1',hue='year',err='ci',ylim=(0,.5)), plot.save('ex','cdm_m.a1.t')
plot.line(M.X,y='ptr_n3m',x='age_5',hue='year',err='ci',ylim=(0,3)),  plot.save('ex','ptr_m.a5.t')
plot.line(M.X,y='vio_a3m',x='age_5',hue='year',err='ci',ylim=(0,1)),  plot.save('ex','vio_m.a5.t')
plot.line(M.X,y='dep',    x='age_5',hue='year',err='ci',ylim=(0,.5)), plot.save('ex','dep_m.a5.t')
plot.line(M.X,y='cdm',    x='age_5',hue='year',err='ci',ylim=(0,.5)), plot.save('ex','cdm_m.a5.t')
# plot: attribs vs attribs
plot.gist(M.X,y='dep',    b=range(2), hue='vio_a3m',type='violin'), plot.save('ex','dep.h.vio')
plot.gist(M.X,y='ptr_n3m',b=range(11),hue='year',   err='pi'),      plot.save('ex','ptr.h.t')
plot.gist(M.X,y='ptr_n3m',b=range(11),hue='dep',    err='pi'),      plot.save('ex','ptr.h.dep')
plot.gist(M.X,y='ptr_n3m',b=range(11),hue='vio_a3m',err='pi'),      plot.save('ex','ptr.h.vio')
plot.gist(M.X,y='cdm',    b=range(2), hue='dep',    type='violin'), plot.save('ex','cdm.h.dep')
plot.gist(M.X,y='cdm',    b=range(2), hue='vio_a3m',type='violin'), plot.save('ex','cdm.h.vio')
