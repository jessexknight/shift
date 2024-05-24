import numpy as np
from sim import net,par,srv
from sim import adur,z3m,z1y,k3m
from utils import pmap
import seaborn as sb

def run(P,zs):
  N = net.Network(P)
  Q = srv.Survey(N,srv.sfun,srv.qfun,z0=z1y,zr=z1y)
  N.run(zs)
  return Q

zs = np.arange(0,1+7*z1y)
Ps = par.get_n_all(range(7),n=256)
Qs = pmap(run,Ps,zs=zs)
X  = srv.concat(Qs)

ax = sb.lineplot(X,y='ptr_n3m',x='z',hue='vio_a3m')
# H = srv.gist(X,'ptr_n3m',range(15),['z','seed','vio_a3m'])
# ax = sb.lineplot(H,y='ptr_n3m',x='b',hue='vio_a3m')
ax.figure.savefig('pyplots.pdf')
