import numpy as np
from model import net,par,out
from model import adur,z1y,k3m
import seaborn as sb

zs = np.arange(0,3*z1y)
seeds = range(7)
Ps = par.get_n_all(seeds,n=100)
Ns = net.run_n(Ps,zs)

M = out.Meta(Ns,['seed'])
M['vio'] = 0 < out.n_vio(M.I,      **k3m(zs[-1]))
M['dep'] = 0 < out.n_dep(M.I,'tot',**k3m(zs[-1]))
M['ptr'] =     out.n_ptr(M.I,'tot',**k3m(zs[-1]))

H = out.gist(M.X,'ptr',range(30),['vio','seed'])
ax = sb.lineplot(H,y='ptr',x='b',hue='vio')
ax.figure.savefig('pyplots.pdf')
