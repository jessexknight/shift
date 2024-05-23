import numpy as np
from sim import net,par,out
from sim import adur,z1y,k3m
import seaborn as sb

zs = np.arange(0,3*z1y)
seeds = range(7)
Ps = par.get_n_all(seeds,n=1000)
Ns = net.run_n(Ps,zs)

M = out.Meta(Ns,['seed'])
M['vio'] = 0 < out.n_vio(M.I,      **k3m(zs[-1]))
M['ptr'] = 0 + out.n_ptr(M.I,'tot',**k3m(zs[-1]))
M['dep'] = 0 + out.attr (M.I,'dep')
M['cdm'] = 0 + out.attr (M.I,'cdm')

# ax = sb.violinplot(M.X,y='ptr',x='vio')
H = out.gist(M.X,'cdm',[0,1],['dep','seed'])
ax = sb.boxplot(H,y='cdm',x='b',hue='dep')
ax.figure.savefig('pyplots.pdf')
