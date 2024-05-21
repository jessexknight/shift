import numpy as np
from model import net,par,out
from model import adur,z1y,z3mf
import seaborn as sb

zs = np.arange(0,3*adur)
seeds = range(7)
Ps = par.get_n_all(seeds,n=100,case='high',ptr_r0_lm=-3,ptr_dur_lm=2)+\
     par.get_n_all(seeds,n=100,case='low', ptr_r0_lm=-5,ptr_dur_lm=6)
Ns = net.run_n(Ps,zs)
M = out.Meta(Ns,['seed','case'])
M['p3m_ptr'] = out.n_ptr(M.I,'tot',**z3mf(zs[-1]))

H = out.gist(M.X,'p3m_ptr',range(15),['case','seed'])
ax = sb.boxplot(H,y='p3m_ptr',x='b',hue='case',fill=False,dodge=False)
ax.set(xlabel='Partners P3M',ylabel='Individuals')
ax.figure.savefig('pyplots.pdf')
