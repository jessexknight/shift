import numpy as np
import matplotlib.pyplot as plt
from model import net,par,out
from model import adur,z1y,z3mf

zs = np.arange(0,z1y*adur)
Ps = par.get_n_sample(range(1))
Ns = net.run_n(Ps,zs)
Y = out.Survey(Ns[0].I)
Y['ptr_tot'] = out.n_ptr(Y.I,'tot')
plt.plot(Y['age'],Y['ptr_tot'],'.')
plt.savefig('pyplots.pdf')
