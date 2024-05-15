import numpy as np
import matplotlib.pyplot as plt
from model import net,par,out,z1y,z3mf

zs = np.arange(0,z1y*2)
Ps = par.get_n_sample(range(7))
Ns = net.run_n(Ps,zs)
for N in Ns:
  plt.hist(out.n_ptr(N.I,'new',**z3mf(zs[-1])),
    bins=range(25),alpha=.2,color='black')
plt.savefig('pyplots.pdf')
