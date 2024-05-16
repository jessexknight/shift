import numpy as np
from model import net,par
from model import adur,z1y

zs = np.arange(0,z1y*adur)
Ps = par.get_n_sample(range(1),n=100)
Ns = net.run_n(Ps,zs,para=False)
