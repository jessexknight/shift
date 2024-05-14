import numpy as np
import matplotlib.pyplot as plt
from model import net,par
import model

n = 1000
zs = np.arange(0,52*model.adur)
D = par.distrs(seed=0)
N = net.Network(D)
N.add_inds(n=n,ages=par.init_ages(D,n))
N.run(zs)

# print(N.I[0].logs)
# print(np.mean([len(I.P)/I.ptr_max for I in N.I]))
# print(np.mean([len(I.logs['vio']) for I in N.I]))
# print(np.mean([I.depressed for I in N.I]))
# plt.hist([len(I.P) for I in N.I])
# plt.hist([P.cdm for I in N.I for P in I.P])
# plt.hist([I.age for I in N.I])
# p3m_ptr = lambda I: len(I.P) + sum(z+model.z3m > zs[-1] for z in I.logs['end_ptr'])
# plt.hist([p3m_ptr(I) for I in N.I])
# plt.savefig('pyplots.pdf')
