import numpy as np
import matplotlib.pyplot as plt
from model import net,par
import model

n = 1000
ks = np.arange(0,52*20)
D = par.distrs()
N = net.Network(D)
N.add_inds(n=n,ages=par.init_ages(D,n))
N.run(ks)

print(np.mean([len(I.P)/I.ptr_max for I in N.I]))
# plt.hist([len(I.P) for I in N.I])
# plt.hist([P.cdm for I in N.I for P in I.P])
plt.hist([I.age for I in N.I])
plt.savefig('pyplots.pdf')



