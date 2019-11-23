import numpy as np
import matplotlib.pyplot as plt

q = 1
nu = 0.3
b = 1

z = np.linspace(0,10,100)
sz = q*(1 - z**3/((b**2 + z**2)**(3/2)))
sr = (q/2)*(1 +2*nu - 2*(1+nu)*z/np.sqrt(b**2 + z**2) + z**3/((b**2 + z**2)**(3/2)))


fig, ax = plt.subplots(nrows=1, ncols=2)
ax[0].plot(z, sz)
ax[0].set_title(r'$\sigma_z$')

ax[1].plot(z, sr)
ax[1].set_title(r'$\sigma_r$')
plt.show()
