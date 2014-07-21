import numpy as np
import pylab as plt

c = 1.
d = 0.001
x = np.linspace(0.00001,0.01,200)
y = d*np.exp(-c*x)/x

plt.plot(x,y,'k-',linewidth=2)
plt.show()
