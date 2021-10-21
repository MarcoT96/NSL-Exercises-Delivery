from math import pi, sqrt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

x=np.linspace(0, 1, 1000)
f1=(pi/2) * np.cos(pi/2 * x)
f2=2*(1-x)

plt.figure(figsize=[10.4, 7.8]) #deafult (6.4,4.8)
plt.xlabel('$x$', fontsize=22)
plt.ylabel('$f(x)$', fontsize=22)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.plot(x, f1, label='$\\frac{\pi}{2} \\ \cos(\\frac{\pi}{2} x)$')
plt.plot(x, f2, label='$2(1-x)$')
plt.legend(fontsize=22, loc='best')
plt.savefig('d(x).png')
