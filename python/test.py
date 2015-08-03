from __future__ import division

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D

import math

x = np.linspace(1., 10., 100000)
y = 1 / x

z = []
for i in range(0, len(x)):
  z.append(math.log(x[i], 10))



plt.xscale('log')

plt.plot(z, y)



plt.show()