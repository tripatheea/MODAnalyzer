
import numpy as np
from rootpy.plotting import Hist, HistStack, Legend, Canvas

import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt





np.random.seed(42)


signal_obs = 126 + 10 * np.random.randn(100)

# create histograms
h1 = Hist(30, 40, 200, title='Background', markersize=0)
h3 = h1.Clone(title='Data')
h3.markersize = 1.2


h3.FillRandom('landau', 1000)
map(h3.Fill, signal_obs)




axes = plt.axes()


rplt.errorbar(h3, xerr=False, emptybins=False, axes=axes)

plt.show()
