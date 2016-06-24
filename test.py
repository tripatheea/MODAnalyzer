import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 100)
y1 = 2 * x
y2 = 3 * x

xthresh = 4.5
diff = np.abs(y1 - y2)
below = diff < xthresh
above = diff >= xthresh


print above

print x[below]

# Plot lines below threshold as dotted...
plt.plot(x[below], y1[below], 'b--')
plt.plot(x[below], y2[below], 'g--')

# Plot lines above threshold as solid...
plt.plot(x[above], y1[above], 'b-')
plt.plot(x[above], y2[above], 'g-')
