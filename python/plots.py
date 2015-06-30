# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat

import sys
from collections import defaultdict

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

input_analysis_file = sys.argv[1]



def parse_file(input_file):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  # Hardest_pT Corr_Hardest_pT Prescale Trigger_Name zg_05 zg_1 zg_2
  
  properties = defaultdict(list)
  '''
  hardest_pts = []
  hardest_corrected_pts = []
  prescales = []
  trigger_names = []
  zg_05s = []
  zg_1s = []
  zg_2s = []
  '''

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        properties['hardest_pts'].append( float( numbers[3] ) )
        properties['prescales'].append( int( numbers[5] ) )

    except:
      pass

  return properties
    



properties = parse_file(input_analysis_file)


print properties['hardest_pts']


x = np.array(properties['hardest_pts'])


# the histogram of the data
n, bins, patches = plt.hist(x, 50, normed=1, weights=properties['prescales'], facecolor='red', alpha=0.75)

# plt.axis([0, 1500, 0, 0.03])
plt.autoscale(True)

plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.grid(True)
plt.yscale('log', nonposy='clip')

plt.show()