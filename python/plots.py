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

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        properties['uncorrected_hardest_pts'].append( float( numbers[3] ) )
        properties['corrected_hardest_pts'].append( float( numbers[4] ) )
        properties['prescales'].append( int( numbers[5] ) )
        properties['trigger_names'].append(  numbers[6] )
        properties['zg_05'].append( float( numbers[7] ) )
        properties['zg_1'].append( float( numbers[8] ) )
        properties['zg_2'].append( float( numbers[9] ) )

    except:
      pass

  return properties
    


def plot_pts(pTs, prescales):

  plt.hist(np.array(pTs), 200, normed=1, weights=prescales, facecolor='red', alpha=1.0)

  plt.autoscale(True)
  # plt.xlim(0, 1000)

  plt.xlabel('$p_T$ (GeV)')
  plt.grid(True)
  plt.yscale('log', nonposy='clip')

  plt.show()


def plot_zg(zgs, prescales):

  plt.hist(np.array(zgs), 50, normed=1, weights=prescales, facecolor='red', alpha=1.0)

  plt.autoscale(True)
  plt.xlim(0, 0.6)

  plt.xlabel('$z_g$')
  plt.grid(True)

  plt.show()


properties = parse_file(input_analysis_file)

plot_pts(properties['uncorrected_hardest_pts'], properties['prescales'])

# plot_zg(properties['zg_05'], properties['prescales'])