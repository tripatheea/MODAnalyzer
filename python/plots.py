# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat

import sys
from collections import defaultdict

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

input_analysis_file = sys.argv[1]



def parse_file(input_file, pT_lower_cut = 0.00):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  # Hardest_pT Corr_Hardest_pT Prescale Trigger_Name zg_05 dr_05 mu_05 zg_1 dr_1 mu_1 zg_2 dr_2 mu_2
  
  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        if float(numbers[3]) > pT_lower_cut:
          properties['uncorrected_hardest_pts'].append( float( numbers[3] ) )
          properties['corrected_hardest_pts'].append( float( numbers[4] ) )
          properties['prescales'].append( int( numbers[5] ) )
          properties['trigger_names'].append(  numbers[6] )
          properties['zg_05'].append( float( numbers[7] ) )
          properties['dr_05'].append( float( numbers[8] ) )
          properties['mu_05'].append( float( numbers[9] ) )
          properties['zg_1'].append( float( numbers[10] ) )
          properties['dr_1'].append( float( numbers[11] ) )
          properties['mu_1'].append( float( numbers[12] ) )
          properties['zg_2'].append( float( numbers[13] ) )
          properties['dr_2'].append( float( numbers[14] ) )
          properties['mu_2'].append( float( numbers[15] ) )

    except:
      pass

  return properties
    


def plot_pts():
  properties = parse_file(input_analysis_file)

  pTs = properties['uncorrected_hardest_pts']
  corrected_pTs = properties['corrected_hardest_pts']
  prescales = properties['prescales']

  plt.hist(np.array(pTs), 500, normed=1, label="Uncorrected", weights=prescales, log=1, histtype='step', color='blue') 
  plt.hist(np.array(corrected_pTs), 500, normed=1, label="Corrected", weights=prescales, log=1, histtype='step', color='red') 

  plt.autoscale(True)
  plt.xlim(0, 2000)

  plt.legend()
  plt.xlabel('$p_T$ (GeV)')
  plt.grid(True)

  plt.show()


def plot_zg():
  properties = parse_file(input_analysis_file, pT_lower_cut = 153)

  zgs = [properties['zg_05'], properties['zg_1'], properties['zg_2']]
  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['0.05', '0.1', '0.2']

  i = 0
  for zg in zgs:
    plt.hist(np.array(zg), 100, normed=1, label=labels[i], weights=prescales, facecolor=colors[i], histtype='step')
    i += 1

  plt.autoscale(True)
  plt.xlim(0, 0.6)

  plt.legend()
  plt.xlabel('Symmetry Measure(z)')
  plt.grid(True)

  plt.show()

def plot_dr():
  properties = parse_file(input_analysis_file, pT_lower_cut = 153)

  drs = [properties['dr_05'], properties['dr_1'], properties['dr_2']]
  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['0.05', '0.1', '0.2']

  i = 0
  for dr in drs:
    plt.hist(np.array(dr), 200, normed=1, label=labels[i], weights=prescales, facecolor=colors[i], histtype='step')
    i += 1

  plt.autoscale(True)
  plt.xlim(0, 0.5)

  plt.legend()
  plt.xlabel('$\Delta$R between Subjets')
  plt.grid(True)

  plt.show()

def plot_mu():
  properties = parse_file(input_analysis_file, pT_lower_cut = 153)

  mus = [properties['mu_05'], properties['mu_1'], properties['mu_2']]
  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['0.05', '0.1', '0.2']

  i = 0
  for mu in mus:
    plt.hist(np.array(mu), 200, normed=1, label=labels[i], weights=prescales, facecolor=colors[i], histtype='step')
    i += 1

  plt.autoscale(True)
  plt.xlim(-0.01, 1)

  plt.legend()
  plt.xlabel('Mass Drop($\mu$)')
  plt.grid(True)

  plt.show()

# plot_pts()
plot_zg()
plot_dr()
plot_mu()