# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat
from __future__ import division

import sys
import math
from collections import defaultdict

# matplotlib
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# RootPy
from rootpy.plotting import Hist, HistStack, Legend
import rootpy.plotting.root2matplotlib as rplt

# Stuff for calculating areas.
from scipy.integrate import simps
from numpy import trapz


from matplotlib import gridspec

from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage, AnnotationBbox



from scipy.stats import binned_statistic

import rootpy.plotting.views

input_analysis_file = sys.argv[1]


mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
plt.rc('font', family='serif', size=43)


def parse_file(input_file, pT_lower_cut = 0.00, pfc_pT_cut = 0.00):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  # Hardest_pT Corr_Hardest_pT Prescale Trigger_Name zg_05 dr_05 mu_05 zg_1 dr_1 mu_1 zg_2 dr_2 mu_2
  
  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        if (float(numbers[3]) > pT_lower_cut) and (float(numbers[17]) > pfc_pT_cut):
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
          properties['hardest_pfc_pdgid'].append( float( numbers[16] ) )
          properties['hardest_pfc_pt'].append( float( numbers[17] ) )

          properties['zg_05_pt_1'].append( float( numbers[18] ) )
          properties['zg_1_pt_1'].append( float( numbers[19] ) )
          properties['zg_2_pt_1'].append( float( numbers[20] ) )

          properties['zg_05_pt_2'].append( float( numbers[21] ) )
          properties['zg_1_pt_2'].append( float( numbers[22] ) )
          properties['zg_2_pt_2'].append( float( numbers[23] ) )

          properties['zg_05_pt_3'].append( float( numbers[24] ) )
          properties['zg_1_pt_3'].append( float( numbers[25] ) )
          properties['zg_2_pt_3'].append( float( numbers[26] ) )

          properties['zg_05_pt_5'].append( float( numbers[27] ) )
          properties['zg_1_pt_5'].append( float( numbers[28] ) )
          properties['zg_2_pt_5'].append( float( numbers[29] ) )

          properties['zg_05_pt_10'].append( float( numbers[30] ) )
          properties['zg_1_pt_10'].append( float( numbers[31] ) )
          properties['zg_2_pt_10'].append( float( numbers[32] ) )

          properties['zg_charged_05'].append( float( numbers[33] ) )
          properties['zg_charged_1'].append( float( numbers[34] ) )
          properties['zg_charged_2'].append( float( numbers[35] ) )

    except:
      pass

  return properties
    
def parse_mc_file(input_file, pT_lower_cut = 0.00, pfc_pT_cut = 0.00):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      properties['zg_05'].append( float( numbers[0] ) )
      properties['zg_1'].append( float( numbers[1] ) )
      properties['zg_2'].append( float( numbers[2] ) )
    except:
      pass
  
  return properties






def plot_zg_th_mc_data(zg_cut, zg_filename, ratio_denominator="theory", data=True, mc=True, theory=True):
  pT_lower_cut = 150
  pfc_pT_cut = 0

  properties_pythia = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_pythia_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)

  zg_pythias = properties_pythia[zg_filename]

  # zg_pythia_hist = Hist(50, 0, 0.6, markersize=5.0, color='blue', linewidth=5)
  # map(zg_pythia_hist.Fill, zg_pythias)
  # pythia_plot = rplt.hist(zg_pythia_hist)

  plt.hist(zg_pythias, normed=1, bins=50, histtype='step')
    
  
  plt.gca().xaxis.set_tick_params(width=5, length=10)
  plt.gca().yaxis.set_tick_params(width=5, length=10)


  plt.autoscale(True)

  plt.show()
  




def parse_theory_file(input_file):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  properties = defaultdict(list)

  points = defaultdict(list)
  for line in lines:
    try:
      numbers = line.split()
      try:
        # Areas range from 4 to 1, 5 to 2 and 6 to 3.
        # When returned, will be 1 to 2, 2 to 3 and 3 to 4.
        
        points[float(numbers[0])].append(float(numbers[4]))
        points[float(numbers[0])].append(float(numbers[1]))
        
        points[float(numbers[0])].append(float(numbers[5]))
        points[float(numbers[0])].append(float(numbers[2]))

        points[float(numbers[0])].append(float(numbers[6]))
        points[float(numbers[0])].append(float(numbers[3]))

      except ValueError:
        pass
    except:
      pass

  return points




  








# plot_zg_th_mc_data('0.1', 'zg_1', 'theory', )
# plot_zg_th_mc_data('0.1', 'zg_1', 'data', theory=1, mc=0, data=0)
plot_zg_th_mc_data('0.1', 'zg_1', 'data', theory=1, mc=1, data=0)
# plot_zg_th_mc_data('0.1', 'zg_1', 'data', theory=1, mc=1, data=1)


