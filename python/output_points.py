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




from rootpy.io import root_open, DoesNotExist


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
          properties['zg_1'].append( float( numbers[10] ) )
          properties['zg_2'].append( float( numbers[13] ) )
          
    except:
      pass

  return properties
    




def plot_zg_th_mc_data(zg_cut, zg_filename):
  pT_lower_cut = 150
  pfc_pT_cut = 0

  properties = parse_file(input_analysis_file, pT_lower_cut, pfc_pT_cut)
  
  zg_data = properties[zg_filename]
  prescales = properties['prescales']

  zg_data_hist = Hist(60, 0.0, 0.6, markersize=2.5, color='black')
  bin_width_data = (zg_data_hist.upperbound() - zg_data_hist.lowerbound()) / zg_data_hist.nbins()

  map(zg_data_hist.Fill, zg_data, prescales)
  
  zg_data_hist.Scale(1.0 / ( zg_data_hist.GetSumOfWeights() * bin_width_data ))


  plotline, caplines, barlinecols = rplt.errorbar(zg_data_hist, xerr=1, yerr=1, emptybins=False, linewidth=5, alpha=1.0)

  x = plotline.get_xdata()
  y = plotline.get_ydata()

  x_errors = []
  y_errors = []


  for x_segment in barlinecols[0].get_segments():
    x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)


  for y_segment in barlinecols[1].get_segments():
    y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)


  f = open('data_' + zg_filename + '.txt','w')
  
  for i in range(0, len(x)):
    f.write(str(x[i]) + '\t' + str(y[i]) + '\t' + str(x_errors[i]) + '\t' + str(y_errors[i]) + '\n') 
  
  f.close() # you can omit in most cases as the destructor will call it

  plt.autoscale()

  plt.show()


  









plot_zg_th_mc_data('0.05', 'zg_05')
plot_zg_th_mc_data('0.1', 'zg_1')
plot_zg_th_mc_data('0.2', 'zg_2')



