# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat
from __future__ import division

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator

from sets import Set

import sys
import math
from collections import defaultdict

# matplotlib
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm


# RootPy
from rootpy.plotting import Hist, HistStack, Legend
import rootpy.plotting.root2matplotlib as rplt
from rootpy.plotting import Hist2D


# Stuff for calculating areas.
from scipy.integrate import simps
from scipy import interpolate
from scipy import optimize

from numpy import trapz


from matplotlib import gridspec

from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png

from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from scipy.stats import norm
from scipy.stats import gamma


from scipy.stats import binned_statistic

import rootpy.plotting.views

input_analysis_file = sys.argv[1]


mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)



def parse_file(input_file, pT_lower_cut = 0.00, pfc_pT_cut = 0.00, pT_upper_cut = 20000.00):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  # Hardest_pT Corr_Hardest_pT Prescale Trigger_Name zg_05 dr_05 mu_05 zg_1 dr_1 mu_1 zg_2 dr_2 mu_2
  
  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        properties['prescales'].append( float( numbers[1] ) )

        properties['rho_05'].append( float( numbers[2] ) )
        properties['rho_10'].append( float( numbers[3] ) )
        properties['rho_15'].append( float( numbers[4] ) )
        properties['rho_20'].append( float( numbers[5] ) )
        properties['rho_25'].append( float( numbers[6] ) )
        properties['rho_30'].append( float( numbers[7] ) )
        properties['rho_35'].append( float( numbers[8] ) )
        properties['rho_40'].append( float( numbers[9] ) )
        properties['rho_45'].append( float( numbers[10] ) )
        properties['rho_50'].append( float( numbers[11] ) )
        

    except:
      if len(numbers) != 0:
        # print "Some kind of error occured while parsing the given file!"
        # print numbers
        # print
        pass

  return properties


def calculate_distribution(y_cut, f):
  C_A = 3
  C_f = 4 / 3
  n_f = 5

  beta_o = (11 * C_A - 2 * n_f) / 12

  try:
    distribution = beta_o - f * C_f * math.log( math.exp( - 3/4) / y_cut, math.e ) - (1 - f) * C_A * math.log( math.exp( - beta_o / C_A) / y_cut, math.e )
  except ZeroDivisionError:
    distribution = None

  return distribution




def plot_rho():
  pT_lower_cut = 100
  properties = parse_file(input_analysis_file, pT_lower_cut)

  prescales = properties['prescales']
  
  rho = {}

  rho['rho_05'] = properties['rho_05']
  rho['rho_10'] = properties['rho_10']
  rho['rho_15'] = properties['rho_15']
  rho['rho_20'] = properties['rho_20']
  rho['rho_25'] = properties['rho_25']
  rho['rho_30'] = properties['rho_30']
  rho['rho_35'] = properties['rho_35']
  rho['rho_40'] = properties['rho_40']
  rho['rho_45'] = properties['rho_45']
  rho['rho_50'] = properties['rho_50']

  for r in rho:
    rho[r] = np.log(rho[r])
    plt.hist(rho[r], label="z\_cut=0." + str(r)[4:], normed=1, histtype='step', lw=5, bins=200)

  plt.autoscale(True)
  # plt.gca().set_ylim(0, 10)

  plt.legend(loc=2)

  plt.gca().set_xlabel("$\\rho = m^2 / (p_T^2~R^2)$", fontsize=75, labelpad=50)


  plt.gcf().set_size_inches(30, 30, forward=1)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

 
  # plt.gca().xaxis.set_minor_locator(MultipleLocator(5))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.05))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)


  # x = np.linspace(math.log(0.000000000001, math.e), math.log(0.5, math.e), 10)
  # labels = [str(round(math.exp(i), 4)) for i in x]
  # plt.xticks(x, labels)

  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  plt.savefig("plots/rho.pdf")



plot_rho()