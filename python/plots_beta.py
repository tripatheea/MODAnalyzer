# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat
from __future__ import division

from subprocess import call

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

from collections import defaultdict

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

def get_version(input_file):

  f = open(input_file, 'r')
  lines = f.read().split("\n")

  properties = defaultdict(list)

  for line in lines:
    numbers = line.split()
    if numbers[0] == "%":
      return numbers[1] + " " + numbers[2] 


def parse_file(input_file, pT_lower_cut = 0.00, jet_quality_level=1):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  # Hardest_pT Corr_Hardest_pT Prescale Trigger_Name zg_05 dr_05 mu_05 zg_1 dr_1 mu_1 zg_2 dr_2 mu_2
  
  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#" and float(numbers[1]) > pT_lower_cut and int(numbers[3]) == 1 and int(numbers[4]) >= jet_quality_level:
        properties['corrected_hardest_pT'].append( float( numbers[1] ) )

        properties['prescales'].append( float( numbers[2] ) )

        properties['trig_jet_matched'].append( int( numbers[3] ) )
        properties['jet_quality'].append( int( numbers[4] ) )

        properties['rho_10'].append( float( numbers[5] ) )
        properties['zg_10'].append( float( numbers[6] ) )
        properties['rho_11'].append( float( numbers[7] ) )
        properties['zg_11'].append( float( numbers[8] ) )
        properties['rho_12'].append( float( numbers[9] ) )
        properties['zg_12'].append( float( numbers[10] ) )
        properties['rho_13'].append( float( numbers[11] ) )
        properties['zg_13'].append( float( numbers[12] ) )
        properties['rho_14'].append( float( numbers[13] ) )
        properties['zg_14'].append( float( numbers[14] ) )
        properties['rho_15'].append( float( numbers[15] ) )
        properties['zg_15'].append( float( numbers[16] ) )
        properties['rho_16'].append( float( numbers[17] ) )
        properties['zg_16'].append( float( numbers[18] ) )
        properties['rho_17'].append( float( numbers[19] ) )
        properties['zg_17'].append( float( numbers[20] ) )
        properties['rho_18'].append( float( numbers[21] ) )
        properties['zg_18'].append( float( numbers[22] ) )
        properties['rho_19'].append( float( numbers[23] ) )
        properties['zg_19'].append( float( numbers[24] ) )
        properties['rho_20'].append( float( numbers[25] ) )
        properties['zg_20'].append( float( numbers[26] ) )
        

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




def plot_rho(pT_lower_cut=150):
  properties = parse_file(input_analysis_file, pT_lower_cut)

  # FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3
  jet_quality_level = 1

  prescales = properties['prescales']
  trig_jet_matched = properties['trig_jet_matched']
  jet_quality = properties['jet_quality']
  
  rho, weights = defaultdict(list), defaultdict(list)
  # cutoffs = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
  # cutoffs = ['0.10', '0.11', '0.12', '0.13', '0.14', '0.15', '0.16', '0.17', '0.18', '0.19', '0.20']
  cutoffs = ['0.13', '0.14', '0.15', '0.16']

  for cutoff in cutoffs:
    for i in range(0, len(properties['zg_' + str(cutoff).replace('.', '')[1:]])):
      if properties[ 'zg_' + str(cutoff).replace('.', '')[1:] ][i] > float(cutoff):  
        rho[str(cutoff)].append(properties['rho_' + str(cutoff).replace('.', '')[1:]][i])
        weights[str(cutoff)].append(prescales[i])

  
  logged_rho = defaultdict(list)
  for r in rho:
    logged_rho[r] = np.log(rho[r])
    
    plt.hist(logged_rho[r], weights=weights[r], label="z\_cut=" + str(r), normed=1, histtype='step', lw=5, bins=50)

  plt.autoscale(True)
  # plt.gca().set_ylim(0, 10)

  legend = plt.legend(loc=2, frameon=0)
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  plt.gca().legend([extra], [r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$"], loc=6, frameon=0, fontsize=60)




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
  
  # plt.savefig("plots/" + get_version(input_analysis_file) + "/rho/rho_" + str(pT_lower_cut) + ".pdf")
  plt.savefig("plots/" + "Version 4" + "/rho/rho_" + str(pT_lower_cut) + ".pdf")
  plt.clf()



plot_rho(pT_lower_cut=100)
plot_rho(pT_lower_cut=150)
plot_rho(pT_lower_cut=200)
plot_rho(pT_lower_cut=250)
plot_rho(pT_lower_cut=300)


call(["python", "/home/aashish/root/macros/MODAnalyzer/utilities/sync_plots.py"])