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
        properties['no_of_const'].append( int( numbers[5] ) )

        properties['rho_10'].append( float( numbers[6] ) )
        properties['zg_10'].append( float( numbers[7] ) )
        properties['rho_11'].append( float( numbers[8] ) )
        properties['zg_11'].append( float( numbers[9] ) )
        properties['rho_12'].append( float( numbers[10] ) )
        properties['zg_12'].append( float( numbers[11] ) )
        properties['rho_13'].append( float( numbers[12] ) )
        properties['zg_13'].append( float( numbers[13] ) )
        properties['rho_14'].append( float( numbers[14] ) )
        properties['zg_14'].append( float( numbers[15] ) )
        properties['rho_15'].append( float( numbers[16] ) )
        properties['zg_15'].append( float( numbers[17] ) )
        properties['rho_16'].append( float( numbers[18] ) )
        properties['zg_16'].append( float( numbers[19] ) )
        properties['rho_17'].append( float( numbers[20] ) )
        properties['zg_17'].append( float( numbers[21] ) )
        properties['rho_18'].append( float( numbers[22] ) )
        properties['zg_18'].append( float( numbers[23] ) )
        properties['rho_19'].append( float( numbers[24] ) )
        properties['zg_19'].append( float( numbers[25] ) )
        properties['rho_20'].append( float( numbers[26] ) )
        properties['zg_20'].append( float( numbers[27] ) )
        
        properties['chrg_rho_10'].append( float( numbers[28] ) )
        properties['chrg_zg_10'].append( float( numbers[29] ) )
        properties['chrg_rho_11'].append( float( numbers[30] ) )
        properties['chrg_zg_11'].append( float( numbers[31] ) )
        properties['chrg_rho_12'].append( float( numbers[32] ) )
        properties['chrg_zg_12'].append( float( numbers[33] ) )
        properties['chrg_rho_13'].append( float( numbers[34] ) )
        properties['chrg_zg_13'].append( float( numbers[35] ) )
        properties['chrg_rho_14'].append( float( numbers[36] ) )
        properties['chrg_zg_14'].append( float( numbers[37] ) )
        properties['chrg_rho_15'].append( float( numbers[38] ) )
        properties['chrg_zg_15'].append( float( numbers[39] ) )
        properties['chrg_rho_16'].append( float( numbers[40] ) )
        properties['chrg_zg_16'].append( float( numbers[41] ) )
        properties['chrg_rho_17'].append( float( numbers[42] ) )
        properties['chrg_zg_17'].append( float( numbers[43] ) )
        properties['chrg_rho_18'].append( float( numbers[44] ) )
        properties['chrg_zg_18'].append( float( numbers[45] ) )
        properties['chrg_rho_19'].append( float( numbers[46] ) )
        properties['chrg_zg_19'].append( float( numbers[47] ) )
        properties['chrg_rho_20'].append( float( numbers[48] ) )
        properties['chrg_zg_20'].append( float( numbers[49] ) )

        properties['chrg_multip'].append( int( numbers[51] ) )
        properties['neu_had_frac'].append( float( numbers[52] ) )
        properties['neu_em_frac'].append( float( numbers[53] ) )
        properties['chrg_had_frac'].append( float( numbers[54] ) )
        properties['chrg_em_frac'].append( float( numbers[55] ) )

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






def plot_charged_rho(pT_lower_cut=150):
  properties = parse_file(input_analysis_file, pT_lower_cut)

  # FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3
  jet_quality_level = 1

  prescales = properties['prescales']
 
  rho, weights = defaultdict(list), defaultdict(list)
  # cutoffs = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
  # cutoffs = ['0.10', '0.11', '0.12', '0.13', '0.14', '0.15', '0.16', '0.17', '0.18', '0.19', '0.20']
  cutoffs = ['0.13', '0.14', '0.15', '0.16']

  for cutoff in cutoffs:
    for i in range(0, len(properties['chrg_zg_' + str(cutoff).replace('.', '')[1:]])):
      if properties[ 'chrg_zg_' + str(cutoff).replace('.', '')[1:] ][i] > float(cutoff):  
        rho[str(cutoff)].append(properties['chrg_rho_' + str(cutoff).replace('.', '')[1:]][i])
        weights[str(cutoff)].append(prescales[i])

  
  logged_rho = defaultdict(list)
  for r in rho:
    logged_rho[r] = np.log(rho[r])
    
    plt.hist(logged_rho[r], weights=weights[r], label="Charged z\_cut=" + str(r), normed=1, histtype='step', lw=5, bins=50)

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

  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  
  # plt.savefig("plots/" + get_version(input_analysis_file) + "/rho/rho_" + str(pT_lower_cut) + ".pdf")
  plt.savefig("plots/" + "Version 4" + "/charged_rho/rho_" + str(pT_lower_cut) + ".pdf")
  plt.clf()


def plot_2d_rho_no_of_constituents(pT_lower_cut=100):

  properties = parse_file(input_analysis_file, pT_lower_cut)

  rhos = properties['rho_13']
  no_of_consts = properties['no_of_const']
  prescales = properties['prescales']

  logged_rho = np.log(rhos)

  rhos = logged_rho

  H, xedges, yedges = np.histogram2d(no_of_consts, rhos, bins=25, weights=prescales, normed=1, range=[[1, 120], [min(rhos), max(rhos)]] )

  H_normalized = []
  for i in range(0, 25):
    current_row = []
    factor = sum(H[i])
    for j in range(0, 25):
      current_row.append(H[i][j] / factor)
    H_normalized.append(current_row)

  H_normalized = np.array(H_normalized)
  H = H_normalized

  H = np.rot90(H)
  H = np.flipud(H)

  for i in range(0, len(H)):
    for j in range(0, len(H[j])):
      if str(H[i][j]) == "nan":
        H[i][j] = 0.

  Hmasked = np.ma.masked_where(H == 0, H) # Mask pixels with a value of zero

  plt.pcolormesh(xedges,yedges, Hmasked)

  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Counts')

  plt.xlabel('No. of Constituents')
  plt.ylabel('$\\rho$', rotation=0, labelpad=30)

  
  plt.gcf().set_size_inches(30, 30, forward=1)
  plt.gcf().set_snap(True)

  
  plt.savefig("plots/" + get_version(input_analysis_file) + "/rho_against_no_of_constituents/pT_lower_cut_" + str(pT_lower_cut) + ".pdf")

  plt.clf()



def plot_rho_jet_quality_parameters(pT_lower_cut=150):
  
  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut)

  parameters = ['no_of_const', 'chrg_multip', 'neu_had_frac', 'neu_em_frac', 'chrg_had_frac', 'chrg_em_frac']
  parameter_label = ['No. of Constituents', 'Charged Multiplicity', 'Neutral Hadron Fraction', 'Neutral EM Fraction', 'Charged Hadron Fraction', 'Charged EM Fraction']

  for index in range(0, len(parameters)):
    parameter = parameters[index]

    quality_parameter = properties[parameter]
    
    rhos = properties['rho_13']
    prescales = properties['prescales']  

    logged_rho = np.log(rhos)

    H, xedges, yedges = np.histogram2d( quality_parameter, logged_rho, bins=25, weights=prescales, normed=1, range=[ [min(quality_parameter), max(quality_parameter)], [ min(logged_rho), max(logged_rho) ] ] )

    H_normalized = []
    for i in range(0, 25):
      current_row = []
      factor = sum(H[i])
      for j in range(0, 25):
        current_row.append(H[i][j] / factor)

      H_normalized.append(current_row)

    H_normalized = np.array(H_normalized)
    H = H_normalized
    H = np.rot90(H)
    H = np.flipud(H)
    
    for a in range(0, len(H)):
      for b in range(0, len(H[j])):
        if str(H[a][b]) == "nan":
          H[a][b] = 0.


    Hmasked = np.ma.masked_where(H == 0, H) # Mask pixels with a value of zero
    plt.pcolormesh(xedges,yedges, Hmasked)

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')

    plt.xlabel(parameter_label[index], fontsize=75)
    plt.ylabel('$ \\rho $', rotation=0, fontsize=75, labelpad=30)

    plt.gcf().set_size_inches(30, 30, forward=1)
    plt.gcf().set_snap(True)


    plt.savefig("plots/" + get_version(input_analysis_file) + "/rho_against_jet_quality_parameters/pT_" + str(pT_lower_cut) + "/rho_13_vs_" + str(parameter) + ".pdf")

    plt.clf()



plot_rho_jet_quality_parameters(100)
plot_rho_jet_quality_parameters(200)
plot_rho_jet_quality_parameters(300)
plot_rho_jet_quality_parameters(400)
plot_rho_jet_quality_parameters(500)
plot_rho_jet_quality_parameters(600)


# plot_rho(pT_lower_cut=100)
# plot_rho(pT_lower_cut=150)
# plot_rho(pT_lower_cut=200)
# plot_rho(pT_lower_cut=250)
# plot_rho(pT_lower_cut=300)

# plot_charged_rho(pT_lower_cut=100)
# plot_charged_rho(pT_lower_cut=150)
# plot_charged_rho(pT_lower_cut=200)
# plot_charged_rho(pT_lower_cut=250)
# plot_charged_rho(pT_lower_cut=300)

# plot_2d_rho_no_of_constituents(100)
# plot_2d_rho_no_of_constituents(200)
# plot_2d_rho_no_of_constituents(300)
# plot_2d_rho_no_of_constituents(400)
# plot_2d_rho_no_of_constituents(500)
# plot_2d_rho_no_of_constituents(600)

call(["python", "/home/aashish/root/macros/MODAnalyzer/utilities/sync_plots.py"])