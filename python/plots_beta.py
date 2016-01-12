# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat
from __future__ import division

from subprocess import call

import time

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



def parse_mc(input_file, pT_lower_cut=0.00, pT_upper_cut = 20000.00, uncorrected_pT_lower_cut = 0.00, softdrop_unc_pT_lower_cut = 0.00, softdrop_cor_pT_lower_cut = 0.00, jet_quality_level=1):
  f = open(input_file, 'r')
  lines = f.read().split("\n")
  
  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        if (float(numbers[1]) > pT_lower_cut) and (float(numbers[1]) < pT_upper_cut):
         
          properties['hardest_pT'].append( float( numbers[1] ) )

          properties['Rg_10'].append( float( numbers[2] ) )
          properties['zg_10'].append( float( numbers[3] ) )
          properties['Rg_11'].append( float( numbers[4] ) )
          properties['zg_11'].append( float( numbers[5] ) )
          properties['Rg_12'].append( float( numbers[6] ) )
          properties['zg_12'].append( float( numbers[7] ) )
          properties['Rg_13'].append( float( numbers[8] ) )
          properties['zg_13'].append( float( numbers[9] ) )
          properties['Rg_14'].append( float( numbers[10] ) )
          properties['zg_14'].append( float( numbers[11] ) )
          properties['Rg_15'].append( float( numbers[12] ) )
          properties['zg_15'].append( float( numbers[13] ) )
          properties['Rg_16'].append( float( numbers[14] ) )
          properties['zg_16'].append( float( numbers[15] ) )
          properties['Rg_17'].append( float( numbers[16] ) )
          properties['zg_17'].append( float( numbers[17] ) )
          properties['Rg_18'].append( float( numbers[18] ) )
          properties['zg_18'].append( float( numbers[19] ) )
          properties['Rg_19'].append( float( numbers[20] ) )
          properties['zg_19'].append( float( numbers[21] ) )
          properties['Rg_20'].append( float( numbers[22] ) )
          properties['zg_20'].append( float( numbers[23] ) )
          
          properties['chrg_Rg_10'].append( float( numbers[24] ) )
          properties['chrg_zg_10'].append( float( numbers[25] ) )
          properties['chrg_Rg_11'].append( float( numbers[26] ) )
          properties['chrg_zg_11'].append( float( numbers[27] ) )
          properties['chrg_Rg_12'].append( float( numbers[28] ) )
          properties['chrg_zg_12'].append( float( numbers[29] ) )
          properties['chrg_Rg_13'].append( float( numbers[30] ) )
          properties['chrg_zg_13'].append( float( numbers[31] ) )
          properties['chrg_Rg_14'].append( float( numbers[32] ) )
          properties['chrg_zg_14'].append( float( numbers[33] ) )
          properties['chrg_Rg_15'].append( float( numbers[34] ) )
          properties['chrg_zg_15'].append( float( numbers[35] ) )
          properties['chrg_Rg_16'].append( float( numbers[36] ) )
          properties['chrg_zg_16'].append( float( numbers[37] ) )
          properties['chrg_Rg_17'].append( float( numbers[38] ) )
          properties['chrg_zg_17'].append( float( numbers[39] ) )
          properties['chrg_Rg_18'].append( float( numbers[40] ) )
          properties['chrg_zg_18'].append( float( numbers[41] ) )
          properties['chrg_Rg_19'].append( float( numbers[42] ) )
          properties['chrg_zg_19'].append( float( numbers[43] ) )
          properties['chrg_Rg_20'].append( float( numbers[44] ) )
          properties['chrg_zg_20'].append( float( numbers[45] ) )
          
    except:
      if len(numbers) != 0:
        # print "Some kind of error occured while parsing the given file!"
        # print numbers
        # print
        pass

  f.close()

  return properties

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

        properties['Rg_10'].append( float( numbers[6] ) )
        properties['zg_10'].append( float( numbers[7] ) )
        properties['Rg_11'].append( float( numbers[8] ) )
        properties['zg_11'].append( float( numbers[9] ) )
        properties['Rg_12'].append( float( numbers[10] ) )
        properties['zg_12'].append( float( numbers[11] ) )
        properties['Rg_13'].append( float( numbers[12] ) )
        properties['zg_13'].append( float( numbers[13] ) )
        properties['Rg_14'].append( float( numbers[14] ) )
        properties['zg_14'].append( float( numbers[15] ) )
        properties['Rg_15'].append( float( numbers[16] ) )
        properties['zg_15'].append( float( numbers[17] ) )
        properties['Rg_16'].append( float( numbers[18] ) )
        properties['zg_16'].append( float( numbers[19] ) )
        properties['Rg_17'].append( float( numbers[20] ) )
        properties['zg_17'].append( float( numbers[21] ) )
        properties['Rg_18'].append( float( numbers[22] ) )
        properties['zg_18'].append( float( numbers[23] ) )
        properties['Rg_19'].append( float( numbers[24] ) )
        properties['zg_19'].append( float( numbers[25] ) )
        properties['Rg_20'].append( float( numbers[26] ) )
        properties['zg_20'].append( float( numbers[27] ) )
        
        properties['chrg_Rg_10'].append( float( numbers[28] ) )
        properties['chrg_zg_10'].append( float( numbers[29] ) )
        properties['chrg_Rg_11'].append( float( numbers[30] ) )
        properties['chrg_zg_11'].append( float( numbers[31] ) )
        properties['chrg_Rg_12'].append( float( numbers[32] ) )
        properties['chrg_zg_12'].append( float( numbers[33] ) )
        properties['chrg_Rg_13'].append( float( numbers[34] ) )
        properties['chrg_zg_13'].append( float( numbers[35] ) )
        properties['chrg_Rg_14'].append( float( numbers[36] ) )
        properties['chrg_zg_14'].append( float( numbers[37] ) )
        properties['chrg_Rg_15'].append( float( numbers[38] ) )
        properties['chrg_zg_15'].append( float( numbers[39] ) )
        properties['chrg_Rg_16'].append( float( numbers[40] ) )
        properties['chrg_zg_16'].append( float( numbers[41] ) )
        properties['chrg_Rg_17'].append( float( numbers[42] ) )
        properties['chrg_zg_17'].append( float( numbers[43] ) )
        properties['chrg_Rg_18'].append( float( numbers[44] ) )
        properties['chrg_zg_18'].append( float( numbers[45] ) )
        properties['chrg_Rg_19'].append( float( numbers[46] ) )
        properties['chrg_zg_19'].append( float( numbers[47] ) )
        properties['chrg_Rg_20'].append( float( numbers[48] ) )
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
  pythia_properties = parse_mc("/home/aashish/pythia_qcd_truth.dat", pT_lower_cut=pT_lower_cut)
  herwig_properties = parse_mc("/home/aashish/herwig_qcd_truth.dat", pT_lower_cut=pT_lower_cut)
  sherpa_properties = parse_mc("/home/aashish/sherpa_qcd_truth.dat", pT_lower_cut=pT_lower_cut)

  # FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3
  jet_quality_level = 1

  prescales = properties['prescales']
  trig_jet_matched = properties['trig_jet_matched']
  jet_quality = properties['jet_quality']
  
  rho, weights = defaultdict(list), defaultdict(list)
  pythia_rho = defaultdict(list)

  # cutoffs = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
  # cutoffs = ['0.10', '0.11', '0.12', '0.13', '0.14', '0.15', '0.16', '0.17', '0.18', '0.19', '0.20']
  cutoffs = ['0.13', '0.14', '0.15', '0.16']

  for cutoff in cutoffs:
    for i in range(0, len(properties['zg_' + str(cutoff).replace('.', '')[1:]])):
      if properties[ 'zg_' + str(cutoff).replace('.', '')[1:] ][i] > float(cutoff):  
        rho[str(cutoff)].append(properties['rho_' + str(cutoff).replace('.', '')[1:]][i])
        weights[str(cutoff)].append(prescales[i])

  for cutoff in cutoffs:
    for i in range(0, len(pythia_properties['zg_' + str(cutoff).replace('.', '')[1:]])):
      # if pythia_properties[ 'zg_' + str(cutoff).replace('.', '')[1:] ][i] > float(cutoff):  
        # pythia_rho[str(cutoff)].append(pythia_properties['rho_' + str(cutoff).replace('.', '')[1:]][i])
        # pass

      pass


      

  
  logged_rho = defaultdict(list)
  pythia_logged_rho = defaultdict(list)

  for r in rho:
    logged_rho[r] = np.log(rho[r])
    pythia_logged_rho[r] = np.log(pythia_rho[r])
    
    plt.hist(logged_rho[r], weights=weights[r], label="z\_cut=" + str(r), normed=1, histtype='step', lw=5, bins=50)
    plt.hist(pythia_logged_rho[r], label="(Pythia) z\_cut=" + str(r), normed=1, histtype='step', lw=5, bins=50)

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






def plot_m_and_pT(pT_lower_cut=150):
  properties = parse_file(input_analysis_file, pT_lower_cut)

  m = properties['jet_13_m']
  pT = properties['jet_13_pT']
  prescales = properties['prescales']



  plt.hist(m, weights=prescales, histtype='step', bins=50, lw=5)
  plt.autoscale(True)
  plt.gca().set_xlabel("$m$", fontsize=75, labelpad=50)
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  plt.savefig("plots/" + "Version 4" + "/rho_experiments/m_pT_" + str(pT_lower_cut) + ".pdf")
  plt.clf()


  plt.hist(np.square(m), weights=prescales, histtype='step', bins=50, lw=5)
  plt.autoscale(True)
  plt.gca().set_xlabel("$m^2$", fontsize=75, labelpad=50)
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  plt.savefig("plots/" + "Version 4" + "/rho_experiments/m_2_pT_" + str(pT_lower_cut) + ".pdf")
  plt.clf()

  
  plt.hist(pT, weights=prescales, histtype='step', bins=50, lw=5)
  plt.autoscale(True)
  plt.gca().set_xlabel("$p_T$", fontsize=75, labelpad=50)
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  plt.savefig("plots/" + "Version 4" + "/rho_experiments/pT_" + str(pT_lower_cut) + ".pdf")
  plt.clf()


  plt.hist(np.square(pT), weights=prescales, histtype='step', bins=50, lw=5)
  plt.autoscale(True)
  plt.gca().set_xlabel("$p_T^2$", fontsize=75, labelpad=50)
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  plt.savefig("plots/" + "Version 4" + "/rho_experiments/pT_2_" + str(pT_lower_cut) + ".pdf")
  plt.clf()


  plt.hist(np.divide(np.square(m), np.square(pT)), weights=prescales, histtype='step', bins=50, lw=5, normed=1)
  plt.autoscale(True)
  plt.gca().set_xlabel("$\\frac{m^2}{p_T^2}$", fontsize=75, labelpad=50)
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  plt.savefig("plots/" + "Version 4" + "/rho_experiments/m_2_over_pT_2_" + str(pT_lower_cut) + ".pdf")
  plt.clf()



def plot_m_and_pT_decomposition(pT_lower_cut=150):
  properties = parse_file(input_analysis_file, pT_lower_cut)

  m = properties['jet_13_m']
  pT = properties['jet_13_pT']
  prescales = properties['prescales']


  def fix_nan(observable, prescales):
    x = []
    y = []
    for i in range(0, len(observable)):
      if not np.isnan(observable[i]):
        x.append(observable[i])
        y.append(prescales[i])

    return x, y

  # plt.hist(np.log(np.divide(np.square(m), np.square(pT))), weights=prescales, histtype='step', bins=50, lw=5, normed=1)
  
  a1 = 2 * np.log(m)
  a2 = 2 * np.log(pT)

  plt.hist(fix_nan(a1, prescales)[0], weights=fix_nan(a1, prescales)[1], label="$2 \\log{m}$", histtype='step', bins=50, lw=5, normed=1)
  plt.hist(fix_nan(a2, prescales)[0], weights=fix_nan(a2, prescales)[1], label="$2 \\log{p_T}$", histtype='step', bins=50, lw=5, normed=1)

  plt.hist(fix_nan(np.subtract(a1, a2), prescales)[0], weights=fix_nan(np.subtract(a1, a2), prescales)[1], label="$2 \\log{m} - 2 \\log{p_T}$", histtype='step', bins=50, lw=5, normed=1)


  plt.plot((0.00, 0.00), (0.00, 0.18), 'k--', lw=5)

  plt.legend(loc=2)
  plt.autoscale(True)
  # plt.gca().set_xlabel("$\\log{\\frac{m^2}{p_T^2}}$", fontsize=75, labelpad=50)
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  plt.savefig("plots/" + "Version 4" + "/rho_experiments/log_m_2_over_pT_2.pdf")
  plt.clf()



def plot_m_and_rho(pT_lower_cut=100):

  properties = parse_file(input_analysis_file, pT_lower_cut)

  rhos = properties['rho_13']
  mass = properties['jet_13_m']
  prescales = properties['prescales']



  H, xedges, yedges = np.histogram2d(mass, rhos, bins=25, weights=prescales, normed=1, range=[[min(mass), max(mass)], [min(rhos), max(rhos)]] )

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

  plt.xlabel('$m$')
  plt.ylabel('$\\rho$', rotation=0, labelpad=30)

  
  plt.gcf().set_size_inches(30, 30, forward=1)
  plt.gcf().set_snap(True)

  
  plt.savefig("plots/" + get_version(input_analysis_file) + "/rho/rho_vs_m.pdf")

  plt.clf()




def plot_pT_and_rho(pT_lower_cut=100):

  properties = parse_file(input_analysis_file, pT_lower_cut)

  rhos = properties['rho_13']
  pT = properties['jet_13_pT']
  prescales = properties['prescales']



  H, xedges, yedges = np.histogram2d(pT, rhos, bins=25, weights=prescales, normed=1, range=[[min(pT), max(pT)], [min(rhos), max(rhos)]] )

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

  plt.xlabel('$p_T$')
  plt.ylabel('$\\rho$', rotation=0, labelpad=30)

  
  plt.gcf().set_size_inches(30, 30, forward=1)
  plt.gcf().set_snap(True)

  
  plt.savefig("plots/" + get_version(input_analysis_file) + "/rho/rho_vs_pT.pdf")

  plt.clf()







def plot_theta_g_plots(pT_lower_cut=150, zg_cut='0.15', zg_filename='zg_15'):


  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut, jet_quality_level=1)
  
  pythia_properties = parse_mc("/home/aashish/pythia_truth_qcd.dat", pT_lower_cut=pT_lower_cut)
  herwig_properties = parse_mc("/home/aashish/herwig_truth_qcd.dat", pT_lower_cut=pT_lower_cut)
  sherpa_properties = parse_mc("/home/aashish/sherpa_truth_qcd.dat", pT_lower_cut=pT_lower_cut)

  # pythia_properties = parse_mc("/home/aashish/pythia_truth_qcd.dat", pT_lower_cut=pT_lower_cut)
  # herwig_properties = parse_mc("/home/aashish/herwig_truth_qcd.dat", pT_lower_cut=pT_lower_cut)
  # sherpa_properties = parse_mc("/home/aashish/sherpa_truth_qcd.dat", pT_lower_cut=pT_lower_cut)

  prescales = properties['prescales']

  z_g = properties[zg_filename]
  pythia_z_g = pythia_properties[zg_filename]
  herwig_z_g = herwig_properties[zg_filename]
  sherpa_z_g = sherpa_properties[zg_filename]

  R_g = properties[zg_filename.replace("zg", "Rg")]
  pythia_R_g = pythia_properties[zg_filename.replace("zg", "Rg")]
  herwig_R_g = herwig_properties[zg_filename.replace("zg", "Rg")]
  sherpa_R_g = sherpa_properties[zg_filename.replace("zg", "Rg")]


  R_g_with_cuts = []
  z_g_with_cuts_on_R_g = []
  prescales_for_R_g_with_cuts = []
  for i in range(0, len(R_g)):
    if float(R_g[i]) > 0.1 and float(R_g[i]) < 0.4:
      R_g_with_cuts.append(R_g[i])
      z_g_with_cuts_on_R_g.append(z_g[i])
      prescales_for_R_g_with_cuts.append(prescales[i])


  def get_mc_with_cuts(z_g, R_g):
    mc_Rg_with_cuts = []
    mc_zg_with_cuts = []
    
    for i in range(0, len(R_g)):
      if float(R_g[i]) > 0.1 and float(R_g[i]) < 0.4:
        mc_Rg_with_cuts.append(R_g[i])
        mc_zg_with_cuts.append(z_g[i])

    return (mc_zg_with_cuts, mc_Rg_with_cuts)

  pythia_R_g_with_cuts, pythia_z_g_with_cuts_on_R_g = get_mc_with_cuts(pythia_z_g, pythia_R_g)
  herwig_R_g_with_cuts, herwig_z_g_with_cuts_on_R_g = get_mc_with_cuts(herwig_z_g, herwig_R_g)
  sherpa_R_g_with_cuts, sherpa_z_g_with_cuts_on_R_g = get_mc_with_cuts(sherpa_z_g, sherpa_R_g)

  theta_g = np.divide(R_g, 0.5)
  theta_g_with_R_g_cuts = np.divide(R_g_with_cuts, 0.5)

  pythia_theta_g = np.divide(pythia_R_g, 0.5)
  pythia_theta_g_with_R_g_cuts = np.divide(pythia_R_g_with_cuts, 0.5)

  herwig_theta_g = np.divide(herwig_R_g, 0.5)
  herwig_theta_g_with_R_g_cuts = np.divide(herwig_R_g_with_cuts, 0.5)

  sherpa_theta_g = np.divide(sherpa_R_g, 0.5)
  sherpa_theta_g_with_R_g_cuts = np.divide(sherpa_R_g_with_cuts, 0.5)

  # ============================================================================================= z_g PLOT OF BEGINS ===========================================================================================================================

  bins_linear_log = np.linspace(math.log(float(zg_cut), math.e), math.log(0.5, math.e), 30)

  # Data.



  theta_g_hist = Hist(bins_linear_log, title='Data', markersize=3.0, color='black')
  bin_width = (theta_g_hist.upperbound() - theta_g_hist.lowerbound()) / theta_g_hist.nbins()
  map(theta_g_hist.Fill, np.log(z_g), prescales)
  theta_g_hist.Scale(1.0 / (theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # bins_linear_log = np.linspace(math.log(float(zg_cut), math.e), math.log(0.5, math.e), 30)
  # theta_g_with_cuts_hist = Hist(bins_linear_log, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (theta_g_with_cuts_hist.upperbound() - theta_g_with_cuts_hist.lowerbound()) / theta_g_with_cuts_hist.nbins()
  # map(theta_g_with_cuts_hist.Fill, np.log(z_g_with_cuts_on_R_g), prescales_for_R_g_with_cuts)
  # theta_g_with_cuts_hist.Scale(1.0 / (theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Data Ends.

  # Monte Carlo.

  # Pythia.

  pythia_theta_g_hist = Hist(bins_linear_log, title='Pythia', markersize=3.0, color='blue')
  bin_width = (pythia_theta_g_hist.upperbound() - pythia_theta_g_hist.lowerbound()) / pythia_theta_g_hist.nbins()
  map(pythia_theta_g_hist.Fill, np.log(pythia_z_g))
  pythia_theta_g_hist.Scale(1.0 / (pythia_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(pythia_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # bins_linear_log = np.linspace(math.log(float(zg_cut), math.e), math.log(0.5, math.e), 30)
  # pythia_theta_g_with_cuts_hist = Hist(bins_linear_log, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (pythia_theta_g_with_cuts_hist.upperbound() - pythia_theta_g_with_cuts_hist.lowerbound()) / pythia_theta_g_with_cuts_hist.nbins()
  # map(pythia_theta_g_with_cuts_hist.Fill, np.log(pythia_z_g_with_cuts_on_R_g))
  # pythia_theta_g_with_cuts_hist.Scale(1.0 / (pythia_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(pythia_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Pythia Ends.

  # Herwig.

  herwig_theta_g_hist = Hist(bins_linear_log, title='Herwig', markersize=3.0, color='green')
  bin_width = (herwig_theta_g_hist.upperbound() - herwig_theta_g_hist.lowerbound()) / herwig_theta_g_hist.nbins()
  map(herwig_theta_g_hist.Fill, np.log(herwig_z_g))
  herwig_theta_g_hist.Scale(1.0 / (herwig_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(herwig_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # bins_linear_log = np.linspace(math.log(float(zg_cut), math.e), math.log(0.5, math.e), 30)
  # herwig_theta_g_with_cuts_hist = Hist(bins_linear_log, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (herwig_theta_g_with_cuts_hist.upperbound() - herwig_theta_g_with_cuts_hist.lowerbound()) / herwig_theta_g_with_cuts_hist.nbins()
  # map(herwig_theta_g_with_cuts_hist.Fill, np.log(herwig_z_g_with_cuts_on_R_g))
  # herwig_theta_g_with_cuts_hist.Scale(1.0 / (herwig_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(herwig_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Herwig Ends.

  # Sherpa.

  sherpa_theta_g_hist = Hist(bins_linear_log, title='Sherpa', markersize=3.0, color='orange')
  bin_width = (sherpa_theta_g_hist.upperbound() - sherpa_theta_g_hist.lowerbound()) / sherpa_theta_g_hist.nbins()
  map(sherpa_theta_g_hist.Fill, np.log(sherpa_z_g))
  sherpa_theta_g_hist.Scale(1.0 / (sherpa_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(sherpa_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # bins_linear_log = np.linspace(math.log(float(zg_cut), math.e), math.log(0.5, math.e), 30)
  # sherpa_theta_g_with_cuts_hist = Hist(bins_linear_log, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (sherpa_theta_g_with_cuts_hist.upperbound() - sherpa_theta_g_with_cuts_hist.lowerbound()) / sherpa_theta_g_with_cuts_hist.nbins()
  # map(sherpa_theta_g_with_cuts_hist.Fill, np.log(sherpa_z_g_with_cuts_on_R_g))
  # sherpa_theta_g_with_cuts_hist.Scale(1.0 / (sherpa_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(sherpa_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Sherpa Ends.


  # Monte Carlo Ends.


  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[1.0, 1.0])
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  labels = ["$ \\textrm{Anti--}k_{t}\\textrm{:}~R = 0.5;\eta<2.4$", "$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  plt.gca().legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.50, 0.72])

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root-6.04.06/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.xlabel('$ z_g $', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=80.)
  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  x = np.linspace(math.log(float(zg_cut), math.e), math.log(0.5, math.e), 6)
  labels = [str(round(math.exp(i), 2)) for i in x]
  plt.xticks(x, labels)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.05))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.gca().autoscale(True)
  plt.gca().set_ylim(0., plt.gca().get_ylim()[1]*1.5)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/" + get_version(input_analysis_file) + "/theta_g/z_g/" + zg_filename + "_pT_lower_" + str(pT_lower_cut) + ".pdf")
  plt.clf()

  # =============================================================================================== z_g PLOT ENDS ===========================================================================================================================

  # ============================================================================================= theta_g PLOT BEGINS ===========================================================================================================================

  bins_linear_log = np.linspace(math.log(0.1, math.e), math.log(1.5, math.e), int((1.5 - 0.1) / 0.04))
  bins_linear_log_cuts = np.linspace(math.log(0.2, math.e), math.log(0.8, math.e), int((0.8 - 0.2) / 0.04))

  # Data.
  theta_g_hist = Hist(bins_linear_log, title='Data', markersize=3.0, color='black')
  bin_width = (theta_g_hist.upperbound() - theta_g_hist.lowerbound()) / theta_g_hist.nbins()
  map(theta_g_hist.Fill, np.log(theta_g), prescales)
  theta_g_hist.Scale(1.0 / (theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (theta_g_with_cuts_hist.upperbound() - theta_g_with_cuts_hist.lowerbound()) / theta_g_with_cuts_hist.nbins()
  # map(theta_g_with_cuts_hist.Fill, np.log(theta_g_with_R_g_cuts), prescales_for_R_g_with_cuts)
  # theta_g_with_cuts_hist.Scale(1.0 / (theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Data Ends.

  # Monte Carlo Begins.

  # Pythia.

  pythia_theta_g_hist = Hist(bins_linear_log, title='Pythia 8.212', markersize=3.0, color='blue')
  bin_width = (pythia_theta_g_hist.upperbound() - pythia_theta_g_hist.lowerbound()) / pythia_theta_g_hist.nbins()
  map(pythia_theta_g_hist.Fill, np.log(pythia_theta_g))
  pythia_theta_g_hist.Scale(1.0 / (pythia_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(pythia_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # pythia_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='purple')
  # bin_width = (pythia_theta_g_with_cuts_hist.upperbound() - pythia_theta_g_with_cuts_hist.lowerbound()) / pythia_theta_g_with_cuts_hist.nbins()
  # map(pythia_theta_g_with_cuts_hist.Fill, np.log(pythia_theta_g_with_R_g_cuts))
  # pythia_theta_g_with_cuts_hist.Scale(1.0 / (pythia_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(pythia_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Pythia Ends.

  # Herwig.

  herwig_theta_g_hist = Hist(bins_linear_log, title='Herwig++ 2.7.1', markersize=3.0, color='green')
  bin_width = (herwig_theta_g_hist.upperbound() - herwig_theta_g_hist.lowerbound()) / herwig_theta_g_hist.nbins()
  map(herwig_theta_g_hist.Fill, np.log(herwig_theta_g))
  herwig_theta_g_hist.Scale(1.0 / (herwig_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(herwig_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # herwig_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='purple')
  # bin_width = (herwig_theta_g_with_cuts_hist.upperbound() - herwig_theta_g_with_cuts_hist.lowerbound()) / herwig_theta_g_with_cuts_hist.nbins()
  # map(herwig_theta_g_with_cuts_hist.Fill, np.log(herwig_theta_g_with_R_g_cuts))
  # herwig_theta_g_with_cuts_hist.Scale(1.0 / (herwig_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(herwig_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Herwig Ends.

  # Sherpa.

  sherpa_theta_g_hist = Hist(bins_linear_log, title='Sherpa 2.2.0', markersize=3.0, color='orange')
  bin_width = (sherpa_theta_g_hist.upperbound() - sherpa_theta_g_hist.lowerbound()) / sherpa_theta_g_hist.nbins()
  map(sherpa_theta_g_hist.Fill, np.log(sherpa_theta_g))
  sherpa_theta_g_hist.Scale(1.0 / (sherpa_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(sherpa_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # sherpa_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='purple')
  # bin_width = (sherpa_theta_g_with_cuts_hist.upperbound() - sherpa_theta_g_with_cuts_hist.lowerbound()) / sherpa_theta_g_with_cuts_hist.nbins()
  # map(sherpa_theta_g_with_cuts_hist.Fill, np.log(sherpa_theta_g_with_R_g_cuts))
  # sherpa_theta_g_with_cuts_hist.Scale(1.0 / (sherpa_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(sherpa_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Sherpa Ends.
  

  # Monte Carlo Ends.



  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[1.0, 1.0])
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  labels = ["$ \\textrm{Anti--}k_{t}\\textrm{:}~R = 0.5;\eta<2.4$", "$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  plt.gca().legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.50, 0.72])

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root-6.04.06/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.xlabel('$ \\theta_g $', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=80.)
  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  x = np.linspace(math.log(0.1, math.e), math.log(1.5, math.e), 6)
  labels = [str(round(math.exp(i), 2)) for i in x]
  plt.xticks(x, labels)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.05))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.gca().autoscale(True)
  plt.gca().set_ylim(0., plt.gca().get_ylim()[1]*1.5)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/" + get_version(input_analysis_file) + "/theta_g/theta_g/" + zg_filename + "_pT_lower_" + str(pT_lower_cut) + ".pdf")
  plt.clf()

  # =============================================================================================== theta_g PLOT ENDS ===========================================================================================================================


  # ============================================================================================= theta_g * z_g PLOT BEGINS ===========================================================================================================================

  bins_linear_log = np.linspace(math.log(0.1 * float(zg_cut), math.e), math.log(0.6, math.e), int((0.6 - 0.1*float(zg_cut)) / 0.02))
  bins_linear_log_cuts = np.linspace(math.log(0.2*0.1, math.e), math.log(0.8*0.5, math.e), int((0.8*0.5 - 0.2*0.1) / 0.02))


  # Data Begins.

  theta_g_hist = Hist(bins_linear_log, title='Data', markersize=3.0, color='black')
  bin_width = (theta_g_hist.upperbound() - theta_g_hist.lowerbound()) / theta_g_hist.nbins()
  map(theta_g_hist.Fill, np.log(np.multiply(theta_g, z_g)), prescales)
  theta_g_hist.Scale(1.0 / (theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (theta_g_with_cuts_hist.upperbound() - theta_g_with_cuts_hist.lowerbound()) / theta_g_with_cuts_hist.nbins()
  # map(theta_g_with_cuts_hist.Fill, np.log(np.multiply(theta_g_with_R_g_cuts, z_g_with_cuts_on_R_g)), prescales_for_R_g_with_cuts)
  # theta_g_with_cuts_hist.Scale(1.0 / (theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Data Ends.

  # Monte Carlo.

  # Pythia.
  
  pythia_theta_g_hist = Hist(bins_linear_log, title='Pythia', markersize=3.0, color='blue')
  bin_width = (pythia_theta_g_hist.upperbound() - pythia_theta_g_hist.lowerbound()) / pythia_theta_g_hist.nbins()
  map(pythia_theta_g_hist.Fill, np.log(np.multiply(pythia_theta_g, pythia_z_g)))
  pythia_theta_g_hist.Scale(1.0 / (pythia_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(pythia_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # pythia_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (pythia_theta_g_with_cuts_hist.upperbound() - pythia_theta_g_with_cuts_hist.lowerbound()) / pythia_theta_g_with_cuts_hist.nbins()
  # map(pythia_theta_g_with_cuts_hist.Fill, np.log(np.multiply(pythia_theta_g_with_R_g_cuts, z_g_with_cuts_on_R_g)))
  # pythia_theta_g_with_cuts_hist.Scale(1.0 / (pythia_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(pythia_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Pythia Ends.

  # Herwig.

  herwig_theta_g_hist = Hist(bins_linear_log, title='Herwig', markersize=3.0, color='green')
  bin_width = (herwig_theta_g_hist.upperbound() - herwig_theta_g_hist.lowerbound()) / herwig_theta_g_hist.nbins()
  map(herwig_theta_g_hist.Fill, np.log(np.multiply(herwig_theta_g, herwig_z_g)))
  herwig_theta_g_hist.Scale(1.0 / (herwig_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(herwig_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # herwig_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (herwig_theta_g_with_cuts_hist.upperbound() - herwig_theta_g_with_cuts_hist.lowerbound()) / herwig_theta_g_with_cuts_hist.nbins()
  # map(herwig_theta_g_with_cuts_hist.Fill, np.log(np.multiply(herwig_theta_g_with_R_g_cuts, z_g_with_cuts_on_R_g)))
  # herwig_theta_g_with_cuts_hist.Scale(1.0 / (herwig_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(herwig_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Herwig Ends.

  # Sherpa.

  sherpa_theta_g_hist = Hist(bins_linear_log, title='Sherpa', markersize=3.0, color='orange')
  bin_width = (sherpa_theta_g_hist.upperbound() - sherpa_theta_g_hist.lowerbound()) / sherpa_theta_g_hist.nbins()
  map(sherpa_theta_g_hist.Fill, np.log(np.multiply(sherpa_theta_g, sherpa_z_g)))
  sherpa_theta_g_hist.Scale(1.0 / (sherpa_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(sherpa_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # sherpa_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (sherpa_theta_g_with_cuts_hist.upperbound() - sherpa_theta_g_with_cuts_hist.lowerbound()) / sherpa_theta_g_with_cuts_hist.nbins()
  # map(sherpa_theta_g_with_cuts_hist.Fill, np.log(np.multiply(sherpa_theta_g_with_R_g_cuts, z_g_with_cuts_on_R_g)))
  # sherpa_theta_g_with_cuts_hist.Scale(1.0 / (sherpa_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(sherpa_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Sherpa Ends.

  # Monte Carlo Ends. 


  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[1.0, 1.0])
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  labels = ["$ \\textrm{Anti--}k_{t}\\textrm{:}~R = 0.5;\eta<2.4$", "$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  plt.gca().legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.50, 0.72])

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root-6.04.06/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.xlabel('$ z_g \\theta_g $', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=80.)
  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  x = np.linspace(math.log(0.1 * float(zg_cut), math.e), math.log(0.6, math.e), 6)
  labels = [str(round(math.exp(i), 2)) for i in x]
  plt.xticks(x, labels)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.05))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.gca().autoscale(True)
  plt.gca().set_ylim(0., plt.gca().get_ylim()[1]*1.5)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/" + get_version(input_analysis_file) + "/theta_g/theta_g_times_zg/" + zg_filename + "_pT_lower_" + str(pT_lower_cut) + ".pdf")
  plt.clf()

  # =============================================================================================== theta_g * z_g PLOT ENDS ===========================================================================================================================

  # ============================================================================================= z_g * theta_g^2 PLOT BEGINS ===========================================================================================================================

  bins_linear_log = np.linspace(math.log(0.1*0.1*float(zg_cut), math.e), math.log(0.6*1*1, math.e), int( (0.5*0.8*0.8 - 0.2*0.2*float(zg_cut)) / 0.01))
  bins_linear_log_cuts = np.linspace(math.log(float(zg_cut)*0.2*0.2, math.e), math.log(0.5*0.8*0.8, math.e), int( (0.5*0.8*0.8 - 0.2*0.2*float(zg_cut)) / 0.01))

  # Data Begins.

  theta_g_hist = Hist(bins_linear_log, title='Data', markersize=3.0, color='black')
  bin_width = (theta_g_hist.upperbound() - theta_g_hist.lowerbound()) / theta_g_hist.nbins()
  map(theta_g_hist.Fill, np.log( np.multiply( z_g, np.square(theta_g) )), prescales)
  theta_g_hist.Scale(1.0 / (theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (theta_g_with_cuts_hist.upperbound() - theta_g_with_cuts_hist.lowerbound()) / theta_g_with_cuts_hist.nbins()
  # map(theta_g_with_cuts_hist.Fill, np.log( np.multiply( z_g_with_cuts_on_R_g, np.square(theta_g_with_R_g_cuts) )), prescales_for_R_g_with_cuts)
  # theta_g_with_cuts_hist.Scale(1.0 / (theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Data Ends.

  # Monte Carlo Begins.

  # Pythia Begins.

  pythia_theta_g_hist = Hist(bins_linear_log, title='Pythia', markersize=3.0, color='blue')
  bin_width = (pythia_theta_g_hist.upperbound() - pythia_theta_g_hist.lowerbound()) / pythia_theta_g_hist.nbins()
  map(pythia_theta_g_hist.Fill, np.log( np.multiply( pythia_z_g, np.square(pythia_theta_g) )))
  pythia_theta_g_hist.Scale(1.0 / (pythia_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(pythia_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # pythia_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (pythia_theta_g_with_cuts_hist.upperbound() - pythia_theta_g_with_cuts_hist.lowerbound()) / pythia_theta_g_with_cuts_hist.nbins()
  # map(pythia_theta_g_with_cuts_hist.Fill, np.log( np.multiply( pythia_z_g_with_cuts_on_R_g, np.square(pythia_theta_g_with_R_g_cuts) )))
  # pythia_theta_g_with_cuts_hist.Scale(1.0 / (pythia_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(pythia_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Pythia Ends.

  # Herwig Begins.

  herwig_theta_g_hist = Hist(bins_linear_log, title='Herwig', markersize=3.0, color='green')
  bin_width = (herwig_theta_g_hist.upperbound() - herwig_theta_g_hist.lowerbound()) / herwig_theta_g_hist.nbins()
  map(herwig_theta_g_hist.Fill, np.log( np.multiply( herwig_z_g, np.square(herwig_theta_g) )))
  herwig_theta_g_hist.Scale(1.0 / (herwig_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(herwig_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # herwig_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (herwig_theta_g_with_cuts_hist.upperbound() - herwig_theta_g_with_cuts_hist.lowerbound()) / herwig_theta_g_with_cuts_hist.nbins()
  # map(herwig_theta_g_with_cuts_hist.Fill, np.log( np.multiply( herwig_z_g_with_cuts_on_R_g, np.square(herwig_theta_g_with_R_g_cuts) )))
  # herwig_theta_g_with_cuts_hist.Scale(1.0 / (herwig_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(herwig_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Herwig Ends.

  # Sherpa Begins.

  sherpa_theta_g_hist = Hist(bins_linear_log, title='Sherpa', markersize=3.0, color='orange')
  bin_width = (sherpa_theta_g_hist.upperbound() - sherpa_theta_g_hist.lowerbound()) / sherpa_theta_g_hist.nbins()
  map(sherpa_theta_g_hist.Fill, np.log( np.multiply( sherpa_z_g, np.square(sherpa_theta_g) )))
  sherpa_theta_g_hist.Scale(1.0 / (sherpa_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(sherpa_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # sherpa_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (sherpa_theta_g_with_cuts_hist.upperbound() - sherpa_theta_g_with_cuts_hist.lowerbound()) / sherpa_theta_g_with_cuts_hist.nbins()
  # map(sherpa_theta_g_with_cuts_hist.Fill, np.log( np.multiply( sherpa_z_g_with_cuts_on_R_g, np.square(sherpa_theta_g_with_R_g_cuts) )))
  # sherpa_theta_g_with_cuts_hist.Scale(1.0 / (sherpa_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(sherpa_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Sherpa Ends.

  # Monte Carlo Ends.


  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[1.0, 1.0])
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  labels = ["$ \\textrm{Anti--}k_{t}\\textrm{:}~R = 0.5;\eta<2.4$", "$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  plt.gca().legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.50, 0.72])

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root-6.04.06/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.xlabel('$ z_g \\theta_g^2 $', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=80.)
  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  x = np.linspace(math.log(0.1*0.1*float(zg_cut), math.e), math.log(0.6*1*1, math.e), 6)
  labels = [str(round(math.exp(i), 3)) for i in x]
  plt.xticks(x, labels)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.05))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.gca().autoscale(True)
  plt.gca().set_ylim(0., plt.gca().get_ylim()[1]*1.5)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/" + get_version(input_analysis_file) + "/theta_g/z_g_times_theta_g^2/" + zg_filename + "_pT_lower_" + str(pT_lower_cut) + ".pdf")
  plt.clf()

  # =============================================================================================== z_g * theta_g^2 PLOT ENDS ===========================================================================================================================

  # ============================================================================================= z_g * theta_g^(0.5) PLOT BEGINS ===========================================================================================================================

  bins_linear_log = np.linspace(math.log(float(zg_cut)*math.sqrt(0.1), math.e), math.log(0.6*1*1, math.e), int( (0.6*1*1 - float(zg_cut)*math.sqrt(0.1)) / 0.02) )
  bins_linear_log_cuts = np.linspace(math.log(float(zg_cut)*math.sqrt(0.2), math.e), math.log(0.5*math.sqrt(0.8), math.e), int( (0.5*math.sqrt(0.8) - float(zg_cut)*math.sqrt(0.2)) / 0.02) )

  # Data.

  theta_g_hist = Hist(bins_linear_log, title='Data', markersize=3.0, color='black')
  bin_width = (theta_g_hist.upperbound() - theta_g_hist.lowerbound()) / theta_g_hist.nbins()
  map(theta_g_hist.Fill, np.log( np.multiply( z_g, np.sqrt(theta_g) )), prescales)
  theta_g_hist.Scale(1.0 / (theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (theta_g_with_cuts_hist.upperbound() - theta_g_with_cuts_hist.lowerbound()) / theta_g_with_cuts_hist.nbins()
  # map(theta_g_with_cuts_hist.Fill, np.log( np.multiply( z_g_with_cuts_on_R_g, np.sqrt(theta_g_with_R_g_cuts) )), prescales_for_R_g_with_cuts)
  # theta_g_with_cuts_hist.Scale(1.0 / (theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Data Ends.

  # Monte Carlo.

  # Pythia.

  pythia_theta_g_hist = Hist(bins_linear_log, title='Pythia', markersize=3.0, color='blue')
  bin_width = (pythia_theta_g_hist.upperbound() - pythia_theta_g_hist.lowerbound()) / pythia_theta_g_hist.nbins()
  map(pythia_theta_g_hist.Fill, np.log( np.multiply( pythia_z_g, np.sqrt(pythia_theta_g) )))
  pythia_theta_g_hist.Scale(1.0 / (pythia_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(pythia_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # pythia_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (pythia_theta_g_with_cuts_hist.upperbound() - pythia_theta_g_with_cuts_hist.lowerbound()) / pythia_theta_g_with_cuts_hist.nbins()
  # map(pythia_theta_g_with_cuts_hist.Fill, np.log( np.multiply( pythia_z_g_with_cuts_on_R_g, np.sqrt(pythia_theta_g_with_R_g_cuts) )))
  # pythia_theta_g_with_cuts_hist.Scale(1.0 / (pythia_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(pythia_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Pythia Ends.

  # Herwig.

  herwig_theta_g_hist = Hist(bins_linear_log, title='Herwig', markersize=3.0, color='green')
  bin_width = (herwig_theta_g_hist.upperbound() - herwig_theta_g_hist.lowerbound()) / herwig_theta_g_hist.nbins()
  map(herwig_theta_g_hist.Fill, np.log( np.multiply( herwig_z_g, np.sqrt(herwig_theta_g) )))
  herwig_theta_g_hist.Scale(1.0 / (herwig_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(herwig_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # herwig_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (herwig_theta_g_with_cuts_hist.upperbound() - herwig_theta_g_with_cuts_hist.lowerbound()) / herwig_theta_g_with_cuts_hist.nbins()
  # map(herwig_theta_g_with_cuts_hist.Fill, np.log( np.multiply( herwig_z_g_with_cuts_on_R_g, np.sqrt(herwig_theta_g_with_R_g_cuts) )))
  # herwig_theta_g_with_cuts_hist.Scale(1.0 / (herwig_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(herwig_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Herwig Ends.

  # Sherpa.

  sherpa_theta_g_hist = Hist(bins_linear_log, title='Sherpa', markersize=3.0, color='orange')
  bin_width = (sherpa_theta_g_hist.upperbound() - sherpa_theta_g_hist.lowerbound()) / sherpa_theta_g_hist.nbins()
  map(sherpa_theta_g_hist.Fill, np.log( np.multiply( sherpa_z_g, np.sqrt(sherpa_theta_g) )))
  sherpa_theta_g_hist.Scale(1.0 / (sherpa_theta_g_hist.GetSumOfWeights() * bin_width))
  rplt.errorbar(sherpa_theta_g_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # sherpa_theta_g_with_cuts_hist = Hist(bins_linear_log_cuts, title='Jets with $0.1 < R_g < 0.4$', markersize=3.0, color='red')
  # bin_width = (sherpa_theta_g_with_cuts_hist.upperbound() - sherpa_theta_g_with_cuts_hist.lowerbound()) / sherpa_theta_g_with_cuts_hist.nbins()
  # map(sherpa_theta_g_with_cuts_hist.Fill, np.log( np.multiply( sherpa_z_g_with_cuts_on_R_g, np.sqrt(sherpa_theta_g_with_R_g_cuts) )))
  # sherpa_theta_g_with_cuts_hist.Scale(1.0 / (sherpa_theta_g_with_cuts_hist.GetSumOfWeights() * bin_width))
  # rplt.errorbar(sherpa_theta_g_with_cuts_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Sherpa Ends.

  # Monte Carlo Ends.


  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[1.0, 1.0])
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  labels = ["$ \\textrm{Anti--}k_{t}\\textrm{:}~R = 0.5;\eta<2.4$", "$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$", "$ \\textrm{Soft~Drop:}~\\boldsymbol{\\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  plt.gca().legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.50, 0.72])

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root-6.04.06/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.xlabel('$ z_g \\theta_g^{1/2} $', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=80.)
  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  x = np.linspace(math.log(float(zg_cut)*math.sqrt(0.1), math.e), math.log(0.6*1*1, math.e), 6)
  labels = [str(round(math.exp(i), 3)) for i in x]
  plt.xticks(x, labels)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.05))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.gca().autoscale(True)
  plt.gca().set_ylim(0., plt.gca().get_ylim()[1]*1.5)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/" + get_version(input_analysis_file) + "/theta_g/z_g_times_sqrt_theta_g/" + zg_filename + "_pT_lower_" + str(pT_lower_cut) + ".pdf")
  plt.clf()

  # =============================================================================================== z_g * theta_g^(0.5) PLOT ENDS ===========================================================================================================================


plot_theta_g_plots(pT_lower_cut=150, zg_cut='0.10', zg_filename='zg_10')
plot_theta_g_plots(pT_lower_cut=300, zg_cut='0.10', zg_filename='zg_10')
# plot_theta_g_plots(pT_lower_cut=500, zg_cut='0.10', zg_filename='zg_10')
# plot_theta_g_plots(pT_lower_cut=600, zg_cut='0.10', zg_filename='zg_10')



plot_theta_g_plots(pT_lower_cut=150, zg_cut='0.15', zg_filename='zg_15')
plot_theta_g_plots(pT_lower_cut=300, zg_cut='0.15', zg_filename='zg_15')
# plot_theta_g_plots(pT_lower_cut=500, zg_cut='0.15', zg_filename='zg_15')
# plot_theta_g_plots(pT_lower_cut=600, zg_cut='0.15', zg_filename='zg_15')



plot_theta_g_plots(pT_lower_cut=150, zg_cut='0.20', zg_filename='zg_20')
plot_theta_g_plots(pT_lower_cut=300, zg_cut='0.20', zg_filename='zg_20')
# plot_theta_g_plots(pT_lower_cut=500, zg_cut='0.20', zg_filename='zg_20')
# plot_theta_g_plots(pT_lower_cut=600, zg_cut='0.20', zg_filename='zg_20')






# plot_m_and_rho()

# plot_pT_and_rho()

# plot_m_and_pT(100)
# plot_m_and_pT(200)
# plot_m_and_pT(300)
# plot_m_and_pT(400)
# plot_m_and_pT(500)


# plot_m_and_pT_decomposition(250)

# plot_rho_jet_quality_parameters(100)
# plot_rho_jet_quality_parameters(200)
# plot_rho_jet_quality_parameters(300)
# plot_rho_jet_quality_parameters(400)
# plot_rho_jet_quality_parameters(500)
# plot_rho_jet_quality_parameters(600)


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