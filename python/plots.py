# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat
from __future__ import division

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

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
        if (float(numbers[4]) > pT_lower_cut) and (float(numbers[17]) > pfc_pT_cut) and (float(numbers[4]) < pT_upper_cut):
          properties['event_number'].append( float( numbers[1] ) )
          properties['run_number'].append( float( numbers[2] ) )

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

          properties['pTs_after_SD'].append( float( numbers[36] ) )

          properties['multiplicity_before_SD'].append( float( numbers[37] ) )
          properties['multiplicity_after_SD'].append( float( numbers[38] ) )

          properties['JEC'].append( float( numbers[39] ) )
          
          properties['jet_area'].append( float( numbers[40] ) )

          properties['jet_mass_before_SD'].append( float( numbers[41] ) )
          properties['jet_mass_after_SD'].append( float( numbers[42] ) )

          properties['fractional_energy_loss'].append( float( numbers[43] ) )

    except:
      if len(numbers) != 0:
        # print "Some kind of error occured while parsing the given file!"
        # print numbers
        # print
        pass

  return properties



def normalize_hist(hist, variable_bin=False):

  if variable_bin:
    hist.Scale( 1.0 / hist.GetSumOfWeights() )

    for i in range(0, hist.GetSize()):
      bin_width = hist.GetXaxis().GetBinWidth(i)
      old_bin_height =  hist.GetBinContent(i)
      
      new_height = old_bin_height / bin_width
      hist.SetBinContent(i, new_height)

      old_error = hist.GetBinError(i)
      new_error = old_error / bin_width
      hist.SetBinError(i, new_error)
  else:
    bin_width = (hist.upperbound() - hist.lowerbound()) / hist.nbins()
    hist.Scale( 1.0 / ( hist.GetSumOfWeights() * bin_width ) )

  return hist


def log_bins(data, weights, number_of_bins=50):
  
  def drop_zeros(a_list):
    return [i for i in a_list if i > 0]

  min_value = math.log( min( drop_zeros(data) ), 10 )
  max_value = math.log( max(data), 10 )

  return np.histogram(data, weights=weights, bins=np.logspace(min_value, max_value, number_of_bins))


def lin_bins(data, weights, number_of_bins=50):
  
  def drop_zeros(a_list):
    return [i for i in a_list if i > 0]

  min_value = min( drop_zeros(data) )
  max_value = max(data)

  return np.histogram(data, weights=weights, bins=np.linspace(min_value, max_value, number_of_bins))


def parse_file_turn_on(input_file, pT_lower_cut = 0.00):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        if (float(numbers[1]) > pT_lower_cut):
          properties['corrected_hardest_pts'].append( float( numbers[3] ) )
          properties['prescales'].append( int( numbers[4] ) )
          properties['trigger_names'].append(  numbers[5] )

    except:
      if len(numbers) != 0:
        print "Some kind of error occured while parsing the given file!"
        print numbers
        print


  return properties



def parse_mc_pt_file(input_file, pT_lower_cut=100., pT_upper_cut=20000.):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  pTs = []

  for line in lines:
    try:
      numbers = line.split()
      if (float(numbers[0]) > pT_lower_cut and float(numbers[0]) < pT_upper_cut):
        pTs.append( float( numbers[0] ) )
    except:
      if len(numbers) != 0:
        # print "Some kind of error occured while parsing the given file!"
        # print numbers
        # print
        pass
  
  return pTs

    
def parse_mc_file(input_file, pT_lower_cut = 0.00, pfc_pT_cut = 0.00, pT_upper_cut=20000):
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
      if len(numbers) != 0:
        # print "Some kind of error occured while parsing the given file!"
        # print numbers
        # print
        pass
  
  return properties




def plot_dr():
  pT_lower_cut = 150
  properties = parse_file(input_analysis_file, pT_lower_cut)

  drs = [properties['dr_05'], properties['dr_1'], properties['dr_2']]
  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['$z_{cut}$ = 0.05', '$z_{cut}$ = 0.1', '$z_{cut}$ = 0.2']

  for i in range(0, len(drs)):
    # no. of bins, xlower, xhigher
    dr_hist = Hist(50, 0.0, 0.5, title=labels[i], markersize=1.0, color=colors[i])

    map(dr_hist.Fill, drs[i], prescales)
    
    dr_hist.Scale(1.0 / dr_hist.GetSumOfWeights())
    
    rplt.errorbar(dr_hist, xerr=False, emptybins=False, markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.autoscale(True)
  # plt.xlim(0.0, 0.5)

  plt.legend()
  plt.xlabel("$\Delta$R between Subjets")
  plt.suptitle("$\Delta$R between Subjets with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV")

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

  plt.savefig("plots/delta_r_distribution_pt_cut_" + str(pT_lower_cut) + ".pdf")
  plt.show()

def plot_mu():
  pT_lower_cut = 150
  properties = parse_file(input_analysis_file, pT_lower_cut)

  mus = [properties['mu_05'], properties['mu_1'], properties['mu_2']]
  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['$z_{cut}$ = 0.05', '$z_{cut}$ = 0.1', '$z_{cut}$ = 0.2']

  for i in range(0, len(mus)):
    # no. of bins, xlower, xhigher
    mu_hist = Hist(50, 0, 1.0, title=labels[i], markersize=1.0, color=colors[i])

    map(mu_hist.Fill, mus[i], prescales)
    
    mu_hist.Scale(1.0 / mu_hist.GetSumOfWeights())
    
    rplt.errorbar(mu_hist, xerr=False, emptybins=False, markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.autoscale(True)
  # plt.xlim(0.00, 1)

  plt.legend()
  plt.xlabel("Mass Drop($\mu$)")
  plt.suptitle("Mass Drop($\mu$) with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV")

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

  plt.savefig("plots/mass_drop_distribution_pt_cut_" + str(pT_lower_cut) + ".pdf")
  plt.show()




def plot_charged_pt():
  properties = parse_file(input_analysis_file)
  pdgid_map = { 1: "d", 130: "$K^0_L$ Meson", 11: "$e^-$", -211: "$\pi^-$", 13: "$\mu^-$", 211: "$\pi^+$", -11: "$e^+$", 22: "$\gamma$", 2: "u", -13: "$\mu^+$" }

  pdgids = properties['hardest_pfc_pdgid']
  pTs = properties['hardest_pfc_pt']
  prescales = properties['prescales']

  pdgid_pts = []
  pdgid_prescales = []

  for i in range(0, len(pdgids)):
    if (abs(pdgids[i]) == 211) and (abs(pdgids[i]) == 11) and (abs(pdgids[i]) == 13):
      pdgid_pts.append(pTs[i])
      pdgid_prescales.append(prescales[i])


  pdgid_hist = Hist(25, 0, 800, markersize=1.0, color='blue')
  
  map(pdgid_hist.Fill, pTs, prescales)

  if pdgid_hist.GetSumOfWeights() != 0:
    pdgid_hist.Scale(1.0 / pdgid_hist.GetSumOfWeights())

  rplt.errorbar(pdgid_hist, xerr=False, emptybins=False, markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.yscale('log')

  plt.autoscale(True)

  plt.xlabel('$p_{T}$ GeV')
  plt.suptitle("Hardest Charged pT Distribution")

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

  plt.savefig("plots/hardest_charged_pt_distribution.pdf")
  plt.show()


def plot_zg_pfc_pt_cut(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', n_bins=10, y_max_limit=20):

  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut)

  zgs = properties[zg_filename]
  zg_pfc_cut_2 = properties[zg_filename + '_pt_2']
  zg_pfc_cut_5 = properties[zg_filename + '_pt_5']

  prescales = properties['prescales']

  zg_hist = Hist(6 * n_bins, 0.0, 0.6, title="All PF Candidates", markersize=1.0, color='black')
  zg_pfc_cut_2_hist = Hist(6 * n_bins, 0.0, 0.6, title="PFCs with $p_T >$ 2 GeV", markersize=1.0, color='red')
  zg_pfc_cut_5_hist = Hist(6 * n_bins, 0.0, 0.6, title="PFCs with $p_T >$ 5 GeV", markersize=1.0, color='blue')

  bin_width_zg = (zg_hist.upperbound() - zg_hist.lowerbound()) / zg_hist.nbins()
  bin_width_zg_pfc_cut_2 = (zg_pfc_cut_2_hist.upperbound() - zg_pfc_cut_2_hist.lowerbound()) / zg_pfc_cut_2_hist.nbins()
  bin_width_zg_pfc_cut_5 = (zg_pfc_cut_5_hist.upperbound() - zg_pfc_cut_5_hist.lowerbound()) / zg_pfc_cut_5_hist.nbins()


  map(zg_hist.Fill, zgs, prescales)
  map(zg_pfc_cut_2_hist.Fill, zg_pfc_cut_2, prescales)
  map(zg_pfc_cut_5_hist.Fill, zg_pfc_cut_5, prescales)

  if zg_hist.GetSumOfWeights() != 0:
    zg_hist.Scale(1.0 / (zg_hist.GetSumOfWeights() * bin_width_zg))

  if zg_pfc_cut_2_hist.GetSumOfWeights() != 0:
    zg_pfc_cut_2_hist.Scale(1.0 / (zg_pfc_cut_2_hist.GetSumOfWeights() * bin_width_zg_pfc_cut_2))

  if zg_pfc_cut_5_hist.GetSumOfWeights() != 0:
    zg_pfc_cut_5_hist.Scale(1.0 / (zg_pfc_cut_5_hist.GetSumOfWeights() * bin_width_zg_pfc_cut_5))
  
  rplt.errorbar(zg_hist, emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.errorbar(zg_pfc_cut_2_hist, emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.errorbar(zg_pfc_cut_5_hist, emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  
  plt.autoscale(True)
  plt.gca().set_ylim(0, y_max_limit)

  legend = plt.legend(frameon=0, fontsize=60, bbox_to_anchor=[1.0, 0.98])
  plt.gca().add_artist(legend)

  
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra, extra, extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.97, 0.60])


  plt.xlabel("$z_g$", fontsize=95)
  plt.ylabel("$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", fontsize=80, rotation=0, labelpad=115, y=0.39)
  

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gca().set_ylim(0, y_max_limit)
  plt.gca().set_xlim(0.0, 0.6)

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "softcut/" + str(zg_filename) + "_pt_cut_" + str(pT_lower_cut) + ".pdf"

  plt.savefig("plots/softcut/" + str(zg_filename) + "_pt_cut_" + str(pT_lower_cut) + ".pdf")
  # plt.show()
  plt.clf()



def plot_turn_on_curves():
  # properties = parse_file_turn_on(input_analysis_file)
  
  properties = parse_file_turn_on('./turn_on.dat')

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']

  expected_trigger_names = ["HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U", "HLT\_Jet15U\_HcalNoiseFiltered" ]
  labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U", "Jet15U-HNF" ]

  colors = ['orange', 'red', 'green', 'blue', 'magenta', 'black']
  colors = colors[::-1]

  pt_hists = []
  for i in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist(60, 0, 300, title=labels[i], markersize=1.0, color=colors[i], linewidth=5))

  for i in range(0, len(pTs)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j].replace("\\", "") in trigger_names[i]:
        pt_hists[j].Fill(pTs[i], prescales[i])

  for k in range(len(pt_hists) - 1, -1, -1):
    rplt.errorbar(pt_hists[k], emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.autoscale(True)
  plt.yscale('log')

  plt.gca().set_ylim(1, 10e10)

  
  # Legend.
  handles, labels = plt.gca().get_legend_handles_labels()
  first_legend = plt.gca().legend(handles, labels, fontsize=60,  frameon=0, borderpad=0.1, bbox_to_anchor=[0.97, 0.98])
  ax = plt.gca().add_artist(first_legend)

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5} $"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.935, 0.65])


  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=95)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)

  fn = get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.9249985), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 30, forward=1)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)


  plt.savefig("plots/turn_on_curves.pdf")
  # plt.show()
  plt.clf()


def plot_hardest_pt_corresponding_triggers():
  properties = parse_file(input_analysis_file)

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']

  expected_trigger_names = ["HLT_Jet100U", "HLT_Jet70U", "HLT_Jet50U", "HLT_Jet30U", "HLT_Jet15U"]
  labels = ["Jet100U", "Jet70U", "Jet50U", "Jet30U", "Jet15U"]

  colors = ['brown', 'red', 'blue', 'magenta', 'green']

  pt_hists = []
  for i in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist(50, 0, 1000, title=labels[i], markersize=1.0, color=colors[i], linewidth=5))

  for i in range(0, len(pTs)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j] in trigger_names[i]:
        pt_hists[j].Fill(pTs[i], prescales[i])

  for k in range(0, len(pt_hists)):
    rplt.errorbar(pt_hists[k], marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.yscale('log')
  plt.autoscale(True)

  plt.gca().set_ylim(1, 10e8)

  plt.legend(loc=0, frameon=0)

  plt.xlabel('$p_T \mathrm{(GeV)}$')

  fn = get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)
  ab = AnnotationBbox(OffsetImage(read_png(fn), zoom=0.15, resample=1, dpi_cor=1), (125, 30e7), frameon=0)
  plt.gca().add_artist(ab)

  preliminary_text = "Prelim. \n(20\%)"
  plt.gcf().text(0.31, 0.83, preliminary_text, fontsize=40, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.savefig("plots/hardest_pt_corresponding_triggers.pdf")
  # plt.show()
  plt.clf()



def plot_2d_hist():

  pT_lower_cut = 150
  properties = parse_file(input_analysis_file, pT_lower_cut)
  zgs = [properties['zg_05'], properties['zg_1'], properties['zg_2']]
  charged_zgs = [properties['zg_charged_05'], properties['zg_charged_1'], properties['zg_charged_2']]
  prescales = properties['prescales']

  
  H, xedges, yedges = np.histogram2d(charged_zgs[0], zgs[0], normed=1, range=[[0, 0.5], [0, 0.5]], weights=prescales, bins=25)
  Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
  plt.pcolormesh(xedges,yedges,Hmasked)
  plt.xlabel('Charged zg\_05')
  plt.ylabel('zg\_05')
  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Counts')  
  plt.gcf().set_size_inches(30, 30, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  plt.savefig("plots/zg_05_vs_charged_zg_05.pdf")
  plt.clf()


  H, xedges, yedges = np.histogram2d(charged_zgs[1], zgs[1], normed=1, range=[[0, 0.5], [0, 0.5]], weights=prescales, bins=25)
  Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
  plt.pcolormesh(xedges,yedges,Hmasked)
  plt.xlabel('Charged zg\_1')
  plt.ylabel('zg\_1')
  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Counts')
  plt.gcf().set_size_inches(30, 30, forward=1)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
  plt.savefig("plots/zg_1_vs_charged_zg_1.pdf")
  plt.clf()


  H, xedges, yedges = np.histogram2d(charged_zgs[2], zgs[2], normed=1, range=[[0, 0.5], [0, 0.5]], weights=prescales, bins=25)
  Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
  plt.pcolormesh(xedges,yedges,Hmasked)
  plt.xlabel('Charged zg\_2')
  plt.ylabel('zg\_2')
  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Counts')
  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)
  plt.savefig("plots/zg_2_vs_charged_zg_2.pdf")
  plt.clf()

 





def plot_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', ratio_denominator="theory", data=True, mc=True, theory=True, n_bins=10, y_max_limit=20, y_limit_ratio_plot=0.5):

  zg_cut = float(zg_cut)

  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut)
  properties_pythia = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_pythia_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)
  properties_herwig = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_herwig_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)

  zg_data = properties[zg_filename]
  
  zg_pythias = properties_pythia[zg_filename]
  zg_herwigs = properties_herwig[zg_filename]

  prescales = properties['prescales']

  data_label = "CMS 2010 Open Data"
  pythia_label = "Pythia 8.205" if mc else ""
  herwig_label = "Herwig++ 2.6.3" if mc else ""
  theory_label = "Theory (MLL)" if theory else ""

  
  gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 

 
  ax0 = plt.subplot(gs[0])
  ax1 = plt.subplot(gs[1])


  # Theory Plots Begin.
  
  points_th_gluon = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Simone/results_7_24_15/band_gluon_pt" + str(pT_lower_cut) + "_zc" + str(zg_cut).replace(".", "") + ".dat")
  points_th_quark = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Simone/results_7_24_15/band_quark_pt" + str(pT_lower_cut) + "_zc" + str(zg_cut).replace(".", "") + ".dat")

  points = defaultdict(list)

  for x in points_th_gluon:
    points[x] = [ points_th_gluon[x][0], points_th_gluon[x][1], points_th_gluon[x][2], points_th_gluon[x][3], points_th_gluon[x][4], points_th_gluon[x][5] ]
    points[x].extend([ points_th_quark[x][0], points_th_quark[x][1], points_th_quark[x][2], points_th_quark[x][3], points_th_quark[x][4], points_th_quark[x][5] ])

  keys = points.keys()
  keys.sort()

  theory_x = keys

  y = []
  for j in range(0, 6):
    y.append([points[i][j] for i in keys])

  # For each x, record three y's viz. max_y, min_y, line_y (i.e. e11 xmu=1).

  theory_y_max = []
  theory_y_min = []
  theory_y_line = []
  for i in range(0, len(theory_x)):
    y_for_current_x = []
    for j in range(0, 6):
      y_for_current_x.append(y[j][i])

    theory_y_min.append(min(y_for_current_x))
    theory_y_line.append(y_for_current_x[1])
    theory_y_max.append(max(y_for_current_x))
    
    
  if theory:
    area_theory_y_max = simps(theory_y_max, theory_x)
    # weighted_theory_y_max = map(lambda x: x / area_theory_y_max, theory_y_max)
    weighted_theory_y_max = theory_y_max
    ax0.plot(theory_x, weighted_theory_y_max, alpha=0.0, color='red')
    
    area_theory_y_line = simps(theory_y_line, theory_x)
    # weighted_theory_y_line = map(lambda x: x / area_theory_y_line, theory_y_line)
    weighted_theory_y_line = theory_y_line
    ax0.plot(theory_x, weighted_theory_y_line, label=theory_label, alpha=1.0, color='red', linewidth=5)

    area_theory_y_min = simps(theory_y_min, theory_x)
    # weighted_theory_y_min = map(lambda x: x / area_theory_y_min, theory_y_min)
    weighted_theory_y_min = theory_y_min
    ax0.plot(theory_x, weighted_theory_y_min, alpha=0.0, color='red')


    ax0.fill_between(theory_x, theory_y_max, theory_y_min, norm=1, where=np.less_equal(theory_y_min, theory_y_max), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)


  # Theory Plot Ends.

  def convert_hist_to_line_plot(hist, n_bins):
    a = []
    b = {}
    bin_width = 0.6 / (6 * n_bins)
    for i in range(0, len(list(hist.x()))):
      a.append(round(list(hist.x())[i] - bin_width / 2., 4))
      a.append(round(list(hist.x())[i], 4))
      a.append(round(list(hist.x())[i] + bin_width / 2., 4))

      if round(list(hist.x())[i] - bin_width / 2., 4) not in b.keys():
        b[round(list(hist.x())[i] - bin_width / 2., 4)] = [ list(hist.y())[i] ]
      else:
        b[round(list(hist.x())[i] - bin_width / 2., 4)].append( list(hist.y())[i] )
      
      if round(list(hist.x())[i], 4) not in b.keys():
        b[round(list(hist.x())[i], 4)] = [ list(hist.y())[i] ]
      else:
        b[round(list(hist.x())[i], 4)].append( list(hist.y())[i] )

      if round(list(hist.x())[i] + bin_width / 2., 4) not in b.keys():
        b[round(list(hist.x())[i] + bin_width / 2., 4)] = [ list(hist.y())[i] ]
      else:
        b[round(list(hist.x())[i] + bin_width / 2., 4)].append( list(hist.y())[i] )

    x = sorted(list(Set(a)))
    a.sort()
    
    c = [b[x[i]] for i in range(0, len(x))]

    y = [item for sublist in c for item in sublist]
    
    a_zero_removed = []
    y_zero_removed = []
    for i in range(0, len(a)):
      if a[i] >= zg_cut and a[i] <= 0.5 and y[i] != 0.0:
        a_zero_removed.append(a[i])
        y_zero_removed.append(y[i])

    return a_zero_removed, y_zero_removed

  
  
  # Data Plot Begins.
  
  zg_data_hist = Hist(6 * n_bins, 0.0, 0.6, title=data_label, markersize=2.5, color='black')
  bin_width_data = (zg_data_hist.upperbound() - zg_data_hist.lowerbound()) / zg_data_hist.nbins()

  map(zg_data_hist.Fill, zg_data, prescales)
  
  zg_data_hist.Scale(1.0 / ( zg_data_hist.GetSumOfWeights() * bin_width_data ))
  
  if data:
    # data_plot, caplines, barlinecols
    data_plot = rplt.errorbar(zg_data_hist, xerr=1, yerr=1, emptybins=False, axes=ax0, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
  else:
    data_plot = rplt.errorbar(zg_data_hist, xerr=1, yerr=1, emptybins=False, axes=ax0, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=0.0)


  data_x_errors, data_y_errors = [], []
  for x_segment in data_plot[2][0].get_segments():
    data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
  for y_segment in data_plot[2][1].get_segments():
    data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

  data_points_x = data_plot[0].get_xdata()
  data_points_y = data_plot[0].get_ydata()

  # print sorted(data_points_x)



  # Data Plots Ends.

  
  # Simulation Plots Begin. 
  
  # Pythia.
  
  zg_pythia_hist = Hist(6 * n_bins, 0, 0.6, title=pythia_label, markersize=5.0, color='blue', linewidth=5)
  bin_width_pythia = (zg_pythia_hist.upperbound() - zg_pythia_hist.lowerbound()) / zg_pythia_hist.nbins()

  map(zg_pythia_hist.Fill, zg_pythias)

  zg_pythia_hist.Scale(1.0 / ( zg_pythia_hist.GetSumOfWeights() * bin_width_pythia ))

  if mc:
    # pythia_plot = rplt.hist(zg_pythia_hist, axes=ax0)
    pythia_plot = ax0.hist(zg_pythias, label=pythia_label, bins=5 * n_bins, normed=1, histtype='step', color='blue', linewidth=5)
  else:
    pythia_plot = ax0.hist(zg_pythias, bins=5 * n_bins, normed=1, histtype='step', color='blue', linewidth=0)

  
  # Pythia Ends.
  
  # Herwig. 
  
  zg_herwig_hist = Hist(6 * n_bins, 0, 0.6, title=herwig_label, markersize=5.0, color='green', linewidth=5)
  bin_width_herwig = (zg_herwig_hist.upperbound() - zg_herwig_hist.lowerbound()) / zg_herwig_hist.nbins()

  map(zg_herwig_hist.Fill, zg_herwigs)

  zg_herwig_hist.Scale(1.0 / ( zg_herwig_hist.GetSumOfWeights() * bin_width_herwig ))

  if mc:
    # herwig_plot = rplt.hist(zg_herwig_hist, axes=ax0)
    herwig_plot = ax0.hist(zg_herwigs, label=herwig_label, bins=5 * n_bins, normed=1, histtype='step', color='green', linewidth=5)
  else:
    herwig_plot = ax0.hist(zg_herwigs, bins=5 * n_bins, normed=1, histtype='step', color='green', linewidth=0)
  
  # Herwig Ends.

  # Simulation Plots End.

  
  # Ratio-Over Plot Begins.


  # Theory-Over-Data Plot.
  
  data_plot_points_x = []
  data_plot_points_y = []

  for i in range(0, len(data_points_x)):
    if float(data_points_x[i]) >= float(zg_cut):
      data_plot_points_x.append(data_points_x[i])
      data_plot_points_y.append(data_points_y[i])


  theory_min_interpolate_function = interpolate.interp1d(theory_x, theory_y_min)
  theory_line_interpolate_function = interpolate.interp1d(theory_x, theory_y_line)
  theory_max_interpolate_function = interpolate.interp1d(theory_x, theory_y_max)

  theory_extrapolated_min = theory_min_interpolate_function(data_plot_points_x)
  theory_extrapolated_line = theory_line_interpolate_function(data_plot_points_x)
  theory_extrapolated_max = theory_max_interpolate_function(data_plot_points_x)

  if ratio_denominator == "data":

    if mc:
      zg_herwig_hist.Divide(zg_data_hist)
      zg_herwig_line_plot = convert_hist_to_line_plot(zg_herwig_hist, n_bins)
      plt.plot(zg_herwig_line_plot[0], zg_herwig_line_plot[1], linewidth=5, color='green')

      zg_pythia_hist.Divide(zg_data_hist)
      zg_pythia_line_plot = convert_hist_to_line_plot(zg_pythia_hist, n_bins)
      plt.plot(zg_pythia_line_plot[0], zg_pythia_line_plot[1], linewidth=5, color='blue')

    if data:
      ratio_data_to_data = [None if n == 0 else m / n for m, n in zip(data_plot_points_y, data_plot_points_y)]
      data_to_data_y_err = [(b / m) for b, m in zip(data_y_errors, data_plot_points_y)]
      data_to_data_x_err = [(b / m) for b, m in zip(data_x_errors, [1] * len(data_plot_points_y))]
      
      plt.errorbar(data_plot_points_x, ratio_data_to_data, xerr=data_to_data_x_err, yerr=data_to_data_y_err, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, color='black')

    if theory:
      ratio_theory_line_to_data = [m / n for m, n in zip(theory_extrapolated_line, data_plot_points_y)]
      ratio_theory_min_to_data = [m / n for m, n in zip(theory_extrapolated_min, data_plot_points_y)]
      ratio_theory_max_to_data = [m / n for m, n in zip(theory_extrapolated_max, data_plot_points_y)]

      zg_theory_line_to_data_hist = Hist(6 * n_bins, 0.0, 0.6)
      map(zg_theory_line_to_data_hist.Fill, data_plot_points_x, ratio_theory_line_to_data)
      zg_theory_line_to_data_plot = convert_hist_to_line_plot(zg_theory_line_to_data_hist, n_bins)
      plt.plot(zg_theory_line_to_data_plot[0], zg_theory_line_to_data_plot[1], linewidth=5, color='red')

      zg_theory_min_to_data_hist = Hist(6 * n_bins, 0.0, 0.6)
      map(zg_theory_min_to_data_hist.Fill, data_plot_points_x, ratio_theory_min_to_data)
      zg_theory_min_to_data_plot = convert_hist_to_line_plot(zg_theory_min_to_data_hist, n_bins)

      zg_theory_max_to_data_hist = Hist(6 * n_bins, 0.0, 0.6)
      map(zg_theory_max_to_data_hist.Fill, data_plot_points_x, ratio_theory_max_to_data)
      zg_theory_max_to_data_plot = convert_hist_to_line_plot(zg_theory_max_to_data_hist, n_bins)

      ax1.fill_between(zg_theory_max_to_data_plot[0], zg_theory_max_to_data_plot[1], zg_theory_min_to_data_plot[1], norm=1, where=np.less_equal(zg_theory_min_to_data_plot[1], zg_theory_max_to_data_plot[1]), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)
      
  elif ratio_denominator == "theory":

    zg_theory_line_hist = Hist(6 * n_bins, 0.0, 0.6, color='red')
    map(zg_theory_line_hist.Fill, data_plot_points_x, theory_extrapolated_line)

    zg_theory_min_hist = Hist(6 * n_bins, 0.0, 0.6, color='pink')
    map(zg_theory_min_hist.Fill, data_plot_points_x, theory_extrapolated_min)

    zg_theory_max_hist = Hist(6 * n_bins, 0.0, 0.6, color='red')
    map(zg_theory_max_hist.Fill, data_plot_points_x, theory_extrapolated_max)


    if mc:
      zg_herwig_hist.Divide(zg_theory_line_hist)
      zg_herwig_line_plot = convert_hist_to_line_plot(zg_herwig_hist, n_bins)
      plt.plot(zg_herwig_line_plot[0], zg_herwig_line_plot[1], linewidth=5, color='green')

      zg_pythia_hist.Divide(zg_theory_line_hist)
      zg_pythia_line_plot = convert_hist_to_line_plot(zg_pythia_hist, n_bins)
      plt.plot(zg_pythia_line_plot[0], zg_pythia_line_plot[1], linewidth=5, color='blue')

    if data:
      zg_data_to_th_y = [b / m for b, m in zip(data_plot_points_y, theory_extrapolated_line)]
      zg_data_to_th_y_err = [b / m for b, m in zip(data_y_errors, theory_extrapolated_line)]
      data_to_th_x_err = [(b / m) for b, m in zip(data_x_errors, [1] * len(zg_data_to_th_y_err))]

      plt.errorbar(data_plot_points_x, zg_data_to_th_y, xerr=data_to_th_x_err, yerr=zg_data_to_th_y_err, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, color='black')
   
    if theory:
      
      zg_theory_min_hist.Divide(zg_theory_line_hist)
      zg_theory_max_hist.Divide(zg_theory_line_hist)

      zg_theory_min_plot = convert_hist_to_line_plot(zg_theory_min_hist, n_bins)
      zg_theory_max_plot = convert_hist_to_line_plot(zg_theory_max_hist, n_bins)

      zg_theory_min_line, = plt.plot(zg_theory_min_plot[0], zg_theory_min_plot[1], linewidth=0)
      zg_theory_max_line, = plt.plot(zg_theory_max_plot[0], zg_theory_max_plot[1], linewidth=0)
      
      x_min, y_min = zg_theory_min_line.get_xdata(), zg_theory_min_line.get_ydata()
      x_max, y_max = zg_theory_max_line.get_xdata(), zg_theory_max_line.get_ydata()

      ax1.fill_between(x_max, y_max, y_min, norm=1, where=np.less_equal(y_min, y_max), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)

      zg_theory_line_hist.Divide(zg_theory_line_hist)
      zg_theory_line_plot = convert_hist_to_line_plot(zg_theory_line_hist, n_bins)
      plt.plot(zg_theory_line_plot[0], zg_theory_line_plot[1], linewidth=5, color='red')

  else:
    raise ValueError("Only 'theory' or 'data' are valid options for calculating ratios!")


  # Normalized-Over-Data Plot Ends.

  ax0.set_xlabel("$z_g$", fontsize=95)
  ax0.set_ylabel("$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", fontsize=80, rotation=0, labelpad=115, y=0.39)
  
  ax1.set_xlabel("$z_g$", fontsize=95)
  
  if ratio_denominator == "data":
    label_pad = 135
  else:
    label_pad = 115

  # ax1.set_ylabel("Ratio           \nto           \n" + ratio_denominator.capitalize() + "           ", fontsize=55, rotation=0, labelpad=250, y=0.31)
  plt.ylabel("Ratio           \nto           \n" + ratio_denominator.capitalize() + "           ", fontsize=55, rotation=0, labelpad=label_pad, y=0.31, axes=ax1)


  # Legend.

  th_line, = ax0.plot(range(1), linewidth=5, color='red')
  th_patch = mpatches.Patch(facecolor='pink', alpha=1.0, linewidth=0)

  if mc:
    pythia_line, = ax0.plot(range(1), linewidth=5, color=zg_pythia_hist.GetLineColor())
    herwig_line, = ax0.plot(range(1), linewidth=5, color=zg_herwig_hist.GetLineColor())
  else:
    pythia_line, = ax0.plot(range(1), linewidth=5, color=zg_pythia_hist.GetLineColor(), alpha=0)
    herwig_line, = ax0.plot(range(1), linewidth=5, color=zg_herwig_hist.GetLineColor(), alpha=0)

  handles = [data_plot, (th_patch, th_line), pythia_line, herwig_line]
  labels = [data_label, theory_label, pythia_label, herwig_label]

  first_legend = ax0.legend(handles, labels, fontsize=60, handler_map = {th_line : HandlerLine2D(marker_pad = 0), pythia_line : HandlerLine2D(marker_pad = 0), herwig_line : HandlerLine2D(marker_pad = 0)}, frameon=0, borderpad=0.1, bbox_to_anchor=[0.97, 0.98])
  ax = ax0.add_artist(first_legend)

  for txt in first_legend.get_texts():
    if ( not data) and txt.get_text() == data_label:
      txt.set_color("white") 

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  
  if pT_upper_cut != 20000:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$" + str(pT_lower_cut) + " < p_{T} < " + str(pT_upper_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]
  else:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]

  # labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5;~p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]
  
  ax0.legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.97, 0.52])


  # Legend Ends.

  ax0.autoscale(True)
  ax1.autoscale(True)
  
  ax0.set_ylim(0, y_max_limit)
  ax1.set_ylim(1.0 - y_limit_ratio_plot, 1.0 + y_limit_ratio_plot)

  ax0.set_xlim(0.0, 0.6)
  ax1.set_xlim(0.0, 0.6)
  

  fig = plt.gcf()

  # 1 - ((1 - 0.895) * 21.429)/30
  if data:
    ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.9249985), xycoords='figure fraction', frameon=0)
    plt.gca().add_artist(ab)
    preliminary_text = "Prelim. (20\%)"
    plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')
  else:
    preliminary_text = ""
    plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')



  fig = plt.gcf()
  fig.set_size_inches(30, 30, forward=1)

  plt.sca(ax0)
  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.02))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.5))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.sca(ax1)
  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.02))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)


  fig.set_snap(True)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "Writint out zg_cut_" + str(zg_filename) + "_pT_lower_" + str(pT_lower_cut) + "_pT_upper_" + str(pT_upper_cut) + "_ratio_over_" + ratio_denominator + "_th_" + str(theory) + "_mc_" + str(mc) + "_data_" + str(data) + ".pdf"
  filename = "plots/zg/zg_cut_" + str(zg_filename) + "_pT_lower_" + str(pT_lower_cut) + "_pT_upper_" + str(pT_upper_cut) + "_ratio_over_" + ratio_denominator + "_th_" + str(theory) + "_mc_" + str(mc) + "_data_" + str(data) + ".pdf"
  
  plt.savefig(filename)
  # plt.show()
  plt.clf()




def parse_theory_file(input_file):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  properties = defaultdict(list)

  points = defaultdict(list)
  for line in lines:
    try:
      numbers = line.split()
      try:
        for i in range(1, len(numbers)):
          points[float(numbers[0])].append(float(numbers[i]))
      except ValueError:
        pass
    except:
      pass

  return points



# 5 GeV bins.

def plot_trigger_efficiency_curves(trigger_1, trigger_2, pT_upper_limit=800):
  
  properties = parse_file_turn_on('./turn_on.dat')
  # properties = parse_file_turn_on('./trigger_proper_turn_on.dat')

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']

  colors = ['magenta', 'blue', 'orange', 'green', 'black', 'red']
  expected_trigger_names = ["HLT\_Jet180U", "HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U" ]

  color = colors[expected_trigger_names.index(trigger_1.replace("_", "\_"))]

  pt_hist_trigger_1 = Hist(500, 0, pT_upper_limit, title=trigger_1[4:], color=color, markersize=1.0, linewidth=5)
  pt_hist_trigger_2 = Hist(500, 0, pT_upper_limit, title=trigger_2[4:], color=color, markersize=1.0, linewidth=5)


  for i in range(0, len(pTs)):
    if trigger_1 in trigger_names[i]:
      pt_hist_trigger_1.Fill(pTs[i], prescales[i])

    # The len thingy is to make sure trigger names like HLT_Jet15U_HcalNoiseFiltered_v3 are excluded.
    # if trigger_2 in trigger_names[i] and len(trigger_names[i]) > (len(trigger_2) + 3):
    if trigger_2 in trigger_names[i]:
      pt_hist_trigger_2.Fill(pTs[i], prescales[i])


  pt_hist_trigger_1.Divide(pt_hist_trigger_2)


  rplt.errorbar(pt_hist_trigger_1, color=color,  markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  plt.gcf().text(0.29, 0.885, "Prelim. (20\%)", fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_snap(1)

  # Horizontal Line.
  plt.plot(list(pt_hist_trigger_1.x()), [1] * len(list(pt_hist_trigger_1.x())), color="black", linewidth=5, linestyle="dashed")

 
  if trigger_1 == "HLT_L1Jet6U":
    lower_pT = 37
  elif trigger_1 == "HLT_Jet15U":
    lower_pT = 56
  elif trigger_1 == "HLT_Jet30U":
    lower_pT = 84
  elif trigger_1 == "HLT_Jet50U":
    lower_pT = 114
  elif trigger_1 == "HLT_Jet70U":
    lower_pT = 153
  elif trigger_1 == "HLT_Jet100U":
    lower_pT = 196
  else:
    lower_pT = 0


  if lower_pT != 0:
    # CMS Vertical line.
    plt.plot([lower_pT, lower_pT], [plt.gca().get_ylim()[0], 1.], color=color, linewidth=3, linestyle="dashed")



  # # MOD Vertical line.
  # efficient_pt_x = 0.0
  # efficient_pt_y = 0.0

  # efficient_pt_x_s = []
  # efficient_pt_y_s = []

  # distance_to_previous_point_close_to_one = 0

  # for i in range(0, len(list(pt_hist_trigger_1.x()))):
  #   if abs(list(pt_hist_trigger_1.y())[i] - 1.00) < 0.1:

  #     if distance_to_previous_point_close_to_one > 25:
  #       # Last point close to one too far away.
  #       # Empty the lists.
  #       efficient_pt_x_s, efficient_pt_y_s = [], []

  #     efficient_pt_x_s.append(list(pt_hist_trigger_1.x())[i])
  #     efficient_pt_y_s.append(list(pt_hist_trigger_1.y())[i])
  #     distance_to_previous_point_close_to_one = 0

  #   else:
  #     distance_to_previous_point_close_to_one += 1


  #   if len(efficient_pt_x_s) > 75:
  #     efficient_pt_x = efficient_pt_x_s[0]
  #     efficient_pt_y = efficient_pt_y_s[0]
  #     break



  mod_efficient_pTs = [325, 260, 196, 153, 114, 84, 50, 32]
  mod_efficient_pT = mod_efficient_pTs[expected_trigger_names.index(trigger_1.replace("_", "\_"))]




  plt.plot([mod_efficient_pT, mod_efficient_pT], [plt.gca().get_ylim()[0], 1.], color="purple", linewidth=3, linestyle="dashed")
  

  if lower_pT != 0:
    plt.gca().annotate("CMS\n" + str(lower_pT) + " GeV", xy=(lower_pT, 1.), xycoords='data', xytext=(-100, 250),  textcoords='offset points', color=color, size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor=color, zorder=9999, connectionstyle="angle3,angleA=0,angleB=90") )
  
  plt.gca().annotate("MOD\n" + str(int(mod_efficient_pT)) + " GeV", xy=(mod_efficient_pT, 1.), xycoords='data', xytext=(250, 200), textcoords='offset points', color="purple", size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor="purple", zorder=9999, connectionstyle="angle3,angleA=45,angleB=-90") )



  plt.yscale('log')
  plt.gca().set_ylim(plt.gca().get_ylim()[0], 100)

  plt.legend(frameon=0)

  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=55, rotation=0)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=50.)

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.gcf().set_snap(True)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/trigger_efficiency/efficiency_curves_" + trigger_1 + ".pdf")
  # plt.show()
  plt.clf()









def plot_all_trigger_efficiency_curves():
  
  properties = parse_file_turn_on('./turn_on.dat')
  # properties = parse_file_turn_on('./trigger_proper_turn_on.dat')

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']


  # colors = ['orange', 'red', 'green', 'blue', 'magenta', 'black']
  # colors = colors[::-1]


  colors = ['black', 'magenta', 'blue', 'green', 'red', 'orange']
  expected_trigger_names = ["HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U", "HLT\_Jet15U\_HcalNoiseFiltered" ]
  labels = ["Jet140U / 100U", "Jet100U / 70U", "Jet70U / 50U", "Jet50U / 30U", "Jet30U / 15U\_HNF", "" ]

  cms_turn_on_pTs = [0, 0, 153, 114, 84]

  pt_hists = []
  for j in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist(60, 0, 300, color=colors[j], title=labels[j], markersize=1.0, linewidth=5))


  for i in range(0, len(expected_trigger_names)):
    for j in range(0, len(pTs)):
      # The len thingy is to make sure trigger names like HLT_Jet15U_HcalNoiseFiltered_v3 are excluded.
      if expected_trigger_names[i].replace("\\", "") in trigger_names[j]:
        pt_hists[i].Fill(pTs[j], prescales[j])

      
  ratio_hists = []
  for i in range(0, len(pt_hists) - 1):
    ratio_hists.append(pt_hists[i] / pt_hists[i + 1])


  for i in range(len(ratio_hists) - 1, -1, -1):
    rplt.errorbar(ratio_hists[i], emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
    
    if cms_turn_on_pTs[i] != 0:
      # plt.plot([cms_turn_on_pTs[i], cms_turn_on_pTs[i]], [plt.gca().get_ylim()[0], 1.], color=colors[i], linewidth=5, linestyle="dashed")
      plt.gca().annotate("CMS\n" + str(cms_turn_on_pTs[i]) + " GeV", xy=(cms_turn_on_pTs[i], 1.), xycoords='data', xytext=(-100, 250),  textcoords='offset points', color=colors[i], size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor=colors[i], zorder=99, connectionstyle="angle3,angleA=0,angleB=90") )


  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.9249985), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_snap(1)

  # Horizontal Line.
  plt.plot(list(pt_hists[0].x()), [1] * len(list(pt_hists[0].x())), color="black", linewidth=5, linestyle="dashed")


  plt.yscale('log')
  plt.gca().set_ylim(plt.gca().get_ylim()[0], 1000)

  # Legend.
  handles, labels = plt.gca().get_legend_handles_labels()
  first_legend = plt.gca().legend(handles, labels, fontsize=60,  frameon=0, borderpad=0.1, bbox_to_anchor=[0.90, 0.98])
  ax = plt.gca().add_artist(first_legend)

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5} $"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.77, 0.70])





  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=95, rotation=0)
  plt.ylabel('$\mathrm{Ratio}$', fontsize=75, rotation=0, labelpad=50.)


  plt.gcf().set_size_inches(30, 30, forward=1)
  plt.gcf().set_snap(True)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/all_efficiency_curves.pdf")
  # plt.show()
  plt.clf()



def plot_2d():

  pT_lower_cut = 150
  properties = parse_file(input_analysis_file, pT_lower_cut)
  zgs = [properties['zg_05'], properties['zg_1'], properties['zg_2']]
  charged_zgs = [properties['zg_charged_05'], properties['zg_charged_1'], properties['zg_charged_2']]
  prescales = properties['prescales']

  zgs = zgs[0]
  charged_zgs = charged_zgs[0]

  H, xedges, yedges = np.histogram2d(zgs, charged_zgs, bins=25, weights=prescales, normed=1, range=[[0.05, 0.5], [0.05, 0.5]] )


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
  
  Hmasked = np.ma.masked_where(H == 0, H) # Mask pixels with a value of zero

  plt.pcolormesh(xedges,yedges, Hmasked)

  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Counts')

  plt.xlabel('Charged $z_g$')
  plt.ylabel('$z_g$')

  plt.gcf().set_size_inches(30, 30, forward=1)
  plt.gcf().set_snap(True)

  

  plt.savefig("test.pdf")

  # plt.show()


def plot_JEC():
  pT_lower_cut = 100
  properties = parse_file(input_analysis_file, pT_lower_cut)

  JEC = properties['JEC']
  
  prescales = properties['prescales']

  plt.hist(JEC, weights=prescales, bins=100, label="JEC Factor", normed=1, histtype='step', linewidth=5)


  plt.xlabel('JEC Factor', fontsize=75)
  plt.ylabel('A.U.', fontsize=75, rotation=0, labelpad=100.)

  plt.autoscale(True)


  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60)
  plt.gca().add_artist(legend)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/JEC.pdf")

  plt.clf()




def plot_jet_area():
  pT_lower_cut = 100
  properties = parse_file(input_analysis_file, pT_lower_cut)

  jet_area = properties['jet_area']
  prescales = properties['prescales']


  plt.hist(jet_area, weights=prescales, bins=100, label="Jet Area", normed=1, histtype='step', linewidth=5)

  jet_area_expanded = []
  for i in range(0, len(jet_area)):
    for j in range(0, prescales[i]):
      jet_area_expanded.append(jet_area[i])


  mu, std = norm.fit(jet_area_expanded)
  xmin, xmax = plt.xlim()
  x = np.linspace(xmin, xmax, 500)
  p = norm.pdf(x, mu, std)
  plt.plot(x, p, 'k', linewidth=5)

  plt.autoscale(True)
  plt.gca().set_ylim(0., 1.1 * plt.gca().get_ylim()[1])

  plt.plot([ math.pi*(0.5**2), math.pi*(0.5**2) ], [ plt.gca().get_ylim()[0], plt.gca().get_ylim()[1] ], color='red', linewidth=5, linestyle="dashed")
  # plt.gca().annotate("$\pi \cdot 0.5^2$", xy=(math.pi*(0.5**2), 20.), xycoords='data', color='red', xytext=(0.9, 20), size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor='red', zorder=99) )

  plt.xlabel('Jet Area', fontsize=75)
  plt.ylabel('A.U.', fontsize=75, rotation=0, labelpad=100.)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60)
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  labels = ["Fit Statistics\n" + "$\mu=" + str(mu.round(3)) + "$\n$" + "\sigma=" + str(std.round(3)) + "$"]
  plt.gca().legend([extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.96, 0.80])


  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.02))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(1))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/jet_area.pdf")

  plt.clf()





def plot_charged_and_all_zgs(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', n_bins=10, y_max_limit=20):

  properties = parse_file(input_analysis_file, pT_lower_cut)

  zgs = properties[zg_filename]
  charged_zgs = properties[zg_filename[0:2] + "_charged_" + zg_filename[3: len(zg_filename)]]
  prescales = properties['prescales']
 
  zg_hist = Hist(6 * n_bins, 0.0, 0.6, title="All PF Candidates", markersize=1.0, color='black')
  zg_charged_hist = Hist(6 * n_bins, 0.0, 0.6, title="Charged PFCs", markersize=1.0, color='red')

  bin_width_zg = (zg_hist.upperbound() - zg_hist.lowerbound()) / zg_hist.nbins()
  bin_width_zg_charged = (zg_charged_hist.upperbound() - zg_charged_hist.lowerbound()) / zg_charged_hist.nbins()


  map(zg_hist.Fill, zgs, prescales)
  map(zg_charged_hist.Fill, charged_zgs, prescales)

  if zg_hist.GetSumOfWeights() != 0:
    zg_hist.Scale(1.0 / (zg_hist.GetSumOfWeights() * bin_width_zg))

  if zg_charged_hist.GetSumOfWeights() != 0:
    zg_charged_hist.Scale(1.0 / (zg_charged_hist.GetSumOfWeights() * bin_width_zg_charged))
  
  rplt.errorbar(zg_hist, emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.errorbar(zg_charged_hist, emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  
  
  legend = plt.legend(frameon=0, fontsize=60, bbox_to_anchor=[0.95, 0.98])
  plt.gca().add_artist(legend)

  
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra, extra, extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]

  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.65])


  plt.xlabel("$z_g$", fontsize=95)
  plt.ylabel("$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", fontsize=80, rotation=0, labelpad=115, y=0.39)
  

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gca().set_ylim(0, y_max_limit)
  plt.gca().set_xlim(0.0, 0.6)

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.02))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.5))

  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "Writing out charged_" + str(zg_filename) + "_pt_cut_" + str(pT_lower_cut) + ".pdf"

  plt.savefig("plots/zg_charged/charged_" + str(zg_filename) + "_pt_cut_" + str(pT_lower_cut) + ".pdf")
  # plt.show()
  plt.clf()



def plot_pts(pT_lower_cut=150, pT_upper_cut=10000):
  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)

  pTs = properties['uncorrected_hardest_pts']
  corrected_pTs = properties['corrected_hardest_pts']
  prescales = properties['prescales']

  herwig_pTs = parse_mc_pt_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_pt_herwig_pp2jj_150pTcut_7TeV.dat", pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)
  pythia_pTs = parse_mc_pt_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_pt_pythia_pp2jj_150pTcut_7TeV.dat", pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)


  pythia_pt_hist = Hist(300, pT_lower_cut, 1500, title="Pythia 8.205", linewidth=5, markersize=5.0, color="blue")
  bin_width_pythia = (pythia_pt_hist.upperbound() - pythia_pt_hist.lowerbound()) / pythia_pt_hist.nbins()

  herwig_pt_hist = Hist(300, pT_lower_cut, 1500, title="Herwig++ 2.6.3", linewidth=5, markersize=5.0, color="green")
  bin_width_herwig = (herwig_pt_hist.upperbound() - herwig_pt_hist.lowerbound()) / herwig_pt_hist.nbins()

  corrected_pt_hist = Hist(300, pT_lower_cut, 1500, title='Corrected', markersize=3.0, color='black')
  bin_width_corrected = (corrected_pt_hist.upperbound() - corrected_pt_hist.lowerbound()) / corrected_pt_hist.nbins()

  uncorrected_pt_hist = Hist(300, pT_lower_cut, 1500, title='Uncorrected', markersize=3.0, color='orange')
  bin_width_uncorrected = (uncorrected_pt_hist.upperbound() - uncorrected_pt_hist.lowerbound()) / uncorrected_pt_hist.nbins()

  map(uncorrected_pt_hist.Fill, pTs, prescales)
  map(corrected_pt_hist.Fill, corrected_pTs, prescales)
  
  map(pythia_pt_hist.Fill, pythia_pTs)
  map(herwig_pt_hist.Fill, herwig_pTs)

  corrected_pt_hist.Scale(1.0 / (corrected_pt_hist.GetSumOfWeights() * bin_width_corrected))
  uncorrected_pt_hist.Scale(1.0 / (uncorrected_pt_hist.GetSumOfWeights() * bin_width_uncorrected))
  pythia_pt_hist.Scale(1.0 / (pythia_pt_hist.GetSumOfWeights() * bin_width_pythia))
  herwig_pt_hist.Scale(1.0 / (herwig_pt_hist.GetSumOfWeights() * bin_width_herwig))

  
  gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 

  ax0 = plt.subplot(gs[0])
  ax1 = plt.subplot(gs[1])




  data_plot = rplt.errorbar(corrected_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  uncorrected_data_plot = rplt.errorbar(uncorrected_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.hist(pythia_pt_hist, axes=ax0, emptybins=False, marker='o',  markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.hist(herwig_pt_hist, axes=ax0, emptybins=False, marker='o',  markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  

  data_x_errors, data_y_errors = [], []
  for x_segment in data_plot[2][0].get_segments():
    data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
  for y_segment in data_plot[2][1].get_segments():
    data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

  data_points_x = data_plot[0].get_xdata()
  data_points_y = data_plot[0].get_ydata()

  data_plot_points_x = []
  data_plot_points_y = []
  for i in range(0, len(data_points_x)):
    data_plot_points_x.append(data_points_x[i])
    data_plot_points_y.append(data_points_y[i])



  uncorrected_data_x_errors, uncorrected_data_y_errors = [], []
  for x_segment in uncorrected_data_plot[2][0].get_segments():
    uncorrected_data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
  for y_segment in uncorrected_data_plot[2][1].get_segments():
    uncorrected_data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

  uncorrected_data_points_x = uncorrected_data_plot[0].get_xdata()
  uncorrected_data_points_y = uncorrected_data_plot[0].get_ydata()

  uncorrected_data_plot_points_x = []
  uncorrected_data_plot_points_y = []
  for i in range(0, len(uncorrected_data_points_x)):
    uncorrected_data_plot_points_x.append(uncorrected_data_points_x[i])
    uncorrected_data_plot_points_y.append(uncorrected_data_points_y[i])




  data_to_data_y_err = [(b / m) for b, m in zip(data_y_errors, data_plot_points_y)]
  data_to_data_x_err = [(b / m) for b, m in zip(data_x_errors, [1] * len(data_plot_points_y))]

  uncorrected_to_corrected_y_err = [(b / m) for b, m in zip(uncorrected_data_y_errors, data_plot_points_y)]
  uncorrected_to_corrected_x_err = [(b / m) for b, m in zip(uncorrected_data_x_errors, [1] * len(data_plot_points_y))]


  # Legends Begin.

  legend = ax0.legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[0.95, 1.0])
  ax0.add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  if pT_upper_cut != 10000:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$" + str(pT_lower_cut) + " < p_{T} < " + str(pT_upper_cut) + "~\mathrm{GeV}$"]
  else:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$"]
  ax0.legend([extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.57])

  # Legends End.



  ax0.set_xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75)
  ax1.set_xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75)
  ax0.set_ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)
  ax1.set_ylabel("Ratio           \nto           \n" + "Data" + "           ", fontsize=55, rotation=0, labelpad=115, y=0.31)


  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.9249985), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  # Ratio Plot.
  pythia_pt_hist.Divide(corrected_pt_hist)
  herwig_pt_hist.Divide(corrected_pt_hist)
  uncorrected_pt_hist.Divide(corrected_pt_hist)
  corrected_pt_hist.Divide(corrected_pt_hist)

  rplt.hist(pythia_pt_hist, axes=ax1, linewidth=5)
  rplt.hist(herwig_pt_hist, axes=ax1, linewidth=5)
  
  rplt.errorbar(corrected_pt_hist, xerr=data_to_data_x_err, yerr=data_to_data_y_err, axes=ax1, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.errorbar(uncorrected_pt_hist, xerr=uncorrected_to_corrected_x_err, yerr=uncorrected_to_corrected_y_err, axes=ax1, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  ax0.set_yscale('log')

  ax0.autoscale(True)
  ax1.autoscale(True)
  
  ax1.set_ylim(0., 2.)



  if ((pT_lower_cut == 150 and pT_upper_cut == 250)):
    ax0.set_xlim(150, 250)
    ax1.set_xlim(150, 250)
    ax0.set_ylim(10e-4, 10e-2)
    pT_minor_ticks = 5
  elif ((pT_lower_cut == 250 and pT_upper_cut == 500)):
    ax0.set_xlim(250, 500)
    ax1.set_xlim(250, 500)
    ax0.set_ylim(10e-5, 10e-2)
    pT_minor_ticks = 10
  elif ((pT_lower_cut == 500 and pT_upper_cut == 10000)):
    ax0.set_xlim(500, 1500)
    ax1.set_xlim(500, 1500)
    ax0.set_ylim(0.0005, 0.05)
    pT_minor_ticks = 50
  elif ((pT_lower_cut == 150 and pT_upper_cut == 10000)):
    ax0.set_xlim(150, 1500)
    ax1.set_xlim(150, 1500)
    ax0.set_ylim(10e-8, 10e-2)
    pT_minor_ticks = 50


  plt.gcf().set_size_inches(30, 30, forward=1)

  plt.sca(ax0)
  plt.gca().xaxis.set_minor_locator(MultipleLocator(pT_minor_ticks))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.sca(ax1)
  plt.gca().xaxis.set_minor_locator(MultipleLocator(pT_minor_ticks))
  # plt.gca().yaxis.set_minor_locator(MultipleLocator(50))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "Printing fractional energy loss with pT > " + str(pT_lower_cut) + " and pT < " + str(pT_upper_cut)

  plt.savefig("plots/pT_distribution/pT_lower_" + str(pT_lower_cut) + "_pT_upper_" + str(pT_upper_cut) + ".pdf")
  # plt.show()
  plt.clf()




def plot_hardest_pt_softdrop(pT_lower_cut=100, pT_upper_cut=20000):
  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)

  pTs_before_SD = properties['corrected_hardest_pts']
  pTs_after_SD = properties['pTs_after_SD']  
  prescales = properties['prescales']

  pT_before_SD_hist = Hist(150, 0, 1500, title='Before SoftDrop', markersize=3.0, color='black')
  bin_width_before = (pT_before_SD_hist.upperbound() - pT_before_SD_hist.lowerbound()) / pT_before_SD_hist.nbins()

  pT_after_SD_hist = Hist(150, 0, 1500, title='After SoftDrop', markersize=3.0, color='red')
  bin_width_after = (pT_after_SD_hist.upperbound() - pT_after_SD_hist.lowerbound()) / pT_after_SD_hist.nbins()

  map(pT_before_SD_hist.Fill, pTs_before_SD, prescales)
  map(pT_after_SD_hist.Fill, pTs_after_SD, prescales)
  
  pT_before_SD_hist.Scale(1.0 / (pT_before_SD_hist.GetSumOfWeights() * bin_width_before))
  pT_after_SD_hist.Scale(1.0 / (pT_after_SD_hist.GetSumOfWeights() * bin_width_after))
  
  pT_before_SD_plot = rplt.errorbar(pT_before_SD_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  pT_after_SD_plot = rplt.errorbar(pT_after_SD_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  
  plt.yscale('log')

  # Legends Begin.

  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[0.89, 1.0])
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  if pT_upper_cut != 20000:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$" + str(pT_lower_cut) + " < p_{T} < " + str(pT_upper_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  else:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  plt.gca().legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.70])

  # Legends End.

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=80.)
  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gca().autoscale(True)
  # plt.ylim(10e-9, 10e-1)
  # plt.gca().set_ylim(0., 1.5 * plt.gca().get_ylim()[1])

  if ((pT_lower_cut == 100 and pT_upper_cut == 200)):
    plt.ylim(10e-5, 10e-1)
  elif ((pT_lower_cut == 200 and pT_upper_cut == 400)):
    plt.ylim(10e-6, 10e-1)
    # plt.xlim(100, 600)
  elif ((pT_lower_cut == 400 and pT_upper_cut == 20000)):
    plt.ylim(10e-6, 10e-2)
  elif ((pT_lower_cut == 100 and pT_upper_cut == 20000)):
    plt.ylim(10e-9, 10e-1)


  plt.gca().xaxis.set_minor_locator(MultipleLocator(50))
  # plt.gca().yaxis.set_minor_locator(MultipleLocator(50))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "Printing fractional energy loss with pT > " + str(pT_lower_cut) + " and pT < " + str(pT_upper_cut)

  plt.savefig("plots/hardest_pt_softdrop/pT_lower_" + str(pT_lower_cut) + "_pT_upper_" + str(pT_upper_cut) + ".pdf")
  # plt.show()
  plt.clf()



def plot_constituent_multiplicity_softdrop(pT_lower_cut=100, pT_upper_cut=20000):
  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)

  multi_before_SD = properties['multiplicity_before_SD']
  multi_after_SD = properties['multiplicity_after_SD']  
  prescales = properties['prescales']

  multi_before_SD_hist = Hist(150, 0, 150, title='Before SoftDrop', markersize=3.0, color='black')
  bin_width_before = (multi_before_SD_hist.upperbound() - multi_before_SD_hist.lowerbound()) / multi_before_SD_hist.nbins()

  multi_after_SD_hist = Hist(150, 0, 150, title='After SoftDrop', markersize=3.0, color='red')
  bin_width_after = (multi_after_SD_hist.upperbound() - multi_after_SD_hist.lowerbound()) / multi_after_SD_hist.nbins()

  map(multi_before_SD_hist.Fill, multi_before_SD, prescales)
  map(multi_after_SD_hist.Fill, multi_after_SD, prescales)
  
  multi_before_SD_hist.Scale(1.0 / (multi_before_SD_hist.GetSumOfWeights() * bin_width_before))
  multi_after_SD_hist.Scale(1.0 / (multi_after_SD_hist.GetSumOfWeights() * bin_width_after))
  
  pT_before_SD_plot = rplt.errorbar(multi_before_SD_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  pT_after_SD_plot = rplt.errorbar(multi_after_SD_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  
  # plt.yscale('log')

  # Legends Begin.

  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[0.89, 1.0])
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  if pT_upper_cut != 20000:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$" + str(pT_lower_cut) + " < p_{T} < " + str(pT_upper_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  else:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  plt.gca().legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.70])

  # Legends End.

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.xlabel('Constituent Multiplicity', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)
  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gca().autoscale(True)
  plt.gca().set_ylim(0., 1.1 * plt.gca().get_ylim()[1])

  plt.gca().xaxis.set_minor_locator(MultipleLocator(5))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.001))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "Printing fractional energy loss with pT > " + str(pT_lower_cut) + " and pT < " + str(pT_upper_cut)

  plt.savefig("plots/constituent_multiplicity_softdrop/pT_lower_" + str(pT_lower_cut) + "_pT_upper_" + str(pT_upper_cut) + ".pdf")
  # plt.show()
  plt.clf()





def plot_fractional_energy_loss(pT_lower_cut=100, pT_upper_cut=20000):
  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)

  fractional_energy_loss = properties['fractional_energy_loss']
  prescales = properties['prescales']

  fractional_energy_loss_hist = Hist(150, 0, 0.5, title='Fractional Jet Energy Loss', markersize=3.0, color='black')
  bin_width_before = (fractional_energy_loss_hist.upperbound() - fractional_energy_loss_hist.lowerbound()) / fractional_energy_loss_hist.nbins()

  map(fractional_energy_loss_hist.Fill, fractional_energy_loss, prescales)
 
  fractional_energy_loss_hist.Scale(1.0 / (fractional_energy_loss_hist.GetSumOfWeights() * bin_width_before))

  rplt.errorbar(fractional_energy_loss_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  
  # plt.yscale('log')

  # Legends Begin.

  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60)
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  if pT_upper_cut != 20000:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$" + str(pT_lower_cut) + " < p_{T} < " + str(pT_upper_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  else:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  plt.gca().legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.93, 0.75])

  # Legends End.

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.xlabel('Fractional Energy Loss', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)
  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gca().autoscale(True)
  plt.gca().set_ylim(0., plt.gca().get_ylim()[1])

  plt.gca().xaxis.set_minor_locator(MultipleLocator(0.01))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(2))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "Printing fractional energy loss with pT > " + str(pT_lower_cut) + " and pT < " + str(pT_upper_cut)

  plt.savefig("plots/fractional_energy_loss/pT_lower_" + str(pT_lower_cut) + "_pT_upper_" + str(pT_upper_cut) + ".pdf")
  # plt.show()
  plt.clf()







def plot_jet_mass_spectrum(pT_lower_cut=100, pT_upper_cut=20000):
  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)

  jet_mass_before_SD = properties['jet_mass_before_SD']
  jet_mass_after_SD = properties['jet_mass_after_SD']  
  prescales = properties['prescales']

  jet_mass_before_SD_hist = Hist(150, 0, 150, title='Before SoftDrop', markersize=3.0, color='black')
  bin_width_before = (jet_mass_before_SD_hist.upperbound() - jet_mass_before_SD_hist.lowerbound()) / jet_mass_before_SD_hist.nbins()

  jet_mass_after_SD_hist = Hist(150, 0, 150, title='After SoftDrop', markersize=3.0, color='red')
  bin_width_after = (jet_mass_after_SD_hist.upperbound() - jet_mass_after_SD_hist.lowerbound()) / jet_mass_after_SD_hist.nbins()

  map(jet_mass_before_SD_hist.Fill, jet_mass_before_SD, prescales)
  map(jet_mass_after_SD_hist.Fill, jet_mass_after_SD, prescales)
  
  jet_mass_before_SD_hist.Scale(1.0 / (jet_mass_before_SD_hist.GetSumOfWeights() * bin_width_before))
  jet_mass_after_SD_hist.Scale(1.0 / (jet_mass_after_SD_hist.GetSumOfWeights() * bin_width_after))
  
  rplt.errorbar(jet_mass_before_SD_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.errorbar(jet_mass_after_SD_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  # Legends Begin.

  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[0.89, 1.0])
  plt.gca().add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  
  if pT_upper_cut != 20000:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$" + str(pT_lower_cut) + " < p_{T} < " + str(pT_upper_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  else:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = 0.05}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = 0.05}$"]
  

  plt.gca().legend([extra, extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.69])

  # # Legends End.

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.xlabel('Jet Mass', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=80.)
  
  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(5))
  plt.gca().yaxis.set_minor_locator(MultipleLocator(0.005))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.gca().autoscale(True)
  plt.gca().set_ylim(0., plt.gca().get_ylim()[1] * 1.2)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "Printing jet mass spectrum with pT > " + str(pT_lower_cut) + " and pT < " + str(pT_upper_cut)

  plt.savefig("plots/jet_mass_spectrum/pT_lower_" + str(pT_lower_cut) + "_pT_upper_" + str(pT_upper_cut) + ".pdf")
  # plt.show()
  plt.clf()



def plot_zg():
  properties = parse_file(input_analysis_file, 150)

  zgs = properties['zg_05'] 
  prescales = properties['prescales']



  # plt.gca().set_xscale('log')

  n_bins = 100

  b, a = log_bins(zgs, prescales, n_bins)

  plt.hist(a[:-1], weights=b, bins=n_bins, normed=1)

  plt.show()




def logged_bin(data, weights, number_of_bins=50):
  
  def drop_zeros(a_list):
    return [i for i in a_list if i > 0]

  # min_value = min( drop_zeros(data) )
  # max_value = max(data)

  min_value = math.log( min( drop_zeros(data) ), 10 )
  max_value = math.log( max(data), 10 )

  return np.histogram(data, weights=weights, bins=np.logspace(min_value, max_value, number_of_bins))
  # return np.histogram(data, weights=weights, bins=np.linspace(min_value, max_value, number_of_bins))



def linear_bin(data, weights, number_of_bins=50):
  
  def drop_zeros(a_list):
    return [i for i in a_list if i > 0]

  min_value = min( drop_zeros(data) )
  max_value = max(data)

  return np.histogram(data, weights=weights, bins=np.linspace(min_value, max_value, number_of_bins))




def test():

  n_bins = 50

  # x = np.arange(5, 10, 0.0001)
  x = np.linspace(0, 2, 10000)
  # y = np.reciprocal(x)

  # x = np.log(x)

  # hist, bins = logged_bin(x, y, 25)
  # hist, bins = linear_bin(x, y, n_bins)

  bins = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1, 2]

  # hist, bins = np.histogram(x, bins=bins, density=True)


  # plt.hist(bins[:-1], weights=hist)

  plt.hist(x, bins=bins, normed=1)


  # plt.hist(bins[:-1], weights=hist, bins=n_bins, color='orange', alpha=0.75, lw=5)
  # plt.errorbar(bins[:-1], hist, lw=0, xerr=True, yerr=True, elinewidth=3, capsize=5, marker="o", markersize=5)


  # plt.gca().set_xscale('log')

  # plt.ylim(0, 8)

  plt.autoscale()

  plt.show()


# test()

def plot_zg_test():
  properties = parse_file(input_analysis_file, 150)

  zgs = properties['zg_02']
  prescales = properties['prescales']

  x = zgs
  y = prescales

  hist, bins = logged_bin(x, y, 25)

  width = 0.7 * (bins[1] - bins[0])
  center = (bins[:-1] + bins[1:]) / 2
  # plt.bar(center, hist, align='center', width=width)

  plt.hist(bins[:-1], weights=hist, color='orange', alpha=0.75, lw=5, bins=200)
  
  # plt.gca().set_xscale('log')

  plt.autoscale()

  plt.show()



def plot_log_zg_th_mc_data(pT_lower_cut, pT_upper_cut, zg_cut, zg_filename, ratio_denominator="theory", data=True, mc=True, theory=True, n_bins=10, y_max_limit=20, y_limit_ratio_plot=0.5):
  pfc_pT_cut = 0

  zg_cut = float(zg_cut)

  properties = parse_file(input_analysis_file, pT_lower_cut, pfc_pT_cut, pT_upper_cut)
  properties_pythia = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_pythia_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut, pT_upper_cut)
  properties_herwig = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_herwig_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut, pT_upper_cut)

  zg_data = properties[zg_filename]
  
  zg_pythias = properties_pythia[zg_filename]
  zg_herwigs = properties_herwig[zg_filename]

  prescales = properties['prescales']

  data_label = "CMS 2010 Open Data"
  pythia_label = "Pythia 8.205" if mc else ""
  herwig_label = "Herwig++ 2.6.3" if mc else ""
  theory_label = "Theory (MLL)" if theory else ""
  
  gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 

  ax0 = plt.subplot(gs[0])
  ax1 = plt.subplot(gs[1])

  def pyplot_hist_to_plot(x, y):
    a = []
    b = []
    for i in range(0, len(x[:-1])):
      if y[i] != 0:
        a.append(x[i])
        b.append(y[i])
    return a, b

  # Theory Plots Begin.
  
  points_th_gluon = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Simone/results_7_24_15/band_gluon_pt" + str(pT_lower_cut) + "_zc" + str(zg_cut).replace(".", "") + ".dat")
  points_th_quark = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Simone/results_7_24_15/band_quark_pt" + str(pT_lower_cut) + "_zc" + str(zg_cut).replace(".", "") + ".dat")

  points = defaultdict(list)

  for x in points_th_gluon:
    points[x] = [ points_th_gluon[x][0], points_th_gluon[x][1], points_th_gluon[x][2], points_th_gluon[x][3], points_th_gluon[x][4], points_th_gluon[x][5] ]
    points[x].extend([ points_th_quark[x][0], points_th_quark[x][1], points_th_quark[x][2], points_th_quark[x][3], points_th_quark[x][4], points_th_quark[x][5] ])

  keys = points.keys()
  keys.sort()

  theory_x = keys

  y = []
  for j in range(0, 6):
    y.append([points[i][j] for i in keys])

  # For each x, record three y's viz. max_y, min_y, line_y (i.e. e11 xmu=1).
  weighted_ys = []
  for i in range(0, len(y)):
    area = simps(y[i], theory_x)
    weighted = []
    for j in range(0, len(y[i])):
      weighted.append( y[i][j] / area )    
    weighted_ys.append(weighted)

  y = weighted_ys

  theory_y_max = []
  theory_y_min = []
  theory_y_line = []
  for i in range(0, len(theory_x)):
    y_for_current_x = []
    for j in range(0, 6):
      y_for_current_x.append(y[j][i])

    theory_y_min.append(theory_x[i] * min(y_for_current_x))
    theory_y_line.append(theory_x[i] * y_for_current_x[1])
    theory_y_max.append(theory_x[i] * max(y_for_current_x))
  
  bins_linear_log = np.linspace(math.log(zg_cut, math.e), math.log(0.5, math.e), 48)
  log_theory_x = np.log(theory_x)

  if theory:
    theory_x_logged = np.log(theory_x)

    y, x = np.histogram(theory_x_logged, weights=theory_y_max, bins=bins_linear_log)
    a, b = pyplot_hist_to_plot(x, y)
    b_max = [x / simps(b, a) for x in b]
    ax0.plot(a, b_max, label=theory_label, color='magenta', lw=5)

    
    y, x = np.histogram(theory_x_logged, weights=theory_y_line, bins=bins_linear_log)
    a, b = pyplot_hist_to_plot(x, y)
    b_line = [x / simps(b, a) for x in b]
    ax0.plot(a, b_line, label=theory_label, color='red', lw=5)

    y, x = np.histogram(theory_x_logged, weights=theory_y_min, bins=bins_linear_log)
    a, b = pyplot_hist_to_plot(x, y)
    b_min = [x / simps(b, a) for x in b]
    ax0.plot(a, b_min, label=theory_label, color='brown', lw=5)

    ax0.fill_between(a, b_max, b_min, norm=1, where=np.less_equal(b_min, b_max), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)
  
  # Theory Plot Ends.

  def convert_hist_to_line_plot(hist, n_bins):
    a = []
    b = {}
    # bin_width = 0.6 / (6 * n_bins)
    bin_width = (hist.upperbound() - hist.lowerbound()) / hist.nbins()

    for i in range(0, len(list(hist.x()))):
      a.append(round(list(hist.x())[i] - bin_width / 2., 4))
      a.append(round(list(hist.x())[i], 4))
      a.append(round(list(hist.x())[i] + bin_width / 2., 4))

      if round(list(hist.x())[i] - bin_width / 2., 4) not in b.keys():
        b[round(list(hist.x())[i] - bin_width / 2., 4)] = [ list(hist.y())[i] ]
      else:
        b[round(list(hist.x())[i] - bin_width / 2., 4)].append( list(hist.y())[i] )
      
      if round(list(hist.x())[i], 4) not in b.keys():
        b[round(list(hist.x())[i], 4)] = [ list(hist.y())[i] ]
      else:
        b[round(list(hist.x())[i], 4)].append( list(hist.y())[i] )

      if round(list(hist.x())[i] + bin_width / 2., 4) not in b.keys():
        b[round(list(hist.x())[i] + bin_width / 2., 4)] = [ list(hist.y())[i] ]
      else:
        b[round(list(hist.x())[i] + bin_width / 2., 4)].append( list(hist.y())[i] )

    x = sorted(list(Set(a)))
    a.sort()
    
    c = [b[x[i]] for i in range(0, len(x))]

    y = [item for sublist in c for item in sublist]
    
    a_zero_removed = []
    y_zero_removed = []
    for i in range(0, len(a)):
      if a[i] >= math.log(zg_cut, math.e) and a[i] <= math.log(0.5, math.e) and y[i] != 0.0:
      # if True:
        a_zero_removed.append(a[i])
        y_zero_removed.append(y[i])

    return a_zero_removed, y_zero_removed

  # Data Plot Begins.
  log_zg_data = np.log(zg_data)
  updated_prescales = []
  for i in range(0, len(zg_data)):
  	updated_prescales.append(prescales[i] * zg_data[i])

  # prescales = updated_prescales

  zg_data_hist = Hist(bins_linear_log, title=data_label, markersize=2.5, color='black')
  bin_width_data = (zg_data_hist.upperbound() - zg_data_hist.lowerbound()) / zg_data_hist.nbins()
  map(zg_data_hist.Fill, log_zg_data, prescales)
  zg_data_hist.Scale(1.0 / ( zg_data_hist.GetSumOfWeights() * bin_width_data ))

  if data:
    data_plot = rplt.errorbar(zg_data_hist, xerr=1, yerr=1, emptybins=False, axes=ax0, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
  else:
    data_plot = rplt.errorbar(zg_data_hist, xerr=1, yerr=1, emptybins=False, axes=ax0, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=0.0)
  


  data_x_errors, data_y_errors = [], []
  for x_segment in data_plot[2][0].get_segments():
    data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
  for y_segment in data_plot[2][1].get_segments():
    data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

  data_points_x = data_plot[0].get_xdata()
  data_points_y = data_plot[0].get_ydata()

  # Data Plots Ends.

  # Simulation Plots Begin. 
  
  # Pythia.
  log_zg_pythias = np.log(zg_pythias)
  y, x = np.histogram(log_zg_pythias, bins=bins_linear_log, normed=True)
  a, b = pyplot_hist_to_plot(x, y)
  log_zg_pythia_hist = Hist(bins_linear_log)
  map(log_zg_pythia_hist.Fill, a, b)

  if mc:
    pythia_plot = ax0.hist(a, histtype='step', normed=True, weights=b, bins=48, label=pythia_label, lw=5, color='blue')
  else:
    pythia_plot = ax0.hist(a, histtype='step', normed=True, weights=b, bins=48, label=pythia_label, lw=0, color='blue')

  # Pythia Ends.
  
  # Herwig. 
  log_zg_herwigs = np.log(zg_herwigs)
  y, x = np.histogram(log_zg_herwigs, bins=bins_linear_log, normed=True)
  a, b = pyplot_hist_to_plot(x, y)
  log_zg_herwig_hist = Hist(bins_linear_log)
  map(log_zg_herwig_hist.Fill, a, b)

  if mc:
    herwig_plot = ax0.hist(a, histtype='step', normed=True, weights=b, bins=48, label=herwig_label, lw=5, color='green')
  else:
    herwig_plot = ax0.hist(a, histtype='step', normed=True, weights=b, bins=48, label=herwig_label, lw=0, color='green')
  
  # Herwig Ends.

  # Simulation Plots End.
  
  # Ratio-Over Plot Begins.

  # Theory-Over-Data Plot.
  
  data_plot_points_x = []
  data_plot_points_y = []

  for i in range(0, len(data_points_x)):
    if float(data_points_x[i]) >= float(math.log(zg_cut, math.e)):
    # if True:
      data_plot_points_x.append(data_points_x[i])
      data_plot_points_y.append(data_points_y[i])




  theory_min_interpolate_function = interpolate.interp1d(log_theory_x, theory_y_min)
  theory_line_interpolate_function = interpolate.interp1d(log_theory_x, theory_y_line)
  theory_max_interpolate_function = interpolate.interp1d(log_theory_x, theory_y_max)


  # data_plot_points_x = sorted(data_plot_points_x, reverse=True)

  theory_extrapolated_min = theory_min_interpolate_function(data_plot_points_x)
  theory_extrapolated_line = theory_line_interpolate_function(data_plot_points_x)
  theory_extrapolated_max = theory_max_interpolate_function(data_plot_points_x)

  if ratio_denominator == "data":

    if mc:
      log_zg_herwig_hist.Divide(zg_data_hist)
      zg_herwig_line_plot = convert_hist_to_line_plot(log_zg_herwig_hist, len(bins_linear_log))
      plt.plot(zg_herwig_line_plot[0], zg_herwig_line_plot[1], linewidth=5, color='green')

      log_zg_pythia_hist.Divide(zg_data_hist)
      zg_pythia_line_plot = convert_hist_to_line_plot(log_zg_pythia_hist, n_bins)
      plt.plot(zg_pythia_line_plot[0], zg_pythia_line_plot[1], linewidth=5, color='blue')

    if data:
      ratio_data_to_data = [None if n == 0 else m / n for m, n in zip(data_plot_points_y, data_plot_points_y)]
      data_to_data_y_err = [(b / m) for b, m in zip(data_y_errors, data_plot_points_y)]
      data_to_data_x_err = [(b / m) for b, m in zip(data_x_errors, [1] * len(data_plot_points_y))]
      
      plt.errorbar(data_plot_points_x, ratio_data_to_data, xerr=data_to_data_x_err, yerr=data_to_data_y_err, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, color='black')

    if theory:
      ratio_theory_line_to_data = [m / n for m, n in zip(theory_extrapolated_line, data_plot_points_y)]
      ratio_theory_min_to_data = [m / n for m, n in zip(theory_extrapolated_min, data_plot_points_y)]
      ratio_theory_max_to_data = [m / n for m, n in zip(theory_extrapolated_max, data_plot_points_y)]


      zg_theory_line_to_data_hist = Hist(bins_linear_log)
      normalized_y_line = [x / simps(ratio_theory_line_to_data, data_plot_points_x) for x in ratio_theory_line_to_data]
      map(zg_theory_line_to_data_hist.Fill, data_plot_points_x, normalized_y_line)
      zg_theory_line_to_data_plot = convert_hist_to_line_plot(zg_theory_line_to_data_hist, n_bins)
      plt.plot(zg_theory_line_to_data_plot[0], zg_theory_line_to_data_plot[1], linewidth=5, color='red')

      zg_theory_min_to_data_hist = Hist(bins_linear_log)
      normalized_y_min = [x / simps(ratio_theory_min_to_data, data_plot_points_x) for x in ratio_theory_min_to_data]
      map(zg_theory_min_to_data_hist.Fill, data_plot_points_x, normalized_y_min)
      zg_theory_min_to_data_plot = convert_hist_to_line_plot(zg_theory_min_to_data_hist, n_bins)
      # plt.plot(zg_theory_min_to_data_plot[0], zg_theory_min_to_data_plot[1], linewidth=5, color='orange')

      zg_theory_max_to_data_hist = Hist(bins_linear_log)
      normalized_y_max = [x / simps(ratio_theory_max_to_data, data_plot_points_x) for x in ratio_theory_max_to_data]
      map(zg_theory_max_to_data_hist.Fill, data_plot_points_x, normalized_y_max)
      zg_theory_max_to_data_plot = convert_hist_to_line_plot(zg_theory_max_to_data_hist, n_bins)
      # plt.plot(zg_theory_max_to_data_plot[0], zg_theory_max_to_data_plot[1], linewidth=5, color='magenta')

      ax1.fill_between(zg_theory_max_to_data_plot[0], zg_theory_max_to_data_plot[1], zg_theory_min_to_data_plot[1], norm=1, where=np.less_equal(zg_theory_min_to_data_plot[1], zg_theory_max_to_data_plot[1]), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)
  
  elif ratio_denominator == "theory":

    zg_theory_line_hist = Hist(bins_linear_log, color='red')
    map(zg_theory_line_hist.Fill, data_plot_points_x, theory_extrapolated_line)

    zg_theory_min_hist = Hist(bins_linear_log, color='pink')
    map(zg_theory_min_hist.Fill, data_plot_points_x, theory_extrapolated_min)

    zg_theory_max_hist = Hist(bins_linear_log, color='red')
    map(zg_theory_max_hist.Fill, data_plot_points_x, theory_extrapolated_max)

  
    if mc:
      log_zg_herwig_hist.Divide(zg_theory_line_hist)
      zg_herwig_line_plot = convert_hist_to_line_plot(log_zg_herwig_hist, n_bins)
      plt.plot(zg_herwig_line_plot[0], zg_herwig_line_plot[1], linewidth=5, color='green')

      log_zg_pythia_hist.Divide(zg_theory_line_hist)
      zg_pythia_line_plot = convert_hist_to_line_plot(log_zg_pythia_hist, n_bins)
      plt.plot(zg_pythia_line_plot[0], zg_pythia_line_plot[1], linewidth=5, color='blue')


    if data:
      zg_data_to_th_y = [b / m for b, m in zip(data_plot_points_y, theory_extrapolated_line)]
      zg_data_to_th_y_err = [b / m for b, m in zip(data_y_errors, theory_extrapolated_line)]
      data_to_th_x_err = [(b / m) for b, m in zip(data_x_errors, [1] * len(zg_data_to_th_y_err))]

      plt.errorbar(data_plot_points_x, zg_data_to_th_y, xerr=data_to_th_x_err, yerr=zg_data_to_th_y_err, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, color='black')
  

    if theory:
      
      zg_theory_min_hist.Divide(zg_theory_line_hist)
      zg_theory_max_hist.Divide(zg_theory_line_hist)

      zg_theory_min_plot = convert_hist_to_line_plot(zg_theory_min_hist, n_bins)
      zg_theory_max_plot = convert_hist_to_line_plot(zg_theory_max_hist, n_bins)

      zg_theory_min_line, = plt.plot(zg_theory_min_plot[0], zg_theory_min_plot[1], linewidth=0)
      zg_theory_max_line, = plt.plot(zg_theory_max_plot[0], zg_theory_max_plot[1], linewidth=0)
      
      x_min, y_min = zg_theory_min_line.get_xdata(), zg_theory_min_line.get_ydata()
      x_max, y_max = zg_theory_max_line.get_xdata(), zg_theory_max_line.get_ydata()

      ax1.fill_between(x_max, y_max, y_min, norm=1, where=np.less_equal(y_min, y_max), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)

      zg_theory_line_hist.Divide(zg_theory_line_hist)
      zg_theory_line_plot = convert_hist_to_line_plot(zg_theory_line_hist, n_bins)
      plt.plot(zg_theory_line_plot[0], zg_theory_line_plot[1], linewidth=5, color='red')

  else:
    raise ValueError("Only 'theory' or 'data' are valid options for calculating ratios!")
  

  # Normalized-Over-Data Plot Ends.

  ax0.set_xlabel("$z_g$", fontsize=95)
  ax0.set_ylabel("$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", fontsize=80, rotation=0, labelpad=115, y=0.39)
  
  ax1.set_xlabel("$z_g$", fontsize=95)
  
  if ratio_denominator == "data":
    label_pad = 135
  else:
    label_pad = 115

  plt.ylabel("Ratio           \nto           \n" + ratio_denominator.capitalize() + "           ", fontsize=55, rotation=0, labelpad=label_pad, y=0.31, axes=ax1)

  # Legend.

  th_line, = ax0.plot(range(1), linewidth=5, color='red')
  th_patch = mpatches.Patch(facecolor='pink', alpha=1.0, linewidth=0)

  if mc:
    pythia_line, = ax0.plot(range(1), linewidth=5, color='blue')
    herwig_line, = ax0.plot(range(1), linewidth=5, color='green')
  else:
    pythia_line, = ax0.plot(range(1), linewidth=5, color='blue', alpha=0)
    herwig_line, = ax0.plot(range(1), linewidth=5, color='green', alpha=0)

  handles = [data_plot, (th_patch, th_line), pythia_line, herwig_line]
  labels = [data_label, theory_label, pythia_label, herwig_label]

  first_legend = ax0.legend(handles, labels, fontsize=60, handler_map = {th_line : HandlerLine2D(marker_pad = 0), pythia_line : HandlerLine2D(marker_pad = 0), herwig_line : HandlerLine2D(marker_pad = 0)}, frameon=0, borderpad=0.1, bbox_to_anchor=[0.90, 0.98])
  ax = ax0.add_artist(first_legend)

  for txt in first_legend.get_texts():
    if ( not data) and txt.get_text() == data_label:
      txt.set_color("white") 

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra, extra]
  # labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5;~p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV};" + "\\abs{ \eta } < 3" + "}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5;~p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]
  ax0.legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.58])

  # Legend Ends.

  ax0.autoscale(True)
  ax1.autoscale(True)
  
  ax0.set_xlim(math.log(float(zg_cut), math.e), math.log(0.6, math.e))
  ax1.set_xlim(math.log(float(zg_cut), math.e), math.log(0.6, math.e))

  if data:
    ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.9249985), xycoords='figure fraction', frameon=0)
    plt.gca().add_artist(ab)
    preliminary_text = "Prelim. (20\%)"
    plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')
  else:
    preliminary_text = ""
    plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')


  x = np.linspace(math.log(0.05, math.e), math.log(0.5, math.e), 10)
  labels = [str(round(math.exp(i), 2)) for i in x]

  plt.sca(ax0)
  plt.xticks(x, labels)

  plt.sca(ax1)
  plt.xticks(x, labels)

  ax0.set_ylim(0, 1.6)
  
  plt.gcf().set_size_inches(30, 30, forward=1)

  ax0.xaxis.set_tick_params(width=5, length=20, labelsize=70)
  ax0.yaxis.set_tick_params(width=5, length=20, labelsize=70)

  ax1.xaxis.set_tick_params(width=5, length=20, labelsize=70)
  ax1.yaxis.set_tick_params(width=5, length=20, labelsize=70)
  

  plt.gcf().set_snap(True)
  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "Writing log_zg_cut_" + str(zg_filename) + "_pt_cut_" + str(pT_lower_cut) + "_ratio_over_" + ratio_denominator + "_th_" + str(theory) + "_mc_" + str(mc) + "_data_" + str(data) + ".pdf"
  filename = "plots/log_zg/log_zg_cut_" + str(zg_filename) + "_pt_cut_" + str(pT_lower_cut) + "_ratio_over_" + ratio_denominator + "_th_" + str(theory) + "_mc_" + str(mc) + "_data_" + str(data) + ".pdf"
  
  plt.savefig(filename)
  # plt.show()
  plt.clf()



def test2():
  properties = parse_file(input_analysis_file, 150)

  zgs = properties['zg_05'] 
  prescales = properties['prescales']

  zgs = np.linspace(0.05, 0.5, 50000)
  prescales = np.reciprocal(zgs)


  x, y = [], []
  for i in range(0, len(zgs)):
    if zgs[i] != 0:
      x.append(zgs[i])
      y.append(prescales[i])

  
  x_logged = np.log(x)

 


  ################################################ MAIN CODE BELOW THIS LINE ###################################################

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  bins_log = np.logspace(math.log(0.05, math.e), math.log(0.5, math.e), 10, base=math.e)
  # bins_linear = np.linspace(0.05, 0.5, 5)
  bins_linear_log = np.linspace(math.log(0.05, math.e), math.log(0.5, math.e), 50)

  
  # plt.hist(x, weights=y, bins=bins_log, histtype='step', lw=5, normed=True)
  # plt.hist(x, weights=y, bins=bins_linear, histtype='step', lw=5, normed=True)
  # plt.hist(x, weights=y, bins=bins_linear_log, histtype='step', lw=5, normed=True)

  plt.hist(x_logged, weights=y, bins=bins_linear_log, histtype='step', lw=5, normed=True)

  

  # This is an attempt to use numpy binning before using pyplot to plot.
  # But that was attempted because of a misunderstanding that numpy would deal with "variable bins" correctly.
  # But in actuality, we don't have to deal with variable bins for log(zg) anyway.

  # hist, bins = np.histogram(x_logged, weights=y, bins=bins_linear_log, normed=True)

  # plt.hist(bins[:-1], weights=hist, histtype='step', bins=50, lw=5)

  # plt.errorbar(bins[:-1], hist, xerr=True, yerr=True, marker="o", markersize=10, lw=0)
  



  plt.show()

  plt.clf()


  '''

  my_hist = Hist([0.05, 0.1, 0.2, 0.5], title="something", markersize=3, color='red')

  my_hist_corrected = Hist([0.05, 0.1, 0.2, 0.5], title="something", markersize=3, color='blue')
  my_hist_crazy_binning = Hist([0.05, 0.1, 0.2, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5], title="something", markersize=3, color='green')
  my_hist_uniform_binning = Hist([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5], title="something", markersize=3, color='magenta')

  map(my_hist.Fill, x, prescales)
  map(my_hist_corrected.Fill, x, prescales)
  
  map(my_hist_crazy_binning.Fill, x, prescales)
  map(my_hist_uniform_binning.Fill, x, prescales)


  def norm_hist(hist):
    bin_width = (hist.upperbound() - hist.lowerbound()) / hist.nbins()
    hist.Scale(1.0 / ( hist.GetSumOfWeights() ))

    for i in range(0, hist.GetSize()):
      bin_width = hist.GetXaxis().GetBinWidth(i)
      old_bin_height =  hist.GetBinContent(i)
      
      new_height = old_bin_height / bin_width
      hist.SetBinContent(i, new_height)

      old_error = hist.GetBinError(i)
      new_error = old_error / bin_width
      hist.SetBinError(i, new_error)

    return hist

  

  bin_width = (my_hist.upperbound() - my_hist.lowerbound()) / my_hist.nbins()
  my_hist.Scale(1.0 / ( my_hist.GetSumOfWeights() * bin_width ))

  my_hist_corrected = norm_hist(my_hist_corrected)
  my_hist_crazy_binning = norm_hist(my_hist_crazy_binning)
  my_hist_uniform_binning = norm_hist(my_hist_uniform_binning)

  # rplt.errorbar(my_hist, xerr=1, yerr=1, emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=0.5)
  rplt.errorbar(my_hist_corrected, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
  rplt.errorbar(my_hist_crazy_binning, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
  rplt.errorbar(my_hist_uniform_binning, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)



  plt.legend()


  plt.locator_params(nbins=5)

  plt.draw()  # This is important, because without this, pyplot won't "draw" the ticks, thereby making the xticks list empty.
  xticks = [ tick.get_text().replace("$", "") for tick in plt.gca().get_xticklabels()]  
  # xticks = [math.log(0.01, math.e), math.log(0.05, math.e), math.log(0.1, math.e), math.log(0.2, math.e), math.log(0.5, math.e), math.log(1, math.e)]

  # print xticks
  # v = [round(math.e**float(i), 2) if i != "" else "" for i in xticks]
  # print v
  v = [i if i != "" else "" for i in xticks]
  plt.gca().xaxis.set_ticklabels( v )


  # plt.gca().set_xscale('log')

  # plt.gca().tick_params(axis='x', which='minor', bottom='off')
  # plt.minorticks_on()

  # minorLocator   = AutoMinorLocator()
  minorLocator = MultipleLocator(0.01)
  plt.gca().xaxis.set_minor_locator(minorLocator)


  plt.tick_params(which='both', width=5)
  plt.tick_params(which='major', length=30)
  plt.tick_params(which='minor', length=10)



  plt.gca().autoscale(True)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)


  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)


  # plt.savefig("test_zg.pdf")


  plt.show()  

  '''


# test2()

def plot_pts_variable_bin():
  pT_lower_cut = 150
  properties = parse_file(input_analysis_file, pT_lower_cut)

  pTs = properties['uncorrected_hardest_pts']
  corrected_pTs = properties['corrected_hardest_pts']
  prescales = properties['prescales']

  event_numbers = properties['event_number']
  run_numbers = properties['run_number']

  print max(pTs)

  # for i in range(0, len(pTs)):
    # if pTs[i] > 10000:
  #     print int(event_numbers[i]), int(run_numbers[i])


  herwig_pTs = parse_mc_pt_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_pt_herwig_pp2jj_150pTcut_7TeV.dat")
  pythia_pTs = parse_mc_pt_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_pt_pythia_pp2jj_150pTcut_7TeV.dat")

  bins = [80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500, 600, 700, 800, 900, 1000, 1250, 1500, 2000]

  pythia_pt_hist = Hist(bins, title="Pythia 8.205", linewidth=5, markersize=5.0, color="blue")
  herwig_pt_hist = Hist(bins, title="Herwig++ 2.6.3", linewidth=5, markersize=5.0, color="green")
  corrected_pt_hist = Hist(bins, title='Corrected', markersize=3.0, color='black')
  uncorrected_pt_hist = Hist(bins, title='Uncorrected', markersize=3.0, color='orange')

  map(uncorrected_pt_hist.Fill, pTs, prescales)
  map(corrected_pt_hist.Fill, corrected_pTs, prescales)
  
  map(pythia_pt_hist.Fill, pythia_pTs)
  map(herwig_pt_hist.Fill, herwig_pTs)


  corrected_pt_hist = normalize_hist(corrected_pt_hist)
  uncorrected_pt_hist = normalize_hist(uncorrected_pt_hist)
  pythia_pt_hist = normalize_hist(pythia_pt_hist)
  herwig_pt_hist = normalize_hist(herwig_pt_hist)

  
  gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 

  ax0 = plt.subplot(gs[0])
  ax1 = plt.subplot(gs[1])


  data_plot = rplt.errorbar(corrected_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  uncorrected_data_plot = rplt.errorbar(uncorrected_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.hist(pythia_pt_hist, axes=ax0, emptybins=False, marker='o',  markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.hist(herwig_pt_hist, axes=ax0, emptybins=False, marker='o',  markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  

  data_x_errors, data_y_errors = [], []
  for x_segment in data_plot[2][0].get_segments():
    data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
  for y_segment in data_plot[2][1].get_segments():
    data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

  data_points_x = data_plot[0].get_xdata()
  data_points_y = data_plot[0].get_ydata()

  data_plot_points_x = []
  data_plot_points_y = []
  for i in range(0, len(data_points_x)):
    data_plot_points_x.append(data_points_x[i])
    data_plot_points_y.append(data_points_y[i])



  uncorrected_data_x_errors, uncorrected_data_y_errors = [], []
  for x_segment in uncorrected_data_plot[2][0].get_segments():
    uncorrected_data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
  for y_segment in uncorrected_data_plot[2][1].get_segments():
    uncorrected_data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

  uncorrected_data_points_x = uncorrected_data_plot[0].get_xdata()
  uncorrected_data_points_y = uncorrected_data_plot[0].get_ydata()

  uncorrected_data_plot_points_x = []
  uncorrected_data_plot_points_y = []
  for i in range(0, len(uncorrected_data_points_x)):
    uncorrected_data_plot_points_x.append(uncorrected_data_points_x[i])
    uncorrected_data_plot_points_y.append(uncorrected_data_points_y[i])




  data_to_data_y_err = [(b / m) for b, m in zip(data_y_errors, data_plot_points_y)]
  data_to_data_x_err = [(b / m) for b, m in zip(data_x_errors, [1] * len(data_plot_points_y))]

  uncorrected_to_corrected_y_err = [(b / m) for b, m in zip(uncorrected_data_y_errors, data_plot_points_y)]
  uncorrected_to_corrected_x_err = [(b / m) for b, m in zip(uncorrected_data_x_errors, [1] * len(data_plot_points_y))]


  ax0.autoscale(True)
  ax0.set_yscale('log')


  # Legends Begin.

  legend = ax0.legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[0.98, 1.0])
  ax0.add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<3$", r"$p_{T} > 150~\mathrm{GeV}$"]
  ax0.legend([extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[1.0, 0.62])

  # Legends End.



  ax0.set_xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75)
  ax1.set_xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75)
  ax0.set_ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)
  ax1.set_ylabel("Ratio           \nto           \n" + "Data" + "           ", fontsize=55, rotation=0, labelpad=115, y=0.31)


  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.9249985), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  # Ratio Plot.
  pythia_pt_hist.Divide(corrected_pt_hist)
  herwig_pt_hist.Divide(corrected_pt_hist)
  uncorrected_pt_hist.Divide(corrected_pt_hist)
  corrected_pt_hist.Divide(corrected_pt_hist)

  rplt.hist(pythia_pt_hist, axes=ax1, linewidth=5)
  rplt.hist(herwig_pt_hist, axes=ax1, linewidth=5)
  
  rplt.errorbar(corrected_pt_hist, xerr=data_to_data_x_err, yerr=data_to_data_y_err, axes=ax1, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.errorbar(uncorrected_pt_hist, xerr=uncorrected_to_corrected_x_err, yerr=uncorrected_to_corrected_y_err, axes=ax1, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  ax1.autoscale(True)
  
  ax0.set_ylim(10e-8, 10e-1)
  ax1.set_ylim(0., 2.)


  plt.gcf().set_size_inches(30, 30, forward=1)


  plt.sca(ax0)
  plt.gca().xaxis.set_minor_locator(MultipleLocator(100))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)
  
  plt.sca(ax1)
  plt.gca().xaxis.set_minor_locator(MultipleLocator(100))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)


  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/pT_distribution_var_bin.pdf")
  # plt.show()
  plt.clf()




def test3():
  properties = parse_file(input_analysis_file, 150)

  pTs = properties['corrected_hardest_pts'] 
  prescales = properties['prescales']

  print max(properties['uncorrected_hardest_pts'])
  print max(properties['corrected_hardest_pts'])


  x, y = [], []
  for i in range(0, len(pTs)):
    if pTs[i] != 0:
      x.append(pTs[i])
      y.append(prescales[i])

  
  x_logged = np.log(x)

  
  plt.hist(x, weights=y, bins=50, histtype='step', lw=5, normed=1)
  plt.savefig("test.pdf")

  plt.clf()

  plt.hist(x_logged, weights=y, bins=50, histtype='step', lw=5, normed=1)
  plt.savefig("test_log.pdf")

  plt.clf()

  # plt.hist(x, weights=y, bins=50, histtype='step', lw=5)

  # plt.draw()  # This is important, because without this, pyplot won't "draw" the ticks, thereby making the xticks list empty.
  # xticks = [ tick.get_text().replace("$", "") for tick in plt.gca().get_xticklabels() ]

  plt.clf()


  ################################################ MAIN CODE BELOW THIS LINE ###################################################

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)



  x_hist = Hist(50, math.log(150, math.e), math.log(10000, math.e), title='CMS Open Data', markersize=3.0, color='black')
  # x_hist = Hist(50, 150, 10000, title='CMS Open Data', markersize=5.0, color='black')
  
  bin_width = (x_hist.upperbound() - x_hist.lowerbound()) / x_hist.nbins()
  
  map(x_hist.Fill, x_logged, y)
  # map(x_hist.Fill, x, y)

  x_hist.Scale(1.0 / (sum([abs(i) for i in x_logged]) * bin_width))
  

  # plt.x_hist(x_logged, weights=y, label="CMS Open Data", bins=50, histtype='step', lw=5, normed=1)
  rplt.errorbar(x_hist, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  plt.locator_params(nbins=5)

  plt.draw()  # This is important, because without this, pyplot won't "draw" the ticks, thereby making the xticks list empty.
  xticks = [ tick.get_text().replace("$", "") for tick in plt.gca().get_xticklabels()]  
  xticks

  v = [round(math.e**float(i), 2) if i != "" else "" for i in xticks]
  # v = [i if i != "" else "" for i in xticks]
  plt.gca().xaxis.set_ticklabels( v )




  # Irrelevant stuff begin.

  legend = plt.gca().legend(loc=1, frameon=0, fontsize=60)
  plt.gca().add_artist(legend)

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5;~p_{T} > " + str(150) + "~\mathrm{GeV}}$"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.58])

  # Legends End.

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')



  
  plt.xlabel('$p_T~(GeV)$', fontsize=75)
  plt.ylabel('A.U.', fontsize=75, rotation=0, labelpad=50.)
  
  

  plt.gca().autoscale(True)
  # plt.gca().set_ylim(-0.001, 0.015)
  # plt.gca().set_ylim(-0.1, 0.015)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)


  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)


  plt.savefig("test_matplotlib.pdf")











# plot_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=250, zg_cut='0.05', zg_filename='zg_05', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=18, y_limit_ratio_plot=0.5)



# Version 3 Begins Here.

# plot_jet_mass_spectrum()
# plot_jet_mass_spectrum(pT_lower_cut=100, pT_upper_cut=200)
# plot_jet_mass_spectrum(pT_lower_cut=200, pT_upper_cut=400)
# plot_jet_mass_spectrum(pT_lower_cut=400)


# plot_fractional_energy_loss()
# plot_fractional_energy_loss(pT_lower_cut=100, pT_upper_cut=200)
# plot_fractional_energy_loss(pT_lower_cut=200, pT_upper_cut=400)
# plot_fractional_energy_loss(pT_lower_cut=400)

# plot_constituent_multiplicity_softdrop()
# plot_constituent_multiplicity_softdrop(pT_lower_cut=100, pT_upper_cut=200)
# plot_constituent_multiplicity_softdrop(pT_lower_cut=200, pT_upper_cut=400)
# plot_constituent_multiplicity_softdrop(pT_lower_cut=400)

# plot_jet_area()

# plot_hardest_pt_softdrop()
# plot_hardest_pt_softdrop(pT_lower_cut=100, pT_upper_cut=200)
# plot_hardest_pt_softdrop(pT_lower_cut=200, pT_upper_cut=400)
# plot_hardest_pt_softdrop(pT_lower_cut=400)


# plot_pts()
# plot_pts(pT_lower_cut=150, pT_upper_cut=250)
# plot_pts(pT_lower_cut=250, pT_upper_cut=500)
# plot_pts(pT_lower_cut=500)


# plot_pts_variable_bin()



# plot_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=18, y_limit_ratio_plot=0.5)

# plot_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=0, data=0, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)
# plot_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=1, data=0, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)
# plot_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)
# plot_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='data', theory=1, mc=1, data=1, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)

# plot_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)



# plot_zg_th_mc_data(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=4, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_zg_th_mc_data(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=4, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_zg_th_mc_data(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=4, y_max_limit=15, y_limit_ratio_plot=1.0)

# plot_zg_th_mc_data(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=2, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_zg_th_mc_data(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=2, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_zg_th_mc_data(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=2, y_max_limit=15, y_limit_ratio_plot=1.0)






# plot_log_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=18, y_limit_ratio_plot=0.5)

# plot_log_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=0, data=0, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)
# plot_log_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)
# plot_log_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=1, data=0, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)
plot_log_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='data', theory=1, mc=1, data=1, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)

# plot_log_zg_th_mc_data(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)


# plot_log_zg_th_mc_data(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=4, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_log_zg_th_mc_data(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=4, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_log_zg_th_mc_data(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=4, y_max_limit=15, y_limit_ratio_plot=1.0)

# plot_log_zg_th_mc_data(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=2, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_log_zg_th_mc_data(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=2, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_log_zg_th_mc_data(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', ratio_denominator='theory', theory=1, mc=1, data=1, n_bins=2, y_max_limit=15, y_limit_ratio_plot=1.0)







# plot_charged_and_all_zgs(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', n_bins=8, y_max_limit=14)
# plot_charged_and_all_zgs(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', n_bins=8, y_max_limit=10)
# plot_charged_and_all_zgs(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', n_bins=8, y_max_limit=10)

# plot_charged_and_all_zgs(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', n_bins=4, y_max_limit=15)
# plot_charged_and_all_zgs(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', n_bins=4, y_max_limit=15)
# plot_charged_and_all_zgs(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', n_bins=4, y_max_limit=15)

# plot_charged_and_all_zgs(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', n_bins=2, y_max_limit=15)
# plot_charged_and_all_zgs(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', n_bins=2, y_max_limit=15)
# plot_charged_and_all_zgs(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', n_bins=2, y_max_limit=15)



# plot_zg_pfc_pt_cut(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', n_bins=8, y_max_limit=18)
# plot_zg_pfc_pt_cut(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', n_bins=8, y_max_limit=10)
# plot_zg_pfc_pt_cut(pT_lower_cut=150, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', n_bins=8, y_max_limit=10)

# plot_zg_pfc_pt_cut(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', n_bins=4, y_max_limit=15)
# plot_zg_pfc_pt_cut(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', n_bins=4, y_max_limit=15)
# plot_zg_pfc_pt_cut(pT_lower_cut=300, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', n_bins=4, y_max_limit=15)

# plot_zg_pfc_pt_cut(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.05', zg_filename='zg_05', n_bins=2, y_max_limit=15)
# plot_zg_pfc_pt_cut(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.1', zg_filename='zg_1', n_bins=2, y_max_limit=15)
# plot_zg_pfc_pt_cut(pT_lower_cut=600, pT_upper_cut=10000, zg_cut='0.2', zg_filename='zg_2', n_bins=2, y_max_limit=15)



# plot_JEC()




# Version 3 Ends Here.





















# plot_fractional_energy_loss()


# plot_JEC()

# plot_jet_area()



# plot_constituent_multiplicity_softdrop()

# plot_hardest_pt_softdrop()



# plot_hardest_pt_corresponding_triggers()




# plot_2d()



# plot_2d_hist()



# # Trigger Efficiency Curves Begin.


# plot_trigger_efficiency_curves("HLT_Jet30U", "HLT_Jet15U", pT_upper_limit=200)
# plot_trigger_efficiency_curves("HLT_Jet50U", "HLT_Jet30U", pT_upper_limit=300)
# plot_trigger_efficiency_curves("HLT_Jet70U", "HLT_Jet50U", pT_upper_limit=350)
# plot_trigger_efficiency_curves("HLT_Jet100U", "HLT_Jet70U", pT_upper_limit=800)
# plot_trigger_efficiency_curves("HLT_Jet140U", "HLT_Jet100U", pT_upper_limit=800)
# plot_trigger_efficiency_curves("HLT_Jet180U", "HLT_Jet140U", pT_upper_limit=1200)


# plot_all_trigger_efficiency_curves()

# # Trigger Efficiency Curves End.








# Triggers Turn-On Curve Begins.

# plot_turn_on_curves()

# Triggers Turn-On Curve Ends.








# # AK5 Distribution Begins.

# plot_pts()

# plot_pts_variable_bin()

# # AK5 Distribution Ends.











# # zg_distribution Begins.

# plot_zg_th_mc_data(150, '0.05', 'zg_05', 'theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=18, y_limit_ratio_plot=0.5)

# plot_zg_th_mc_data(150, '0.1', 'zg_1', 'theory', theory=1, mc=0, data=0, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)
# plot_zg_th_mc_data(150, '0.1', 'zg_1', 'theory', theory=1, mc=1, data=0, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)
# plot_zg_th_mc_data(150, '0.1', 'zg_1', 'theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)
# plot_zg_th_mc_data(150, '0.1', 'zg_1', 'data', theory=1, mc=1, data=1, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)


# plot_zg_th_mc_data(150, '0.2', 'zg_2', 'theory', theory=1, mc=1, data=1, n_bins=8, y_max_limit=10, y_limit_ratio_plot=0.5)

# plot_zg_th_mc_data(300, '0.05', 'zg_05', 'theory', theory=1, mc=1, data=1, n_bins=4, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_zg_th_mc_data(300, '0.1', 'zg_1', 'theory', theory=1, mc=1, data=1, n_bins=4, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_zg_th_mc_data(300, '0.2', 'zg_2', 'theory', theory=1, mc=1, data=1, n_bins=4, y_max_limit=15, y_limit_ratio_plot=1.0)

# plot_zg_th_mc_data(600, '0.05', 'zg_05', 'theory', theory=1, mc=1, data=1, n_bins=2, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_zg_th_mc_data(600, '0.1', 'zg_1', 'theory', theory=1, mc=1, data=1, n_bins=2, y_max_limit=15, y_limit_ratio_plot=1.0)
# plot_zg_th_mc_data(600, '0.2', 'zg_2', 'theory', theory=1, mc=1, data=1, n_bins=2, y_max_limit=15, y_limit_ratio_plot=1.0)

# # zg_distribution Ends.








# # Charged zg Begins.

# plot_charged_and_all_zgs(150, '0.05', 'zg_05', n_bins=8, y_max_limit=18)
# plot_charged_and_all_zgs(150, '0.1', 'zg_1', n_bins=8, y_max_limit=10)
# plot_charged_and_all_zgs(150, '0.2', 'zg_2', n_bins=8, y_max_limit=10)

# plot_charged_and_all_zgs(300, '0.05', 'zg_05', n_bins=4, y_max_limit=15)
# plot_charged_and_all_zgs(300, '0.1', 'zg_1', n_bins=4, y_max_limit=15)
# plot_charged_and_all_zgs(300, '0.2', 'zg_2', n_bins=4, y_max_limit=15)

# plot_charged_and_all_zgs(600, '0.05', 'zg_05', n_bins=2, y_max_limit=15)
# plot_charged_and_all_zgs(600, '0.1', 'zg_1', n_bins=2, y_max_limit=15)
# plot_charged_and_all_zgs(600, '0.2', 'zg_2', n_bins=2, y_max_limit=15)

# # Charged zg Ends.











# # zg with PFC pT_cut Begins.


# plot_zg_pfc_pt_cut(150, '0.05', 'zg_05', n_bins=8, y_max_limit=18)
# plot_zg_pfc_pt_cut(150, '0.1', 'zg_1', n_bins=8, y_max_limit=10)
# plot_zg_pfc_pt_cut(150, '0.2', 'zg_2', n_bins=8, y_max_limit=10)

# plot_zg_pfc_pt_cut(300, '0.05', 'zg_05', n_bins=4, y_max_limit=15)
# plot_zg_pfc_pt_cut(300, '0.1', 'zg_1', n_bins=4, y_max_limit=15)
# plot_zg_pfc_pt_cut(300, '0.2', 'zg_2', n_bins=4, y_max_limit=15)

# plot_zg_pfc_pt_cut(600, '0.05', 'zg_05', n_bins=2, y_max_limit=15)
# plot_zg_pfc_pt_cut(600, '0.1', 'zg_1', n_bins=2, y_max_limit=15)
# plot_zg_pfc_pt_cut(600, '0.2', 'zg_2', n_bins=2, y_max_limit=15)



# # zg with PFC pT_cut Ends.


