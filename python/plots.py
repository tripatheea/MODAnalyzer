# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat
from __future__ import division

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

# RootPy
from rootpy.plotting import Hist, HistStack, Legend
import rootpy.plotting.root2matplotlib as rplt

# Stuff for calculating areas.
from scipy.integrate import simps
from scipy import interpolate
from scipy import optimize

from numpy import trapz


from matplotlib import gridspec

from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png

from matplotlib.offsetbox import OffsetImage, AnnotationBbox



from scipy.stats import binned_statistic

import rootpy.plotting.views

input_analysis_file = sys.argv[1]


mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

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
        if (float(numbers[4]) > pT_lower_cut) and (float(numbers[17]) > pfc_pT_cut):
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
      if len(numbers) != 0:
        # print "Some kind of error occured while parsing the given file!"
        # print numbers
        # print
        pass

  return properties




def parse_file_turn_on(input_file, pT_lower_cut = 0.00):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        if (float(numbers[1]) > pT_lower_cut):
          properties['corrected_hardest_pts'].append( float( numbers[1] ) )
          properties['prescales'].append( int( numbers[2] ) )
          properties['trigger_names'].append(  numbers[3] )

    except:
      if len(numbers) != 0:
        print "Some kind of error occured while parsing the given file!"
        print numbers
        print


  return properties



def parse_mc_pt_file(input_file, pT_lower_cut = 0.00):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  pTs = []

  for line in lines:
    try:
      numbers = line.split()
      pTs.append( float( numbers[0] ) )
    except:
      if len(numbers) != 0:
        # print "Some kind of error occured while parsing the given file!"
        # print numbers
        # print
        pass
  
  return pTs

    
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
    
    rplt.errorbar(dr_hist, xerr=False, emptybins=False)

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
    
    rplt.errorbar(mu_hist, xerr=False, emptybins=False)

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

  rplt.errorbar(pdgid_hist, xerr=False, emptybins=False)

  plt.yscale('log')

  plt.autoscale(True)

  plt.xlabel('$p_{T}$ GeV')
  plt.suptitle("Hardest Charged pT Distribution")

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

  plt.savefig("plots/hardest_charged_pt_distribution.pdf")
  plt.show()


def plot_zg_pfc_pt_cut(pT_lower_cut, zg_cut, zg_filename, n_bins=10, y_max_limit=20):

  properties = parse_file(input_analysis_file, pT_lower_cut)

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
  
  rplt.errorbar(zg_hist, emptybins=False, ls='None', marker='o', markersize=8, capthick=3, pickradius=3, elinewidth=3)
  rplt.errorbar(zg_pfc_cut_2_hist, emptybins=False, ls='None', marker='o', markersize=8, capthick=3, pickradius=3, elinewidth=3)
  rplt.errorbar(zg_pfc_cut_5_hist, emptybins=False, ls='None', marker='o', markersize=8, capthick=3, pickradius=3, elinewidth=3)
  
  plt.autoscale(True)
  plt.gca().set_ylim(0, y_max_limit)

  legend = plt.legend(frameon=0, fontsize=60, bbox_to_anchor=[0.935, 0.98])
  plt.gca().add_artist(legend)

  
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra, extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5;~p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.60])


  plt.xlabel("$z_g$", fontsize=95)
  plt.ylabel("$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", fontsize=80, rotation=0, labelpad=115, y=0.39)
  

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gca().set_ylim(0, y_max_limit)
  plt.gca().set_xlim(0.0, 0.6)

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/soft_cut/" + str(zg_filename) + "_pt_cut_" + str(pT_lower_cut) + ".pdf")
  # plt.show()
  plt.clf()


def plot_charged_and_all_zgs(pT_lower_cut, zg_cut, zg_filename, n_bins=10, y_max_limit=20):

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
  
  rplt.errorbar(zg_hist, emptybins=False, ls='None', marker='o', markersize=8, capthick=3, pickradius=3, elinewidth=3)
  rplt.errorbar(zg_charged_hist, emptybins=False, ls='None', marker='o', markersize=8, capthick=3, pickradius=3, elinewidth=3)
  
  plt.autoscale(True)
  plt.gca().set_ylim(0, y_max_limit)

  legend = plt.legend(frameon=0, fontsize=60, bbox_to_anchor=[0.86, 0.98])
  plt.gca().add_artist(legend)

  
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra, extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5;~p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.70])


  plt.xlabel("$z_g$", fontsize=95)
  plt.ylabel("$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", fontsize=80, rotation=0, labelpad=115, y=0.39)
  

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gca().set_ylim(0, y_max_limit)
  plt.gca().set_xlim(0.0, 0.6)

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/charged_zg/charged_" + str(zg_filename) + "_pt_cut_" + str(pT_lower_cut) + ".pdf")
  # plt.show()
  plt.clf()




def plot_turn_on_curves():
  # properties = parse_file_turn_on(input_analysis_file)
  
  properties = parse_file_turn_on('./turn_on.dat')

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']

  expected_trigger_names = ["HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U", "HLT\_Jet15U\_HcalNoiseFiltered" ]
  labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U", "Jet15U\_HNF" ]

  colors = ['red', 'blue', 'green', 'magenta', 'black', 'orange']
  colors = colors[::-1]

  pt_hists = []
  for i in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist(60, 0, 300, title=labels[i], markersize=1.0, color=colors[i], linewidth=5))

  for i in range(0, len(pTs)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j].replace("\\", "") in trigger_names[i]:
        pt_hists[j].Fill(pTs[i], prescales[i])

  for k in range(0, len(pt_hists)):
    rplt.errorbar(pt_hists[k])

  plt.autoscale(True)
  plt.yscale('log')

  plt.gca().set_ylim(1, 10e10)

  
  # Legend.
  handles, labels = plt.gca().get_legend_handles_labels()
  first_legend = plt.gca().legend(handles[::-1], labels[::-1], fontsize=60,  frameon=0, borderpad=0.1, bbox_to_anchor=[0.97, 0.98])
  ax = plt.gca().add_artist(first_legend)

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5} $"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.93, 0.50])


  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=95)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)

  fn = get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.885, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)


  plt.savefig("plots/turn_on_curves.pdf")
  # plt.show()
  plt.clf()

def plot_pts():
  properties = parse_file(input_analysis_file, 150)

  pTs = properties['uncorrected_hardest_pts']
  corrected_pTs = properties['corrected_hardest_pts']
  prescales = properties['prescales']

  herwig_pTs = parse_mc_pt_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_pt_herwig_pp2jj_150pTcut_7TeV.dat")
  pythia_pTs = parse_mc_pt_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_pt_pythia_pp2jj_150pTcut_7TeV.dat")

  pythia_pt_hist = Hist(100, 0, 1000, title="Pythia 8.205", linewidth=5, markersize=5.0, color="blue")
  bin_width_pythia = (pythia_pt_hist.upperbound() - pythia_pt_hist.lowerbound()) / pythia_pt_hist.nbins()

  herwig_pt_hist = Hist(100, 0, 1000, title="Herwig++ 2.6.3", linewidth=5, markersize=5.0, color="green")
  bin_width_herwig = (herwig_pt_hist.upperbound() - herwig_pt_hist.lowerbound()) / herwig_pt_hist.nbins()

  corrected_pt_hist = Hist(100, 0, 1000, title='Corrected', markersize=3.0, color='black')
  bin_width_corrected = (corrected_pt_hist.upperbound() - corrected_pt_hist.lowerbound()) / corrected_pt_hist.nbins()

  uncorrected_pt_hist = Hist(100, 0, 1000, title='Uncorrected', markersize=3.0, color='orange')
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

  data_plot = rplt.errorbar(corrected_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=8, capthick=3, linewidth=3, pickradius=5, elinewidth=3)
  rplt.errorbar(uncorrected_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=8, capthick=3, linewidth=3, pickradius=5, elinewidth=3)
  rplt.hist(pythia_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=8, capthick=3, linewidth=5, pickradius=5, elinewidth=3)
  rplt.hist(herwig_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=8, capthick=3, linewidth=5, pickradius=5, elinewidth=3)
  

  data_x_errors, data_y_errors = [], []
  for x_segment in data_plot[2][0].get_segments():
    data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
  for y_segment in data_plot[2][1].get_segments():
    data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

  data_points_x = data_plot[0].get_xdata()
  data_points_y = data_plot[0].get_ydata()



  ax0.autoscale(True)
  ax0.set_yscale('log')


  # Legends Begin.

  legend = ax0.legend(loc=1, frameon=0, fontsize=60)
  ax0.add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5}$"]
  ax0.legend([extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.91, 0.62])

  # Legends End.




  ax0.set_xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75)
  ax1.set_xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75)
  ax0.set_ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.9249985), xycoords='figure fraction', frameon=0)
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
  rplt.errorbar(corrected_pt_hist, axes=ax1, emptybins=False, marker='o', markersize=8, capthick=3, linewidth=5, pickradius=5, elinewidth=3)
  rplt.errorbar(uncorrected_pt_hist, axes=ax1, emptybins=False, marker='o', markersize=8, capthick=3, linewidth=5, pickradius=5, elinewidth=3)

  ax1.autoscale(True)
  
  ax0.set_ylim(10e-8, 10e-1)
  ax1.set_ylim(0., 2.)


  plt.gcf().set_size_inches(30, 30, forward=1)


  ax0.xaxis.set_tick_params(width=5, length=20, labelsize=70)
  ax0.yaxis.set_tick_params(width=5, length=20, labelsize=70)

  ax1.xaxis.set_tick_params(width=5, length=20, labelsize=70)
  ax1.yaxis.set_tick_params(width=5, length=20, labelsize=70)

  plt.ylabel("Ratio           \nto           \n" + "Data" + "           ", fontsize=55, rotation=0, labelpad=135, y=0.31, axes=ax1)


  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/ak5_pt_distribution.pdf")
  # plt.show()
  plt.clf()


def plot_hardest_pt_corresponding_triggers():
  properties = parse_file(input_analysis_file)

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']

  expected_trigger_names = ["HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U", "HLT\_Jet15U"]

  colors = ['red', 'blue', 'pink', 'green']

  pt_hists = []
  for i in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist(50, 0, 1000, title=expected_trigger_names[i], markersize=1.0, color=colors[i], linewidth=5))

  for i in range(0, len(pTs)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j].replace("\\", "") in trigger_names[i]:
        pt_hists[j].Fill(pTs[i], prescales[i])

  for k in range(0, len(pt_hists)):
    rplt.errorbar(pt_hists[k], marker='o', markersize=8, linewidth=3, capthick=3, pickradius=5, elinewidth=3)

  plt.yscale('log')
  plt.autoscale(True)

  plt.gca().set_ylim(1, 10e8)

  plt.legend(loc=0, frameon=0)

  plt.xlabel('$p_T \mathrm{(GeV)}$')

  fn = get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)
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

  
  H, xedges, yedges = np.histogram2d(zgs[0], charged_zgs[0], normed=1, range=[[0, 0.5], [0, 0.5]], weights=prescales, bins=100)
   
  # Mask zeros
  Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
   
  # Plot 2D histogram using pcolor
  
  plt.pcolormesh(xedges,yedges,Hmasked)
  plt.xlabel('Charged zg_05')
  plt.ylabel('zg_05')
  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Counts')

  plt.savefig("plots/zg_05_vs_charged_zg_05.pdf")
  plt.show()

  H, xedges, yedges = np.histogram2d(zgs[1], charged_zgs[1], normed=1, range=[[0, 0.5], [0, 0.5]], weights=prescales, bins=100)
   
  # Mask zeros
  Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
   
  # Plot 2D histogram using pcolor
  
  plt.pcolormesh(xedges,yedges,Hmasked)
  plt.xlabel('Charged zg_1')
  plt.ylabel('zg_1')
  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Counts')

  plt.savefig("plots/zg_1_vs_charged_zg_1.pdf")
  plt.show()

  H, xedges, yedges = np.histogram2d(zgs[2], charged_zgs[2], normed=1, range=[[0, 0.5], [0, 0.5]], weights=prescales, bins=100)
   
  # Mask zeros
  Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
   
  # Plot 2D histogram using pcolor
  
  plt.pcolormesh(xedges,yedges,Hmasked)
  plt.xlabel('Charged zg_2')
  plt.ylabel('zg_2')
  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Counts')

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

  plt.savefig("plots/zg_2_vs_charged_zg_2.pdf")
  plt.show()







def plot_zg_th_mc_data(pT_lower_cut, zg_cut, zg_filename, ratio_denominator="theory", data=True, mc=True, theory=True, n_bins=10, y_max_limit=20, y_limit_ratio_plot=0.5):
  pfc_pT_cut = 0

  zg_cut = float(zg_cut)

  properties = parse_file(input_analysis_file, pT_lower_cut, pfc_pT_cut)
  properties_pythia = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_pythia_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)
  properties_herwig = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_herwig_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)

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

  norm_data_prescales = map(lambda x: x / ( zg_data_hist.GetSumOfWeights() * bin_width_data ), prescales)
  
  if data:
    # data_plot, caplines, barlinecols
    data_plot = rplt.errorbar(zg_data_hist, xerr=1, yerr=1, emptybins=False, axes=ax0, ls='None', marker='o', markersize=8, pickradius=3, capthick=3, elinewidth=3, alpha=1.0)
  else:
    data_plot = rplt.errorbar(zg_data_hist, xerr=1, yerr=1, emptybins=False, axes=ax0, ls='None', marker='o', markersize=8, pickradius=3, capthick=3, elinewidth=3, alpha=0.0)


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
      
      plt.errorbar(data_plot_points_x, ratio_data_to_data, xerr=data_to_data_x_err, yerr=data_to_data_y_err, ls='None', marker='o', markersize=8, linewidth=3, capthick=3, pickradius=5, elinewidth=3, color='black')

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

      plt.errorbar(data_plot_points_x, zg_data_to_th_y, xerr=data_to_th_x_err, yerr=zg_data_to_th_y_err, ls='None', capthick=3, marker='o', markersize=8, linewidth=3, pickradius=5, elinewidth=3, color='black')
   
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

  first_legend = ax0.legend(handles, labels, fontsize=60, handler_map = {th_line : HandlerLine2D(marker_pad = 0), pythia_line : HandlerLine2D(marker_pad = 0), herwig_line : HandlerLine2D(marker_pad = 0)}, frameon=0, borderpad=0.1, bbox_to_anchor=[0.90, 0.98])
  ax = ax0.add_artist(first_legend)

  for txt in first_legend.get_texts():
    if ( not data) and txt.get_text() == data_label:
      txt.set_color("white") 

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra, extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5;~p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}}$", r"$ \textrm{Soft~Drop:}~\boldsymbol{\beta = 0;~z_{\mathrm{cut}} = " + str(zg_cut) + "}$"]
  # labels = ["", ""]
  ax0.legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.58])


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
    ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.9249985), xycoords='figure fraction', frameon=0)
    plt.gca().add_artist(ab)
    preliminary_text = "Prelim. (20\%)"
    plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')
  else:
    preliminary_text = ""
    plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')



  fig = plt.gcf()
  fig.set_size_inches(30, 30, forward=1)

  ax0.xaxis.set_tick_params(width=5, length=20, labelsize=70)
  ax0.yaxis.set_tick_params(width=5, length=20, labelsize=70)

  ax1.xaxis.set_tick_params(width=5, length=20, labelsize=70)
  ax1.yaxis.set_tick_params(width=5, length=20, labelsize=70)

  fig.set_snap(True)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  filename = "plots/zg/zg.cut." + str(zg_filename) + ".pt.cut." + str(pT_lower_cut) + ".ratio.over." + ratio_denominator + ".th." + str(theory) + ".mc." + str(mc) + ".data." + str(data) + ".pdf"
  
  print filename

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

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']

  colors = ['magenta', 'blue', 'orange', 'green', 'black', 'red']
  expected_trigger_names = ["HLT\_Jet180U", "HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U" ]

  color = colors[expected_trigger_names.index(trigger_1.replace("_", "\_"))]

  pt_hist_trigger_1 = Hist(2000, 0, pT_upper_limit, title=trigger_1.replace("_", "\_"), color=color, markersize=1.0, linewidth=5)
  pt_hist_trigger_2 = Hist(2000, 0, pT_upper_limit, title=trigger_2.replace("_", "\_"), color=color, markersize=1.0, linewidth=5)


  for i in range(0, len(pTs)):
    if trigger_1 in trigger_names[i]:
      pt_hist_trigger_1.Fill(pTs[i], prescales[i])

    # The len thingy is to make sure trigger names like HLT_Jet15U_HcalNoiseFiltered_v3 are excluded.
    # if trigger_2 in trigger_names[i] and len(trigger_names[i]) > (len(trigger_2) + 3):
    if trigger_2 in trigger_names[i]:
      pt_hist_trigger_2.Fill(pTs[i], prescales[i])


  pt_hist_trigger_1.Divide(pt_hist_trigger_2)


  rplt.errorbar(pt_hist_trigger_1, color=color)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
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
  else:
    lower_pT = 0


  if lower_pT != 0:
    # CMS Vertical line.
    plt.plot([lower_pT, lower_pT], [plt.gca().get_ylim()[0], 1.], color="purple", linewidth=3, linestyle="dashed")



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


  # plt.plot([efficient_pt_x, efficient_pt_x], [plt.gca().get_ylim()[0], 1000], color="gray", linewidth=3, linestyle="dashed")
  

  if lower_pT != 0:
    plt.gca().annotate("CMS: " + str(lower_pT) + " GeV", xy=(lower_pT, 1.), xycoords='data', xytext=(-100, 250),  textcoords='offset points', color="purple", size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor="purple", zorder=99, connectionstyle="angle3,angleA=0,angleB=90") )
  
  # plt.gca().annotate("MOD: " + str(int(efficient_pt_x)) + " GeV", xy=(efficient_pt_x, efficient_pt_y), xycoords='data', xytext=(250, 200), textcoords='offset points', color="gray", size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor="gray", zorder=99, connectionstyle="angle3,angleA=45,angleB=-90") )



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

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']

  colors = ['magenta', 'blue', 'green', 'red', 'black']
  expected_trigger_names = ["HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U" ]
  labels = ["Jet140U / 100U", "Jet100U / 70U", "Jet70U / 50U", "Jet50U / 30U", "Jet30U / 15U" ]


  cms_turn_on_pTs = [0, 0, 153, 114]


  pt_hists = []
  for j in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist(60, 0, 300, color=colors[j], title=labels[j], markersize=1.0, linewidth=5))

  for i in range(0, len(expected_trigger_names)):
    for j in range(0, len(pTs)):
      # The len thingy is to make sure trigger names like HLT_Jet15U_HcalNoiseFiltered_v3 are excluded.
      if expected_trigger_names[i].replace("\\", "") in trigger_names[j]:
        pt_hists[i].Fill(pTs[j], prescales[j])

      
  for i in range(0, len(pt_hists) - 1):
    ratio_hist = pt_hists[i] / pt_hists[i + 1]
    rplt.errorbar(ratio_hist)

    # CMS Vertical line.
    if cms_turn_on_pTs[i] != 0:
      plt.plot([cms_turn_on_pTs[i], cms_turn_on_pTs[i]], [plt.gca().get_ylim()[0], 1.], color=colors[i], linewidth=5, linestyle="dashed")
      plt.gca().annotate("CMS\n" + str(cms_turn_on_pTs[i]) + " GeV", xy=(cms_turn_on_pTs[i], 1.), xycoords='data', xytext=(-100, 250),  textcoords='offset points', color=colors[i], size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor=colors[i], zorder=99, connectionstyle="angle3,angleA=0,angleB=90") )


  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  plt.gcf().text(0.29, 0.885, "Prelim. (20\%)", fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_snap(1)

  # Horizontal Line.
  plt.plot(list(pt_hists[0].x()), [1] * len(list(pt_hists[0].x())), color="black", linewidth=5, linestyle="dashed")



  plt.yscale('log')
  plt.gca().set_ylim(plt.gca().get_ylim()[0], 1000)

  # Legend.
  handles, labels = plt.gca().get_legend_handles_labels()
  first_legend = plt.gca().legend(handles[::-1], labels[::-1], fontsize=60,  frameon=0, borderpad=0.1, bbox_to_anchor=[0.90, 0.98])
  ax = plt.gca().add_artist(first_legend)

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5} $"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.825, 0.63])





  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=95, rotation=0)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=50.)


  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.gcf().set_snap(True)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/all_efficiency_curves.pdf")
  # plt.show()
  plt.clf()



# # Trigger Efficiency Curves Begin.

# All Plots in the same one: 50U to 140U. Upto pT = 300 GeV. 

# For turn-on curves, remove 180U. Also upper limit 0 to 300 GeV.

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

plot_charged_and_all_zgs(150, '0.05', 'zg_05', n_bins=8, y_max_limit=18)
plot_charged_and_all_zgs(150, '0.1', 'zg_1', n_bins=8, y_max_limit=10)
plot_charged_and_all_zgs(150, '0.2', 'zg_2', n_bins=8, y_max_limit=10)

plot_charged_and_all_zgs(300, '0.05', 'zg_05', n_bins=4, y_max_limit=15)
plot_charged_and_all_zgs(300, '0.1', 'zg_1', n_bins=4, y_max_limit=15)
plot_charged_and_all_zgs(300, '0.2', 'zg_2', n_bins=4, y_max_limit=15)

plot_charged_and_all_zgs(600, '0.05', 'zg_05', n_bins=2, y_max_limit=15)
plot_charged_and_all_zgs(600, '0.1', 'zg_1', n_bins=2, y_max_limit=15)
plot_charged_and_all_zgs(600, '0.2', 'zg_2', n_bins=2, y_max_limit=15)

# # Charged zg Ends.











# # zg with PFC pT_cut Begins.

# Everything, 2 GeV, 5 GeV. Same plots as charged. 


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


