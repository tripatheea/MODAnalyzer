# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat
from __future__ import division

import sys
import math
from collections import defaultdict

# matplotlib
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

plt.rc('font', family='serif', size=30)


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

def plot_pts():
  properties = parse_file(input_analysis_file)

  pTs = properties['uncorrected_hardest_pts']
  corrected_pTs = properties['corrected_hardest_pts']
  prescales = properties['prescales']

  # no. of bins, xlower, xhigher
  uncorrected_pt_hist = Hist(50, 0, 800, title='Uncorrected', markersize=1.0, color='green')
  corrected_pt_hist = Hist(50, 0, 800, title='Corrected', markersize=1.0, color='red')

  map(uncorrected_pt_hist.Fill, pTs, prescales)
  map(corrected_pt_hist.Fill, corrected_pTs, prescales)

  uncorrected_pt_hist.Scale(1.0 / uncorrected_pt_hist.GetSumOfWeights())
  corrected_pt_hist.Scale(1.0 / corrected_pt_hist.GetSumOfWeights())

  rplt.errorbar(uncorrected_pt_hist, xerr=False, emptybins=False)
  rplt.errorbar(corrected_pt_hist, xerr=False, emptybins=False)

  plt.autoscale(True)
   
  plt.yscale('log')

  plt.legend()
  plt.xlabel('$p_T$ (GeV)')
  plt.suptitle("$p_T$ (GeV) Spectrum of anti-kT Jets (R = 0.5)")
  plt.grid(True)

  plt.savefig("plots/ak5_pt_distribution.pdf")
  plt.show()


def plot_zg():
  pT_lower_cut = 150
  pfc_pT_cut = 0
  properties = parse_file(input_analysis_file, pT_lower_cut, pfc_pT_cut)
  properties_pythia = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_pythia_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)
  properties_herwig = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_herwig_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)

  zgs = [properties['zg_05'], properties['zg_1'], properties['zg_2']]
  
  zg_pythias = [properties_pythia['zg_05'], properties_pythia['zg_1'], properties_pythia['zg_2']]
  zg_herwigs = [properties_herwig['zg_05'], properties_herwig['zg_1'], properties_herwig['zg_2']]

  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  colors_2 = ['gray', 'orange', 'magenta']

  labels = ['CMS $z_{cut}$ = 0.05', 'CMS $z_{cut}$ = 0.1', 'CMS $z_{cut}$ = 0.2']
  pythia_labels = ['Pythia 8 $z_{cut}$ = 0.05', 'Pythia 8 $z_{cut}$ = 0.1', 'Pythia 8 $z_{cut}$ = 0.2']
  herwig_labels = ['Herwig 2 $z_{cut}$ = 0.05', 'Herwig 2 $z_{cut}$ = 0.1', 'Herwig 2 $z_{cut}$ = 0.2']
  
  for i in range(0, len(zgs)):
    # no. of bins, xlower, xhigher
    zg_hist = Hist(50, 0.0, 0.5, title=labels[i], markersize=1.0, color=colors[i])

    map(zg_hist.Fill, zgs[i], prescales)
    
    zg_hist.Scale(1.0 / zg_hist.GetSumOfWeights())
    
    rplt.errorbar(zg_hist, xerr=False, emptybins=False)

  # for j in range(0, len(zg_pythias)):
  #   zg_pythia_hist = Hist(50, 0, 0.5, title=pythia_labels[j], markersize=1.0, color=colors[j])

  #   map(zg_pythia_hist.Fill, zg_pythias[j])

  #   zg_pythia_hist.Scale(1.0 / zg_pythia_hist.GetSumOfWeights())

  #   rplt.hist(zg_pythia_hist)

  # for j in range(0, len(zg_herwigs)):
  #   zg_herwig_hist = Hist(50, 0, 0.5, title=herwig_labels[j], markersize=1.0, color=colors_2[j], linestyle="2") # 1=solid, 2=dash, 3=dot, 4=dash-dot

  #   map(zg_herwig_hist.Fill, zg_herwigs[j])

  #   zg_herwig_hist.Scale(1.0 / zg_herwig_hist.GetSumOfWeights())

  #   rplt.hist(zg_herwig_hist)


  plt.autoscale(True)

  plt.legend()
  plt.xlabel("Symmetry Measure(z)")
  plt.suptitle("Symmetry Measure(z) with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV")
  plt.grid(True)

  plt.savefig("plots/zg_distribution_pt_cut_" + str(pT_lower_cut) + ".pdf")
  plt.show()


def plot_zg_pfc_pt_cut(pfc_pT_cut):
  pT_lower_cut = 150
  properties = parse_file(input_analysis_file, pT_lower_cut, pfc_pT_cut)

  zgs = [properties['zg_05_pt_' + str(int(pfc_pT_cut))], properties['zg_1_pt_' + str(int(pfc_pT_cut))], properties['zg_2_pt_' + str(int(pfc_pT_cut))]]

  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['CMS $z_{cut}$ = 0.05', 'CMS $z_{cut}$ = 0.1', 'CMS $z_{cut}$ = 0.2']
  
  for i in range(0, len(zgs)):
    # no. of bins, xlower, xhigher
    zg_hist = Hist(50, 0.0, 0.5, title=labels[i], markersize=1.0, color=colors[i])

    map(zg_hist.Fill, zgs[i], prescales)
    
    zg_hist.Scale(1.0 / zg_hist.GetSumOfWeights())
    
    rplt.errorbar(zg_hist, xerr=False, emptybins=False)

  plt.autoscale(True)

  plt.legend()
  plt.xlabel("Symmetry Measure(z)")
  plt.suptitle("Symmetry Measure(z) with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV & PFC $p_{T cut}$ = " + str(pfc_pT_cut) + " GeV")
  plt.grid(True)

  plt.savefig("plots/zg_distribution_pt_cut_" + str(pT_lower_cut) + "_pfc_pt_cut_" + str(pfc_pT_cut) + ".pdf")
  plt.show()


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
  plt.grid(True)

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
  plt.grid(True)

  plt.savefig("plots/mass_drop_distribution_pt_cut_" + str(pT_lower_cut) + ".pdf")
  plt.show()



def plot_pdgid_pt():
  properties = parse_file(input_analysis_file)
  pdgid_map = { 1: "d", 130: "$K^0_L$ Meson", 11: "$e^-$", -211: "$\pi^-$", 13: "$\mu^-$", 211: "$\pi^+$", -11: "$e^+$", 22: "$\gamma$", 2: "u", -13: "$\mu^+$" }

  pdgids = properties['hardest_pfc_pdgid']
  pTs = properties['hardest_pfc_pt']
  prescales = properties['prescales']

  pdgid_pts = defaultdict(list)
  pdgid_prescales = defaultdict(list)

  for i in range(0, len(pdgids)):
    pdgid_pts[pdgids[i]].append(pTs[i])
    pdgid_prescales[pdgids[i]].append(prescales[i])


  for pdgid in pdgid_pts:
    pdgid_hist = Hist(100, 0, 400, title="pdgId = " + str(pdgid), markersize=1.0, color='blue')

    map(pdgid_hist.Fill, pdgid_pts[pdgid], pdgid_prescales[pdgid])

    if pdgid_hist.GetSumOfWeights() != 0:
      pdgid_hist.Scale(1.0 / pdgid_hist.GetSumOfWeights())

    rplt.errorbar(pdgid_hist, xerr=False, emptybins=False)

    plt.yscale('log')

    plt.autoscale(True)

    plt.xlabel('$p_{T}$ GeV')
    plt.suptitle("Hardest " + pdgid_map[pdgid] + " s (pdgid=" + str(int(pdgid)) + ") pT Distribution")
    plt.grid(True)

    plt.savefig("plots/hardest_pdgid_" + str(int(pdgid)) + "_pt_distribution.pdf")
    plt.show()


def plot_pdgid_zg():
  pT_lower_cut = 153
  properties = parse_file(input_analysis_file, pT_lower_cut)

  pdgid_map = { 1: "d", 130: "$K^0_L$ Meson", 11: "$e^-$", -211: "$\pi^-$", 13: "$\mu^-$", 211: "$\pi^+$", -11: "$e^+$", 22: "$\gamma$", 2: "u", -13: "$\mu^+$" }

  labels = ['$z_{cut}$ = 0.05', '$z_{cut}$ = 0.1', '$z_{cut}$ = 0.2']

  pdgids = properties['hardest_pfc_pdgid']
  zg_05s = properties['zg_05']
  zg_1s = properties['zg_1']
  zg_2s = properties['zg_2']
  prescales = properties['prescales']

  pdgid_zg_05s = defaultdict(list)
  pdgid_zg_1s = defaultdict(list)
  pdgid_zg_2s = defaultdict(list)

  pdgid_prescales = defaultdict(list)

  for i in range(0, len(pdgids)):
    pdgid_zg_05s[pdgids[i]].append(zg_05s[i])
    pdgid_zg_1s[pdgids[i]].append(zg_1s[i])
    pdgid_zg_2s[pdgids[i]].append(zg_2s[i])

    pdgid_prescales[pdgids[i]].append(prescales[i])


  for pdgid in pdgid_zg_05s:

    pdgid_zg_05_hist = Hist(25, 0, 0.5, title=labels[0], markersize=1.0, color='blue')
    pdgid_zg_1_hist = Hist(25, 0, 0.5, title=labels[1], markersize=1.0, color='red')
    pdgid_zg_2_hist = Hist(25, 0, 0.5, title=labels[2], markersize=1.0, color='green')

    map(pdgid_zg_05_hist.Fill, pdgid_zg_05s[pdgid], pdgid_prescales[pdgid])
    map(pdgid_zg_1_hist.Fill, pdgid_zg_1s[pdgid], pdgid_prescales[pdgid])
    map(pdgid_zg_2_hist.Fill, pdgid_zg_2s[pdgid], pdgid_prescales[pdgid])
    
    if pdgid_zg_05_hist.GetSumOfWeights() != 0:
      pdgid_zg_05_hist.Scale(1.0 / pdgid_zg_05_hist.GetSumOfWeights())

    if pdgid_zg_1_hist.GetSumOfWeights() != 0:
      pdgid_zg_1_hist.Scale(1.0 / pdgid_zg_1_hist.GetSumOfWeights())

    if pdgid_zg_2_hist.GetSumOfWeights() != 0:
      pdgid_zg_2_hist.Scale(1.0 / pdgid_zg_2_hist.GetSumOfWeights())

    rplt.errorbar(pdgid_zg_05_hist, xerr=False, emptybins=False)
    rplt.errorbar(pdgid_zg_1_hist, xerr=False, emptybins=False)
    rplt.errorbar(pdgid_zg_2_hist, xerr=False, emptybins=False)

    plt.autoscale(True)

    plt.legend()
    plt.xlabel('Symmetery Measure (z)')
    plt.suptitle("Hardest " + pdgid_map[pdgid] + " s (pdgid=" + str(int(pdgid)) + ") Symmetry Measure(z) with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV")
    plt.grid(True)

    plt.savefig("plots/hardest_pdgid_" + str(int(pdgid)) + "_zg_distribution.pdf")
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
  plt.grid(True)

  plt.savefig("plots/hardest_charged_pt_distribution.pdf")
  plt.show()


def plot_charged_zgs():
  pT_lower_cut = 150
  properties = parse_file(input_analysis_file, pT_lower_cut)

  charged_zgs = [properties['zg_charged_05'], properties['zg_charged_1'], properties['zg_charged_2']]
  prescales = properties['prescales']



  colors = ['red', 'blue', 'green']
  labels = ['$z_{cut}$ = 0.05', '$z_{cut}$ = 0.1', '$z_{cut}$ = 0.2']

  for i in range(0, len(charged_zgs)):
    # no. of bins, xlower, xhigher
    zg_hist = Hist(50, 0.0, 0.5, title=labels[i], markersize=1.0, color=colors[i])

    map(zg_hist.Fill, charged_zgs[i], prescales)
    
    if zg_hist.GetSumOfWeights() != 0:
      zg_hist.Scale(1.0 / zg_hist.GetSumOfWeights())
    
    rplt.errorbar(zg_hist, xerr=False, emptybins=False)

  plt.autoscale(True)

  plt.legend()
  plt.xlabel("Symmetry Measure(z)")
  plt.suptitle("Symmetry Measure(z) with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV for charged PFCs")
  plt.grid(True)

  plt.savefig("plots/zg_distribution_charged_pt_cut_" + str(pT_lower_cut) + ".pdf")
  plt.show()


def plot_charged_zg_pfc_pt_cut(pfc_pT_cut):
  pT_lower_cut = 150
  properties = parse_file(input_analysis_file, pT_lower_cut, pfc_pT_cut)

  pdgids = properties['hardest_pfc_pdgid']
  zgs = [properties['zg_05_pt_' + str(int(pfc_pT_cut))], properties['zg_1_pt_' + str(int(pfc_pT_cut))], properties['zg_2_pt_' + str(int(pfc_pT_cut))]]

  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['CMS $z_{cut}$ = 0.05', 'CMS $z_{cut}$ = 0.1', 'CMS $z_{cut}$ = 0.2']
  
  for i in range(0, len(zgs)):
    if (abs(pdgids[i]) == 211) and (abs(pdgids[i]) == 11) and (abs(pdgids[i]) == 13):
      # no. of bins, xlower, xhigher
      zg_hist = Hist(50, 0.0, 0.5, title=labels[i], markersize=1.0, color=colors[i])

      map(zg_hist.Fill, zgs[i], prescales)
      
      zg_hist.Scale(1.0 / zg_hist.GetSumOfWeights())
      
      rplt.errorbar(zg_hist, xerr=False, emptybins=False)

  plt.autoscale(True)

  plt.legend()
  plt.xlabel("Symmetry Measure(z)")
  plt.suptitle("Symmetry Measure(z) with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV & PFC $p_{T cut}$ = " + str(pfc_pT_cut) + " GeV")
  plt.grid(True)

  plt.savefig("plots/zg_charged_distribution_pt_cut_" + str(pT_lower_cut) + "_pfc_pt_cut_" + str(pfc_pT_cut) + ".pdf")
  plt.show()



def plot_hardest_pt_corresponding_triggers():
  properties = parse_file(input_analysis_file)

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']

  expected_trigger_names = ["HLT_Jet70U", "HLT_Jet50U", "HLT_Jet30U", "HLT_Jet15U", "HLT_L1Jet6U"]

  colors = ['red', 'blue', 'pink', 'green', 'orange']

  pt_hists = []
  for i in range(0, len(expected_trigger_names)):
    # no. of bins, xlower, xhigher
    pt_hists.append(Hist(50, 0, 500, title=expected_trigger_names[i], markersize=1.0, color=colors[i]))


  for i in range(0, len(pTs)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j] in trigger_names[i]:
        pt_hists[j].Fill(pTs[i], prescales[i])

  
  for k in range(0, len(pt_hists)):
    # if pt_hists[k].GetSumOfWeights() != 0:
      # pt_hists[k].Scale(1.0 / pt_hists[k].GetSumOfWeights())
    rplt.errorbar(pt_hists[k], xerr=False, emptybins=False)


  plt.yscale('log')
  
  plt.autoscale(True)
   
  plt.yscale('log')

  plt.legend()
  plt.xlabel('$p_T$ (GeV)')
  plt.suptitle("$p_T$ (GeV) Spectrum of anti-kT Jets (R = 0.5) & Corresponding Triggers")
  plt.grid(True)

  plt.savefig("plots/hardest_pt_corresponding_triggers.pdf")
  plt.show()

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

  plt.savefig("plots/zg_2_vs_charged_zg_2.pdf")
  plt.show()






def plot_zg_th_mc_data(zg_cut, zg_filename):
  pT_lower_cut = 150
  pfc_pT_cut = 0

  properties = parse_file(input_analysis_file, pT_lower_cut, pfc_pT_cut)
  properties_pythia = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_pythia_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)
  properties_herwig = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_herwig_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)

  zg_data = properties[zg_filename]
  
  zg_pythias = properties_pythia[zg_filename]
  zg_herwigs = properties_herwig[zg_filename]

  prescales = properties['prescales']

  data_label = 'CMS 2010 Open Data'
  pythia_label = 'Pythia 8'
  herwig_label = 'Herwig++'
  theory_label = 'Theory (MLL)'

  
  gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1]) 

 
  ax0 = plt.subplot(gs[0])
  ax1 = plt.subplot(gs[1])


  # Theory Plots Begin.
  
  points_th_gluon = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Simone/results/band_gluon_pt" + str(pT_lower_cut) + "_zc01.dat")
  points_th_quark = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Simone/results/band_quark_pt" + str(pT_lower_cut) + "_zc01.dat")

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
    
    

  area_theory_y_max = simps(theory_y_max, theory_x)
  # weighted_theory_y_max = map(lambda x: x / area_theory_y_max, theory_y_max)
  weighted_theory_y_max = theory_y_max
  ax0.plot(theory_x, weighted_theory_y_max, alpha=0.0, color='red')
  
  area_theory_y_line = simps(theory_y_line, theory_x)
  # weighted_theory_y_line = map(lambda x: x / area_theory_y_line, theory_y_line)
  weighted_theory_y_line = theory_y_line
  ax0.plot(theory_x, weighted_theory_y_line, label=theory_label, alpha=1.0, color='red')

  area_theory_y_min = simps(theory_y_min, theory_x)
  # weighted_theory_y_min = map(lambda x: x / area_theory_y_min, theory_y_min)
  weighted_theory_y_min = theory_y_min
  ax0.plot(theory_x, weighted_theory_y_min, alpha=0.0, color='red')


  ax0.fill_between(theory_x, theory_y_max, theory_y_min, norm=1, where=np.less_equal(theory_y_min, theory_y_max), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)


  # Theory Plot Ends.
  
  
  # Data Plot Begins.
  zg_data_hist = Hist(75, 0.0, 0.6, title=data_label, markersize=0.75, color='black')
  bin_width_data = (zg_data_hist.upperbound() - zg_data_hist.lowerbound()) / zg_data_hist.nbins()

  map(zg_data_hist.Fill, zg_data, prescales)
  
  zg_data_hist.Scale(1.0 / ( zg_data_hist.GetSumOfWeights() * bin_width_data ))

  norm_data_prescales = map(lambda x: x / ( zg_data_hist.GetSumOfWeights() * bin_width_data ), prescales)
  
  data_plot = rplt.errorbar(zg_data_hist, emptybins=False, axes=ax0)



  # Data Plots Ends.

  
  # Simulation Plots Begin. 
  
  # Pythia.
  
  zg_pythia_hist = Hist(75, 0, 0.6, title=pythia_label, markersize=1.0, color='blue')
  bin_width_pythia = (zg_pythia_hist.upperbound() - zg_pythia_hist.lowerbound()) / zg_pythia_hist.nbins()

  map(zg_pythia_hist.Fill, zg_pythias)

  zg_pythia_hist.Scale(1.0 / ( zg_pythia_hist.GetSumOfWeights() * bin_width_pythia ))

  pythia_plot = rplt.hist(zg_pythia_hist, axes=ax0)
  
  # Pythia Ends.
  
  # Herwig. 
  
  zg_herwig_hist = Hist(75, 0, 0.6, title=herwig_label, markersize=1.0, color='green')
  bin_width_herwig = (zg_herwig_hist.upperbound() - zg_herwig_hist.lowerbound()) / zg_herwig_hist.nbins()

  map(zg_herwig_hist.Fill, zg_herwigs)

  zg_herwig_hist.Scale(1.0 / ( zg_herwig_hist.GetSumOfWeights() * bin_width_herwig ))

  herwig_plot = rplt.hist(zg_herwig_hist, axes=ax0)
  
  # Herwig Ends.

  # Simulation Plots End.

  
  # Normalized-Over-Data Plot Begins.

  

  zg_herwig_hist.Divide(zg_data_hist)
  rplt.hist(zg_herwig_hist, axes=ax1)

  zg_pythia_hist.Divide(zg_data_hist)
  rplt.hist(zg_pythia_hist, axes=ax1)

  zg_data_hist.Divide(zg_data_hist)
  rplt.errorbar(zg_data_hist, axes=ax1)

  # Theory-Over-Data Plot.
  

  data_points_x = data_plot[0].get_xdata()
  data_points_y = data_plot[0].get_ydata()

  data_plot_points_x = []
  data_plot_points_y = []

  for i in range(0, len(data_points_x)):
    if float(data_points_x[i]) >= float(zg_cut):
      data_plot_points_x.append(data_points_x[i])
      data_plot_points_y.append(data_points_y[i])




  theory_extrapolated_min = []
  theory_extrapolated_line = []
  theory_extrapolated_max = []


  j = 0
  for i in range(0, len(data_plot_points_x)):
    x = data_plot_points_x[i]

    
    if x >= theory_x[j] and x <= theory_x[j + 1]:
      x1, x2 = theory_x[j], theory_x[j + 1]
      
      y1, y2 = theory_y_min[j], theory_y_min[j + 1]
      y_min = ( ( (x - x1) / (x2 - x1) ) * (y2 - y1) ) + y1  

      y1, y2 = theory_y_line[j], theory_y_line[j + 1]
      y_line = ( ( (x - x1) / (x2 - x1) ) * (y2 - y1) ) + y1  

      y1, y2 = theory_y_max[j], theory_y_max[j + 1]
      y_max = ( ( (x - x1) / (x2 - x1) ) * (y2 - y1) ) + y1  

    elif x > theory_x[j + 1]:
      j += 1

      x1, x2 = theory_x[j], theory_x[j + 1]
      
      y1, y2 = theory_y_min[j], theory_y_min[j + 1]
      y_min = ( ( (x - x1) / (x2 - x1) ) * (y2 - y1) ) + y1  

      y1, y2 = theory_y_line[j], theory_y_line[j + 1]
      y_line = ( ( (x - x1) / (x2 - x1) ) * (y2 - y1) ) + y1  

      y1, y2 = theory_y_max[j], theory_y_max[j + 1]
      y_max = ( ( (x - x1) / (x2 - x1) ) * (y2 - y1) ) + y1  
    elif x == 0:
      y_min, y_line, y_max = 0, 0, 0
    else:
      raise ValueError("Some very weird zg value found!", x)

    theory_extrapolated_min.append(y_min)
    theory_extrapolated_line.append(y_line)
    theory_extrapolated_max.append(y_max)


  ratio_theory_line_to_data = [m / n for m, n in zip(theory_extrapolated_line, data_plot_points_y)]
  ax1.plot(data_plot_points_x, ratio_theory_line_to_data, alpha=1.0, color='red')

  ratio_theory_min_to_data = [m / n for m, n in zip(theory_extrapolated_min, data_plot_points_y)]
  ax1.plot(data_plot_points_x, ratio_theory_min_to_data, alpha=0.0, color='red')

  ratio_theory_max_to_data = [m / n for m, n in zip(theory_extrapolated_max, data_plot_points_y)]
  ax1.plot(data_plot_points_x, ratio_theory_max_to_data, alpha=0.0, color='red')

  ax1.fill_between(data_plot_points_x, ratio_theory_max_to_data, ratio_theory_min_to_data, norm=1, where=np.less_equal(ratio_theory_min_to_data, ratio_theory_max_to_data), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)
  
  # Legend.

  handles, labels = ax0.get_legend_handles_labels()
  
  handles = handles[::-1]
  labels = labels[::-1]

  for i in range(0, len(labels)):
    if labels[i] == theory_label:
      handles[i] = mpatches.Patch(facecolor='red', edgecolor='red', alpha=1.0, linewidth=0, label=theory_label, hatch='-')

      a = handles[i]
      a.set_facecolor('pink')
      handles[i] = a

  first_legend = ax0.legend(handles, labels, loc=0, frameon=0, borderpad=0.1)
  ax = ax0.add_artist(first_legend)
  
  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra, extra]
  labels = ["$p_{T cut}$ = " + str(pT_lower_cut) + "; R = $0.5$", "$z_{cut}$ = " + zg_cut + "; $\\beta$ = 0"]
  ax0.legend(handles, labels, loc=2, frameon=0, borderpad=0.1)


  # Legend Ends.


  ax0.autoscale(True)
  ax1.autoscale(True)
  
  ax0.set_ylim(0, 10)
  ax1.set_ylim(0.5, 1.8)


  fn = get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)
  ab = AnnotationBbox(OffsetImage(read_png(fn), zoom=1.0), (0.55, 0.80), boxcoords="offset points")
  ax0.add_artist(ab)

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

  plt.savefig("plots/zg_distribution_data_mc_th_pt_cut_" + str(pT_lower_cut) + ".pdf")
  

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


def plot_th():
  pT_lower_cut = 300
  pfc_pT_cut = 0
  
  points_th = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Simone/results/band_gluon_pt" + str(pT_lower_cut) + "_zc005.dat")

  
  points = defaultdict(list)

  for x in points_th:
    points[x] = [ points_th[x][0], points_th[x][1], points_th[x][2], points_th[x][3], points_th[x][4], points_th[x][5] ]

  keys = points.keys()
  keys.sort()

  x = keys

  y = []
  for j in range(0, 6):
    y.append([points[i][j] for i in keys])

  plt.plot(x, y[0], x, y[1], x, y[2], x, y[3], x, y[4], x, y[5], alpha=0.5, color='green')

  for i in range(0, 5):
    plt.fill_between(x, y[i], y[i + 1], where=np.less_equal(y[i], y[i + 1]), facecolor='green', interpolate=True, alpha=0.5)
  

  







# plot_th()

# plot_zg_th_mc_data('0.1', 'zg_1')


# plot_pts()

# plot_zg()

# plot_dr()
# plot_mu()

# plot_pdgid_pt()

# plot_pdgid_zg()

# plot_charged_pt()

# plot_charged_zgs()


# plot_zg_pfc_pt_cut(pfc_pT_cut=1)
# plot_zg_pfc_pt_cut(pfc_pT_cut=2)
# plot_zg_pfc_pt_cut(pfc_pT_cut=3)
# plot_zg_pfc_pt_cut(pfc_pT_cut=5)
# plot_zg_pfc_pt_cut(pfc_pT_cut=10)


# plot_charged_zg_pfc_pt_cut(pfc_pT_cut=1)
# plot_charged_zg_pfc_pt_cut(pfc_pT_cut=2)
# plot_charged_zg_pfc_pt_cut(pfc_pT_cut=3)
# plot_charged_zg_pfc_pt_cut(pfc_pT_cut=5)
# plot_charged_zg_pfc_pt_cut(pfc_pT_cut=10)

# plot_hardest_pt_corresponding_triggers()


plot_2d_hist()
