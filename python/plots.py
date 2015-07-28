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

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

  plt.savefig("plots/ak5_pt_distribution.pdf")
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

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

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

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

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

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

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

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

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

  fig = plt.gcf()
  fig.set_size_inches(20, 20, forward=1)

  plt.savefig("plots/zg_2_vs_charged_zg_2.pdf")
  plt.show()





def plot_turn_on_curves():
  # properties = parse_file_turn_on(input_analysis_file)
  
  properties = parse_file_turn_on('./turn_on.dat')

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescales']

  # expected_trigger_names = ["HLT_DiJetAve15U", "HLT_DiJetAve30U", "HLT_DiJetAve50U", "HLT_DiJetAve70U", "HLT_EcalOnly_SumEt160", "HLT_HT100U", "HLT_HT120U", "HLT_HT140U", "HLT_Jet100U", "HLT_Jet15U_HcalNoiseFiltered", "HLT_QuadJet20U", "HLT_QuadJet25U", "HLT_Jet70U", "HLT_Jet50U", "HLT_Jet30U", "HLT_Jet15U", "HLT_L1Jet6U"]
  expected_trigger_names = ["HLT_Jet180U", "HLT_Jet140U", "HLT_Jet100U", "HLT_Jet70U", "HLT_Jet50U", "HLT_Jet30U" ]

  colors = ['red', 'blue', 'orange', 'green', 'black', 'pink']





  pt_hists = []
  for i in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist(1000, 0, 500, title=expected_trigger_names[i], markersize=1.0, color=colors[i], linewidth=5))

  for i in range(0, len(pTs)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j] in trigger_names[i]:
        pt_hists[j].Fill(pTs[i], prescales[i])

  for k in range(0, len(pt_hists)):
    rplt.errorbar(pt_hists[k], xerr=0, yerr=0)
    # rplt.hist(pt_hists[k])

  # plt.axvspan(153, pt_hists[0].upperbound(), fc="red", linewidth=0, alpha=0.5)
  # plt.axvspan(114, 153, fc="blue", linewidth=0, alpha=0.5)
  # plt.axvspan(84, 114, fc="orange", linewidth=0, alpha=0.5)
  # plt.axvspan(56, 84, fc="green", linewidth=0, alpha=0.5)
  # plt.axvspan(37, 56, fc="black", linewidth=0, alpha=0.5)
  # plt.axvspan(37, 56, fc="pink", linewidth=0, alpha=0.5)


  fn = get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)
  ab = AnnotationBbox(OffsetImage(read_png(fn), zoom=1.0), (305, 1375), frameon=0, boxcoords='figure points')
  plt.gca().add_artist(ab)

  plt.gcf().text(0.300, 0.93, "Preliminary \n(25% sample)", fontsize=40, weight='bold', color='#444444', multialignment='center')

  plt.autoscale(True)
  plt.yscale('log')

  plt.legend()
  plt.xlabel('$p_T$ (GeV)')


  plt.gcf().set_size_inches(20, 20, forward=1)


  plt.savefig("plots/turn_on_curves.pdf")
  plt.show()




def plot_zg_th_mc_data(zg_cut, zg_filename, ratio_denominator="theory", data=True, mc=True, theory=True):
  pT_lower_cut = 150
  pfc_pT_cut = 0

  zg_cut = float(zg_cut)

  properties = parse_file(input_analysis_file, pT_lower_cut, pfc_pT_cut)
  properties_pythia = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_pythia_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)
  properties_herwig = parse_mc_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_herwig_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)

  zg_data = properties[zg_filename]
  
  zg_pythias = properties_pythia[zg_filename]
  zg_herwigs = properties_herwig[zg_filename]

  prescales = properties['prescales']

  data_label = "CMS 2010 Open Data  " if data else "                    "
  pythia_label = "Pythia 8.205        " if mc else "                    "
  herwig_label = "Herwig++ 2.6.3      " if mc else "                    "
  theory_label = "Theory (MLL)        " if theory else "                    "

  
  gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1]) 

 
  ax0 = plt.subplot(gs[0])
  ax1 = plt.subplot(gs[1])


  # Theory Plots Begin.
  
  points_th_gluon = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Simone/results_7_24_15/band_gluon_pt" + str(pT_lower_cut) + "_zc01.dat")
  points_th_quark = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Simone/results_7_24_15/band_quark_pt" + str(pT_lower_cut) + "_zc01.dat")

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


  
  
  # Data Plot Begins.
  
  zg_data_hist = Hist(60, 0.0, 0.6, title=data_label, markersize=2.5, color='black')
  bin_width_data = (zg_data_hist.upperbound() - zg_data_hist.lowerbound()) / zg_data_hist.nbins()

  map(zg_data_hist.Fill, zg_data, prescales)
  
  zg_data_hist.Scale(1.0 / ( zg_data_hist.GetSumOfWeights() * bin_width_data ))

  norm_data_prescales = map(lambda x: x / ( zg_data_hist.GetSumOfWeights() * bin_width_data ), prescales)
  
  if data:
    data_plot, caplines, barlinecols = rplt.errorbar(zg_data_hist, xerr=1, yerr=1, emptybins=False, axes=ax0, ls='None', marker='o', markersize=8, pickradius=3, elinewidth=3, alpha=1.0)
  else:
    zg_data_hist.SetTitle("")
    data_plot, caplines, barlinecols = rplt.errorbar(zg_data_hist, xerr=1, yerr=1, emptybins=False, axes=ax0, ls='None', marker='o', markersize=8, pickradius=3, elinewidth=3, alpha=0.0)


  data_x_errors, data_y_errors = [], []
  for x_segment in barlinecols[0].get_segments():
    data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 10.)
  for y_segment in barlinecols[1].get_segments():
    data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 10.)

  data_points_x = data_plot.get_xdata()
  data_points_y = data_plot.get_ydata()

  # print sorted(data_points_x)



  # Data Plots Ends.

  
  # Simulation Plots Begin. 
  
  # Pythia.
  
  zg_pythia_hist = Hist(60, 0, 0.6, title=pythia_label, markersize=5.0, color='blue', linewidth=5)
  bin_width_pythia = (zg_pythia_hist.upperbound() - zg_pythia_hist.lowerbound()) / zg_pythia_hist.nbins()

  map(zg_pythia_hist.Fill, zg_pythias)

  zg_pythia_hist.Scale(1.0 / ( zg_pythia_hist.GetSumOfWeights() * bin_width_pythia ))

  if mc:
    # pythia_plot = rplt.hist(zg_pythia_hist, axes=ax0)
    pythia_plot = ax0.hist(zg_pythias, label=pythia_label, bins=60, normed=1, histtype='step', color='blue', linewidth=5)
  else:
    pythia_plot = ax0.hist(zg_pythias, bins=60, normed=1, histtype='step', color='blue', linewidth=0)

  
  # Pythia Ends.
  
  # Herwig. 
  
  zg_herwig_hist = Hist(60, 0, 0.6, title=herwig_label, markersize=5.0, color='green', linewidth=5)
  bin_width_herwig = (zg_herwig_hist.upperbound() - zg_herwig_hist.lowerbound()) / zg_herwig_hist.nbins()

  map(zg_herwig_hist.Fill, zg_herwigs)

  zg_herwig_hist.Scale(1.0 / ( zg_herwig_hist.GetSumOfWeights() * bin_width_herwig ))

  if mc:
    # herwig_plot = rplt.hist(zg_herwig_hist, axes=ax0)
    herwig_plot = ax0.hist(zg_herwigs, label=herwig_label, bins=60, normed=1, histtype='step', color='green', linewidth=5)
  else:
    herwig_plot = ax0.hist(zg_herwigs, bins=60, normed=1, histtype='step', color='green', linewidth=0)
  
  # Herwig Ends.

  # Simulation Plots End.

  
  # Ratio-Over Plot Begins.


  # Theory-Over-Data Plot.
  

  

  data_plot_points_x = []
  data_plot_points_y = []

  for i in range(0, len(data_points_x)):
    if float(data_points_x[i]) > float(zg_cut):
      data_plot_points_x.append(data_points_x[i])
      data_plot_points_y.append(data_points_y[i])

  


  theory_extrapolated_min = []
  theory_extrapolated_line = []
  theory_extrapolated_max = []


  j = 0
  for i in range(0, len(data_plot_points_x)):
    x = float(data_plot_points_x[i])
    
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
      raise ValueError("Some very weird zg value found!", x, theory_x[j], theory_x[j + 1])

    theory_extrapolated_min.append(y_min)
    theory_extrapolated_line.append(y_line)
    theory_extrapolated_max.append(y_max)


  if ratio_denominator == "data":

    if mc:
      zg_herwig_hist.Divide(zg_data_hist)

      a = []
      b = {}
      for i in range(0, len(list(zg_herwig_hist.x()))):
        a.append(round(list(zg_herwig_hist.x())[i] - 0.01 / 2., 4))
        a.append(round(list(zg_herwig_hist.x())[i], 4))
        a.append(round(list(zg_herwig_hist.x())[i] + 0.01 / 2., 4))

        if round(list(zg_herwig_hist.x())[i] - 0.01 / 2., 4) not in b.keys():
          b[round(list(zg_herwig_hist.x())[i] - 0.01 / 2., 4)] = [ float("{0:.3f}".format(list(zg_herwig_hist.y())[i])) ]
        else:
          # print "Bing", round(list(zg_herwig_hist.x())[i] - 0.01 / 2., 4)
          b[round(list(zg_herwig_hist.x())[i] - 0.01 / 2., 4)].append(float("{0:.3f}".format(list(zg_herwig_hist.y())[i])))
        
        if round(list(zg_herwig_hist.x())[i], 4) not in b.keys():
          b[round(list(zg_herwig_hist.x())[i], 4)] = [ float("{0:.3f}".format(list(zg_herwig_hist.y())[i])) ]
        else:
          # print "Bing", round(list(zg_herwig_hist.x())[i], 4)
          b[round(list(zg_herwig_hist.x())[i], 4)].append(float("{0:.3f}".format(list(zg_herwig_hist.y())[i])))

        if round(list(zg_herwig_hist.x())[i] + 0.01 / 2., 4) not in b.keys():
          b[round(list(zg_herwig_hist.x())[i] + 0.01 / 2., 4)] = [ float("{0:.3f}".format(list(zg_herwig_hist.y())[i])) ]
        else:
          # print "Bing", round(list(zg_herwig_hist.x())[i] + 0.01 / 2., 4)
          b[round(list(zg_herwig_hist.x())[i] + 0.01 / 2., 4)].append(float("{0:.3f}".format(list(zg_herwig_hist.y())[i])))

      
      


      
      x = sorted(list(Set(a)))
      a.sort()
      
      c = []

      for i in range(0, len(x)):
        c.append(b[x[i]])
      
      y = [item for sublist in c for item in sublist]

      # a = a[0:len(a) - 50]
      # y = y[0:len(y) - 50]

      a_zero_removed = []
      y_zero_removed = []
      for i in range(0, len(a)):
        if a[i] > zg_cut and a[i] < 0.5:
          a_zero_removed.append(a[i])
          y_zero_removed.append(y[i])




      # plt.hist(list(zg_herwig_hist.x()) , histtype='step', bins=60, weights=list(zg_herwig_hist.y()), axes=ax1, color='green', linewidth=5)

      plt.plot(a_zero_removed, y_zero_removed, linewidth=5, color='green')


      zg_pythia_hist.Divide(zg_data_hist)
      # plt.hist(list(zg_pythia_hist.x()), histtype='step', bins=60, weights=list(zg_pythia_hist.y()), axes=ax1, color='blue', linewidth=5)

    if data:
      ratio_data_to_data = [None if n == 0 else m / n for m, n in zip(data_plot_points_y, data_plot_points_y)]
      data_to_data_y_err = [b / m for b, m in zip(data_y_errors, data_plot_points_y)]
      plt.errorbar(data_plot_points_x, ratio_data_to_data, yerr=data_to_data_y_err, ls='None', marker='o', markersize=8, pickradius=15, elinewidth=15, color='black')
   
    if theory:
      ratio_theory_line_to_data = [m / n for m, n in zip(theory_extrapolated_line, data_plot_points_y)]
      ratio_theory_min_to_data = [m / n for m, n in zip(theory_extrapolated_min, data_plot_points_y)]
      ratio_theory_max_to_data = [m / n for m, n in zip(theory_extrapolated_max, data_plot_points_y)]

      ax1.plot(data_plot_points_x, ratio_theory_line_to_data, alpha=1.0, color='red', linewidth=5)

      ax1.fill_between(data_plot_points_x, ratio_theory_max_to_data, ratio_theory_min_to_data, norm=1, where=np.less_equal(ratio_theory_min_to_data, ratio_theory_max_to_data), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)
      
  elif ratio_denominator == "theory":

    zg_theory_line_hist = Hist(60, 0.0, 0.6, color='red')
    map(zg_theory_line_hist.Fill, data_plot_points_x, theory_extrapolated_line)

    zg_theory_min_hist = Hist(60, 0.0, 0.6, color='pink')
    map(zg_theory_min_hist.Fill, data_plot_points_x, theory_extrapolated_min)

    zg_theory_max_hist = Hist(60, 0.0, 0.6, color='red')
    map(zg_theory_max_hist.Fill, data_plot_points_x, theory_extrapolated_max)


    if mc:
      zg_herwig_hist.Divide(zg_theory_line_hist)
      plt.hist(list(zg_herwig_hist.x()), histtype='step', bins=60, weights=list(zg_herwig_hist.y()), axes=ax1, color='green', linewidth=5)

      zg_pythia_hist.Divide(zg_theory_line_hist)
      plt.hist(list(zg_pythia_hist.x()), histtype='step', bins=60, weights=list(zg_pythia_hist.y()), axes=ax1, color='blue', linewidth=5)

    if data:
      zg_data_to_th_y = [b / m for b, m in zip(data_plot_points_y, theory_extrapolated_line)]
      zg_data_to_th_y_err = [b / m for b, m in zip(data_y_errors, theory_extrapolated_line)]

      plt.errorbar(data_plot_points_x, zg_data_to_th_y, yerr=zg_data_to_th_y_err, ls='None', marker='o', markersize=8, pickradius=15, elinewidth=15, color='black')
   
    if theory:
      
      zg_theory_min_hist.Divide(zg_theory_line_hist)
      zg_theory_max_hist.Divide(zg_theory_line_hist)
      ax1.fill_between(list(zg_theory_max_hist.x()), list(zg_theory_max_hist.y()), list(zg_theory_min_hist.y()), norm=1, where=np.less_equal(list(zg_theory_min_hist.y()), list(zg_theory_max_hist.y())), facecolor='red', interpolate=True, alpha=0.2, linewidth=0.0)

      zg_theory_line_hist.Divide(zg_theory_line_hist)
      plt.hist(list(zg_theory_line_hist.x()), histtype='step', bins=60, weights=list(zg_theory_line_hist.y()), axes=ax1, linewidth=5, color='red')

  else:
    raise ValueError("Only 'theory' or 'data' are valid options for calculating ratios!")


  # Normalized-Over-Data Plot Ends.

  ax0.set_xlabel("$z_g$", fontsize=85)
  ax0.set_ylabel("$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$     ", fontsize=85, rotation=0)
  
  ax1.set_xlabel("$z_g$", fontsize=85)
  ax1.set_ylabel("Ratio           \nto           \n" + ratio_denominator.capitalize() + "           ", fontsize=45, rotation=0)

  ax0.tick_params(axis='x', labelsize=60)
  ax0.tick_params(axis='y', labelsize=60)

  ax1.tick_params(axis='x', labelsize=60)
  ax1.tick_params(axis='y', labelsize=60)

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

  first_legend = ax0.legend(handles, labels, handler_map = {th_line : HandlerLine2D(marker_pad = 0)}, frameon=0, borderpad=0.1)
  ax = ax0.add_artist(first_legend)

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra, extra]
  labels = ["Anti-$k_T$: $R = 0.5$; $p_{T}$ > " + str(pT_lower_cut) + " GeV", "Soft Drop: $\\beta$ = 0; $z_{\mathrm{cut}}$ = " + str(zg_cut)]
  ax0.legend(handles, labels, loc=2, frameon=0, borderpad=0.1)


  # Legend Ends.



  ax0.autoscale(True)
  ax1.autoscale(True)
  
  ax0.set_ylim(0, 10)
  ax1.set_ylim(0.5, 1.5)

  ax0.set_xlim(0.0, 0.6)
  ax1.set_xlim(0.0, 0.6)
  

  fig = plt.gcf()

  if data:
    fn = get_sample_data("/home/aashish/CMS/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)
    ab = AnnotationBbox(OffsetImage(read_png(fn), zoom=1.0), (0.185, 0.97), boxcoords=("figure fraction"), frameon=0)
    ax0.add_artist(ab)

  if data:
    preliminary_text = "Preliminary \n(25% sample)"
  else:
    preliminary_text = "Preliminary"

  fig.text(0.255, 0.95, preliminary_text, fontsize=40, weight='bold', color='#444444', multialignment='center')



  fig = plt.gcf()
  fig.set_size_inches(30, 30, forward=1)

  ax0.xaxis.set_tick_params(width=5, length=10)
  ax0.yaxis.set_tick_params(width=5, length=10)

  ax1.xaxis.set_tick_params(width=5, length=10)
  ax1.yaxis.set_tick_params(width=5, length=10)

  fig.set_snap(True)

  plt.savefig("plots/zg_distribution_data_mc_th_pt_cut_" + str(pT_lower_cut) + "_ratio_over_" + ratio_denominator + "_th_" + str(theory) + "_mc_" + str(mc) + "_data_" + str(data) + ".pdf")



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
        for i in range(1, len(numbers)):
          points[float(numbers[0])].append(float(numbers[i]))
      except ValueError:
        pass
    except:
      pass

  return points




  






# plot_turn_on_curves()


# plot_zg_th_mc_data('0.1', 'zg_1', 'theory', theory=1, mc=0, data=0)
# plot_zg_th_mc_data('0.1', 'zg_1', 'theory', theory=1, mc=1, data=0)
# plot_zg_th_mc_data('0.1', 'zg_1', 'theory', theory=1, mc=1, data=1)

plot_zg_th_mc_data('0.1', 'zg_1', 'data', theory=1, mc=1, data=1)


# plot_pts()



# plot_dr()
# plot_mu()


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


# plot_2d_hist()
