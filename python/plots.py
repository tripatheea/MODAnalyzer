# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat

import sys
import math
from collections import defaultdict

# matplotlib
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# RootPy
from rootpy.plotting import Hist, HistStack, Legend, Canvas
import rootpy.plotting.root2matplotlib as rplt

import rootpy.plotting.views

input_analysis_file = sys.argv[1]



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
    
def parse_theory_file(input_file, pT_lower_cut = 0.00, pfc_pT_cut = 0.00):
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
  properties_pythia = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_pythia_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)
  properties_herwig = parse_theory_file("/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Andrew/fastjet_sudakov_safe_herwig_pp2jj_" + str(pT_lower_cut) + "pTcut_7TeV.dat", pT_lower_cut, pfc_pT_cut)

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

  for j in range(0, len(zg_pythias)):
    zg_pythia_hist = Hist(50, 0, 0.5, title=pythia_labels[j], markersize=1.0, color=colors[j])

    map(zg_pythia_hist.Fill, zg_pythias[j])

    zg_pythia_hist.Scale(1.0 / zg_pythia_hist.GetSumOfWeights())

    rplt.hist(zg_pythia_hist)

  for j in range(0, len(zg_herwigs)):
    zg_herwig_hist = Hist(50, 0, 0.5, title=herwig_labels[j], markersize=1.0, color=colors_2[j], linestyle="2") # 1=solid, 2=dash, 3=dot, 4=dash-dot

    map(zg_herwig_hist.Fill, zg_herwigs[j])

    zg_herwig_hist.Scale(1.0 / zg_herwig_hist.GetSumOfWeights())

    rplt.hist(zg_herwig_hist)


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
  plt.xlim(0.01, 0.5)

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
  plt.xlabel('zg_05')
  plt.ylabel('Charged zg_05')
  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Counts')

  plt.savefig("plots/zg_vs_charged_zg.pdf")
  plt.show()

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