# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat

import sys
from collections import defaultdict

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

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

    except:
      pass

  return properties
    


def plot_pts():
  properties = parse_file(input_analysis_file)

  pTs = properties['uncorrected_hardest_pts']
  corrected_pTs = properties['corrected_hardest_pts']
  prescales = properties['prescales']

  err = np.exp(- np.array(prescales))


  plt.hist(np.array(pTs), 1000, normed=1, label="Uncorrected", weights=prescales, log=1, histtype='step', color='blue') 
  plt.hist(np.array(corrected_pTs), 1000, normed=1, label="Corrected", weights=prescales, log=1, histtype='step', color='red') 


  plt.autoscale(True)
  plt.xlim(0, 1000)

  plt.legend()
  plt.xlabel('$p_T$ (GeV)')
  plt.suptitle("$p_T$ (GeV) Spectrum of anti-kT Jets (R = 0.5)")
  plt.grid(True)

  plt.savefig("plots/ak5_pt_distribution.pdf")
  plt.show()


def plot_zg():
  pT_lower_cut = 150
  pfc_pT_cut = 50
  properties = parse_file(input_analysis_file, pT_lower_cut, pfc_pT_cut)

  zgs = [properties['zg_05'], properties['zg_1'], properties['zg_2']]
  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['$z_{cut}$ = 0.05', '$z_{cut}$ = 0.1', '$z_{cut}$ = 0.2']

  i = 0
  for zg in zgs:
    plt.hist(np.array(zg), 100, normed=1, label=labels[i], weights=prescales, facecolor=colors[i], histtype='step')
    i += 1

  plt.autoscale(True)
  # plt.xlim(0.01, 0.5)

  plt.legend()
  plt.xlabel("Symmetry Measure(z)")
  plt.suptitle("Symmetry Measure(z) with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV")
  plt.grid(True)

  plt.savefig("plots/zg_distribution.pdf")
  plt.show()

def plot_dr():
  pT_lower_cut = 153
  properties = parse_file(input_analysis_file, pT_lower_cut)

  drs = [properties['dr_05'], properties['dr_1'], properties['dr_2']]
  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['$z_{cut}$ = 0.05', '$z_{cut}$ = 0.1', '$z_{cut}$ = 0.2']

  i = 0
  for dr in drs:
    plt.hist(np.array(dr), 200, normed=1, label=labels[i], weights=prescales, facecolor=colors[i], histtype='step')
    i += 1

  plt.autoscale(True)
  # plt.xlim(0.0, 0.5)

  plt.legend()
  plt.xlabel("$\Delta$R between Subjets")
  plt.suptitle("$\Delta$R between Subjets with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV")
  plt.grid(True)

  plt.savefig("plots/delta_r_distribution.pdf")
  plt.show()

def plot_mu():
  pT_lower_cut = 153
  properties = parse_file(input_analysis_file, pT_lower_cut)

  mus = [properties['mu_05'], properties['mu_1'], properties['mu_2']]
  prescales = properties['prescales']

  colors = ['red', 'blue', 'green']
  labels = ['$z_{cut}$ = 0.05', '$z_{cut}$ = 0.1', '$z_{cut}$ = 0.2']

  i = 0
  for mu in mus:
    plt.hist(np.array(mu), 200, normed=1, label=labels[i], weights=prescales, facecolor=colors[i], histtype='step')
    i += 1

  plt.autoscale(True)
  plt.xlim(0.00, 1)

  plt.legend()
  plt.xlabel("Mass Drop($\mu$)")
  plt.suptitle("Mass Drop($\mu$) with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV")
  plt.grid(True)

  plt.savefig("plots/mass_drop_distribution.pdf")
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
    plt.hist(np.array(pdgid_pts[pdgid]), 100, normed=1, label="pdgId = " + str(pdgid), weights=pdgid_prescales[pdgid], log=1, histtype='step')

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
    plt.hist(np.array(pdgid_zg_05s[pdgid]), 100, normed=1, label=labels[0], weights=pdgid_prescales[pdgid], log=1, histtype='step')
    plt.hist(np.array(pdgid_zg_1s[pdgid]), 100, normed=1, label=labels[1], weights=pdgid_prescales[pdgid], log=1, histtype='step')
    plt.hist(np.array(pdgid_zg_2s[pdgid]), 100, normed=1, label=labels[2], weights=pdgid_prescales[pdgid], log=1, histtype='step')

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
    if (pdgids[i] != 1) and (pdgids[i] != 2) and (pdgids[i] != 22):
      pdgid_pts.append(pTs[i])
      pdgid_prescales.append(prescales[i])


  
  plt.hist(np.array(pdgid_pts), 100, normed=1, weights=pdgid_prescales, log=1, histtype='step')

  plt.autoscale(True)

  plt.xlabel('$p_{T}$ GeV')
  plt.suptitle("Hardest Charged pT Distribution")
  plt.grid(True)

  plt.savefig("plots/hardest_charged_pt_distribution.pdf")
  plt.show()


def plot_charged_zgs():
  pT_lower_cut = 153
  properties = parse_file(input_analysis_file, pT_lower_cut)

  pdgids = properties['hardest_pfc_pdgid']
  zg_05s = properties['zg_05']
  zg_1s = properties['zg_1']
  zg_2s = properties['zg_2']
  prescales = properties['prescales']

  pdgid_zgs = [[], [], []]

  pdgid_prescales = []

  for i in range(0, len(pdgids)):
    if (pdgids[i] != 1) and (pdgids[i] != 2) and (pdgids[i] != 22):
      pdgid_zgs[0].append(zg_05s[i])
      pdgid_zgs[1].append(zg_1s[i])
      pdgid_zgs[2].append(zg_2s[i])

      pdgid_prescales.append(prescales[i])


  colors = ['red', 'blue', 'green']
  labels = ['$z_{cut}$ = 0.05', '$z_{cut}$ = 0.1', '$z_{cut}$ = 0.2']

  i = 0
  for zg in pdgid_zgs:
    plt.hist(np.array(zg), 100, normed=1, label=labels[i], weights=pdgid_prescales, facecolor=colors[i], histtype='step')
    i += 1

  plt.autoscale(True)
  plt.xlim(0.01, 0.5)

  plt.legend()
  plt.xlabel("Symmetry Measure(z)")
  plt.suptitle("Symmetry Measure(z) with $p_{T cut}$ = " + str(pT_lower_cut) + " GeV for charged PFCs")
  plt.grid(True)

  plt.savefig("plots/zg_distribution_charged.pdf")
  plt.show()


plot_pts()
# plot_zg()
# plot_dr()
# plot_mu()

# plot_pdgid_pt()

# plot_pdgid_zg()

# plot_charged_pt()

# plot_charged_zgs()