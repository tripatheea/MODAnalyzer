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
from scipy import arange, array, exp

from scipy.stats import binned_statistic

import rootpy.plotting.views



mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)


def parse_mc(input_file, pT_lower_cut=150.00, pT_upper_cut = 20000.00, uncorrected_pT_lower_cut = 0.00, softdrop_unc_pT_lower_cut = 0.00, softdrop_cor_pT_lower_cut = 0.00, jet_quality_level=1):
  f = open(input_file, 'r')
  lines = f.read().split("\n")
  
  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":


        properties['hardest_pts'].append( float( numbers[1] ) )
        properties['zg_05'].append( float( numbers[2] ) )
        properties['dr_05'].append( float( numbers[3] ) )
        properties['mu_05'].append( float( numbers[4] ) )
        properties['zg_1'].append( float( numbers[5] ) )
        properties['dr_1'].append( float( numbers[6] ) )
        properties['mu_1'].append( float( numbers[7] ) )
        properties['zg_2'].append( float( numbers[8] ) )
        properties['dr_2'].append( float( numbers[9] ) )
        properties['mu_2'].append( float( numbers[10] ) )
        
        properties['zg_05_pt_1'].append( float( numbers[11] ) )
        properties['zg_1_pt_1'].append( float( numbers[12] ) )
        properties['zg_2_pt_1'].append( float( numbers[13] ) )

        properties['zg_05_pt_2'].append( float( numbers[14] ) )
        properties['zg_1_pt_2'].append( float( numbers[15] ) )
        properties['zg_2_pt_2'].append( float( numbers[16] ) )

        properties['zg_05_pt_3'].append( float( numbers[17] ) )
        properties['zg_1_pt_3'].append( float( numbers[18] ) )
        properties['zg_2_pt_3'].append( float( numbers[19] ) )

        properties['zg_05_pt_5'].append( float( numbers[20] ) )
        properties['zg_1_pt_5'].append( float( numbers[21] ) )
        properties['zg_2_pt_5'].append( float( numbers[22] ) )

        properties['zg_05_pt_10'].append( float( numbers[23] ) )
        properties['zg_1_pt_10'].append( float( numbers[24] ) )
        properties['zg_2_pt_10'].append( float( numbers[25] ) )

        properties['zg_charged_05'].append( float( numbers[26] ) )
        properties['zg_charged_1'].append( float( numbers[27] ) )
        properties['zg_charged_2'].append( float( numbers[28] ) )

        properties['pTs_after_SD'].append( float( numbers[29] ) )

        properties['multiplicity_before_SD'].append( float( numbers[30] ) )
        properties['multiplicity_after_SD'].append( float( numbers[31] ) )
        properties['jet_mass_before_SD'].append( float( numbers[32] ) )
        properties['jet_mass_after_SD'].append( float( numbers[33] ) )

        properties['charged_multiplicity_before_SD'].append( float( numbers[34] ) )
        properties['charged_multiplicity_after_SD'].append( float( numbers[35] ) )
        properties['charged_jet_mass_before_SD'].append( float( numbers[36] ) )
        properties['charged_jet_mass_after_SD'].append( float( numbers[37] ) )

        properties['chrg_dr_05'].append( float( numbers[38] ) )
        properties['chrg_dr_1'].append( float( numbers[39] ) )
        properties['chrg_dr_2'].append( float( numbers[40] ) )

        properties['fractional_energy_loss'].append( float( numbers[41] ) )
        properties['eta'].append( float( numbers[42] ) )
        # properties['jet_area'].append( float( numbers[43] ) )


          
    except:
      if len(numbers) != 0:
        # print "Some kind of error occured while parsing the given file!"
        # print numbers
        # print
        pass

  f.close()

  return properties


def plot_pt(pT_lower_cut):

	properties_reco = [parse_mc("/home/aashish/pythia_reco.dat"), parse_mc("/home/aashish/herwig_reco.dat"), parse_mc("/home/aashish/sherpa_reco.dat")]
	properties_truth = [parse_mc("/home/aashish/pythia_truth.dat"), parse_mc("/home/aashish/herwig_truth.dat"), parse_mc("/home/aashish/sherpa_truth.dat")]
	labels = ["pythia", "herwig", "sherpa"]

	for prop_reco, prop_truth, label in zip(properties_reco, properties_truth, labels):
	
		x = prop_truth['hardest_pts']
		y = prop_reco['hardest_pts']

		H, xedges, yedges = np.histogram2d(x, y, bins=200, normed=1 )

		H = np.rot90(H)
		H = np.flipud(H)

		Hmasked = np.ma.masked_where(H == 0, H) # Mask pixels with a value of zero

		plt.pcolormesh(xedges,yedges, Hmasked)

		cbar = plt.colorbar()
		cbar.ax.set_ylabel('Counts')

		plt.xlim(0, 200)
		plt.ylim(0, 200)

		plt.xlabel('Truth $p_T / \mathrm{GeV}$', fontsize=50, labelpad=75)
		plt.ylabel('Reco $p_T / \mathrm{GeV}$', fontsize=50, labelpad=75)

		plt.gcf().set_size_inches(30, 30, forward=1)
		plt.gcf().set_snap(True)

		plt.savefig("plots/With MC/2D/pT_" + label + ".pdf")

		plt.clf()


def plot_jet_multiplicity(pT_lower_cut, charged=False):

	properties_reco = [parse_mc("/home/aashish/pythia_reco.dat"), parse_mc("/home/aashish/herwig_reco.dat"), parse_mc("/home/aashish/sherpa_reco.dat")]
	properties_truth = [parse_mc("/home/aashish/pythia_truth.dat"), parse_mc("/home/aashish/herwig_truth.dat"), parse_mc("/home/aashish/sherpa_truth.dat")]
	labels = ["pythia", "herwig", "sherpa"]

	for prop_reco, prop_truth, label in zip(properties_reco, properties_truth, labels):
		
		if charged:
			param = 'charged_multiplicity_before_SD'
		else:
			param = 'multiplicity_before_SD'

		x = prop_truth[param]
		y = prop_reco[param]

		H, xedges, yedges = np.histogram2d(x, y, bins=25, normed=1 )

		H = np.rot90(H)
		H = np.flipud(H)

		Hmasked = np.ma.masked_where(H == 0, H) # Mask pixels with a value of zero

		plt.pcolormesh(xedges,yedges, Hmasked)

		cbar = plt.colorbar()
		cbar.ax.set_ylabel('Counts')

		plt.xlim(0, 35)
		plt.ylim(0, 35)

		plt.xlabel('Truth Jet Multiplicity', fontsize=50, labelpad=75)
		plt.ylabel('Reco Jet Multiplicity', fontsize=50, labelpad=75)

		plt.gcf().set_size_inches(30, 30, forward=1)
		plt.gcf().set_snap(True)

		if charged:
			plt.savefig("plots/With MC/2D/charged_jet_multiplicity_" + label + ".pdf")
		else:
			plt.savefig("plots/With MC/2D/jet_multiplicity_" + label + ".pdf")

		plt.clf()



def plot_jet_mass(pT_lower_cut, charged=False):

	properties_reco = [parse_mc("/home/aashish/pythia_reco.dat"), parse_mc("/home/aashish/herwig_reco.dat"), parse_mc("/home/aashish/sherpa_reco.dat")]
	properties_truth = [parse_mc("/home/aashish/pythia_truth.dat"), parse_mc("/home/aashish/herwig_truth.dat"), parse_mc("/home/aashish/sherpa_truth.dat")]
	labels = ["pythia", "herwig", "sherpa"]

	for prop_reco, prop_truth, label in zip(properties_reco, properties_truth, labels):
		
		if charged:
			param = 'charged_jet_mass_before_SD'
		else:
			param = 'jet_mass_before_SD'

		x = prop_truth[param]
		y = prop_reco[param]

		H, xedges, yedges = np.histogram2d(x, y, bins=100, normed=1 )

		H = np.rot90(H)
		H = np.flipud(H)

		Hmasked = np.ma.masked_where(H == 0, H) # Mask pixels with a value of zero

		plt.pcolormesh(xedges,yedges, Hmasked)

		cbar = plt.colorbar()
		cbar.ax.set_ylabel('Counts')


		plt.xlim(0, 50)
		plt.ylim(0, 50)

		plt.xlabel('Truth Jet Mass / GeV', fontsize=50, labelpad=75)
		plt.ylabel('Reco Jet Mass / GeV', fontsize=50, labelpad=75)

		plt.gcf().set_size_inches(30, 30, forward=1)
		plt.gcf().set_snap(True)

		if charged:
			plt.savefig("plots/With MC/2D/charged_jet_mass_" + label + ".pdf")
		else:
			plt.savefig("plots/With MC/2D/jet_mass_" + label + ".pdf")

		plt.clf()





def plot_eta(pT_lower_cut):

	properties_reco = [parse_mc("/home/aashish/pythia_reco.dat"), parse_mc("/home/aashish/herwig_reco.dat"), parse_mc("/home/aashish/sherpa_reco.dat")]
	properties_truth = [parse_mc("/home/aashish/pythia_truth.dat"), parse_mc("/home/aashish/herwig_truth.dat"), parse_mc("/home/aashish/sherpa_truth.dat")]
	labels = ["pythia", "herwig", "sherpa"]

	for prop_reco, prop_truth, label in zip(properties_reco, properties_truth, labels):
	
		x = prop_truth['eta']
		y = prop_reco['eta']

		H, xedges, yedges = np.histogram2d(x, y, bins=200, normed=1 )

		H = np.rot90(H)
		H = np.flipud(H)

		Hmasked = np.ma.masked_where(H == 0, H) # Mask pixels with a value of zero

		plt.pcolormesh(xedges,yedges, Hmasked)

		cbar = plt.colorbar()
		cbar.ax.set_ylabel('Counts')

		plt.xlim(0, 3)
		plt.ylim(0, 3)

		plt.xlabel('Truth $\eta$', fontsize=50, labelpad=75)
		plt.ylabel('Reco $\eta$', fontsize=50, labelpad=75)

		plt.gcf().set_size_inches(30, 30, forward=1)
		plt.gcf().set_snap(True)

		plt.savefig("plots/With MC/2D/eta_" + label + ".pdf")

		plt.clf()



# plot_pt(150)

plot_jet_multiplicity(150, charged=True)
plot_jet_multiplicity(150, charged=False)

# plot_jet_mass(150, charged=True)
# plot_jet_mass(150, charged=False)

# plot_eta(150)