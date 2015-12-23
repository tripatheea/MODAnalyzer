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



def parse_file(input_file, pT_lower_cut=150., pT_upper_cut=20000., softdrop_pT_lower_cut=0., softdrop_pT_upper_cut=20000.):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  # FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3
  
  properties = defaultdict(list)

  keywords = []
  keywords_set = False
  for line in lines:
    try:
      numbers = line.split()

      if numbers[0] == "#" and (not keywords_set):
        keywords = numbers[2:]
        keywords_set = True
      elif numbers[0] == "Entry":
        pT_index = keywords.index("hardest_pT") + 1
        softdrop_pT_index = keywords.index("pT_after_SD") + 1

        if float(numbers[pT_index]) > pT_lower_cut and float(numbers[pT_index]) < pT_upper_cut and float(numbers[softdrop_pT_index]) > softdrop_pT_lower_cut and float(numbers[softdrop_pT_index]) < softdrop_pT_upper_cut:
          for i in range(len(keywords)):
            keyword = keywords[i]
            properties[keyword].append( float(numbers[i + 1]) ) # + 1 because we ignore the first keyword "Entry".
        else:
        	for i in range(len(keywords)):
						keyword = keywords[i]
						properties[keyword].append( -1.0 )

    except:
      pass

  return properties

def plot_pt(pT_lower_cut):

	properties_reco = [parse_file("/home/aashish/pythia_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_reco.dat", pT_lower_cut=pT_lower_cut)]
	properties_truth = [parse_file("/home/aashish/pythia_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_truth.dat", pT_lower_cut=pT_lower_cut)]
	labels = ["pythia", "herwig", "sherpa"]


	for prop_reco, prop_truth, label in zip(properties_reco, properties_truth, labels):
	
		x = prop_truth['hardest_pT']
		y = prop_reco['hardest_pT']


		H, xedges, yedges = np.histogram2d(x, y, bins=200, normed=1 )

		H = np.rot90(H)
		H = np.flipud(H)

		Hmasked = np.ma.masked_where(H == 0, H) # Mask pixels with a value of zero

		plt.pcolormesh(xedges,yedges, Hmasked)

		cbar = plt.colorbar()
		cbar.ax.set_ylabel('Counts')

		plt.xlim(150, max(max(x), max(y)))
		plt.ylim(150, max(max(x), max(y)))

		plt.xlabel('Truth $p_T / \mathrm{GeV}$', fontsize=50, labelpad=75)
		plt.ylabel('Reco $p_T / \mathrm{GeV}$', fontsize=50, labelpad=75)

		plt.gcf().set_size_inches(30, 30, forward=1)
		plt.gcf().set_snap(True)

		plt.savefig("plots/With MC/2D/pT_" + label + ".pdf")


		plt.clf()


def plot_jet_multiplicity(pT_lower_cut, charged=False):

	properties_reco = [parse_file("/home/aashish/pythia_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_reco.dat", pT_lower_cut=pT_lower_cut)]
	properties_truth = [parse_file("/home/aashish/pythia_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_truth.dat", pT_lower_cut=pT_lower_cut)]
	labels = ["pythia", "herwig", "sherpa"]

	for prop_reco, prop_truth, label in zip(properties_reco, properties_truth, labels):
		
		if charged:
			param = 'chrg_mul_pre_SD'
		else:
			param = 'mul_pre_SD'

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

	properties_reco = [parse_file("/home/aashish/pythia_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_reco.dat", pT_lower_cut=pT_lower_cut)]
	properties_truth = [parse_file("/home/aashish/pythia_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_truth.dat", pT_lower_cut=pT_lower_cut)]
	labels = ["pythia", "herwig", "sherpa"]

	for prop_reco, prop_truth, label in zip(properties_reco, properties_truth, labels):
		
		if charged:
			param = 'chrg_mass_pre_SD'
		else:
			param = 'mass_pre_SD'

		x = prop_truth[param]
		y = prop_reco[param]

		print len(x), len(y)

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

	properties_reco = [parse_file("/home/aashish/pythia_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_reco.dat", pT_lower_cut=pT_lower_cut)]
	properties_truth = [parse_file("/home/aashish/pythia_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_truth.dat", pT_lower_cut=pT_lower_cut)]
	labels = ["pythia", "herwig", "sherpa"]

	for prop_reco, prop_truth, label in zip(properties_reco, properties_truth, labels):
	
		x = prop_truth['hardest_eta']
		y = prop_reco['hardest_eta']

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




def plot_zg(pT_lower_cut):

	properties_reco = [parse_file("/home/aashish/pythia_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_reco.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_reco.dat", pT_lower_cut=pT_lower_cut)]
	properties_truth = [parse_file("/home/aashish/pythia_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/herwig_truth.dat", pT_lower_cut=pT_lower_cut), parse_file("/home/aashish/sherpa_truth.dat", pT_lower_cut=pT_lower_cut)]
	labels = ["pythia", "herwig", "sherpa"]


	for prop_reco, prop_truth, label in zip(properties_reco, properties_truth, labels):

	
		x = prop_truth['zg_05']
		y = prop_reco['zg_05']



		H, xedges, yedges = np.histogram2d(x, y, bins=50, normed=1 )
		

		H = np.rot90(H)
		H = np.flipud(H)

		Hmasked = np.ma.masked_where(H == 0, H) # Mask pixels with a value of zero

		plt.pcolormesh(xedges,yedges, Hmasked)

		cbar = plt.colorbar()
		cbar.ax.set_ylabel('Counts')

		plt.xlim(0, 0.5)
		plt.ylim(0, 0.5)

		plt.xlabel('Truth $z_g$', fontsize=50, labelpad=75)
		plt.ylabel('Reco $z_g$', fontsize=50, labelpad=75)

		plt.gcf().set_size_inches(30, 30, forward=1)
		plt.gcf().set_snap(True)

		plt.savefig("plots/With MC/2D/zg_05_" + str(label) + ".pdf")
		# plt.show()

		plt.clf()


plot_zg(150)


plot_pt(150)

plot_jet_multiplicity(150, charged=True)
plot_jet_multiplicity(150, charged=False)

plot_jet_mass(150, charged=True)
plot_jet_mass(150, charged=False)

plot_eta(150)