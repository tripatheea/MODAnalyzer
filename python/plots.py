from __future__ import division

from subprocess import call

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FormatStrFormatter

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


from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

from scipy.stats import norm
from scipy.stats import gamma
from scipy import arange, array, exp

from scipy.stats import binned_statistic

import rootpy.plotting.views

import hists


input_analysis_file = sys.argv[1]


mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)



plot_labels = { "data": "CMS 2010 Open Data", "pythia": "Pythia 8.215", "herwig": "Herwig 7", "sherpa": "Sherpa 2.2.0", "theory": "Theory (MLL)" }
plot_colors = {"theory": "red", "pythia": "blue", "herwig": "green", "sherpa": "purple", "pythia_post": "red", "data": "black", "data_post": "red"}



def parse_file(input_file, pT_lower_cut=150., pT_upper_cut=20000., softdrop_pT_lower_cut=0., softdrop_pT_upper_cut=20000., eta_cut=2.4):
	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	all_hists = hists.all_hist_templates()

	keywords = []
	keywords_set = False

	with open(input_file) as infile:
		for line in infile:
			try:
				numbers = line.split()

				if numbers[0] == "#" and (not keywords_set):
					keywords = numbers[2:]
					keywords_set = True
				elif numbers[0] == "Entry":
					pT_index = keywords.index("hardest_pT") + 1
					softdrop_pT_index = keywords.index("pT_after_SD") + 1

					eta_index = keywords.index("hardest_eta") + 1

					prescale_index = keywords.index("prescale") + 1

					if abs(float(numbers[eta_index])) < eta_cut and float(numbers[pT_index]) > pT_lower_cut and float(numbers[pT_index]) < pT_upper_cut and float(numbers[softdrop_pT_index]) > softdrop_pT_lower_cut and float(numbers[softdrop_pT_index]) < softdrop_pT_upper_cut:
						
						

						for i in range(len(keywords)):
							keyword = keywords[i]

							# Get first half of the keyword (anything before the first '_').
							key = keyword.split('_')[0]

							


							
							all_hists[key].fill_array( float(numbers[i + 1]), float(numbers[prescale_index]) ) # + 1 because we ignore the first keyword "Entry".

			except:
				pass

	return all_hists







input_analysis_file = sys.argv[1]

print parse_file(input_analysis_file)['eta'].nbins()
