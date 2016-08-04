from __future__ import division



import hists


from MODPlot import *


import trigger_parse





from subprocess import call

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FormatStrFormatter

from sets import Set

import time as time
import copy


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

import matplotlib.ticker as mtick

from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox, AnchoredOffsetbox, HPacker

from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea

from scipy.stats import norm
from scipy.stats import gamma
from scipy import arange, array, exp

from scipy.stats import binned_statistic

import rootpy.plotting.views






default_dir = "plots/Version 5_2/"




parsed_linear = trigger_parse.load_root_files_to_hist()



def normalize_hist(hist):
	bin_width = (hist.upperbound() - hist.lowerbound()) / hist.nbins()
	
	if hist.GetSumOfWeights() != 0.0:
		hist.Scale(1.0 / ( hist.GetSumOfWeights() * bin_width ))

	return hist

def trigger_efficiency_plot():
	mod_hists = parsed_linear[0]


	print mod_hists

	colors = ['green', 'magenta', 'blue', 'red', 'brown', 'orange']
	labels = ["Jet140U / 100U", "Jet100U / 70U", "Jet70U / 50U", "Jet50U / 30U", "Jet30U / 15U\_HNF", "" ]
	hist_labels = [("Jet140U", "Jet100U"), ("Jet100U", "Jet70U"), ("Jet70U", "Jet50U"), ("Jet50U", "Jet30U"), ("Jet30U", "Jet15U_HcalNoiseFiltered") ]
	lower_pTs = [140, 100, 70, 50, 30, 15]

	# rplt.hist(mod_hists[0].hist())

	for i in range(len(hist_labels) - 1):
		
		first_hist, second_hist = mod_hists[hist_labels[i][0]], mod_hists[hist_labels[i][1]]

		print first_hist
		# rplt.errorbar()

		new_hist = normalize_hist( first_hist.hist() / second_hist.hist() )

		new_hist.SetColor(colors[i])

		rplt.errorbar(new_hist)

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.savefig(default_dir + "trigger_efficiency.pdf")

	plt.clf()




start = time.time()

trigger_efficiency_plot()

end = time.time()

print "Finished all plotting in {} seconds.".format(end - start)
