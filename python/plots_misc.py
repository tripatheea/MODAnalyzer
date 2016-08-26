
from __future__ import division

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


def parse_weights(input_filename):
	
	weights = []

	with open(input_filename) as f:
		for line in f:
			
			components = line.split(";")
			
			try:	
				weights.append(float(components[3]))
			except ValueError as e:
				pass
				
	return weights


def plot_weights(input_filename):
	weights = parse_weights(input_filename)

	hist = Hist(np.logspace(math.log(float(1e-8), math.e), math.log(1, math.e), 25, base=np.e))

	hist.fill_array(weights)

	# plt.hist(weights, bins=200, normed=1)
	rplt.hist(hist)

	print max(weights)


	# print "Total number of events = ", len(weights)

	plt.xscale("log")
	plt.yscale("log")

	plt.autoscale()

	plt.xlabel("Weight")
	plt.ylabel("# of Events")

	plt.savefig("plots/weights.pdf")
	plt.clf()


# plot_weights(sys.argv[1])

