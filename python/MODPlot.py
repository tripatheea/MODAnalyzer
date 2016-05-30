from __future__ import division

from subprocess import call

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FormatStrFormatter

from sets import Set

import time as time

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


from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox, AnchoredOffsetbox, HPacker

from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea

from scipy.stats import norm
from scipy.stats import gamma
from scipy import arange, array, exp

from scipy.stats import binned_statistic

import rootpy.plotting.views


logo_location = "/home/aashish/root/macros/MODAnalyzer/mod_logo.png"
logo_text = "Prelim. (20\%)"






mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)





class MODPlot:

	def __init__(self, hists, plot_types, plot_colors, plot_labels, x_label, y_label, x_lims=(0, -1), y_lims=(0, -1)):
		
		self._hists = hists
		self._plot_types = plot_types
		self._plot_colors = plot_colors
		self._plot_labels = plot_labels

		self._x_label = x_label
		self._y_label = y_label

		self._x_lims = x_lims
		self._y_lims = y_lims
		



	def set_logo(self):
		
		logo_offset_image = OffsetImage(read_png(get_sample_data(logo_location, asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
		text_box = TextArea(logo_text, textprops=dict(color='#444444', fontsize=50, weight='bold'))

		logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

		anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.)

		plt.gca().add_artist(anchored_box)
		

	def normalize_hists(self):

		for i in range(len(self._hists)):
			bin_width = (self._hists[i].upperbound() - self._hists[i].lowerbound()) / self._hists[i].nbins()
			self._hists[i].Scale(1.0 / ( self._hists[i].GetSumOfWeights() * bin_width ))

	def set_formatting(self):
		for i in range(len(self._hists)):
			self._hists[i].SetColor(self._plot_colors[i])
			self._hists[i].SetTitle(self._plot_labels[i])
			self._hists[i].SetLineWidth(8)


	def plot(self, filename):
		
		self.set_logo()
		self.set_formatting()

		self.normalize_hists()


		z_indices = range(len(self._hists), 0, -1)
		z_indices[0] *= 10

		for i in range(len(self._hists)):
			hist = self._hists[i]
			plot_type = self._plot_types[i]

			if plot_type == 'hist':
				rplt.hist(hist, zorder=z_indices[i], emptybins=False)
			elif plot_type == 'error':
				rplt.errorbar(hist, zorder=z_indices[i], emptybins=False, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)


		plt.legend(frameon=False)

		plt.autoscale()

		if self._x_lims[1] == -1:
			plt.gca().set_xlim( self._x_lims[0], plt.gca().get_xlim()[1] )
		else:
			plt.gca().set_xlim( self._x_lims[0], self._x_lims[1] )

		if self._y_lims[1] == -1:
			plt.gca().set_ylim( self._y_lims[0], plt.gca().get_ylim()[1] * 1.125 )
		else:
			plt.gca().set_ylim( self._y_lims[0], self._y_lims[1] )


		plt.xlabel(self._x_label)
		plt.ylabel(self._y_label)

		plt.tick_params(which='major', width=5, length=25, labelsize=70)
		plt.tick_params(which='minor', width=3, length=15)

		plt.gcf().set_size_inches(30, 24, forward=1)
		plt.gcf().set_snap(True)

		plt.savefig(filename)

		plt.clf()


