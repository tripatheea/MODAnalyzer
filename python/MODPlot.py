
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


logo_location = "/home/aashish/root/macros/MODAnalyzer/mod_logo.png"
logo_text = "Prelim. (20\%)"






mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)





class MODPlot:

	def __init__(self, hists, plot_types, plot_colors, plot_labels, line_styles, x_scale='linear', y_scale='linear', ratio_plot=False, ratio_to_index=-1, ratio_label="", x_label="", y_label="", x_lims=(0, -1), y_lims=(0, -1)):
		
		self._hists = hists
		self._plot_types = plot_types
		self._plot_colors = plot_colors
		self._plot_labels = plot_labels
		self._line_styles = line_styles

		self._x_scale = x_scale
		self._y_scale = y_scale

		self._ratio_plot = ratio_plot
		self._ratio_to_index = ratio_to_index
		self._ratio_label = ratio_label

		self._x_label = x_label
		self._y_label = y_label

		self._x_lims = x_lims
		self._y_lims = y_lims

		self._plt = plt
		self._rplt = rplt
		



	def logo_box(self):
		
		logo_offset_image = OffsetImage(read_png(get_sample_data(logo_location, asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
		text_box = TextArea(logo_text, textprops=dict(color='#444444', fontsize=50, weight='bold'))

		logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

		anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.)

		return anchored_box
		

	def normalize_hists(self):

		for i in range(len(self._hists)):
			bin_width = (self._hists[i].hist().upperbound() - self._hists[i].hist().lowerbound()) / self._hists[i].hist().nbins()
			if self._hists[i].hist().GetSumOfWeights() != 0.0:
				self._hists[i].hist().Scale(1.0 / ( self._hists[i].hist().GetSumOfWeights() * bin_width ))

	def set_formatting(self):
		for i in range(len(self._hists)):
			self._hists[i].hist().SetColor(self._plot_colors[i])
			self._hists[i].hist().SetTitle(self._plot_labels[i])
			self._hists[i].hist().SetLineWidth(8)

	def convert_hist_to_line_plot(self, hist, x_range):	# x_range = range of x for which it's non-zero. 

		n_bins = hist.nbins()

		a = []
		b = {}
		bin_width = (x_range[1] - x_range[0]) / n_bins

		for i in range(0, len(list(hist.x()))):

			a.append(round(list(hist.x())[i] - bin_width / 2., 4))
			a.append(round(list(hist.x())[i], 4))
			a.append(round(list(hist.x())[i] + bin_width / 2., 4))

			if round(list(hist.x())[i] - bin_width / 2., 4) not in b.keys():
				b[round(list(hist.x())[i] - bin_width / 2., 4)] = [ list(hist.y())[i] ]
			else:
				b[round(list(hist.x())[i] - bin_width / 2., 4)].append( list(hist.y())[i] )
			
			if round(list(hist.x())[i], 4) not in b.keys():
				b[round(list(hist.x())[i], 4)] = [ list(hist.y())[i] ]
			else:
				b[round(list(hist.x())[i], 4)].append( list(hist.y())[i] )

			if round(list(hist.x())[i] + bin_width / 2., 4) not in b.keys():
				b[round(list(hist.x())[i] + bin_width / 2., 4)] = [ list(hist.y())[i] ]
			else:
				b[round(list(hist.x())[i] + bin_width / 2., 4)].append( list(hist.y())[i] )

		x = sorted(list(Set(a)))
		a.sort()
		
		c = [b[x[i]] for i in range(0, len(x))]

		y = [item for sublist in c for item in sublist]
		
		a_zero_removed = []
		y_zero_removed = []
		for i in range(0, len(a)):
			# if a[i] >= x_range[0] and a[i] <= x_range[1] and y[i] != 0.0:
			if a[i] >= x_range[0] and a[i] <= x_range[1]:
				a_zero_removed.append(a[i])
				y_zero_removed.append(y[i])

		return a_zero_removed, y_zero_removed



	def plot(self):
		
		if self._ratio_plot:
			gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
			ax0 = self._plt.subplot(gs[0])
			ax1 = self._plt.subplot(gs[1])
		else:
			ax0 = self._plt.gca()

		# Set the logo.
		ax0.add_artist(self.logo_box())

		# Set basic plot element formattings. 
		self.set_formatting()

		# Normalize all the histograms.
		self.normalize_hists()


		z_indices = range(len(self._hists), 0, -1)
		z_indices[0] *= 10

		# First, draw the regular "non-ratio" plot.


		legend_handles = []
		for i in range(len(self._hists)):
			hist = self._hists[i].hist()
			plot_type = self._plot_types[i]

			if plot_type == 'hist':
				hist.SetLineStyle(self._line_styles[i])
				plot = self._rplt.hist(hist, axes=ax0, zorder=z_indices[i], emptybins=False)
			elif plot_type == 'error':
				plot = self._rplt.errorbar(hist, axes=ax0, zorder=z_indices[i], emptybins=False, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)

			legend_handles.append( plot )


		if self._y_scale == 'log':
			ax0.set_yscale('log')


		# Ratio plot.

		if self._ratio_plot:

			denominator_hist = self._hists[self._ratio_to_index].hist()

			
			for i in range(len(self._hists)):
				ratio_hist = copy.deepcopy( self._hists[i].hist() )
				ratio_hist.Divide(denominator_hist)

				plot_type = self._plot_types[i]

				if plot_type == 'hist':
					# rplt.hist(ratio_hist, axes=ax1, zorder=z_indices[i], emptybins=False)
					


					line_plot = self.convert_hist_to_line_plot(ratio_hist, ratio_hist.bounds())

					self._plt.plot(line_plot[0], line_plot[1], ls=self._line_styles[i], axes=ax1, lw=8, color=ratio_hist.GetColor()[0])

				elif plot_type == 'error':
					self._rplt.errorbar(ratio_hist, axes=ax1, zorder=z_indices[i], emptybins=False, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)


		# Ratio plot ends.


		handles, labels = legend_handles, self._plot_labels
		legend = ax0.legend(handles, labels, loc=1, frameon=0, fontsize=60)
		ax0.add_artist(legend)

		# Any additional texts.
		extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
		# ax0.legend([extra]*len(labels), labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[1.00, 0.53])

		for position, anchor_location, text in self._hists[0].additional_text():
			texts = text.split("\n")
			ax0.legend( [extra] * len(texts), texts, frameon=0, borderpad=0, fontsize=60, bbox_to_anchor=position, loc=anchor_location)
		
		plt.autoscale()

		if self._x_lims[1] == -1:
			ax0.set_xlim( self._x_lims[0], ax0.get_xlim()[1] )
		else:
			ax0.set_xlim( self._x_lims[0], self._x_lims[1] )

		if self._y_lims[1] == -1:
			ax0.set_ylim( self._y_lims[0], ax0.get_ylim()[1] * 1.125 )
		else:
			ax0.set_ylim( self._y_lims[0], self._y_lims[1] )

		if self._ratio_plot:
			ax1.set_xlim( ax0.get_xlim()[0], ax0.get_xlim()[1] )
			ax1.set_ylim(0., 2.)

		# Axes labels.

		ax0.set_xlabel(self._x_label, fontsize=75)
		ax0.set_ylabel(self._y_label, rotation=0, fontsize=75, labelpad=75)

		if self._ratio_plot:
			ax1.set_xlabel(self._x_label, fontsize=75)
			ax1.set_ylabel(self._ratio_label, rotation=0, fontsize=55, labelpad=115, y=0.31)

		# Axes labels end.


		self._plt.sca(ax0)

		self._plt.tick_params(which='major', width=5, length=25, labelsize=70)
		self._plt.tick_params(which='minor', width=3, length=15)

		if self._ratio_plot:
			self._plt.sca(ax1)
			
			self._plt.tick_params(which='major', width=5, length=25, labelsize=70)
			self._plt.tick_params(which='minor', width=3, length=15)


		if self._ratio_plot:
			self._plt.gcf().set_size_inches(30, 30, forward=1)
		else:
			self._plt.gcf().set_size_inches(30, 24, forward=1)

		


		self._plt.gcf().set_snap(True)

		self._plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

	def save_plot(self, filename):
		self._plt.savefig(filename)
		self._plt.clf()

	def get_plt(self):
		return self._plt

	def get_rplt(self):
		return self._rplt





plot_labels = { "data": "CMS 2010 Open Data", "pythia": "Pythia 8.215", "herwig": "Herwig 7", "sherpa": "Sherpa 2.2.0", "theory": "Theory (MLL)" }
plot_colors = {"theory": "red", "pythia": "blue", "herwig": "green", "sherpa": "purple", "pythia_post": "red", "data": "black", "data_post": "red"}

plot_types = ['error', 'hist', 'hist', 'hist']
colors = [ plot_colors['data'], plot_colors['pythia'], plot_colors['herwig'], plot_colors['sherpa'] ]
labels = [ plot_labels['data'], plot_labels['pythia'], plot_labels['herwig'], plot_labels['sherpa'] ]
line_styles = [ "", "solid", "dashed", "dotted" ]


def create_multi_page_plot(filename, hists):
	# mod_hists is a list of MODHist objects.

	with PdfPages(filename) as pdf:

		for mod_hist in hists:
			plot = MODPlot(mod_hist, plot_types=plot_types, plot_colors=colors, plot_labels=labels, line_styles=line_styles, y_scale=mod_hist[0].y_scale(), ratio_plot=True, ratio_to_index=1, ratio_label="Ratio\nto\nPythia", x_label=mod_hist[0].x_label(), y_label=mod_hist[0].y_label(), x_lims=mod_hist[0].x_range(), y_lims=mod_hist[0].y_range())
			plot.plot()
		
			pdf.savefig()
			plot.get_plt().close()