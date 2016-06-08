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


from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox, AnchoredOffsetbox, HPacker

from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea
from matplotlib.backends.backend_pdf import PdfPages

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

	def __init__(self, hists, plot_types, plot_colors, plot_labels, multi_page=False, x_scale='linear', y_scale='linear', ratio_plot=False, ratio_to_index=-1, ratio_label="", x_label="", y_label=""):
		
		self._hists = hists
		self._plot_types = plot_types
		self._plot_colors = plot_colors
		self._plot_labels = plot_labels

		self._multi_page = multi_page

		self._x_scale = x_scale
		self._y_scale = y_scale

		self._ratio_plot = ratio_plot
		self._ratio_to_index = ratio_to_index
		self._ratio_label = ratio_label

		self._x_label = x_label
		self._y_label = y_label

		



	def logo_box(self):
		
		logo_offset_image = OffsetImage(read_png(get_sample_data(logo_location, asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
		text_box = TextArea(logo_text, textprops=dict(color='#444444', fontsize=50, weight='bold'))

		logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

		anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.)

		return anchored_box
		

	def normalize_hists(self):
		for i in range(len(self._hists)):
			for j in range(len(self._hists[i])):
				bin_width = (self._hists[i][j].hist().upperbound() - self._hists[i][j].hist().lowerbound()) / self._hists[i][j].hist().nbins()

				if self._hists[i][j].hist().GetSumOfWeights() != 0.0:
					self._hists[i][j].hist().Scale(1.0 / ( self._hists[i][j].hist().GetSumOfWeights() * bin_width ))
	

	def set_formatting(self):
		for i in range(len(self._hists)):
			for j in range(len(self._hists[i])):
				self._hists[i][j].hist().SetColor(self._plot_colors[j])
				self._hists[i][j].hist().SetTitle(self._plot_labels[j])
				self._hists[i][j].hist().SetLineWidth(8)

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
				
				if y[i] == 0.0:
					y_zero_removed.append( None )
				else:
					y_zero_removed.append(y[i])



		return a_zero_removed, y_zero_removed



	def plot(self, filename):
		
		# Normalize all the histograms.
		self.normalize_hists()

		# Set basic plot element formattings. 
		self.set_formatting()


		with PdfPages(filename) as pdf:


			for k in range(len(self._hists)):	# k is in a sense the page number of our multi-page PDF plot. For each value of k, we produce one complete plot (with a set of 4 histograms viz. data, pythia, herwig, sherpa).

				print "Printing page {} of the plot.".format(k + 1)

				if self._ratio_plot:
					gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
					ax0 = plt.subplot(gs[0])
					ax1 = plt.subplot(gs[1])
				else:
					ax0 = plt.gca()

				# Set the logo.
				ax0.add_artist(self.logo_box())

				z_indices = range(len(self._hists[k]), 0, -1)
				z_indices[0] *= 10

				# First, draw the regular "non-ratio" plot.

				plt.sca(ax0)

				all_plots = []
				
				points_x_s, points_y_s = [], []

				for i in range(len(self._hists[k])):
					hist = self._hists[k][i].hist()
					plot_type = self._plot_types[i]

					if plot_type == 'hist':
						# plot = rplt.hist(hist, ax=ax0, zorder=z_indices[i], emptybins=False)
						plot = rplt.hist(hist, zorder=z_indices[i], emptybins=False)
						points_x, points_y = plot[1].get_xdata(), plot[1].get_ydata()

					elif plot_type == 'error':
						# plot = rplt.errorbar(hist, ax=ax0, zorder=z_indices[i], emptybins=False, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
						plot = rplt.errorbar(hist, zorder=z_indices[i], emptybins=False, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
						points_x, points_y = plot[0].get_xdata(), plot[0].get_ydata()
					
					points_x_s.append(points_x)
					points_y_s.append(points_y)

					

					all_plots.append(plot)

				if self._y_scale == 'log':
					ax0.set_yscale('log')


				# Ratio plot.

				plt.sca(ax1)

				if self._ratio_plot:

					denominator_hist = self._hists[k][self._ratio_to_index].hist()

					
					for i in range(len(self._hists[k])):

						if i != self._ratio_to_index:
							ratio_hist = copy.deepcopy( self._hists[k][i].hist() )
							ratio_hist.Divide(denominator_hist)
						else:
							ratio_hist = self._hists[k][i].hist().empty_clone(color=self._hists[k][i].hist().GetColor()[0])

							print ratio_hist.bounds()[1]

							map( ratio_hist.Fill, np.linspace(self._hists[k][i].hist().bounds()[0], self._hists[k][i].hist().bounds()[1], self._hists[k][i].hist().nbins()), [1.] * self._hists[k][i].hist().nbins() )

							print ratio_hist.bounds()[1]
							
						plot_type = self._plot_types[i]

						if plot_type == 'hist':
							# rplt.hist(ratio_hist, axes=ax1, zorder=z_indices[i], emptybins=False)

							line_plot = self.convert_hist_to_line_plot(ratio_hist, ratio_hist.bounds())

							# plt.plot(line_plot[0], line_plot[1], ax=ax1, lw=8, color=ratio_hist.GetColor()[0])
							plt.plot(line_plot[0], line_plot[1], lw=8, color=ratio_hist.GetColor()[0])

						elif plot_type == 'error':

							# We need to calculate ratio of errors ourselves.
							
							x_errors, y_errors = [], []
							for x_segment in all_plots[i][2][0].get_segments():
								x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
							for y_segment in all_plots[i][2][1].get_segments():
								y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

							ratio_y_err = [ b / m if m != 0 else None for b, m in zip(y_errors, points_y_s[self._ratio_to_index]) ]
							
							# rplt.errorbar(ratio_hist, ax=ax1, zorder=z_indices[i], emptybins=False, xerr=1, yerr=ratio_y_err, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
							rplt.errorbar(ratio_hist, zorder=z_indices[i], emptybins=False, xerr=1, yerr=ratio_y_err, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)


				# Ratio plot ends.


				# Legend.

				plt.sca(ax0)
				legend_handles = all_plots
				handles, labels = legend_handles, self._plot_labels
				legend = ax0.legend(handles, labels, loc=1, frameon=0, fontsize=60)
				ax0.add_artist(legend)


				plt.autoscale()

				print self._hists[k][0].x_lims()


				'''
				if self._x_lims[1] == -1:
					print "hey"
					ax0.set_xlim( self._x_lims[0], ax0.get_xlim()[1] )
				else:
					ax0.set_xlim( self._x_lims[0], self._x_lims[1] )

				'''
				'''
				if self._y_lims[1] == -1:
					ax0.set_ylim( self._y_lims[0], ax0.get_ylim()[1] * 1.125 )
				else:
					ax0.set_ylim( self._y_lims[0], self._y_lims[1] )

				if self._ratio_plot:
					ax1.set_xlim( ax0.get_xlim()[0], ax0.get_xlim()[1] )
					ax1.set_ylim(0., 2.)
				'''


				# Axes labels.

				ax0.set_xlabel(self._x_label, fontsize=75)
				ax0.set_ylabel(self._y_label, rotation=0, fontsize=75, labelpad=75)

				if self._ratio_plot:
					ax1.set_xlabel(self._x_label, fontsize=75)
					ax1.set_ylabel(self._ratio_label, rotation=0, fontsize=55, labelpad=115, y=0.31)

				# Axes labels end.

				
				plt.sca(ax0)
				plt.tick_params(which='major', width=5, length=25, labelsize=70)
				plt.tick_params(which='minor', width=3, length=15)

				if self._ratio_plot:
					plt.sca(ax1)
					
					plt.tick_params(which='major', width=5, length=25, labelsize=70)
					plt.tick_params(which='minor', width=3, length=15)


				if self._ratio_plot:
					plt.gcf().set_size_inches(30, 30, forward=1)
				else:
					plt.gcf().set_size_inches(30, 24, forward=1)

				plt.gcf().set_snap(True)

				plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

				# plt.savefig(filename)
				pdf.savefig()

				# plt.clf()
				plt.close()


