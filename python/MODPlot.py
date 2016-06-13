
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

	def __init__(self, hists, plot_types, plot_colors, plot_labels, line_styles, x_scale='linear', y_scale='linear', mark_regions=[], ratio_plot=False, ratio_to_index=-1, ratio_label="", x_label="", y_label="", x_lims=(0, -1), y_lims=(0, -1)):
		
		self._hists = hists
		self._plot_types = plot_types
		self._plot_colors = plot_colors
		self._plot_labels = plot_labels
		self._line_styles = line_styles

		self._x_scale = x_scale
		self._y_scale = y_scale

		self._mark_regions = mark_regions

		self._ratio_plot = ratio_plot
		self._ratio_to_index = ratio_to_index
		self._ratio_label = ratio_label

		self._x_label = x_label
		self._y_label = y_label

		self._x_lims = x_lims
		self._y_lims = y_lims

		self._plt = plt
		self._rplt = rplt
		

		self._plot_points_x_s = []
		self._plot_points_y_s = []
	
	def extrap1d(self, interpolator):
		xs = interpolator.x
		ys = interpolator.y

		def pointwise(x):
				if x < xs[0]:
						return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
				elif x > xs[-1]:
						return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
				else:
						return interpolator(x)

		def ufunclike(xs):
				return array(map(pointwise, array(xs)))

		return ufunclike




	def logo_box(self):
		
		logo_offset_image = OffsetImage(read_png(get_sample_data(logo_location, asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
		text_box = TextArea(logo_text, textprops=dict(color='#444444', fontsize=50, weight='bold'))

		logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

		anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.)

		return anchored_box
		

	def normalize_hists(self):

		for i in range(len(self._hists)):
			if self._plot_types[i] != "theory":
				bin_width = (self._hists[i].hist().upperbound() - self._hists[i].hist().lowerbound()) / self._hists[i].hist().nbins()
					
				if self._hists[i].hist().GetSumOfWeights() != 0.0:
					self._hists[i].hist().Scale(1.0 / ( self._hists[i].hist().GetSumOfWeights() * bin_width ))

	def set_formatting(self):
		for i in range(len(self._hists)):
			if self._plot_types[i] != "theory":
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
		
		if "error" in self._plot_types:
			z_indices[self._plot_types.index("error")] *= 10

		if "theory" in self._plot_types:
			z_indices[self._plot_types.index("theory")] *= 5

		# First, draw the regular "non-ratio" plot.


		legend_handles = []
		for i in range(len(self._hists)):

			plot_type = self._plot_types[i]

			if plot_type == 'hist':
				self._hists[i].hist().SetLineStyle(self._line_styles[i])
				plot = self._rplt.hist(self._hists[i].hist(), axes=ax0, zorder=z_indices[i], emptybins=False)
				legend_handles.append( plot )
			elif plot_type == 'error':
				plot = self._rplt.errorbar(self._hists[i].hist(), axes=ax0, zorder=z_indices[i], emptybins=True, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
				legend_handles.append( plot )
			elif plot_type == "theory":

				if self._x_scale == "log":
					x_s = np.exp( self._hists[i][0] )
					y_line_s = np.array( self._hists[i][2] ) / simps(y=self._hists[i][2], x=x_s)
					y_min_s = np.array( self._hists[i][1] ) / simps(y=self._hists[i][1], x=x_s)
					y_max_s = np.array( self._hists[i][3] ) / simps(y=self._hists[i][3], x=x_s)
				else:
					x_s = self._hists[i][0]
					y_line_s = self._hists[i][2]
					y_min_s = self._hists[i][1]
					y_max_s = self._hists[i][3]

				
				# Distribution.
				ax0.plot(x_s, y_line_s, lw=5, zorder=z_indices[i], color=self._plot_colors[i])

				ax0.fill_between(x_s, y_max_s, y_min_s, zorder=z_indices[i], where=np.less_equal(y_min_s, y_max_s), facecolor=self._plot_colors[i], color=self._plot_colors[i], interpolate=True, alpha=0.3, linewidth=8.)

				theory_min_interpolate_function = self.extrap1d(interpolate.interp1d(x_s, y_min_s))
				theory_line_interpolate_function = self.extrap1d(interpolate.interp1d(x_s, y_line_s))
				theory_max_interpolate_function = self.extrap1d(interpolate.interp1d(x_s, y_max_s))

				theory_extrapolated_min = theory_min_interpolate_function(self._plot_points_x_s[self._plot_types.index("error")])
				theory_extrapolated_line = theory_line_interpolate_function(self._plot_points_x_s[self._plot_types.index("error")])
				theory_extrapolated_max = theory_max_interpolate_function(self._plot_points_x_s[self._plot_types.index("error")])

				
			if plot_type == "theory":
				self._plot_points_x_s.append( [] )
				self._plot_points_y_s.append( [] )
			elif plot_type == "error":
				data_points_x = plot[0].get_xdata()
				data_points_y = plot[0].get_ydata()

				data_plot_points_x = []
				data_plot_points_y = []

				for i in range(0, len(data_points_x)):
					if float(data_points_x[i]) > 0.0:	# This is just to ignore the points at the 0th bin.
						data_plot_points_x.append(data_points_x[i])
						data_plot_points_y.append(data_points_y[i])

				self._plot_points_x_s.append( data_plot_points_x )
				self._plot_points_y_s.append( data_plot_points_y )
			elif plot_type == "hist":

				
				bin_width = (self._hists[i].hist().upperbound() - self._hists[i].hist().lowerbound()) / self._hists[i].hist().nbins()

				data_points_x = [x + bin_width / 2 for x in plot[1].get_xdata()][:-1]
				data_points_y = plot[1].get_ydata()

				data_plot_points_x = []
				data_plot_points_y = []

				for i in range(0, len(data_points_x)):
					if float(data_points_x[i]) > 0.0:	# This is just to ignore the points at the 0th bin.
						data_plot_points_x.append(data_points_x[i])
						data_plot_points_y.append(data_points_y[i])

				self._plot_points_x_s.append( data_plot_points_x )
				self._plot_points_y_s.append( data_plot_points_y )
				


		if self._y_scale == 'log':
			ax0.set_yscale('log')


		# Ratio plot.

		if self._ratio_plot:

			denominator_hist = self._hists[self._ratio_to_index].hist()

			
			for i in range(len(self._hists)):

				plot_type = self._plot_types[i]

				if plot_type == "theory":
					
					ratio_theory_line = [None if n == 0 else m / n for m, n in zip(theory_extrapolated_line, self._plot_points_y_s[self._ratio_to_index])]
					ratio_theory_min = [None if n == 0 else m / n for m, n in zip(theory_extrapolated_min, self._plot_points_y_s[self._ratio_to_index])]
					ratio_theory_max = [None if n == 0 else m / n for m, n in zip(theory_extrapolated_max, self._plot_points_y_s[self._ratio_to_index])]

					ratio_theory_line_hist = Hist(self._hists[self._plot_types.index("error")].hist().nbins(), self._hists[self._plot_types.index("error")].hist().bounds()[0], self._hists[self._plot_types.index("error")].hist().bounds()[1])
					ratio_theory_min_hist = Hist(self._hists[self._plot_types.index("error")].hist().nbins(), self._hists[self._plot_types.index("error")].hist().bounds()[0], self._hists[self._plot_types.index("error")].hist().bounds()[1])
					ratio_theory_max_hist = Hist(self._hists[self._plot_types.index("error")].hist().nbins(), self._hists[self._plot_types.index("error")].hist().bounds()[0], self._hists[self._plot_types.index("error")].hist().bounds()[1])
					
					ratio_theory_line_hist.fill_array(self._plot_points_x_s[self._ratio_to_index], ratio_theory_line)
					ratio_theory_min_hist.fill_array(self._plot_points_x_s[self._ratio_to_index], ratio_theory_min)
					ratio_theory_max_hist.fill_array(self._plot_points_x_s[self._ratio_to_index], ratio_theory_max)

					line_plot = self.convert_hist_to_line_plot(ratio_theory_line_hist, ratio_theory_line_hist.bounds())
					min_plot = self.convert_hist_to_line_plot(ratio_theory_min_hist, ratio_theory_min_hist.bounds())
					max_plot = self.convert_hist_to_line_plot(ratio_theory_max_hist, ratio_theory_max_hist.bounds())
					
					self._plt.plot(line_plot[0], line_plot[1], zorder=z_indices[i], axes=ax1, lw=8, color=self._plot_colors[self._plot_types.index("theory")])
					ax1.fill_between( max_plot[0], max_plot[1], min_plot[1], zorder=z_indices[i], where=np.less_equal(min_plot[1], max_plot[1]), color=self._plot_colors[i], facecolor=self._plot_colors[i], interpolate=True, alpha=0.3, linewidth=0.0)

					pass
				else:
					ratio_hist = copy.deepcopy( self._hists[i].hist() )
					ratio_hist.Divide(denominator_hist)

				

				if plot_type == 'hist':
					rplt.hist(ratio_hist, axes=ax1, zorder=z_indices[i], emptybins=False)
					
					line_plot = self.convert_hist_to_line_plot(ratio_hist, ratio_hist.bounds())

					# self._plt.plot(line_plot[0], line_plot[1], ls=self._line_styles[i], axes=ax1, lw=8, color=ratio_hist.GetColor()[0])

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
		


		



		# Any possible markers.
		for marker in self._mark_regions:
			ax0.plot([marker, marker], [ax0.get_ylim()[0], ax0.get_ylim()[1] / 2], color='red', linewidth=5, linestyle="dashed")
			

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

		
		if self._hists[0].axes_label_pi():
			plt.sca(ax0)
			plt.gca().set_xticks( [0, round(0.5*np.pi, 3), round(np.pi, 3), round(1.5*np.pi, 3), round(2*np.pi, 3)] )
			plt.gca().set_xticklabels( ["0", "$\pi / 2$", "$\pi$", "$3 \pi / 2$", "$2 \pi$"] )

			if self._ratio_plot:
				plt.sca(ax1)
				plt.gca().set_xticks( [0, round(0.5*np.pi, 3), round(np.pi, 3), round(1.5*np.pi, 3), round(2*np.pi, 3)] )
				plt.gca().set_xticklabels( ["0", "$\pi / 2$", "$\pi$", "$3 \pi / 2$", "$2 \pi$"] )				


		
		self._plt.sca(ax0)
		self._plt.autoscale()
		self._plt.xscale(self._x_scale)
		self._plt.gca().xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

		self._plt.sca(ax1)
		self._plt.autoscale()
		self._plt.xscale(self._x_scale)
		self._plt.gca().xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

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

		


		self._plt.gcf().set_snap(True)

		self._plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

	def save_plot(self, filename):
		self._plt.savefig(filename)
		self._plt.clf()

	def get_plt(self):
		return self._plt

	def get_rplt(self):
		return self._rplt





plot_labels = { "data": "CMS 2010 Open Data", "pythia": "Pythia 8.215", "herwig": "Herwig 7.0.1", "sherpa": "Sherpa 2.2.0", "theory": "Theory (MLL)" }
plot_colors = {"theory": "red", "pythia": "blue", "herwig": "green", "sherpa": "purple", "pythia_post": "red", "data": "black", "data_post": "red"}

global_plot_types = ['error', 'hist', 'hist', 'hist']
global_colors = [ plot_colors['data'], plot_colors['pythia'], plot_colors['herwig'], plot_colors['sherpa'] ]
global_labels = [ plot_labels['data'], plot_labels['pythia'], plot_labels['herwig'], plot_labels['sherpa'] ]
global_line_styles = [ "", "solid", "dashed", "dotted" ]


def create_multi_page_plot(filename, hists, theory=False, x_scale='linear'):
	# mod_hists is a list of MODHist objects.

	labels = copy.deepcopy(global_labels)
	types = copy.deepcopy(global_plot_types)
	colors = copy.deepcopy(global_colors)
	line_styles = copy.deepcopy(global_line_styles)


	if theory and "theory" not in types:
		types.insert(1, "theory")
		colors.insert(1, plot_colors['theory'])
		labels.insert(1, plot_labels['theory'])
		line_styles.insert(1, "")

	with PdfPages(filename) as pdf:

		for mod_hists in hists:	# mod_hists contains a list [ data_mod_hist, pythia_mod_hist, ... ]

			if mod_hists[0].x_scale() == x_scale:
				
				if theory:
					ratio_to_index = 2
				else:
					ratio_to_index = 1

				plot = MODPlot(mod_hists, plot_types=types, plot_colors=colors, plot_labels=labels, line_styles=line_styles, x_scale=mod_hists[0].x_scale(), y_scale=mod_hists[0].y_scale(), ratio_plot=True, ratio_to_index=ratio_to_index, ratio_label="Ratio\nto\nPythia", mark_regions=mod_hists[0].mark_regions(), x_label=mod_hists[0].x_label(), y_label=mod_hists[0].y_label(), x_lims=mod_hists[0].x_range(), y_lims=mod_hists[0].y_range())
				plot.plot()
			
				pdf.savefig()
				plot.get_plt().close()