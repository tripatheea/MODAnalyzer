
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


logo_location = "/home/aashish/root/macros/MODAnalyzer/mod_logo.png"
logo_text = "Preliminary"






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

		anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0., bbox_to_anchor=[0.16, 0.98], bbox_transform = plt.gcf().transFigure)

		return anchored_box
		

	def normalize_hists(self):

		for i in range(len(self._hists)):
			if self._plot_types[i] != "theory":
				if self._x_scale == "log":
					bin_width = (np.log(self._hists[i].hist().upperbound()) - np.log(self._hists[i].hist().lowerbound())) / self._hists[i].hist().nbins()
				else:
					bin_width = (self._hists[i].hist().upperbound() - self._hists[i].hist().lowerbound()) / self._hists[i].hist().nbins()
				
				if self._hists[i].hist().GetSumOfWeights() != 0.0:
					self._hists[i].hist().Scale(1.0 / ( self._hists[i].hist().GetSumOfWeights() * bin_width ))

	def set_formatting(self):

		for i in range(len(self._hists)):
			if self._plot_types[i] != "theory":
				self._hists[i].hist().SetColor(self._plot_colors[i])
				self._hists[i].hist().SetTitle(self._plot_labels[i])
				self._hists[i].hist().SetLineWidth(8)


	def convert_log_hist_to_line_plot(self, hist, helper_x_s):
		n_bins = hist.nbins()

		line_x_s, line_y_s = [], []

		x_s = list(hist.x())
		y_s = list(hist.y())

		corrections = []
		for i in range(len(helper_x_s)):
			corrections.append( abs(helper_x_s[i] - x_s[i] ) )

		for i in range(len(x_s)):

			y = y_s[i]
			x = x_s[i]

			x_left = x - corrections[i]
			x_right = x + corrections[i]

			left_segment = [x - fraction * corrections[i] for fraction in np.linspace(0., 1., 10)]
			right_segment = [x + fraction * corrections[i] for fraction in np.linspace(0., 1., 10)]

			# line_x_s.extend( [x_left, x_right] )
			# line_y_s.extend( [y, y] )

			line_x_s.extend( left_segment )
			line_x_s.extend( right_segment )

			line_y_s.extend( [y] * (2 * len(left_segment)) )

		return line_x_s, line_y_s



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


	def get_np_correction_boundary(self):

		# First, figure out what kind of thing this is. 
		var = self._x_label

		# Next, find out the lowest pT in the bin. 
		# Just use lower boundary of pT for now.


		lambda_value = 1.

		lower_pT = 85
		lower_pT = 0.1
		
		for keyword, cond in self._hists[0].conditions():
			if keyword == "hardest_pT":
				lower_pT = cond[0]
		
		# Next, find z_cut. 
		# TODO: This is too hacky. Fix this + comment it properly.
		z_cut = float(self._hists[0].additional_text()[0][2].split("=")[-1].split("$")[0])

		if "e_g^{0.5}":
			return (lambda_value / lower_pT)**0.5
		elif "e_g^2" in var:
			return (lambda_value**2) / (lower_pT * lower_pT * z_cut)
		elif "e_g" in var:
			return lambda_value / lower_pT
		elif "theta_g" in var:
			min_value = lambda_value / (lower_pT * z_cut)
		else:
			min_value = 0.0
		
		return min_value

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
			
				# self._hists[i].hist().SetLineStyle(self._line_styles[i])

				plot = self._rplt.hist(self._hists[i].hist(), axes=ax0, zorder=z_indices[i], emptybins=False)
				plot[1].set_dashes(self._line_styles[i])
				legend_handles.append( plot )

				# print legend_handles[i][1].get_dashes()
			
			elif plot_type == 'error':
			
				plot = self._rplt.errorbar(self._hists[i].hist(), axes=ax0, zorder=z_indices[i], emptybins=True, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
				legend_handles.append( plot )
			
			elif plot_type == "theory":

				if self._x_scale == "log":
					x_s = np.exp( self._hists[i][0] )
				else:
					
					x_s = self._hists[i][0]

				# There are two portions of x. 
				# Find the index where we need to switch.
				np_correction_index = 0
				for ab in range(len(x_s)):
					if x_s[ab] > self.get_np_correction_boundary():
						np_correction_index = ab
						break


					
				
				y_line_s = self._hists[i][2]
				y_min_s = self._hists[i][1]
				y_max_s = self._hists[i][3]

				# print x_s[ : np_correction_index]
				# Distribution.
				ax0.plot(x_s[ : np_correction_index + 1], y_line_s[ : np_correction_index + 1], lw=8, zorder=z_indices[i], color=self._plot_colors[i], ls="dotted")
				ax0.plot(x_s[np_correction_index : ], y_line_s[np_correction_index : ], lw=8, zorder=z_indices[i], color=self._plot_colors[i])

				ax0.fill_between(x_s[np_correction_index : ], y_max_s[np_correction_index : ], y_min_s[np_correction_index : ], zorder=z_indices[i], where=np.less_equal(y_min_s[np_correction_index : ], y_max_s[np_correction_index : ]), facecolor=self._plot_colors[i], color=self._plot_colors[i], interpolate=True, alpha=0.2, linewidth=0.)

				theory_min_interpolate_function = self.extrap1d(interpolate.interp1d(x_s, y_min_s))
				theory_line_interpolate_function = self.extrap1d(interpolate.interp1d(x_s, y_line_s))
				theory_max_interpolate_function = self.extrap1d(interpolate.interp1d(x_s, y_max_s))

				if self._hists[0].x_range() != (0, -1):
					
					# This is important for log plots because we don't want to consider the bins that's not visibile while figuring out how many bins to use for theory.



					# x_s_to_use = filter( lambda x: x >= self._hists[0].x_range()[0], self._plot_points_x_s[self._plot_types.index("error")] )
					x_s_to_use = self._plot_points_x_s[self._plot_types.index("error")]
				
				else:
					x_s_to_use = self._plot_points_x_s[self._plot_types.index("error")]


				theory_extrapolated_min = theory_min_interpolate_function(x_s_to_use)
				theory_extrapolated_line = theory_line_interpolate_function(x_s_to_use)
				theory_extrapolated_max = theory_max_interpolate_function(x_s_to_use)

				# print "the data x points were", (x_s_to_use)
				
			if plot_type == "theory":
				self._plot_points_x_s.append( x_s )
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

				if self._x_scale != "log":
					data_points_x = [x + bin_width / 2 for x in plot[1].get_xdata()][:-1]
				else:
					data_points_x = plot[1].get_xdata()[:-1]
				
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

					if self._x_scale == "log":

						theory_x_ss = self._plot_points_x_s[self._plot_types.index("error")]
						theory_y_line_ss = theory_line_interpolate_function(theory_x_ss)
						theory_y_min_ss = theory_min_interpolate_function(theory_x_ss)
						theory_y_max_ss = theory_max_interpolate_function(theory_x_ss)

						to_ratio_x_ss = self._plot_points_x_s[self._ratio_to_index]
						to_ratio_y_ss = self._plot_points_y_s[self._ratio_to_index]

						theory_to_ratio_line = [a / b for a, b in zip(theory_y_line_ss, to_ratio_y_ss)]
						theory_to_ratio_min = [a / b for a, b in zip(theory_y_min_ss, to_ratio_y_ss)]
						theory_to_ratio_max = [a / b for a, b in zip(theory_y_max_ss, to_ratio_y_ss)]

						ratio_theory_line_hist = self._hists[self._ratio_to_index].hist().empty_clone()
						ratio_theory_min_hist = self._hists[self._ratio_to_index].hist().empty_clone()
						ratio_theory_max_hist = self._hists[self._ratio_to_index].hist().empty_clone()

						ratio_theory_line_hist.SetLineColor(self._plot_colors[i])
						ratio_theory_line_hist.SetLineWidth(8)

						map(ratio_theory_line_hist.Fill, theory_x_ss, theory_to_ratio_line)
						map(ratio_theory_min_hist.Fill, theory_x_ss, theory_to_ratio_min)
						map(ratio_theory_max_hist.Fill, theory_x_ss, theory_to_ratio_max)


						line_x, line_y = self.convert_log_hist_to_line_plot(ratio_theory_line_hist, to_ratio_x_ss)
						min_x, min_y = self.convert_log_hist_to_line_plot(ratio_theory_min_hist, to_ratio_x_ss)
						max_x, max_y = self.convert_log_hist_to_line_plot(ratio_theory_max_hist, to_ratio_x_ss)

						
						line_sorted_x = sorted(line_x)
						line_sorted_y = [x for (y,x) in sorted(zip(line_x, line_y))]

						min_sorted_x = sorted(min_x)
						min_sorted_y = [x for (y,x) in sorted(zip(min_x, min_y))]

						max_sorted_x = sorted(max_x)
						max_sorted_y = [x for (y,x) in sorted(zip(max_x, max_y))]

						# There are two portions of x. 
						# Find the index where we need to switch.
						np_correction_index = 0
						for ab in range(len(line_sorted_x)):
							if line_sorted_x[ab] > self.get_np_correction_boundary():
								np_correction_index = ab
								break

						self._plt.plot(line_sorted_x[ : np_correction_index + 1], line_sorted_y[ : np_correction_index + 1], zorder=z_indices[i], axes=ax1, lw=8, color=self._plot_colors[self._plot_types.index("theory")], ls="dotted")
						self._plt.plot(line_sorted_x[np_correction_index : ], line_sorted_y[np_correction_index : ], zorder=z_indices[i], axes=ax1, lw=8, color=self._plot_colors[self._plot_types.index("theory")])

						ax1.fill_between( max_sorted_x[np_correction_index : ], max_sorted_y[np_correction_index : ], min_sorted_y[np_correction_index : ], zorder=z_indices[i], where=np.less_equal(min_sorted_y[np_correction_index : ], max_sorted_y[np_correction_index : ]), color=self._plot_colors[i], facecolor=self._plot_colors[i], interpolate=True, alpha=0.2, linewidth=0.0)
						
					else:

						ratio_theory_line = [None if n == 0 else m / n for m, n in zip(theory_extrapolated_line, self._plot_points_y_s[self._ratio_to_index])]
						ratio_theory_min = [None if n == 0 else m / n for m, n in zip(theory_extrapolated_min, self._plot_points_y_s[self._ratio_to_index])]
						ratio_theory_max = [None if n == 0 else m / n for m, n in zip(theory_extrapolated_max, self._plot_points_y_s[self._ratio_to_index])]

						ratio_theory_line_hist = Hist(self._hists[self._plot_types.index("error")].hist().nbins(), self._hists[self._plot_types.index("error")].hist().bounds()[0], self._hists[self._plot_types.index("error")].hist().bounds()[1])
						ratio_theory_min_hist = Hist(self._hists[self._plot_types.index("error")].hist().nbins(), self._hists[self._plot_types.index("error")].hist().bounds()[0], self._hists[self._plot_types.index("error")].hist().bounds()[1])
						ratio_theory_max_hist = Hist(self._hists[self._plot_types.index("error")].hist().nbins(), self._hists[self._plot_types.index("error")].hist().bounds()[0], self._hists[self._plot_types.index("error")].hist().bounds()[1])
						
						ratio_theory_line_hist.fill_array(self._plot_points_x_s[self._ratio_to_index], ratio_theory_line)
						ratio_theory_min_hist.fill_array(self._plot_points_x_s[self._ratio_to_index], ratio_theory_min)
						ratio_theory_max_hist.fill_array(self._plot_points_x_s[self._ratio_to_index], ratio_theory_max)

						line_x, line_y = self.convert_hist_to_line_plot(ratio_theory_line_hist, ratio_theory_line_hist.bounds())
						min_x, min_y = self.convert_hist_to_line_plot(ratio_theory_min_hist, ratio_theory_min_hist.bounds())
						max_x, max_y = self.convert_hist_to_line_plot(ratio_theory_max_hist, ratio_theory_max_hist.bounds())
						
						# There are two portions of x. 
						# Find the index where we need to switch.
						np_correction_index = 0
						for ab in range(len(line_x)):
							if line_x[ab] > self.get_np_correction_boundary():
								np_correction_index = ab
								break

						self._plt.plot(line_x, line_y, zorder=z_indices[i], axes=ax1, lw=8, color=self._plot_colors[self._plot_types.index("theory")], ls="dotted")
						self._plt.plot(line_x[np_correction_index : ], line_y[np_correction_index : ], zorder=z_indices[i], axes=ax1, lw=8, color=self._plot_colors[self._plot_types.index("theory")])

						ax1.fill_between( max_x[np_correction_index : ], max_y[np_correction_index : ], min_y[np_correction_index : ], zorder=z_indices[i], where=np.less_equal(min_y[np_correction_index : ], max_y[np_correction_index : ]), color=self._plot_colors[i], facecolor=self._plot_colors[i], interpolate=True, alpha=0.2, linewidth=0.0)
						
				else:
					ratio_hist = copy.deepcopy( self._hists[i].hist() )
					ratio_hist.Divide(denominator_hist)

				

				if plot_type == 'hist':
					plot = rplt.hist(ratio_hist, axes=ax1, zorder=z_indices[i], emptybins=False)
					plot[1].set_dashes(self._line_styles[i])
					
					line_plot = self.convert_hist_to_line_plot(ratio_hist, ratio_hist.bounds())

					# self._plt.plot(line_plot[0], line_plot[1], ls=self._line_styles[i], axes=ax1, lw=8, color=ratio_hist.GetColor()[0])

				elif plot_type == 'error':
					self._rplt.errorbar(ratio_hist, axes=ax1, zorder=z_indices[i], emptybins=False, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)


		# Ratio plot ends.


		handles, labels = legend_handles, self._plot_labels

		handler_map = {}
		if "theory" in self._plot_types:
			th_line, = ax0.plot(range(1), linewidth=8, color='red')
			th_patch = mpatches.Patch(facecolor='red', alpha=0.2, linewidth=0., edgecolor='red')

			handles.insert( self._plot_types.index("theory"), (th_patch, th_line))
			# labels.insert( self._plot_types.index("theory"), self._plot_labels[self._plot_types.index("theory")])

			handler_map[th_line ] = HandlerLine2D(marker_pad=0)

		for i in range(len(handles)):
			if self._plot_types[i] == "hist":
				line, = ax0.plot(range(1), linewidth=8, color=self._plot_colors[i], dashes=self._line_styles[i])
				handles[i] = line

				handler_map[line] = HandlerLine2D(marker_pad=0)


		legend = ax0.legend(handles, labels, loc=1, frameon=0, fontsize=50, handler_map=handler_map )
		ax0.add_artist(legend)

		# Any additional texts.
		extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
		# ax0.legend([extra]*len(labels), labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[1.00, 0.53])

		for position, anchor_location, text in self._hists[0].additional_text():
			texts = text.split("\n")
			ax0.legend( [extra] * len(texts), texts, frameon=0, borderpad=0, fontsize=50, bbox_to_anchor=position, loc=anchor_location)
		


		



		# Any possible markers.
		for marker in self._mark_regions:
			ax0.plot([marker, marker], [ax0.get_ylim()[0], ax0.get_ylim()[1] / 2], color='red', linewidth=5, linestyle="dashed")
			

		# Axes labels.

		ax0.set_xlabel(self._x_label, fontsize=60)
		ax0.set_ylabel(self._y_label, rotation=0, fontsize=50, labelpad=85)

		if self._ratio_plot:
			ax1.set_xlabel(self._x_label, fontsize=60)
			ax1.set_ylabel(self._ratio_label, rotation=0, fontsize=55, labelpad=115, y=0.31)

		# Axes labels end.


		self._plt.sca(ax0)

		self._plt.tick_params(which='major', width=5, length=25, labelsize=50)
		self._plt.tick_params(which='minor', width=3, length=15)

		if self._ratio_plot:
			self._plt.sca(ax1)
			
			self._plt.tick_params(which='major', width=5, length=25, labelsize=50)
			self._plt.tick_params(which='minor', width=3, length=15)


		if self._ratio_plot:
			self._plt.gcf().set_size_inches(30, 30, forward=1)
		else:
			self._plt.gcf().set_size_inches(30, 24, forward=1)

		
		
		
		self._plt.sca(ax0)
		self._plt.autoscale()
		self._plt.xscale(self._x_scale)
		self._plt.gca().xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

		if self._ratio_plot:
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

		if self._hists[0].axes_label_pi():
			plt.sca(ax0)
			plt.gca().set_xticks( [0, round(0.5*np.pi, 3), round(np.pi, 3), round(1.5*np.pi, 3), round(2*np.pi, 3)] )
			plt.gca().set_xticklabels( ["0", "$\pi / 2$", "$\pi$", "$3 \pi / 2$", "$2 \pi$"] )

			if self._ratio_plot:
				plt.sca(ax1)
				plt.gca().set_xticks( [0, round(0.5*np.pi, 3), round(np.pi, 3), round(1.5*np.pi, 3), round(2*np.pi, 3)] )
				plt.gca().set_xticklabels( ["0", "$\pi / 2$", "$\pi$", "$3 \pi / 2$", "$2 \pi$"] )				

		if self._x_scale == "log":
			plt.sca(ax0)
			ax0.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
			plt.sca(ax1)
			ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))



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
global_line_styles = [ [], [], [21, 7], [7, 7] ]

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



def create_data_only_plot(filename, hists, labels, types, colors, line_styles, ratio_plot=True, ratio_to_label="", ratio_to_index=1):
	# mod_hists is a list of MODHist objects.

	with PdfPages(filename) as pdf:

		for mod_hists in hists:	# mod_hists contains a list [ data_mod_hist, pythia_mod_hist, ... ]

			plot = MODPlot(mod_hists, plot_types=types, plot_colors=colors, plot_labels=labels, line_styles=line_styles, ratio_plot=ratio_plot, ratio_to_index=ratio_to_index, ratio_label=ratio_to_label, x_scale=mod_hists[0].x_scale(), y_scale=mod_hists[0].y_scale(), x_label=mod_hists[0].x_label(), y_label=mod_hists[0].y_label(), x_lims=mod_hists[-1].x_range(), y_lims=mod_hists[-1].y_range())
			plot.plot()
		
			pdf.savefig()
			plot.get_plt().close()