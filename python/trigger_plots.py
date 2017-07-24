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




mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)




default_dir = "plots/Version 6/"


logo_location = "/home/aashish/root/macros/MODAnalyzer/mod_logo.png"
logo_text = ""




parsed_linear = trigger_parse.load_root_files_to_hist()


def logo_box(x, y):
	
	logo_offset_image = OffsetImage(read_png(get_sample_data(logo_location, asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
	text_box = TextArea(logo_text, textprops=dict(color='#444444', fontsize=50, weight='bold'))

	logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

	anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0., bbox_to_anchor=[x, y], bbox_transform = plt.gcf().transFigure)

	return anchored_box



def normalize_hist(hist):
	bin_width = (hist.upperbound() - hist.lowerbound()) / hist.nbins()
	
	if hist.GetSumOfWeights() != 0.0:
		hist.Scale(1.0 / ( hist.GetSumOfWeights() * bin_width ))

	return hist



def trigger_turn_on_curves():

	mod_hists = parsed_linear[0]

	colors = ['green', 'magenta', 'blue', 'red', 'orange', '#cecece']
	labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U", "Jet15\_HNF" ]
	hist_labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U", "Jet15U_HcalNoiseFiltered" ]
	lower_pTs = [140, 100, 70, 50, 30, 15]


	for i in range(len(hist_labels)):
		
		# hist = normalize_hist( mod_hists[hist_labels[i]].hist() )
		hist = mod_hists[hist_labels[i]].hist()

		n_bins = hist.nbins()

		# Loop through those bins.
		for j in range(n_bins):
			current_bin_x = hist.GetBinCenter(j)
			current_bin_y = hist.GetBinContent(j),
			
			if current_bin_x < lower_pTs[i]:
				hist.SetBinContent(j, 0.)


		hist.SetColor(colors[i])
		hist.SetTitle(labels[i])

		rplt.errorbar(hist, zorder=range(len(hist_labels))[len(hist_labels) - i - 1], emptybins=False, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)


	# Info about R, pT_cut, etc.
	extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
	handles = [extra]
	labels = ["AK5; $\left| \eta \\right| < 2.4$"]
	info_legend = plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.28, 0.92])
	plt.gca().add_artist(info_legend)


	plt.gca().set_xlabel("Trigger Jet $p_T$ [GeV]", fontsize=70, labelpad=50)
	plt.gca().set_ylabel("Rescaled Events", fontsize=70, labelpad=25)

	plt.gca().add_artist(logo_box(0.105, 0.99))

	plt.gca().xaxis.set_minor_locator(MultipleLocator(10))
	plt.gca().set_yscale('log')

	handles, labels = plt.gca().get_legend_handles_labels()


	legend = plt.legend(handles[::-1], labels[::-1], frameon=0, fontsize=60, bbox_to_anchor=[0.97, 0.99])

	ax = plt.gca().add_artist(legend)

	outside_text = plt.gca().legend( [extra], ["CMS 2010 Open Data"], frameon=0, borderpad=0, fontsize=50, bbox_to_anchor=(1.0, 1.005), loc='lower right')
	plt.gca().add_artist(outside_text)




	plt.autoscale()
	plt.gca().set_ylim(1e2, 1e10)
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.tight_layout()

	plt.savefig(default_dir + "trigger_turn_on.pdf")

	plt.clf()





def trigger_efficiency_plot():
	mod_hists = parsed_linear[0]

	
	colors = ['green', 'magenta', 'blue', 'red', 'orange', '#cecece']
	labels = ["Jet140U / 100U", "Jet100U / 70U", "Jet70U / 50U", "Jet50U / 30U", "Jet30U / 15U\_HNF", "" ]
	hist_labels = [("Jet140U", "Jet100U"), ("Jet100U", "Jet70U"), ("Jet70U", "Jet50U"), ("Jet50U", "Jet30U"), ("Jet30U", "Jet15U_HcalNoiseFiltered") ]
	lower_pTs = [140, 100, 70, 50, 30, 15]

	# rplt.hist(mod_hists[0].hist())


	plots = []
	for i in range(len(hist_labels) ):
		
		first_hist, second_hist = mod_hists[hist_labels[i][0]], mod_hists[hist_labels[i][1]]

		ratio_hist = first_hist.hist() / second_hist.hist()

		ratio_hist.SetColor(colors[i])
		ratio_hist.SetTitle(labels[i])

		# rplt.errorbar(ratio_hist, emptybins=False, xerr=1, yerr=1, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
		plots.append( rplt.errorbar(ratio_hist, emptybins=True, alpha=0.) )


	for i in range(len(plots)):
	    data_plot = plots[i]

	    data_x_errors, data_y_errors = [], []
	    for x_segment in data_plot[2][0].get_segments():
	      data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
	    for y_segment in data_plot[2][1].get_segments():
	      data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

	    data_points_x = data_plot[0].get_xdata()
	    data_points_y = data_plot[0].get_ydata()

	    filtered_x, filtered_y, filtered_x_err, filtered_y_err = [], [], [], []
	    for x, y, xerr, yerr in zip(data_points_x, data_points_y, data_x_errors, data_y_errors):
	      if x > lower_pTs[i]:
	        filtered_x.append(x)
	        filtered_y.append(y)
	        filtered_x_err.append(xerr)
	        filtered_y_err.append(yerr)

	    plt.errorbar(filtered_x, filtered_y, zorder=range(len(plots))[len(plots) - i - 1], color=colors[i], markeredgecolor=colors[i], label=labels[i], xerr=filtered_x_err, yerr=filtered_y_err, ls='None', alpha=1.0, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)


	extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)

	plt.autoscale()
	plt.gca().set_ylim(1e-3, 7e4)

	cms_turn_on_pTs = [250, 200, 150, 115, 85]
	for i in range(0, len(hist_labels)):

		'''
		if cms_turn_on_pTs[i] != 0:
			source = "MOD"
			plt.gca().annotate(str(cms_turn_on_pTs[i]) + " GeV", xy=(cms_turn_on_pTs[i], 1.), xycoords='data', xytext=(-100, 350),  textcoords='offset points', color=colors[i], size=50, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor=colors[i], zorder=99, connectionstyle="angle3,angleA=0,angleB=90") )
		'''

		plt.plot([ cms_turn_on_pTs[i], cms_turn_on_pTs[i] ], [ 1e0, 2e1 ], lw=10, ls="dashed", color=colors[i])

		plt.gca().text((cms_turn_on_pTs[i]), 3e1, str(cms_turn_on_pTs[i]) + " GeV", color=colors[i], horizontalalignment='center')
		  
	# Horizontal Line.
	plt.plot([0] + list(mod_hists[hist_labels[i][0]].hist().x()), [1] * (1 + len(list(mod_hists[hist_labels[i][0]].hist().x()))), color="black", linewidth=5, linestyle="dashed")





	# Info about R, pT_cut, etc.
	
	handles = [extra]
	labels = ["AK5; $\left| \eta \\right| < 2.4$"]
	info_legend = plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.28, 0.92])
	plt.gca().add_artist(info_legend)


	plt.gca().set_xlabel("Trigger Jet $p_T$ [GeV]", fontsize=70, labelpad=50)
	plt.gca().set_ylabel("Ratio", fontsize=70, labelpad=50)

	plt.gca().add_artist(logo_box(0.114, 0.98))

	plt.gca().xaxis.set_minor_locator(MultipleLocator(10))
	plt.gca().set_yscale('log')

	handles, labels = plt.gca().get_legend_handles_labels()
	legend = plt.legend(handles[::-1][:-5], labels[::-1][:-5], fontsize=60, frameon=0, bbox_to_anchor=[0.99, 0.99])
	ax = plt.gca().add_artist(legend)


	outside_text = plt.gca().legend( [extra], ["CMS 2010 Open Data"], frameon=0, borderpad=0, fontsize=50, bbox_to_anchor=(1.0, 1.005), loc='lower right')
	plt.gca().add_artist(outside_text)
	

	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.tight_layout()

	plt.gcf().set_size_inches(30, 24, forward=1)
	plt.savefig(default_dir + "trigger_efficiency.pdf")

	plt.clf()





def trigger_prescales():
	mod_hists = parsed_linear[0]

	prescale_hists = {}
	for key, value in mod_hists.items():
		if "prescale" in key:
			prescale_hists[key.split("_prescale")[0]] = value



	
	colors = ['green', 'magenta', 'blue', 'red', 'orange'][::-1]
	legend_labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U"][::-1]
	hist_labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U" ][::-1]
	

	plots = []
	for i in range(len(hist_labels) ):
		mod_hist = prescale_hists[hist_labels[i]]
		hist = mod_hist.hist()

		hist.SetColor(colors[i])
		hist.SetTitle(legend_labels[i])
		hist.SetLineWidth(5)

		# plots.append( rplt.hist(hist, linewidth=10) )

		# , markeredgecolor=colors[i], label=labels[i], xerr=filtered_x_err, yerr=filtered_y_err, ls='None', alpha=1.0, marker='o'
		plots.append( rplt.errorbar(hist, linewidth=5, markeredgecolor=colors[i], marker="o", markersize=15, pickradius=20, yerr=False, capthick=5, capsize=8, elinewidth=5) )

	

	plt.autoscale()
	plt.gca().set_ylim(1e1, 5e7)
	plt.gca().set_xlim(0.5, 1e4)

	average_prescales = [1.00, 1.93, 5.36, 100.31, 851.39][::-1]
	for i in range(0, len(hist_labels)):

		plt.plot([ average_prescales[i], average_prescales[i] ], [ 3e1, 1e2 ], lw=10, ls="dashed", color=colors[i])

		plt.gca().text((average_prescales[i]), 2e1, str(average_prescales[i]), color=colors[i], horizontalalignment='center')
		  





	# # Info about R, pT_cut, etc.
	extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
	handles = [extra]
	labels = ["AK5; $\left| \eta \\right| < 2.4$"]
	info_legend = plt.gca().legend(handles, labels, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.30, 0.98])
	plt.gca().add_artist(info_legend)


	plt.gca().set_xlabel("Trigger Prescale", fontsize=70, labelpad=10)
	plt.gca().set_ylabel("Events", fontsize=70, labelpad=50)

	plt.gca().add_artist(logo_box(0.103, 1.00))


	plt.gca().set_xscale('log')
	plt.gca().set_yscale('log')


	handles = plots
	legend = plt.legend(handles, legend_labels, fontsize=60, frameon=0, bbox_to_anchor=[0.99, 0.99])
	ax = plt.gca().add_artist(legend)

	# plt.gca().xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=False))
	
	outside_text = plt.gca().legend( [extra], ["CMS 2010 Open Data"], frameon=0, borderpad=0, fontsize=50, bbox_to_anchor=(1.0, 1.005), loc='lower right')
	plt.gca().add_artist(outside_text)


	plt.tick_params(which='major', width=5, length=25, labelsize=70, pad=10)
	plt.tick_params(which='minor', width=3, length=15)

	# plt.tight_layout()

	plt.gcf().set_size_inches(30, 24, forward=1)
	plt.savefig(default_dir + "trigger_prescales.pdf")

	plt.clf()




start = time.time()


# trigger_turn_on_curves()

# trigger_efficiency_plot()

trigger_prescales()

end = time.time()

print "Finished all plotting in {} seconds.".format(end - start)
