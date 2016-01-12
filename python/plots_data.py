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


from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

from scipy.stats import norm
from scipy.stats import gamma
from scipy import arange, array, exp

from scipy.stats import binned_statistic

import rootpy.plotting.views

input_analysis_file = sys.argv[1]


mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)



plot_labels = { "data": "CMS 2010 Open Data", "pythia": "Pythia 8.212", "herwig": "Herwig++ 2.7.1", "sherpa": "Sherpa 2.2.0", "theory": "Theory (MLL)" }
plot_colors = {"theory": "red", "pythia": "blue", "herwig": "green", "sherpa": "orange", "pythia_post": "magenta", "data": "black"}



def get_version(input_file):

  f = open(input_file, 'r')
  lines = f.read().split("\n")

  properties = defaultdict(list)

  for line in lines:
    numbers = line.split()
    if numbers[0] == "%":
      return numbers[1] + " " + numbers[2] 


def parse_file(input_file, pT_lower_cut=150., pT_upper_cut=20000.):
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

        corrected_pT_index = keywords.index("cor_hardest_pT") + 1

        if float(numbers[corrected_pT_index]) > pT_lower_cut and float(numbers[corrected_pT_index]) < pT_upper_cut:
          for i in range(len(keywords)):
            keyword = keywords[i]

            if keyword == 'trigger_name':
              properties[keyword].append( numbers[i + 1] ) # + 1 because we ignore the first keyword "Entry".
            else:
              properties[keyword].append( float(numbers[i + 1]) ) # + 1 because we ignore the first keyword "Entry".

    except:
      pass


  return properties




def normalize_hist(hist, variable_bin=False):

  if variable_bin:
    hist.Scale( 1.0 / hist.GetSumOfWeights() )

    for i in range(0, hist.GetSize()):
      bin_width = hist.GetXaxis().GetBinWidth(i)
      old_bin_height =  hist.GetBinContent(i)
      
      new_height = old_bin_height / bin_width
      hist.SetBinContent(i, new_height)

      old_error = hist.GetBinError(i)
      new_error = old_error / bin_width
      hist.SetBinError(i, new_error)
  else:
    bin_width = (hist.upperbound() - hist.lowerbound()) / hist.nbins()
    hist.Scale( 1.0 / ( hist.GetSumOfWeights() * bin_width ) )

  return hist


def log_bins(data, weights, number_of_bins=50):
  
  def drop_zeros(a_list):
    return [i for i in a_list if i > 0]

  min_value = math.log( min( drop_zeros(data) ), 10 )
  max_value = math.log( max(data), 10 )

  return np.histogram(data, weights=weights, bins=np.logspace(min_value, max_value, number_of_bins))


def lin_bins(data, weights, number_of_bins=50):
  
  def drop_zeros(a_list):
    return [i for i in a_list if i > 0]

  min_value = min( drop_zeros(data) )
  max_value = max(data)

  return np.histogram(data, weights=weights, bins=np.linspace(min_value, max_value, number_of_bins))





def plot_pTs(pT_lower_cut=100, pT_upper_cut=10000):
  properties = parse_file(input_analysis_file, pT_lower_cut=pT_lower_cut, pT_upper_cut=pT_upper_cut)

  uncorrected_pTs = properties['uncor_hardest_pT']
  corrected_pTs = properties['cor_hardest_pT']
  prescales = properties['prescale']

  corrected_pt_hist = Hist(100, 0, 1000, title='Jet Energy Corrected', markersize=3.0, color='black')
  bin_width_corrected = (corrected_pt_hist.upperbound() - corrected_pt_hist.lowerbound()) / corrected_pt_hist.nbins()

  uncorrected_pt_hist = Hist(100, 0, 1000, title='Jet Energy Uncorrected', markersize=3.0, color='orange')
  bin_width_uncorrected = (uncorrected_pt_hist.upperbound() - uncorrected_pt_hist.lowerbound()) / uncorrected_pt_hist.nbins()

  map(uncorrected_pt_hist.Fill, uncorrected_pTs, prescales)
  map(corrected_pt_hist.Fill, corrected_pTs, prescales)

  corrected_pt_hist.Scale(1.0 / (corrected_pt_hist.GetSumOfWeights() * bin_width_corrected))
  uncorrected_pt_hist.Scale(1.0 / (uncorrected_pt_hist.GetSumOfWeights() * bin_width_uncorrected))
  
  gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 

  ax0 = plt.subplot(gs[0])
  ax1 = plt.subplot(gs[1])

  data_plot = rplt.errorbar(corrected_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  uncorrected_data_plot = rplt.errorbar(uncorrected_pt_hist, axes=ax0, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  

  data_x_errors, data_y_errors = [], []
  for x_segment in data_plot[2][0].get_segments():
    data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
  for y_segment in data_plot[2][1].get_segments():
    data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

  data_points_x = data_plot[0].get_xdata()
  data_points_y = data_plot[0].get_ydata()

  data_plot_points_x = []
  data_plot_points_y = []
  for i in range(0, len(data_points_x)):
    data_plot_points_x.append(data_points_x[i])
    data_plot_points_y.append(data_points_y[i])



  uncorrected_data_x_errors, uncorrected_data_y_errors = [], []
  for x_segment in uncorrected_data_plot[2][0].get_segments():
    uncorrected_data_x_errors.append((x_segment[1][0] - x_segment[0][0]) / 2.)
  for y_segment in uncorrected_data_plot[2][1].get_segments():
    uncorrected_data_y_errors.append((y_segment[1][1] - y_segment[0][1]) / 2.)

  uncorrected_data_points_x = uncorrected_data_plot[0].get_xdata()
  uncorrected_data_points_y = uncorrected_data_plot[0].get_ydata()

  uncorrected_data_plot_points_x = []
  uncorrected_data_plot_points_y = []
  for i in range(0, len(uncorrected_data_points_x)):
    uncorrected_data_plot_points_x.append(uncorrected_data_points_x[i])
    uncorrected_data_plot_points_y.append(uncorrected_data_points_y[i])




  data_to_data_y_err = [(b / m) for b, m in zip(data_y_errors, data_plot_points_y)]
  data_to_data_x_err = [(b / m) for b, m in zip(data_x_errors, [1] * len(data_plot_points_y))]

  uncorrected_to_corrected_y_err = [(b / m) for b, m in zip(uncorrected_data_y_errors, data_plot_points_y)]
  uncorrected_to_corrected_x_err = [(b / m) for b, m in zip(uncorrected_data_x_errors, [1] * len(data_plot_points_y))]


  # Legends Begin.

  legend = ax0.legend(loc=1, frameon=0, fontsize=60, bbox_to_anchor=[1.0, 1.0])
  ax0.add_artist(legend)

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  if pT_upper_cut != 10000:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<2.4$", r"$p_{T} \in [" + str(pT_lower_cut) + ", " + str(pT_upper_cut) + "]~\mathrm{GeV}$"]
  else:
    labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~R = 0.5;\eta<2.4$", r"$p_{T} > " + str(pT_lower_cut) + "~\mathrm{GeV}$"]
  ax0.legend([extra, extra], labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.92, 0.70])

  # Legends End.

  ax0.set_xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75, labelpad=45)
  ax1.set_xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75, labelpad=45)
  ax0.set_ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)
  ax1.set_ylabel("Ratio           \nto           \n" + "Corrected" + "           ", fontsize=55, rotation=0, labelpad=115, y=0.31)


  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root-6.04.06/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.26, 0.93), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.32, 0.9215, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  # Ratio Plot.
  uncorrected_pt_hist.Divide(corrected_pt_hist)
  corrected_pt_hist.Divide(corrected_pt_hist)

  rplt.errorbar(corrected_pt_hist, xerr=data_to_data_x_err, yerr=data_to_data_y_err, axes=ax1, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)
  rplt.errorbar(uncorrected_pt_hist, xerr=uncorrected_to_corrected_x_err, yerr=uncorrected_to_corrected_y_err, axes=ax1, emptybins=False, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  ax0.set_yscale('log')

  ax0.autoscale(True)
  ax1.autoscale(True)
  
  ax0.set_ylim(10e-8, 0.5 * 10e-1)
  ax1.set_ylim(0., 2.)

  ax0.set_xlim(0, 1000)
  ax1.set_xlim(0, 1000)

  plt.gcf().set_size_inches(30, 30, forward=1)

  plt.sca(ax0)
  plt.gca().xaxis.set_minor_locator(MultipleLocator(25))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.sca(ax1)
  plt.gca().xaxis.set_minor_locator(MultipleLocator(25))
  # plt.gca().yaxis.set_minor_locator(MultipleLocator(50))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  print "Printing pT spectrum with pT > " + str(pT_lower_cut) + " and pT < " + str(pT_upper_cut)

  plt.savefig("plots/" + get_version(input_analysis_file) + "/pT_distribution/pT_data_lower_" + str(pT_lower_cut) + "_pT_upper_" + str(pT_upper_cut) + ".pdf")
  # plt.show()
  plt.clf()



def plot_hardest_pt_corresponding_triggers():

  properties = parse_file(input_analysis_file, pT_lower_cut=0)

  pTs = properties['cor_hardest_pT']
  trigger_names = properties['trigger_name']
  prescales = properties['prescale']

  expected_trigger_names = ["HLT_Jet180U", "HLT_Jet140U", "HLT_Jet100U", "HLT_Jet70U", "HLT_Jet50U", "HLT_Jet30U", "HLT_Jet15U"]
  labels = ["Jet180u", "Jet140u", "Jet100u", "Jet70u", "Jet50u", "Jet30u", "Jet15u"]

  colors = ['purple', 'orange', 'brown', 'red', 'blue', 'magenta', 'green']

  pt_hists = []
  for i in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist(50, 0, 1000, title=labels[i].upper(), markersize=1.0, color=colors[i], linewidth=5))

  found = False
  for i in range(0, len(pTs)):
    found = False
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j] in trigger_names[i]:
        pt_hists[j].Fill(pTs[i], prescales[i])
        found = True

    if not found:
      print trigger_names[i]
      pt_hists[ len(pt_hists) - 1].Fill(pTs[i], prescales[i])


  for k in range(0, len(pt_hists)):
    rplt.errorbar(pt_hists[k], marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.yscale('log')
  plt.autoscale(True)

  plt.gca().set_ylim(0.1, 10e8)

  plt.legend(loc=0, frameon=0)

  plt.xlabel('$p_T \mathrm{(GeV)}$', fontsize=65)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(25))
  # plt.gca().yaxis.set_minor_locator(MultipleLocator(50))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)


  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root-6.04.06/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.20, 0.845), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.265, 0.835, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.savefig("plots/" + get_version(input_analysis_file) + "/hardest_pT_corresponding_triggers.pdf")
  # plt.show()
  plt.clf()





# plot_pTs(pT_lower_cut=100)
plot_hardest_pt_corresponding_triggers()