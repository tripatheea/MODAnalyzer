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
            properties[keyword].append( float(numbers[i + 1]) ) # + 1 because we ignore the first keyword "Entry".

    except:
      pass


  return properties



def parse_file_turn_on(input_file, pT_lower_cut = 0.00, jet_quality_level=1):
  f = open(input_file, 'r')
  lines = f.read().split("\n")

  # FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3

  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        if (float(numbers[1]) > pT_lower_cut) and (int(numbers[3]) == 1) and (int(numbers[4]) >= jet_quality_level):
          properties['event_number'].append( float( numbers[1] ) )
          properties['corrected_hardest_pts'].append( float( numbers[5] ) )
          properties['prescale'].append( int( numbers[6] ) )
          properties['trigger_names'].append(  numbers[7] )

    except:
      if len(numbers) != 0:
        print "Some kind of error occured while parsing the given file!"
        print numbers
        print


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



def plot_turn_on_curves():
  # properties = parse_file_turn_on(input_analysis_file)
  
  properties = parse_file_turn_on('/home/aashish/turn_on.dat')

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescale']

  expected_trigger_names = ["HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U", "HLT\_Jet15U\_HcalNoiseFiltered" ]
  labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U", "Jet15U-HNF" ]

  colors = ['orange', 'red', 'green', 'blue', 'magenta', 'black']
  colors = colors[::-1]

  pt_hists = []
  for i in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist(60, 0, 300, title=labels[i], markersize=1.0, color=colors[i], linewidth=5))

  for i in range(0, len(pTs)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j].replace("\\", "") in trigger_names[i]:
        pt_hists[j].Fill(pTs[i], prescales[i])

  for k in range(len(pt_hists) - 1, -1, -1):
    rplt.errorbar(pt_hists[k], emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.autoscale(True)
  plt.yscale('log')

  plt.gca().set_ylim(1, 10e10)

  
  # Legend.
  handles, labels = plt.gca().get_legend_handles_labels()
  first_legend = plt.gca().legend(handles, labels, fontsize=60,  frameon=0, borderpad=0.1, bbox_to_anchor=[0.97, 0.98])
  ax = plt.gca().add_artist(first_legend)

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra]
  labels = [r"$ \textrm{Anti--}k_{t}\textrm{:}~\boldsymbol{R = 0.5} $"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.935, 0.65])


  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=95)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)

  fn = get_sample_data("/home/aashish/root-6.04.06/macros/MODAnalyzer/mod_logo.png", asfileobj=False)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root-6.04.06/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.9249985), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.29, 0.9178555, preliminary_text, fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_size_inches(30, 30, forward=1)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)


  plt.savefig("plots/" + get_version(input_analysis_file) + "/turn_on_curves.pdf")
  # plt.show()
  plt.clf()




# plot_turn_on_curves()

plot_hardest_pt_corresponding_triggers(pT_lower_cut=100)