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


from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox, AnchoredOffsetbox, HPacker

from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea

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


def parse_file(input_file, keywords_to_populate, pT_lower_cut=150., pT_upper_cut=20000., softdrop_pT_lower_cut=0., softdrop_pT_upper_cut=20000., eta_cut=2.4):

  # We'll populate only those fileds that are in the list keywords_to_populate.


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
        # pT_index = keywords.index("cor_hardest_pT") + 1
        pT_index = keywords.index("hardest_pT") + 1

        # eta_index = keywords.index("hardest_eta") + 1

        # if abs(float(numbers[eta_index])) < eta_cut and float(numbers[pT_index]) > pT_lower_cut and float(numbers[pT_index]) < pT_upper_cut:
        if float(numbers[pT_index]) > pT_lower_cut and float(numbers[pT_index]) < pT_upper_cut:
          for i in range(len(keywords)):
            keyword = keywords[i]

            if keyword in keywords_to_populate:
              if keyword == "trigger_name":
                properties[keyword].append( numbers[i + 1] )  # + 1 because we ignore the first keyword "Entry".
              else:
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
        if (float(numbers[5]) > pT_lower_cut) and (int(numbers[3]) == 1) and (int(numbers[4]) >= jet_quality_level):
          properties['event_number'].append( float( numbers[1] ) )
          properties['Cor_Hardest_pT'].append( float( numbers[5] ) )
          properties['prescale'].append( float( numbers[6] ) )
          properties['trigger_names'].append(  numbers[7] )

    except Exception as e:
      if len(numbers) != 0:
        print "Some kind of error occured while parsing the given file!"
        print e
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
  properties = parse_file_turn_on(input_analysis_file, pT_lower_cut=0)
  

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescale']

  expected_trigger_names = ["HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U", "HLT\_Jet15U\_HcalNoiseFiltered" ]
  labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U", "Jet15\_HNF" ]
  lower_pTs = [140, 100, 70, 50, 30, 15]

  colors = ['orange', 'brown', 'red', 'blue', 'magenta', 'green']
  colors = colors[::-1]

  pt_hists = []
  for i in range(0, len(expected_trigger_names)):
    pt_hists.append(Hist( int( (300 - lower_pTs[i]) / 5), lower_pTs[i], 300, title=labels[i], markersize=1.0, color=colors[i], linewidth=5))

  for i in range(0, len(pTs)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j].replace("\\", "") in trigger_names[i]:
        pt_hists[j].Fill(pTs[i], prescales[i])

  for k in range(len(pt_hists) - 1, -1, -1):
    rplt.errorbar(pt_hists[k], emptybins=False, ls='None', marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.autoscale(True)
  plt.yscale('log')

  plt.gca().set_ylim(0.1, 10e9)

  
  # Legend.
  handles, labels = plt.gca().get_legend_handles_labels()
  first_legend = plt.gca().legend(handles, labels, fontsize=60,  frameon=0, borderpad=0.1, bbox_to_anchor=[1.00, 1.00])
  ax = plt.gca().add_artist(first_legend)

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra]
  labels = ["$ \\textrm{Anti--}k_{t}\\textrm{:}~R = 0.5$\n$\eta<2.4$"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.99, 0.62])


  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=75)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=75.)



  logo_offset_image = OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
  text_box = TextArea("Prelim.", textprops=dict(color='#444444', fontsize=50, weight='bold'))

  logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

  anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.)
  plt.gca().add_artist(anchored_box)


  plt.gcf().set_size_inches(30, 30, forward=1)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(10))
  # plt.gca().yaxis.set_minor_locator(MultipleLocator(50))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)


  plt.savefig("plots/Version 5_2/turn_on_curves.pdf")
  # plt.show()
  plt.clf()


def plot_trigger_efficiency_curves(trigger_1, trigger_2, pT_upper_limit=800):
  
  properties = parse_file_turn_on(input_analysis_file, pT_lower_cut=0)

  event_numbers = properties['event_number']

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescale']

  colors = ['magenta', 'blue', 'orange', 'green', 'black', 'red']
  expected_trigger_names = ["HLT\_Jet180U", "HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U" ]

  color = colors[expected_trigger_names.index(trigger_1.replace("_", "\_"))]

  pt_hist_trigger_1 = Hist(500, 0, pT_upper_limit, title=trigger_1[4:], color=color, markersize=1.0, linewidth=5)
  pt_hist_trigger_2 = Hist(500, 0, pT_upper_limit, title=trigger_2[4:], color=color, markersize=1.0, linewidth=5)



  # trigger_names_proper = defaultdict(list)
  # pTs_proper = defaultdict(list)
  # prescales_proper = defaultdict(list)

  # for i in range(0, len(trigger_names)):
  #   trigger_names_proper[event_numbers[i]].append(trigger_names[i])
  #   pTs_proper[event_numbers[i]].append(pTs[i])
  #   prescales_proper[event_numbers[i]].append(prescales[i])

  # for event_number in trigger_names_proper:
  #   if trigger_1 in trigger_names_proper[event_number] and trigger_2 in trigger_names_proper[event_number]:
      
  #     for i in range(0, len(trigger_names_proper[event_number])):
  #       if trigger_1 in trigger_names_proper[event_number][i]:
  #         trigger_1_index = i

  #     for i in range(0, len(trigger_names_proper[event_number])):
  #       if trigger_2 in trigger_names_proper[event_number][i]:
  #         trigger_2_index = i

  #     print pTs[trigger_1_index]

  #     pt_hist_trigger_1.Fill(pTs[trigger_1_index], prescales[trigger_1_index])
  #     pt_hist_trigger_2.Fill(pTs[trigger_2_index], prescales[trigger_2_index])



  for i in range(0, len(pTs)):
    if trigger_1 in trigger_names[i]:
      pt_hist_trigger_1.Fill(pTs[i], prescales[i])

    # The len thingy is to make sure trigger names like HLT_Jet15U_HcalNoiseFiltered_v3 are excluded.
    # if trigger_2 in trigger_names[i] and len(trigger_names[i]) > (len(trigger_2) + 3):
    if trigger_2 in trigger_names[i]:
      pt_hist_trigger_2.Fill(pTs[i], prescales[i])


  pt_hist_trigger_1.Divide(pt_hist_trigger_2)


  rplt.errorbar(pt_hist_trigger_1, color=color,  markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/Users/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.23, 0.895), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  plt.gcf().text(0.29, 0.885, "Prelim.", fontsize=50, weight='bold', color='#444444', multialignment='center')

  plt.gcf().set_snap(1)

  # Horizontal Line.
  plt.plot(list(pt_hist_trigger_1.x()), [1] * len(list(pt_hist_trigger_1.x())), color="black", linewidth=5, linestyle="dashed")

 
  if trigger_1 == "HLT_L1Jet6U":
    lower_pT = 37
  elif trigger_1 == "HLT_Jet15U":
    lower_pT = 56
  elif trigger_1 == "HLT_Jet30U":
    lower_pT = 84
  elif trigger_1 == "HLT_Jet50U":
    lower_pT = 114
  elif trigger_1 == "HLT_Jet70U":
    lower_pT = 153
  elif trigger_1 == "HLT_Jet100U":
    lower_pT = 196
  else:
    lower_pT = 0


  if lower_pT != 0:
    # CMS Vertical line.
    plt.plot([lower_pT, lower_pT], [plt.gca().get_ylim()[0], 1.], color=color, linewidth=3, linestyle="dashed")



  # # MOD Vertical line.
  # efficient_pt_x = 0.0
  # efficient_pt_y = 0.0

  # efficient_pt_x_s = []
  # efficient_pt_y_s = []

  # distance_to_previous_point_close_to_one = 0

  # for i in range(0, len(list(pt_hist_trigger_1.x()))):
  #   if abs(list(pt_hist_trigger_1.y())[i] - 1.00) < 0.1:

  #     if distance_to_previous_point_close_to_one > 25:
  #       # Last point close to one too far away.
  #       # Empty the lists.
  #       efficient_pt_x_s, efficient_pt_y_s = [], []

  #     efficient_pt_x_s.append(list(pt_hist_trigger_1.x())[i])
  #     efficient_pt_y_s.append(list(pt_hist_trigger_1.y())[i])
  #     distance_to_previous_point_close_to_one = 0

  #   else:
  #     distance_to_previous_point_close_to_one += 1


  #   if len(efficient_pt_x_s) > 75:
  #     efficient_pt_x = efficient_pt_x_s[0]
  #     efficient_pt_y = efficient_pt_y_s[0]
  #     break



  # mod_efficient_pTs = [325, 260, 196, 153, 114, 84, 50, 32]
  # mod_efficient_pT = mod_efficient_pTs[expected_trigger_names.index(trigger_1.replace("_", "\_"))]




  # plt.plot([mod_efficient_pT, mod_efficient_pT], [plt.gca().get_ylim()[0], 1.], color="purple", linewidth=3, linestyle="dashed")
  

  # if lower_pT != 0:
  #   plt.gca().annotate("CMS\n" + str(lower_pT) + " GeV", xy=(lower_pT, 1.), xycoords='data', xytext=(-100, 250),  textcoords='offset points', color=color, size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor=color, zorder=9999, connectionstyle="angle3,angleA=0,angleB=90") )
  
  # plt.gca().annotate("MOD\n" + str(int(mod_efficient_pT)) + " GeV", xy=(mod_efficient_pT, 1.), xycoords='data', xytext=(250, 200), textcoords='offset points', color="purple", size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor="purple", zorder=9999, connectionstyle="angle3,angleA=45,angleB=-90") )



  plt.yscale('log')
  plt.gca().set_ylim(plt.gca().get_ylim()[0], 100)

  plt.legend(frameon=0)

  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=55, rotation=0)
  plt.ylabel('$\mathrm{A.U.}$', fontsize=75, rotation=0, labelpad=50.)

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)
  plt.gcf().set_snap(True)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/Version 4/trigger_efficiency/efficiency_curves_" + trigger_1 + ".pdf")
  # plt.show()
  plt.clf()






def plot_all_trigger_efficiency_curves():
  
  properties = parse_file_turn_on(input_analysis_file, pT_lower_cut=0)

  pTs = properties['corrected_hardest_pts']
  trigger_names = properties['trigger_names']
  prescales = properties['prescale']

 
  # colors = ['orange', 'red', 'green', 'blue', 'magenta', 'black']
  # colors = colors[::-1]

  
  colors = ['green', 'magenta', 'blue', 'red', 'brown', 'orange']
  expected_trigger_names = ["HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U", "HLT\_Jet15U\_HcalNoiseFiltered" ]
  labels = ["Jet140U / 100U", "Jet100U / 70U", "Jet70U / 50U", "Jet50U / 30U", "Jet30U / 15U\_HNF", "" ]
  lower_pTs = [140, 100, 70, 50, 30, 15]
  # lower_pTs = lower_pTs[::-1]

  cms_turn_on_pTs = [250, 200, 150, 115, 85]

  pt_hists = []
  for j in range(len(expected_trigger_names) - 1, -1, -1):
    pt_hists.append(Hist(50, 0, 300))


  for i in range(0, len(expected_trigger_names)):
    for j in range(0, len(pTs)):
      # The len thingy is to make sure trigger names like HLT_Jet15U_HcalNoiseFiltered_v3 are excluded.
      if expected_trigger_names[i].replace("\\", "") in trigger_names[j]:
        pt_hists[i].Fill(pTs[j], prescales[j])

      
  ratio_hists = []
  for i in range(0, len(pt_hists) - 1):
    ratio_hists.append(pt_hists[i] / pt_hists[i + 1])


  data_plots = []
  for i in range(0, len(pt_hists) - 1):
    data_plot = rplt.errorbar(ratio_hists[i], emptybins=False, alpha=0.0)
    data_plots.append(data_plot)

    if cms_turn_on_pTs[i] != 0:
      # plt.plot([cms_turn_on_pTs[i], cms_turn_on_pTs[i]], [plt.gca().get_ylim()[0], 1.], color=colors[i], linewidth=5, linestyle="dashed")
      # if cms_turn_on_pTs[i] > 153:
      #   source = "MOD"
      # else:
      #   source = "CMS"

      source = "MOD"
      plt.gca().annotate(str(cms_turn_on_pTs[i]) + " GeV", xy=(cms_turn_on_pTs[i], 1.), xycoords='data', xytext=(-100, 350),  textcoords='offset points', color=colors[i], size=40, va="center", ha="center", arrowprops=dict(arrowstyle="simple", facecolor=colors[i], zorder=99, connectionstyle="angle3,angleA=0,angleB=90") )

  



  for i in range(len(ratio_hists) - 1, -1, -1):
    data_plot = data_plots[i]


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

    plt.errorbar(filtered_x, filtered_y, color=colors[i], markeredgecolor=colors[i], label=labels[i], xerr=filtered_x_err, yerr=filtered_y_err, ls='None', alpha=1.0, marker='o', markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5)

  plt.gca().xaxis.set_tick_params(width=5, length=20, labelsize=70)
  plt.gca().yaxis.set_tick_params(width=5, length=20, labelsize=70)

  logo_offset_image = OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
  text_box = TextArea("Prelim.", textprops=dict(color='#444444', fontsize=50, weight='bold'))

  logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

  anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.)
  plt.gca().add_artist(anchored_box)


  plt.gcf().set_snap(1)

  # Horizontal Line.
  plt.plot([0] + list(pt_hists[0].x()), [1] * (1 + len(list(pt_hists[0].x()))), color="black", linewidth=5, linestyle="dashed")


  plt.yscale('log')
  plt.gca().set_ylim(0.0001, 20000)



  # Legend.
  handles, labels = plt.gca().get_legend_handles_labels()
  first_legend = plt.gca().legend(handles, labels, fontsize=60,  frameon=0, borderpad=0.1, bbox_to_anchor=[1.0, 1.0])
  ax = plt.gca().add_artist(first_legend)

  # Info about R, pT_cut, etc.
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  handles = [extra]
  labels = ["$ \\textrm{Anti--}k_{t}\\textrm{:}~R = 0.5$ \n $\eta<2.4$"]
  plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=60, bbox_to_anchor=[0.32, 0.85])





  plt.xlabel('$p_T~\mathrm{(GeV)}$', fontsize=65, rotation=0)
  plt.ylabel('$\mathrm{Ratio}$', fontsize=65, rotation=0, labelpad=50.)

  plt.gca().xaxis.set_minor_locator(MultipleLocator(10))
  # plt.gca().yaxis.set_minor_locator(MultipleLocator(50))
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)


  plt.gcf().set_size_inches(30, 30, forward=1)
  plt.gcf().set_snap(True)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/Version 5_2/all_efficiency_curves.pdf")
  # plt.show()
  plt.clf()




def plot_prescales():
  keywords = ['cor_hardest_pT', 'prescale', 'trigger_name']
  properties = parse_file(input_analysis_file, pT_lower_cut=0, keywords_to_populate=keywords)

  # properties = parse_file_turn_on(input_analysis_file, pT_lower_cut=0)
  

  
  trigger_names = properties['trigger_name']
  prescales = properties['prescale']

  expected_trigger_names = ["HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U" ]
  labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U" ]
  lower_pTs = [140, 100, 70, 50, 30]

  colors = ['brown', 'red', 'blue', 'magenta', 'green']
  colors = colors[::-1]

  prescale_hists = []
  for i in range(0, len(expected_trigger_names)):
    # prescale_hists.append(Hist(100, 1, 1e5, title=labels[i], markersize=1.0, color=colors[i], linewidth=5))
    prescale_hists.append([])

  for i in range(0, len(prescales)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j].replace("\\", "") in trigger_names[i]:
        prescale_hists[j].append( prescales[i] )

  for k in range(len(prescale_hists) - 1, -1, -1):
    # rplt.hist(prescale_hists[k], emptybins=False, lw=8, histtype='step')
    plt.hist(prescale_hists[k], label=labels[k], color=colors[k], lw=8, histtype='step')



  # handles = []
  # # for trigger, prescales in triggers.items():
  # for trigger_name in expected_trigger_names:
    
  #   zorder = expected_trigger_names.index(trigger_name)
  #   color = colors[expected_trigger_names.index(trigger_name)]
  #   label = trigger_name.replace("_", "\_")

  #   plt.hist(triggers[trigger_name], zorder=zorder, color=color, label=label, lw=8, histtype='step')
    


  logo_offset_image = OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
  text_box = TextArea("Preliminary", textprops=dict(color='#444444', fontsize=50, weight='bold'))

  logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

  anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0., bbox_to_anchor=[0.14, 1.0], bbox_transform = plt.gcf().transFigure)
  plt.gca().add_artist(anchored_box)


  plt.xlabel("Trigger Prescale", fontsize=65, rotation=0)
  plt.ylabel("No. of \n Events", fontsize=65, rotation=0, labelpad=50.)

  plt.legend(frameon=False)

  plt.xscale('log')
  plt.yscale('log')

  plt.autoscale()
  plt.ylim(plt.gca().get_ylim()[0], plt.gca().get_ylim()[1] * 2)

  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)

  plt.gcf().set_size_inches(30, 30, forward=1)
  plt.gcf().set_snap(True)

  plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

  plt.savefig("plots/Version 5_2/trigger_prescales.pdf")





def calculate_average_prescales():
  keywords = ['cor_hardest_pT', 'prescale', 'trigger_name']
  properties = parse_file(input_analysis_file, pT_lower_cut=0, keywords_to_populate=keywords)

  # properties = parse_file_turn_on(input_analysis_file, pT_lower_cut=0)
  

  
  trigger_names = properties['trigger_name']
  prescales = properties['prescale']

  expected_trigger_names = ["HLT\_Jet140U", "HLT\_Jet100U", "HLT\_Jet70U", "HLT\_Jet50U", "HLT\_Jet30U" ]
  labels = ["Jet140U", "Jet100U", "Jet70U", "Jet50U", "Jet30U" ]
  lower_pTs = [140, 100, 70, 50, 30, 15]

  colors = ['orange', 'brown', 'red', 'blue', 'magenta', 'green']
  colors = colors[::-1]

  prescale_hists = []
  for i in range(0, len(expected_trigger_names)):
    # prescale_hists.append(Hist(100, 1, 1e5, title=labels[i], markersize=1.0, color=colors[i], linewidth=5))
    prescale_hists.append([])

  for i in range(0, len(prescales)):
    for j in range(0, len(expected_trigger_names)):
      if expected_trigger_names[j].replace("\\", "") in trigger_names[i]:
        try:
          prescale_hists[j].append( float(prescales[i]) )
        except:
          pass

  for j in range(0, len(expected_trigger_names)):
    print expected_trigger_names[j], np.mean(prescale_hists[j])




plot_turn_on_curves()
# plot_all_trigger_efficiency_curves()

# plot_prescales()


# calculate_average_prescales()