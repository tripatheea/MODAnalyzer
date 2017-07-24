from __future__ import division

from subprocess import call


import time
import datetime


from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FormatStrFormatter

import matplotlib.dates as mdates

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

from scipy.stats import norm
from scipy.stats import gamma
from scipy import arange, array, exp

from scipy.stats import binned_statistic

import rootpy.plotting.views

from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter, DayLocator



input_analysis_file = sys.argv[1]


mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)



def parse_file(input_file):
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
        for i in range(len(keywords)):
          keyword = keywords[i]
          properties[keyword].append( float(numbers[i + 1]) ) # + 1 because we ignore the first keyword "Entry".

    except:
      pass


  return properties

def filter_run_lumi(data):
  lumis = {}
  for lumi_block, run_number, avg_inst_lumi, intg_del_lumi, intg_rec_lumi, time in zip(data['lumi_block'], data['Run_Number'], data['avg_inst_lumi'], data['intg_del_lumi'], data['intg_rec_lumi'], data['time']):
    
    key = (run_number, lumi_block)
    # key = lumi_block
    # key = run_number

    value = {'avg_inst_lumi': avg_inst_lumi, 'intg_del_lumi': intg_del_lumi, 'intg_rec_lumi': intg_rec_lumi, 'time': time}
    lumis[key] = value

  properties = defaultdict(list)

  for key in lumis:
    current_element = lumis[key]
    for prop in current_element:
      properties[prop].append( current_element[prop] )

  return properties


def plot_integrated_recorded_lumi(cumulative=False):
  properties = filter_run_lumi( parse_file(input_analysis_file) )

  timestamps = sorted(properties['time'])
  intg_rec_lumi = [x for (y,x) in sorted(zip(properties['time'], properties['intg_rec_lumi']))]
  intg_del_lumi = [x for (y,x) in sorted(zip(properties['time'], properties['intg_del_lumi']))]

  # Convert from (ub)-1 to (pb)-1
  intg_rec_lumi = [x * 1e-6 for x in intg_rec_lumi]
  intg_del_lumi = [x * 1e-6 for x in intg_del_lumi]



  dates = [datetime.datetime.fromtimestamp( int(time) ) for time in timestamps]


  total_rec_lumi = sum(intg_rec_lumi)
  total_del_lumi = sum(intg_del_lumi)
  max_lumi = max(intg_rec_lumi)

  print "Total Rec. Luminosity: {}".format(total_rec_lumi)
  print "Total Del. Luminosity: {}".format(total_del_lumi)
  # print "Max Luminosity: {}".format(max_lumi)

  # rec = [(36.1 * i) / total_del_lumi for i in intg_rec_lumi]
  # deli = [(36.1 * i) / total_del_lumi for i in intg_del_lumi]
  # intg_rec_lumi, intg_del_lumi = rec, deli

  
  # dates = [datetime.datetime.fromtimestamp( int(time) ).strftime('%m-%d') for time in timestamps]

  if cumulative:
    label = "CMS Recorded: " + str(round(total_rec_lumi, 2)) + " $\mathrm{pb}^{-1}$"
  else:
    label = "CMS Recorded, max " + str(round(max_lumi, 2)) + " $\mathrm{pb}^{-1}$/day"


  print "Min date = ", min(dates)
  print "Max date = ", max(dates)

  plt.hist(mpl.dates.date2num(dates), label="Recorded", weights=intg_rec_lumi, lw=8, bins=50, cumulative=cumulative, histtype='step', color='orange')
  plt.hist(mpl.dates.date2num(dates), label="Delivered", weights=intg_del_lumi, lw=8, bins=50, cumulative=cumulative, histtype='step', color='green')
  

  years = mdates.YearLocator()   # every year
  months = mdates.MonthLocator()  # every month
  yearsFmt = mdates.DateFormatter('%Y-%M')

  # Format the ticks
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
  # plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

  

  plt.xlabel("Date (2010)", labelpad=25, fontsize=70)

  plt.gca().ticklabel_format(axis='y', style='sci')

  if cumulative:
    plt.ylabel("Integrated Luminosity [A.U.]", labelpad=50, fontsize=70)
  else:
    plt.ylabel("Average Luminosity [A.U.]", labelpad=50, fontsize=70)


  plt.gca().get_yaxis().set_ticks([])

  print max(intg_rec_lumi),
  
  handles, labels = plt.gca().get_legend_handles_labels()
  legend = plt.gca().legend(handles[::-1], labels[::-1], bbox_to_anchor=[0.007, 0.99], frameon=False, loc='upper left', fontsize=70)
  plt.gca().add_artist(legend)

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)



  # plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))

  plt.autoscale()

  if cumulative:
    plt.ylim(0, 330)
  else:
    plt.ylim(0, 55)

    
  plt.xlim(datetime.date(2010, 9, 20), datetime.date(2010, 10, 31))

  plt.locator_params(axis='x', nbins=5)

  # months = DayLocator([4, 10, 16, 22, 29])
  # months = DayLocator([1, 8, 15, 22, 29])
  months = WeekdayLocator(mdates.FR)
  

  monthsFmt = DateFormatter("%b '%y")
  plt.gca().xaxis.set_major_locator(months)
  # plt.gca().xaxis.set_major_formatter(monthsFmt)

  # plt.gca().set_xticks(plt.gca().get_xticks()[1:])
  
  plt.tick_params(which='major', width=5, length=25, labelsize=70)
  plt.tick_params(which='minor', width=3, length=15)
  
  # plt.gca().xaxis.set_minor_locator(MultipleLocator(0.02))

  if cumulative:
    plt.gca().yaxis.set_minor_locator(MultipleLocator(10))
  else:
    plt.gca().yaxis.set_minor_locator(MultipleLocator(2))


  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  outside_text = plt.gca().legend( [extra], ["CMS 2010 Open Data"], frameon=0, borderpad=0, fontsize=50, bbox_to_anchor=(1.0, 1.005), loc='lower right')
  plt.gca().add_artist(outside_text)

  # plt.xlim()

  if cumulative:
    logo = [0.062, 0.985]
  else:
    logo = [0.051, 0.978]

  logo_offset_image = OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
  text_box = TextArea("", textprops=dict(color='#444444', fontsize=50, weight='bold'))
  logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)
  anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0., bbox_to_anchor=logo, bbox_transform = plt.gcf().transFigure)

  plt.gca().add_artist(anchored_box)

  plt.tight_layout()
  plt.gcf().set_size_inches(30, 24, forward=1)

  if cumulative:
    plt.savefig("plots/Version 6/lumi_cumulative.pdf")
  else:
    plt.savefig("plots/Version 6/lumi.pdf")

  plt.clf()
   








def plot_inst_lumi():
  properties = parse_file(input_analysis_file)

  timestamps = sorted(properties['time'])
  inst_lumi = [x for (y,x) in sorted(zip(properties['time'], properties['avg_inst_lumi']))]

  max_lumi = max(inst_lumi)

  print "Max Inst Luminosity: {}".format(max_lumi)

  inst_lumi = [x * 1e-3 for x in inst_lumi]

  dates = [datetime.datetime.fromtimestamp( int(time) ) for time in timestamps]
  # dates = [datetime.datetime.fromtimestamp( int(time) ).strftime('%m-%d') for time in timestamps]

  
  label = "Max. Inst. Lumi.: " + str(round(max_lumi, 2)) + " Hz/$\mathrm{\mu b}$"

  plt.hist(mpl.dates.date2num(dates), label=label, weights=inst_lumi, bins=25, color='orange', edgecolor='darkorange')

  years = mdates.YearLocator()   # every year
  months = mdates.MonthLocator()  # every month
  yearsFmt = mdates.DateFormatter('%Y-%M')

  # format the ticks
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
  plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

  plt.xlabel("Date (UTC)", labelpad=40, fontsize=60)

  
  plt.ylabel("Peak Delivered Luminosity (Hz/$\mathrm{\mu b}$)", labelpad=50, fontsize=60)


  
  plt.legend(bbox_to_anchor=[0.007, 0.85], frameon=False, loc='upper left', fontsize=50)

  plt.gcf().set_size_inches(30, 21.4285714, forward=1)

  plt.gcf().autofmt_xdate()

  plt.autoscale()
  plt.xlim(datetime.date(2010, 3, 30), datetime.date(2010, 10, 31))

  # plt.gca().xaxis.set_minor_locator(MultipleLocator(0.02))

  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  outside_text = plt.gca().legend( [extra], ["CMS 2010 Open Data"], frameon=0, borderpad=0, fontsize=50, bbox_to_anchor=(1.0, 1.005), loc='lower right')
  plt.gca().add_artist(outside_text)
  
  plt.tick_params(which='major', width=5, length=15, labelsize=60)

  ab = AnnotationBbox(OffsetImage(read_png(get_sample_data("/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.15, resample=1, dpi_cor=1), (0.205, 0.840), xycoords='figure fraction', frameon=0)
  plt.gca().add_artist(ab)
  preliminary_text = "Prelim. (20\%)"
  plt.gcf().text(0.270, 0.825, preliminary_text, fontsize=60, weight='bold', color='#444444', multialignment='center')


  plt.savefig("plots/Version 6/lumi/inst_lumi_number.pdf")

  plt.clf()




plot_integrated_recorded_lumi(cumulative=True)
plot_integrated_recorded_lumi(cumulative=False)


# plot_inst_lumi(number=False)
# plot_inst_lumi(number=True)





# call(["python", "/home/aashish/root/macros/MODAnalyzer/utilities/sync_plots.py"])