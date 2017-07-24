from __future__ import division


import hists


from MODPlot import *


import two_dim_parse


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


mpl.rcParams['axes.linewidth'] = 5.0  # set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)


default_dir = "plots/Version 6/"


logo_location = "/home/aashish/root/macros/MODAnalyzer/mod_logo.png"
logo_text = ""


parsed_linear = two_dim_parse.load_root_files_to_hist(log=False)
parsed_log = two_dim_parse.load_root_files_to_hist(log=True)


def normalize_hist(hist):
    bin_width = (hist.upperbound() - hist.lowerbound()) / hist.nbins()

    if hist.GetSumOfWeights() != 0.0:
        hist.Scale(1.0 / (hist.GetSumOfWeights() * bin_width))

    return hist


def two_dim_plots(track=False):

    colors = ['Greys', 'Blues', 'Greens', 'purple']
    hatch_colors = ['gray', 'blue', 'green', 'purple']
    sources = ['data', 'pythia', 'herwig', 'sherpa']
    source_labels = ["CMS 2010 Open Data", "Pythia 8.219", "Herwig 7.0.3", "Sherpa 2.2.1"]

    if track:
        sources = [sources[0]]

    # colors = ['Greys']
    # hatch_colors = ['gray']
    # sources = ['data']
    # source_labels = ["CMS 2010 Open Data"]

    startcolor = 'white'  # a dark olive
    endcolor = 'purple'    # medium dark red
    cmap2 = mpl.colors.LinearSegmentedColormap.from_list('purple', [startcolor, endcolor])
    cm.register_cmap(cmap=cmap2)

    lower_boundaries = [85, 115, 150, 200, 85, 150, 250]
    upper_boundaries = [115, 150, 200, 250, 100000., 100000., 100000.]

    lambda_value = 2.
    z_cut = 0.1

    a = 0
    for a in range(len(sources)):
        color, source, source_label, hatch_color = colors[a], sources[a], source_labels[a], hatch_colors[a]

        if track:
            filename = "big5_zg_vs_rg_" + source + "_track_linear.pdf"
        else:
            filename = "big5_zg_vs_rg_" + source + "_linear.pdf"

        with PdfPages("plots/Version 5/zg_against_theta_g/" + filename) as pdf:

            if track:
                var = ('track_zg_10', 'track_rg_10')
            else:
                var = ('zg_10', 'rg_10')

            for b in range(len(parsed_linear[a][var])):

                lower, upper = lower_boundaries[b], upper_boundaries[b]

                np_correction_boundary = lambda_value / (lower * z_cut)
                # print np_correction_boundary

                # print len(parsed_linear),
                # print len(parsed_linear[a][('zg_10', 'rg_10')])

                hist = parsed_linear[a][var][b].hist()

                b += 1

                x_s = []
                y_s = []
                z_s = []

                for i in range(1, hist.nbins(0) + 1):
                    for j in range(1, hist.nbins(1) + 1):

                        z = hist.GetBinContent(i, j)
                        x = hist.GetXaxis().GetBinCenter(i)
                        y = hist.GetYaxis().GetBinCenter(j)

                        x_s.append(x)
                        y_s.append(y)
                        z_s.append(z)

                H, xedges, yedges = np.histogram2d(x_s, y_s, bins=[25, 25], range=[[0.0, 0.5], [0.0, 1.0]], weights=z_s, normed=True)

                H = np.array(H)

                H = np.rot90(H)
                H = np.flipud(H)
                Hmasked = np.ma.masked_where(H == 0, H)  # Mask pixels with a value of zero

                plt.pcolor(xedges, yedges, Hmasked, cmap=color, vmin=0, vmax=20)

                cbar = plt.colorbar(ticks=[4 * i for i in range(6)])
                cbar.ax.tick_params(labelsize=70)
                # cbar.ax.set_ylabel('$\\frac{1}{\sigma} \\frac{\mathrm{d}^2 \sigma}{\mathrm{d} z_g \mathrm{d} \\theta_g}$', labelpad=150, fontsize=105, rotation=0)
                # cbar.ax.set_ylim(0, 24)

                if upper != 100000.:
                    # $p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$; AK5 \n $\left| \eta \\right| < 2.4$; $p_T^{\mathrm{jet}} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$
                    label = "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$\n AK5; $\left| \eta \\right| < 2.4$ \n $p_T^{\mathrm{jet}} \in [" + str(
                        lower) + ", " + str(upper) + "]~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$"
                else:
                    label = "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$\n AK5; $\left| \eta \\right| < 2.4$ \n $p_T^{\mathrm{jet}} >" + \
                        str(lower) + "~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$"

                extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)

                labels = label.split("\n")
                additional_legend = plt.gca().legend([extra] * len(labels), labels, frameon=0, borderpad=0,
                                                     fontsize=50, bbox_to_anchor=[0.95, 0.98], loc="upper right")
                plt.gca().add_artist(additional_legend)

                shift = max([t.get_window_extent(plt.gcf().canvas.get_renderer()).width for t in additional_legend.get_texts()])

                for t in additional_legend.get_texts():
                    t.set_ha('right')  # ha is alias for horizontalalignment
                    t.set_position((shift, 0))

                # Data Source Label.
                label = []
                label.extend([source_label])
                additional_legend = plt.gca().legend([extra] * len(label), label, frameon=0, borderpad=0,
                                                     fontsize=60, bbox_to_anchor=[1.02, 1.08], loc="upper right")
                plt.gca().add_artist(additional_legend)

                if track:
                    plt.xlabel('Track $z_g$', fontsize=90, labelpad=40)
                    plt.ylabel('Track $\\theta_g$', rotation=90, fontsize=90, labelpad=40)
                else:
                    plt.xlabel('$z_g$', fontsize=90, labelpad=40)
                    plt.ylabel('$\\theta_g$', rotation=0, fontsize=90, labelpad=40)

                # plt.gca().xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=False))
                # plt.gca().yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=False))

                plt.tick_params(which='major', width=5, length=25, labelsize=70)
                plt.tick_params(which='minor', width=3, length=15)

                plt.gca().xaxis.set_minor_locator(MultipleLocator(0.02))
                plt.gca().yaxis.set_minor_locator(MultipleLocator(0.02))

                plt.xlim(0.0, 0.5)
                plt.ylim(0.0, 1.0)

                logo_offset_image = OffsetImage(read_png(get_sample_data(
                    "/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
                text_box = TextArea(logo_text, textprops=dict(color='#444444', fontsize=50, weight='bold'))
                logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

                if track:
                    anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.,
                                                     bbox_to_anchor=[0.109, 0.985], bbox_transform=plt.gcf().transFigure)
                else:
                    anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.,
                                                     bbox_to_anchor=[0.086, 0.985], bbox_transform=plt.gcf().transFigure)

                plt.gca().add_artist(anchored_box)

                plt.gcf().set_size_inches(30, 25, forward=1)
                plt.gcf().set_snap(True)
                plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

                # rplt.hist2d(hist)

                # plt.savefig("plots/Version 5/zg_against_theta_g/linear/" + source + "_zg_against_theta_g.png")
                pdf.savefig()
                plt.close()


def two_dim_log_plots(track=False):

    colors = ['Greys', 'Blues', 'Greens', 'purple']
    hatch_colors = ['gray', 'blue', 'green', 'purple']
    sources = ['data', 'pythia', 'herwig', 'sherpa']
    source_labels = ["CMS 2010 Open Data", "Pythia 8.219", "Herwig 7.0.3", "Sherpa 2.2.1"]

    if track:
        sources = [sources[0]]

    startcolor = 'white'  # a dark olive
    endcolor = 'purple'    # medium dark red
    cmap2 = mpl.colors.LinearSegmentedColormap.from_list('purple', [startcolor, endcolor])
    cm.register_cmap(cmap=cmap2)
    lower_boundaries = [85, 115, 150, 200, 85, 150, 250]
    upper_boundaries = [115, 150, 200, 250, 100000., 100000., 100000.]

    lambda_value = 2.
    z_cut = 0.1

    a = 0
    for a in range(len(sources)):
        color, source, source_label, hatch_color = colors[a], sources[a], source_labels[a], hatch_colors[a]

        if track:
            filename = "big5_zg_vs_rg_" + source + "_track_log.pdf"
        else:
            filename = "big5_zg_vs_rg_" + source + "_log.pdf"

        with PdfPages("plots/Version 5/zg_against_theta_g/" + filename) as pdf:

            if track:
                var = ('track_zg_10', 'track_rg_10')
            else:
                var = ('zg_10', 'rg_10')

            for b in range(len(parsed_log[a][var])):

                lower, upper = lower_boundaries[b], upper_boundaries[b]
                np_correction_boundary = lambda_value / (lower * z_cut)

                # print len(parsed_log),
                # print len(parsed_log[a][('zg_10', 'rg_10')])

                hist = parsed_log[a][var][b].hist()

                b += 1

                x_s = []
                y_s = []
                z_s = []

                np_region_x_s, np_region_y_s, np_region_z_s = [], [], []

                # print
                total_integral = hist.integral()

                x_bounds = hist.bounds(axis=0)
                y_bounds = hist.bounds(axis=1)

                area = (math.log(x_bounds[1], math.e) - math.log(x_bounds[0], math.e)) * \
                    (math.log(y_bounds[1], math.e) - math.log(y_bounds[0], math.e))

                sum_of_bin_contents = sum([hist.GetBinContent(i, j) for i in range(1, hist.nbins(0) + 1) for j in range(1, hist.nbins(1) + 1)])

                # print sum_of_bin_contents, total_integral

                # for i in range(1, hist.nbins(0) + 2):
                for i in range(1, hist.nbins(0) + 1):
                    # for j in range(1, hist.nbins(1) + 2):
                    for j in range(1, hist.nbins(1) + 1):

                        z = hist.nbins(0) * hist.nbins(1) * hist.GetBinContent(i, j) / (total_integral * area)
                        # x = hist.GetXaxis().GetBinCenter(i) - (hist.GetXaxis().GetBinWidth(i) / 2)
                        # y = hist.GetYaxis().GetBinCenter(j) - (hist.GetYaxis().GetBinWidth(j) / 2)

                        # print z

                        x = hist.GetXaxis().GetBinCenter(i)
                        y = hist.GetYaxis().GetBinCenter(j)

                        x_s.append(x)
                        y_s.append(y)
                        z_s.append(z)

                        # if i == 25:
                        # 	print x, y

                # print i, j

                # print max(z_s)

                # DO NOT USE NORMED=True HERE.
                H, xedges, yedges = np.histogram2d(x_s, y_s, bins=[np.logspace(math.log(float(0.1), math.e), math.log(
                    0.5, math.e), (25 + 1), base=np.e), np.logspace(math.log(float(0.01), math.e), math.log(1.0, math.e), (25 + 1), base=np.e)], weights=z_s)

                H_normalized = np.array(H)
                H = H_normalized

                H = np.rot90(H)
                H = np.flipud(H)
                Hmasked = np.ma.masked_where(H == 0, H)  # Mask pixels with a value of zero

                plt.pcolor(xedges, yedges, Hmasked, cmap=color, vmin=0, vmax=0.5)

                cbar = plt.colorbar(ticks=[0.1 * i for i in range(6)])
                cbar.ax.tick_params(labelsize=70)
                # cbar.ax.set_ylabel('$\\frac{z_g \\theta_g}{\sigma} \\frac{\mathrm{d}^2 \sigma}{\mathrm{d} z_g \mathrm{d} \\theta_g}$', labelpad=150, fontsize=105, rotation=0)

                label = []
                if upper != 100000.:
                    # $p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$; AK5 \n $\left| \eta \\right| < 2.4$; $p_T^{\mathrm{jet}} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$
                    label = "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $\left| \eta \\right| < 2.4$ \n $p_T^{\mathrm{jet}} \in [" + str(
                        lower) + ", " + str(upper) + "]~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$"
                else:
                    label = "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $\left| \eta \\right| < 2.4$ \n $p_T^{\mathrm{jet}} >" + \
                        str(lower) + "~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$"

                extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)

                labels = label.split("\n")
                additional_legend = plt.gca().legend([extra] * len(labels), labels, frameon=0, borderpad=0,
                                                     fontsize=50, bbox_to_anchor=[1.41, 0.00], loc="lower right")

                plt.gca().add_artist(additional_legend)

                for t in additional_legend.get_texts():
                    t.set_ha('right')  # ha is alias for horizontalalignment

                # Data Source Label.
                label = []
                label.extend([source_label])
                additional_legend = plt.gca().legend([extra] * len(label), label, frameon=0, borderpad=0,
                                                     fontsize=60, bbox_to_anchor=[1.02, 1.08], loc="upper right")
                plt.gca().add_artist(additional_legend)

                if track:
                    plt.xlabel('Track $z_g$', fontsize=90, labelpad=40)
                    plt.ylabel('Track $\\theta_g$', rotation=90, fontsize=90, labelpad=40)
                else:
                    plt.xlabel('$z_g$', fontsize=90, labelpad=40)
                    plt.ylabel('$\\theta_g$', rotation=0, fontsize=90, labelpad=40)

                plt.tick_params(which='major', width=5, length=25, labelsize=70)
                plt.tick_params(which='minor', width=3, length=15)

                plt.xscale('log')
                plt.yscale('log')

                plt.xlim(0.1, 0.5)
                plt.ylim(0.01, 1.0)

                plt.gca().xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=False))
                plt.gca().yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=False))

                # plt.xticks([0.1, 0.2, 0.5, 1.0])
                plt.xticks([0.1, 0.2, 0.3, 0.4, 0.5])
                plt.yticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0])

                logo_offset_image = OffsetImage(read_png(get_sample_data(
                    "/home/aashish/root/macros/MODAnalyzer/mod_logo.png", asfileobj=False)), zoom=0.25, resample=1, dpi_cor=1)
                text_box = TextArea(logo_text, textprops=dict(color='#444444', fontsize=50, weight='bold'))
                logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

                if track:
                    anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.,
                                                     bbox_to_anchor=[0.127, 0.985], bbox_transform=plt.gcf().transFigure)
                else:
                    anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0.,
                                                     bbox_to_anchor=[0.104, 0.985], bbox_transform=plt.gcf().transFigure)

                plt.gca().add_artist(anchored_box)

                plt.gcf().set_size_inches(30, 25, forward=1)
                plt.gcf().set_snap(True)
                plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)

                # rplt.hist2d(hist)

                # plt.savefig("plots/Version 5/zg_against_theta_g/log/" + source + "_zg_against_theta_g.png")
                pdf.savefig()
                plt.close()


start = time.time()


# two_dim_plots(track=False)
# two_dim_plots(track=True)

two_dim_log_plots(track=False)
two_dim_log_plots(track=True)

end = time.time()

print "Finished all plotting in {} seconds.".format(end - start)
