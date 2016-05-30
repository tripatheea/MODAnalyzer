from __future__ import division
import numpy as np



# RootPy
from rootpy.plotting import Hist

def all_hist_templates():
	all_hists = {}

	# all_hists['zg'] = Hist(25, 0.0, 0.5, title="zg", color="black")
	all_hists['hardest_pT'] = Hist(100, 5, 1005, title="pT", color="black")
	all_hists['hardest_eta'] = Hist(25, -5, 5, title="eta", color="black")
	all_hists['hardest_phi'] = Hist(25, 0, 2*np.pi, title="phi", color="black")

	return all_hists

def get_hist_template(var):
	all_hists = all_hist_templates()
	return all_hists[var]