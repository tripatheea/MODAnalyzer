from __future__ import division




# RootPy
from rootpy.plotting import Hist

def all_hist_templates():
	all_hists = {}

	all_hists['zg'] = Hist(25, 0.0, 0.5, markersize=2.5)
	all_hists['pT'] = Hist(25, 0, 1000, markersize=2.5)
	all_hists['eta'] = Hist(25, -5, 5, markersize=2.5)

	return all_hists

def get_hist_template(var):
	all_hists = all_hist_templates()
	return all_hists[var]