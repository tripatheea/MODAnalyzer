from __future__ import division
import numpy as np
import copy


# RootPy
from rootpy.plotting import Hist


class MODHist:

	def __init__(self, hist, conditions=[], use_prescale=True, x_lims=(0, -1), y_lims=(0, -1)):
		# hist should be a RootPy Hist object (TH1). 
		# conditions should be a list of lambda functions that need to be satisfied by 'value' for this hist to be filled. 

		self._hist = hist
		self._conditions = conditions
		self._use_prescale = use_prescale

		self._x_lims = x_lims
		self._y_lims = y_lims

	def hist(self):
		return self._hist

	def conditions(self):
		return self._conditions

	def use_prescale(self):
		return self._use_prescale

	def x_lims(self):
		return self._x_lims

	def y_lims(self):
		return self._y_lims





def all_hist_templates():
	all_hists = {}

	'''
	# all_hists['zg'] = Hist(25, 0.0, 0.5, title="zg", color="black")
	all_hists['hardest_pT'] = Hist(100, 5, 1005, title="pT", color="black")
	all_hists['hardest_eta'] = Hist(25, -5, 5, title="eta", color="black")
	all_hists['hardest_phi'] = Hist(25, 0, 2*np.pi, title="phi", color="black")
	
	all_hists['mul_pre_SD'] = Hist(50, -1, 99, title="constituent_multiplicity", color="black")

	all_hists['zg_10'] = Hist(50, 0., 0.5, title="zg", color="black")
	'''

	hardest_pT_hist = Hist(100, 5, 1005, title="pT")
	hardest_eta_hist = Hist(25, -5, 5, title="eta")

	# pT_boundaries = [85, 115, 150, 200, 250, 325]
	pT_boundaries = [85, 115]

	all_hists['hardest_pT'] = []
	all_hists['hardest_eta'] = []

	for i in range(len(pT_boundaries) - 1):
		all_hists['hardest_pT'].append( copy.deepcopy( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_pT', lambda x: x > pT_boundaries[i] and x < pT_boundaries[i + 1])], use_prescale=False, x_lims=(pT_boundaries[i], pT_boundaries[i + 1])) ) )
		all_hists['hardest_eta'].append( copy.deepcopy( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_pT', lambda x: x > pT_boundaries[i] and x < pT_boundaries[i + 1])], use_prescale=False) ) )


	# all_hists['hardest_pT'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_pT', lambda x: x > 85)], use_prescale=True) )
	all_hists['hardest_eta'].append( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_eta', lambda x: x > 85)], use_prescale=True) )
	
	# all_hists['hardest_pT'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_pT', lambda x: x > 150)], use_prescale=True) )
	all_hists['hardest_eta'].append( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_eta', lambda x: x > 150)], use_prescale=True) )
	 
	
	return all_hists

def get_hist_template(var):
	all_hists = all_hist_templates()
	return all_hists[var]