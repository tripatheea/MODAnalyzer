from __future__ import division
import numpy as np
import copy
import math

from collections import defaultdict
# RootPy
from rootpy.plotting import Hist, Hist2D


class MODHist:

	def __init__(self, hist, conditions=[], use_prescale=True, x_label="", y_label="", x_scale='linear', y_scale='linear', mark_regions=[], x_range=(0, -1), y_range=(0, -1), additional_text=[], legend_location=(1., 1.), axes_label_pi=False):
		# hist should be a RootPy Hist object (TH1). 
		# conditions should be a list of lambda functions that need to be satisfied by 'value' for this hist to be filled. 

		self._hist = hist
		self._conditions = conditions
		self._use_prescale = use_prescale

		self._x_label = x_label
		self._y_label = y_label

		self._x_scale = x_scale
		self._y_scale = y_scale

		self._mark_regions = mark_regions

		self._x_range = x_range
		self._y_range = y_range

		self._additional_text = additional_text

		self._axes_label_pi = axes_label_pi

		self._legend_location = legend_location

	def replace_hist(self, hist):
		self._hist = hist

	def hist(self):
		return self._hist

	def conditions(self):
		return self._conditions

	def use_prescale(self):
		return self._use_prescale

	def x_label(self):
		return self._x_label

	def y_label(self):
		return self._y_label

	def x_scale(self):
		return self._x_scale

	def y_scale(self):
		return self._y_scale

	def x_range(self):
		return self._x_range

	def y_range(self):
		return self._y_range

	def additional_text(self):
		return self._additional_text

	def axes_label_pi(self):
		return self._axes_label_pi

	def mark_regions(self):
		return self._mark_regions

	def legend_location(self):
		return self._legend_location



def multi_page_log_plot_hist_templates():
	all_hists = {}
	

	hardest_pT_hist = Hist(np.logspace(math.log(float(5), math.e), math.log(2005, math.e), 50, base=np.e))
	pT_D_hist = Hist(np.logspace(math.log(float(0.01), math.e), math.log(1, math.e), 50, base=np.e))
	
		
	frac_pT_loss_hist = Hist(np.logspace(math.log(float(0.001), math.e), math.log(1.0, math.e), 50, base=np.e))
	softkill_pT_loss_hist = Hist(np.logspace(math.log(float(0.0001), math.e), math.log(0.1, math.e), 50, base=np.e))


	lha_hist =  Hist(np.logspace(math.log(float(0.1), math.e), math.log(1.0, math.e), 50, base=np.e))
	width_hist =  Hist(np.logspace(math.log(float(0.01), math.e), math.log(1.0, math.e), 50, base=np.e))
	thrust_hist =  Hist(np.logspace(math.log(float(0.001), math.e), math.log(1.0, math.e), 50, base=np.e))


	zg_10_hist = Hist(np.logspace(math.log(float(0.01), math.e), math.log(0.5, math.e), 50, base=np.e))
	rg_10_hist = Hist(np.logspace(math.log(float(0.01), math.e), math.log(1.0, math.e), 50, base=np.e))
	e1_10_hist = Hist(np.logspace(math.log(float(0.001), math.e), math.log(1.0, math.e), 50, base=np.e))
	e05_10_hist = Hist(np.logspace(math.log(float(0.01), math.e), math.log(1.0, math.e), 50, base=np.e))
	e2_10_hist = Hist(np.logspace(math.log(float(0.0001), math.e), math.log(1.0, math.e), 50, base=np.e))





	zg_05_hist = Hist(np.logspace(math.log(float(0.05), math.e), math.log(0.5, math.e), 50, base=np.e))
	rg_05_hist = Hist(np.logspace(math.log(float(0.005), math.e), math.log(1.0, math.e), 50, base=np.e))
	e1_05_hist = Hist(np.logspace(math.log(float(0.0005), math.e), math.log(1.0, math.e), 50, base=np.e))
	e05_05_hist = Hist(np.logspace(math.log(float(0.005), math.e), math.log(1.0, math.e), 50, base=np.e))
	e2_05_hist = Hist(np.logspace(math.log(float(0.00005), math.e), math.log(1.0, math.e), 50, base=np.e))





	zg_20_hist = Hist(np.logspace(math.log(float(0.01), math.e), math.log(0.5, math.e), 50, base=np.e))
	rg_20_hist = Hist(np.logspace(math.log(float(0.01), math.e), math.log(1.0, math.e), 50, base=np.e))
	e1_20_hist = Hist(np.logspace(math.log(float(0.001), math.e), math.log(1.0, math.e), 50, base=np.e))
	e05_20_hist = Hist(np.logspace(math.log(float(0.01), math.e), math.log(1.0, math.e), 50, base=np.e))
	e2_20_hist = Hist(np.logspace(math.log(float(0.0001), math.e), math.log(1.0, math.e), 50, base=np.e))


	pT_boundaries = [85, 115, 150, 200, 250]


	all_hists['frac_pT_loss'] = []
	all_hists['softkill_pT_loss'] = []


	all_hists['LHA_pre_SD'], all_hists['track_LHA_pre_SD'] = [], []
	all_hists['LHA_post_SD'], all_hists['track_LHA_post_SD'] = [], []
	all_hists['width_pre_SD'], all_hists['track_width_pre_SD'] = [], []
	all_hists['width_post_SD'], all_hists['track_width_post_SD'] = [], []
	all_hists['thrust_pre_SD'], all_hists['track_thrust_pre_SD'] = [], []
	all_hists['thrust_post_SD'], all_hists['track_thrust_post_SD'] = [], []


	all_hists['zg_10'], all_hists['track_zg_10'] = [], []
	all_hists['rg_10'], all_hists['track_rg_10'] = [], []
	all_hists['e1_10'], all_hists['track_e1_10'] = [], []
	all_hists['e2_10'], all_hists['track_e2_10'] = [], []
	all_hists['e05_10'], all_hists['track_e05_10'] = [], []

	
	all_hists['zg_05'], all_hists['track_zg_05'] = [], []
	all_hists['rg_05'], all_hists['track_rg_05'] = [], []
	all_hists['e1_05'], all_hists['track_e1_05'] = [], []
	all_hists['e2_05'], all_hists['track_e2_05'] = [], []
	all_hists['e05_05'], all_hists['track_e05_05'] = [], []


	all_hists['zg_20'], all_hists['track_zg_20'] = [], []
	all_hists['rg_20'], all_hists['track_rg_20'] = [], []
	all_hists['e1_20'], all_hists['track_e1_20'] = [], []
	all_hists['e2_20'], all_hists['track_e2_20'] = [], []
	all_hists['e05_20'], all_hists['track_e05_20'] = [], []

	




	all_hists['hardest_pT'] = []
	all_hists['track_pT_D_pre_SD'] = []
	all_hists['track_pT_D_post_SD'] = []
	all_hists['pT_D_pre_SD'] = []
	all_hists['pT_D_post_SD'] = []
	


	for i in range(len(pT_boundaries) - 1):

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['hardest_pT'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$p_T$", y_label="A.U.", y_scale='log', x_range=(pT_boundaries[i], pT_boundaries[i + 1]), y_range=(1e-5, 1e2), additional_text=additional_text ) ) 




		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['frac_pT_loss'].append( MODHist(frac_pT_loss_hist, conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Fractional $p_T$ Loss", y_label="A.U.", y_range=(0., 1.2), additional_text=additional_text ) ) 


		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['softkill_pT_loss'].append( MODHist(softkill_pT_loss_hist, conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, y_range=(0., 0.9), x_label="Softkill $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 



		additional_text = [ ( (-0.08, 1.0), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$p_T^D$", y_label="$\\frac{p_T^D}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", y_range=(0., 1.2), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 1.0), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $p_T^D$", y_label="$\\frac{p_T^D}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", y_range=(0., 1.2), additional_text=additional_text ) ) 



		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$p_T^D$", y_label="$\\frac{p_T^D}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $p_T^D$", y_label="$\\frac{p_T^D}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 






		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, y_range=(0.0, 1.5), x_label="LHA", y_label="$\\frac{e^{(1/2)}}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', y_range=(0.0, 1.5), use_prescale=False, x_label="Track LHA", y_label="$\\frac{e^{(1/2)}}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 



		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', y_range=(0.0, 1.6), use_prescale=False, x_label="LHA", y_label="$\\frac{e^{(1/2)}}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track LHA", y_label="$\\frac{e^{(1/2)}}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 




		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, y_range=(0.0, 1.0), x_label="Width", y_label="$\\frac{e^{(1)}}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, y_range=(0.0, 1.0), x_label="Track Width", y_label="$\\frac{e^{(1)}}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 


		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', y_range=(0.0, 1.0), use_prescale=False, x_label="Width", y_label="$\\frac{e^{(1)}}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track Width", y_label="$\\frac{e^{(1)}}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 




		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, y_range=(0.0, 0.8), x_label="Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, y_range=(0.0, 0.6), x_label="Track Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 


		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', y_range=(0.0, 0.75), use_prescale=False, x_label="Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 



		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['zg_10'].append( MODHist(copy.deepcopy(zg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.01, 1.0), y_range=(0, 2), additional_text=additional_text ) ) 
	
		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['rg_10'].append( MODHist(copy.deepcopy(rg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['e1_10'].append( MODHist(copy.deepcopy(e1_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$e_g^1$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0.0, 0.75), x_range=(0.001, 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['e2_10'].append( MODHist(copy.deepcopy(e2_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0.0, 0.45), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['e05_10'].append( MODHist(copy.deepcopy(e05_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0.0, 1.4), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_zg_10'].append( MODHist(copy.deepcopy(zg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.01, 1.0), y_range=(0, 2), additional_text=additional_text ) ) 
	
		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_rg_10'].append( MODHist(copy.deepcopy(rg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_e1_10'].append( MODHist(copy.deepcopy(e1_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0, 0.75), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_e2_10'].append( MODHist(copy.deepcopy(e2_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), y_range=(0, 0.45), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_e05_10'].append( MODHist(copy.deepcopy(e05_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g a _g})}$", x_range=(0.01, 1.0), y_range=(0, 1.4), additional_text=additional_text ) ) 






		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 1.7), additional_text=additional_text ) ) 
	
		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.3), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.6), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0., 1.3), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 3), additional_text=additional_text ) ) 
	
		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 2.), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g a _g})}$", x_range=(0.01, 1.0), additional_text=additional_text ) ) 






		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['zg_20'].append( MODHist(copy.deepcopy(zg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.2, 1.0), y_range=(0, 2.5), additional_text=additional_text ) ) 
	
		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['rg_20'].append( MODHist(copy.deepcopy(rg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 0.8), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['e1_20'].append( MODHist(copy.deepcopy(e1_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.8), x_range=(0.001, 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['e2_20'].append( MODHist(copy.deepcopy(e2_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.45), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['e05_20'].append( MODHist(copy.deepcopy(e05_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0., 1.4), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_zg_20'].append( MODHist(copy.deepcopy(zg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 3), additional_text=additional_text ) ) 
	
		additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_rg_20'].append( MODHist(copy.deepcopy(rg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 2.), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_e1_20'].append( MODHist(copy.deepcopy(e1_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0, 0.8), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_e2_20'].append( MODHist(copy.deepcopy(e2_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), y_range=(0, 0.45), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_e05_20'].append( MODHist(copy.deepcopy(e05_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], x_scale='log', use_prescale=False, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g ta _g})}$", x_range=(0.01, 1.0), y_range=(0, 1.5), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['hardest_pT'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_scale='log', x_label="$p_T$", y_label="A.U.", y_scale='log', x_range=(1e1, 2e3), y_range=(1e-9, 1e2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['hardest_pT'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_scale='log', x_label="$p_T$", y_label="A.U.", y_scale='log', x_range=(1e1, 2e3), y_range=(1e-9, 1e2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['hardest_pT'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_scale='log', x_label="$p_T$", y_label="A.U.", y_scale='log', x_range=(1e1, 2e3), y_range=(1e-9, 1e2), additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['frac_pT_loss'].append( MODHist(frac_pT_loss_hist, conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_scale='log', x_label="Fractional $p_T$ Loss", y_label="A.U.", y_range=(0., 1.2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['frac_pT_loss'].append( MODHist(frac_pT_loss_hist, conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_scale='log', x_label="Fractional $p_T$ Loss", y_label="A.U.", y_range=(0., 0.8), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['frac_pT_loss'].append( MODHist(frac_pT_loss_hist, conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_scale='log', x_label="Fractional $p_T$ Loss", y_label="A.U.", y_range=(0., 1.2), additional_text=additional_text ) ) 





	additional_text = [ ( (-0.08, 1.0), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >85~\mathrm{GeV}$" ) ]
	all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D}$", y_range=(0., 1.2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 1.0), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >150~\mathrm{GeV}$" ) ]
	all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D}$", y_range=(0., 1.2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 1.0), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >250~\mathrm{GeV}$" ) ]
	all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D}$", y_range=(0., 1.2), additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 1.0), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >85~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D$", y_range=(0., 1.2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 1.0), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >150~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D$", y_range=(0., 1.2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 1.0), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >250~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D$", y_range=(0., 1.2), additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >85~\mathrm{GeV}$" ) ]
	all_hists['pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >150~\mathrm{GeV}$" ) ]
	all_hists['pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >250~\mathrm{GeV}$" ) ]
	all_hists['pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >85~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >150~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D$", additional_text=additional_text ) ) 


	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} >250~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{p_T^D \mathrm{d} \sigma}{\sigma \mathrm{d} p_T^D$", additional_text=additional_text ) ) 







	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['softkill_pT_loss'].append( MODHist(softkill_pT_loss_hist, conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=False, x_label="Softkill $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['softkill_pT_loss'].append( MODHist(softkill_pT_loss_hist, conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=False, x_label="Softkill $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['softkill_pT_loss'].append( MODHist(softkill_pT_loss_hist, conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=False, x_label="Softkill $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="LHA", y_range=(0.0, 1.5), y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="LHA", y_range=(0.0, 1.5), y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="LHA", y_range=(0.0, 1.5), y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track LHA", y_range=(0.0, 1.5), y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track LHA", y_range=(0.0, 1.5), y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track LHA", y_range=(0.0, 1.5), y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, y_range=(0.0, 1.6), x_label="LHA", y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, y_range=(0.0, 1.6), x_label="LHA", y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, y_range=(0.0, 1.6), x_label="LHA", y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track LHA", y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track LHA", y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track LHA", y_label="$\\frac{e^{(1/2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Width", y_range=(0.0, 1.0), y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Width", y_range=(0.0, 1.0), y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Width", y_range=(0.0, 1.0), y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, y_range=(0.0, 1.0), x_label="Track Width", y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, y_range=(0.0, 1.0), x_label="Track Width", y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, y_range=(0.0, 1.0), x_label="Track Width", y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 





	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Width", y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Width", y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Width", y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track Width", y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track Width", y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track Width", y_label="$\\frac{e^{(1)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, y_range=(0.0, 0.8), x_label="Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, y_range=(0.0, 0.6), x_label="Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, y_range=(0.0, 0.8), x_label="Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, y_range=(0.0, 0.6), x_label="Track Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, y_range=(0.0, 0.6), x_label="Track Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, y_range=(0.0, 0.6), x_label="Track Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 





	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, y_range=(0.0, 0.75), x_label="Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, y_range=(0.0, 0.75), x_label="Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, y_range=(0.0, 0.75), x_label="Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 





	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track Thrust", y_label="$\\frac{e^{(2)} \mathrm{d} \sigma}{\sigma \mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 





	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['zg_10'].append( MODHist(copy.deepcopy(zg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.01, 1.0), y_range=(0, 2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['rg_10'].append( MODHist(copy.deepcopy(rg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e1_10'].append( MODHist(copy.deepcopy(e1_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0.0, 0.75), x_range=(0.001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e2_10'].append( MODHist(copy.deepcopy(e2_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0.0, 0.45), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e05_10'].append( MODHist(copy.deepcopy(e05_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0.0, 1.4), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_zg_10'].append( MODHist(copy.deepcopy(zg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.01, 1.0), y_range=(0, 2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_rg_10'].append( MODHist(copy.deepcopy(rg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e1_10'].append( MODHist(copy.deepcopy(e1_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0, 0.75), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e2_10'].append( MODHist(copy.deepcopy(e2_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), y_range=(0, 0.45), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e05_10'].append( MODHist(copy.deepcopy(e05_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g a _g})}$", x_range=(0.01, 1.0), y_range=(0, 1.4), additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 1.7), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.3), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.6), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0., 1.3), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 3), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 2.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g a _g})}$", x_range=(0.01, 1.0), additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['zg_20'].append( MODHist(copy.deepcopy(zg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.2, 1.0), y_range=(0, 2.5), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['rg_20'].append( MODHist(copy.deepcopy(rg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 0.8), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e1_20'].append( MODHist(copy.deepcopy(e1_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.8), x_range=(0.001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e2_20'].append( MODHist(copy.deepcopy(e2_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.45), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e05_20'].append( MODHist(copy.deepcopy(e05_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0., 1.4), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_zg_20'].append( MODHist(copy.deepcopy(zg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 3), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_rg_20'].append( MODHist(copy.deepcopy(rg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 2.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e1_20'].append( MODHist(copy.deepcopy(e1_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0, 0.8), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e2_20'].append( MODHist(copy.deepcopy(e2_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), y_range=(0, 0.45), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e05_20'].append( MODHist(copy.deepcopy(e05_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g ta _g})}$", x_range=(0.01, 1.0), y_range=(0, 1.5), additional_text=additional_text ) ) 





	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['zg_10'].append( MODHist(copy.deepcopy(zg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.01, 1.0), y_range=(0, 2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['rg_10'].append( MODHist(copy.deepcopy(rg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e1_10'].append( MODHist(copy.deepcopy(e1_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0.0, 0.75), x_range=(0.001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e2_10'].append( MODHist(copy.deepcopy(e2_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0.0, 0.45), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e05_10'].append( MODHist(copy.deepcopy(e05_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0.0, 1.4), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_zg_10'].append( MODHist(copy.deepcopy(zg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.01, 1.0), y_range=(0, 2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_rg_10'].append( MODHist(copy.deepcopy(rg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e1_10'].append( MODHist(copy.deepcopy(e1_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0, 0.75), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e2_10'].append( MODHist(copy.deepcopy(e2_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), y_range=(0, 0.45), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e05_10'].append( MODHist(copy.deepcopy(e05_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g a _g})}$", x_range=(0.01, 1.0), y_range=(0, 1.4), additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 1.7), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.3), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.6), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0., 1.3), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 3), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 2.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g a _g})}$", x_range=(0.01, 1.0), additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['zg_20'].append( MODHist(copy.deepcopy(zg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.2, 1.0), y_range=(0, 2.5), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['rg_20'].append( MODHist(copy.deepcopy(rg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 0.8), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e1_20'].append( MODHist(copy.deepcopy(e1_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.8), x_range=(0.001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e2_20'].append( MODHist(copy.deepcopy(e2_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.45), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e05_20'].append( MODHist(copy.deepcopy(e05_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0., 1.4), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_zg_20'].append( MODHist(copy.deepcopy(zg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 3), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_rg_20'].append( MODHist(copy.deepcopy(rg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 2.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e1_20'].append( MODHist(copy.deepcopy(e1_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0, 0.8), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e2_20'].append( MODHist(copy.deepcopy(e2_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), y_range=(0, 0.45), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e05_20'].append( MODHist(copy.deepcopy(e05_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g ta _g})}$", x_range=(0.01, 1.0), y_range=(0, 1.5), additional_text=additional_text ) ) 




	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['zg_10'].append( MODHist(copy.deepcopy(zg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.01, 1.0), y_range=(0, 2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['rg_10'].append( MODHist(copy.deepcopy(rg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e1_10'].append( MODHist(copy.deepcopy(e1_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0.0, 0.75), x_range=(0.001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e2_10'].append( MODHist(copy.deepcopy(e2_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0.0, 0.45), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e05_10'].append( MODHist(copy.deepcopy(e05_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0.0, 1.4), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_zg_10'].append( MODHist(copy.deepcopy(zg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.01, 1.0), y_range=(0, 2), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_rg_10'].append( MODHist(copy.deepcopy(rg_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e1_10'].append( MODHist(copy.deepcopy(e1_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0, 0.75), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e2_10'].append( MODHist(copy.deepcopy(e2_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), y_range=(0, 0.45), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e05_10'].append( MODHist(copy.deepcopy(e05_10_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g a _g})}$", x_range=(0.01, 1.0), y_range=(0, 1.4), additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 1.7), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 1.3), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0., 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.6), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0., 1.3), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 3), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 2.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g a _g})}$", x_range=(0.01, 1.0), additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['zg_20'].append( MODHist(copy.deepcopy(zg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{z_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.2, 1.0), y_range=(0, 2.5), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['rg_20'].append( MODHist(copy.deepcopy(rg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 0.8), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e1_20'].append( MODHist(copy.deepcopy(e1_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.8), x_range=(0.001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e2_20'].append( MODHist(copy.deepcopy(e2_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", y_range=(0., 0.45), x_range=(0.0001, 1.0), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e05_20'].append( MODHist(copy.deepcopy(e05_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{e_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", y_range=(0., 1.4), x_range=(0.01, 1.0), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_zg_20'].append( MODHist(copy.deepcopy(zg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0, 3), additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_rg_20'].append( MODHist(copy.deepcopy(rg_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 2.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e1_20'].append( MODHist(copy.deepcopy(e1_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.001, 1.0), y_range=(0, 0.8), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e2_20'].append( MODHist(copy.deepcopy(e2_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0001, 1.0), y_range=(0, 0.45), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e05_20'].append( MODHist(copy.deepcopy(e05_20_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], x_scale='log', use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{e_g ta _g})}$", x_range=(0.01, 1.0), y_range=(0, 1.5), additional_text=additional_text ) ) 




	return all_hists



def multi_page_plot_hist_templates():
	
	all_hists = {}

	hardest_pT_hist = Hist(200, 5, 2005, title="pT")
	hardest_phi_hist = Hist(25, 0.0, 2*np.pi, title="phi")
	hardest_eta_hist = Hist(50, -5, 5, title="eta")
	constituent_mul_hist, track_constituent_mul_hist = Hist(60, -0.5, 59.5, title="mul"), Hist(40, -0.5, 39.5, title="mul")
	pT_D_hist = Hist(50, 0, 1, title="pT_D")
	mass_hist, track_mass_hist = Hist(40, 0, 80, title="mass"), Hist(40, 0, 40, title="mass")

	lha_hist = Hist(50, 0, 1)
	width_hist = Hist(60, 0, 0.6)
	thrust_hist = Hist(50, 0, 0.5)

	jec_hist = Hist(50, 1.0, 1.15)
	hardest_area_hist = Hist(50, -0.5, 0.5)

	zg_hist = Hist(50, 0.0, 0.5)
	rg_hist = Hist(60, 0.0, 1.2)
	e1_hist = Hist(60, 0.0, 0.6)
	e05_hist = Hist(60, 0.0, 0.6)
	e2_hist = Hist(60, 0.0, 0.6)


	zg_05_hist = Hist(20, 0.0, 0.5)
	rg_05_hist = Hist(30, 0.0, 1.2)
	e1_05_hist = Hist(30, 0.0, 0.6)
	e05_05_hist = Hist(30, 0.0, 0.6)
	e2_05_hist = Hist(30, 0.0, 0.6)

	pT_boundaries = [85, 115, 150, 200, 250]

	all_hists['hardest_pT'] = []
	all_hists['uncor_hardest_pT'] = []
	all_hists['frac_pT_loss'] = []
	all_hists['jec'] = []
	all_hists['hardest_area'] = []
	all_hists['softkill_pT_loss'] = []
	all_hists['pT_after_SD'] = []
	all_hists['hardest_eta'] = []
	all_hists['hardest_phi'] = []
	all_hists['mul_pre_SD'], all_hists['track_mul_pre_SD'] = [], []
	all_hists['mul_post_SD'], all_hists['track_mul_post_SD'] = [], []
	all_hists['pT_D_pre_SD'], all_hists['track_pT_D_pre_SD'] = [], []
	all_hists['pT_D_post_SD'], all_hists['track_pT_D_post_SD'] = [], []
	all_hists['mass_pre_SD'], all_hists['track_mass_pre_SD'] = [], []
	all_hists['mass_post_SD'], all_hists['track_mass_post_SD'] = [], []

	all_hists['LHA_pre_SD'], all_hists['track_LHA_pre_SD'] = [], []
	all_hists['LHA_post_SD'], all_hists['track_LHA_post_SD'] = [], []
	all_hists['width_pre_SD'], all_hists['track_width_pre_SD'] = [], []
	all_hists['width_post_SD'], all_hists['track_width_post_SD'] = [], []
	all_hists['thrust_pre_SD'], all_hists['track_thrust_pre_SD'] = [], []
	all_hists['thrust_post_SD'], all_hists['track_thrust_post_SD'] = [], []


	all_hists['zg_10'], all_hists['track_zg_10'] = [], []
	all_hists['rg_10'], all_hists['track_rg_10'] = [], []
	all_hists['e1_10'], all_hists['track_e1_10'] = [], []
	all_hists['e2_10'], all_hists['track_e2_10'] = [], []
	all_hists['e05_10'], all_hists['track_e05_10'] = [], []

	all_hists['zg_05'], all_hists['track_zg_05'] = [], []
	all_hists['rg_05'], all_hists['track_rg_05'] = [], []
	all_hists['e1_05'], all_hists['track_e1_05'] = [], []
	all_hists['e2_05'], all_hists['track_e2_05'] = [], []
	all_hists['e05_05'], all_hists['track_e05_05'] = [], []

	all_hists['zg_20'], all_hists['track_zg_20'] = [], []
	all_hists['rg_20'], all_hists['track_rg_20'] = [], []
	all_hists['e1_20'], all_hists['track_e1_20'] = [], []
	all_hists['e2_20'], all_hists['track_e2_20'] = [], []
	all_hists['e05_20'], all_hists['track_e05_20'] = [], []
	
	

	for i in range(len(pT_boundaries) - 1):
		
		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['hardest_pT'].append( MODHist(Hist(25, pT_boundaries[i], pT_boundaries[i + 1]), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$p_T$ [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T}$", y_scale='log', x_range=(pT_boundaries[i], pT_boundaries[i + 1]), y_range=(5e-3, 5e-1), additional_text=additional_text ) ) 



		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['uncor_hardest_pT'].append( MODHist(Hist(25, pT_boundaries[i], pT_boundaries[i + 1]), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1])), ('hardest_eta', (-2.4, 2.4))], use_prescale=False, x_label="$p_T$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T}$", x_scale='linear', y_scale='log', x_range=(pT_boundaries[i], pT_boundaries[i + 1]), y_range=(5e-3, 5e-1), additional_text=additional_text ) ) 


		additional_text = [ ( (0.50, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['frac_pT_loss'].append( MODHist(Hist(25, 0., 0.2), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Fractional $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 



		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['pT_after_SD'].append( MODHist(Hist(25, pT_boundaries[i], pT_boundaries[i + 1]), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$p_T$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T}$", y_scale='log', x_range=(pT_boundaries[i], pT_boundaries[i + 1]), y_range=(1e-3, 1e1), additional_text=additional_text ) ) 




		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['hardest_eta'].append( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$\eta$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} \eta}$", x_range=(-5., 5.), y_range=[(0., 0.35), (0., 0.35), (0., 0.38), (0., 0.38)][i], mark_regions=[(-2.4, 0.24, 'right', 0.20), (2.4, 0.24, 'left', 0.20)], additional_text=additional_text ) ) 




		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['hardest_phi'].append( MODHist(copy.deepcopy(hardest_phi_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$\phi$", x_range=(0, 2 * np.pi), y_range=[0.0, 0.40], y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} \phi}$", additional_text=additional_text, axes_label_pi=True ) ) 


		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['jec'].append( MODHist(jec_hist, conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1])), ('hardest_eta', (-2.4, 2.4))], use_prescale=False, x_label="JEC", y_label="A.U.", x_range=(1.0, 1.15), x_scale='linear', additional_text=additional_text ) ) 


		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['hardest_area'].append( MODHist(hardest_area_hist, conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1])), ('hardest_eta', (-2.4, 2.4))], use_prescale=False, x_label="Jet Area", y_label="A.U.", x_range=(-0.5, 0.5), x_scale='linear', additional_text=additional_text ) ) 



		additional_text = [ ( (0.50, 0.68), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['softkill_pT_loss'].append( MODHist(Hist(25, 0., 0.2), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Softkill $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 




		additional_text = [ ( (1.0, 0.68), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['mul_pre_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Constituent Multiplicity", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", x_range=(0, 45), y_range=(0, 0.11), additional_text=additional_text, legend_location=[(0.96, 1.0), (0.94, 1.0), (0.94, 1.0), (0.94, 1.0)][i] ) ) 

		additional_text = [ ( (1.0, 0.68), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_mul_pre_SD'].append( MODHist(copy.deepcopy(track_constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track Constituent Multiplicity", x_range=(0, 45), y_range=(0, 0.11), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text, legend_location=[(0.96, 1.0), (0.94, 1.0), (0.94, 1.0), (0.94, 1.0)][i] ) ) 
		


		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['mul_post_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Constituent Multiplicity", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_mul_post_SD'].append( MODHist(copy.deepcopy(track_constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track Constituent Multiplicity", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text ) ) 
		




		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", y_range=(0, 10), additional_text=additional_text, legend_location=[(0.96, 1.0), (0.94, 1.0), (0.94, 1.0), (0.94, 1.0)][i] ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", y_range=(0, 10), additional_text=additional_text, legend_location=[(0.96, 1.0), (0.94, 1.0), (0.94, 1.0), (0.94, 1.0)][i] ) ) 



		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 





		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['mass_pre_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", x_range=[(0.0, 40.0), (0.0, 40.0), (0.0, 50.0), (0.0, 70.0)][i], y_range=[(0, 0.12), (0, 0.11), (0, 0.08), (0, 0.10)][i], additional_text=additional_text, legend_location=[(0.96, 1.0), (0.94, 1.0), (0.94, 1.0), (0.94, 1.0)][i] ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_mass_pre_SD'].append( MODHist(copy.deepcopy(track_mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", x_range=[(0.0, 30.0), (0.0, 30.0), (0.0, 40.0), (0.0, 40.0)][i], y_range=[(0, 0.12), (0, 0.11), (0, 0.08), (0, 0.10)][i], additional_text=additional_text, legend_location=[(0.96, 1.0), (0.94, 1.0), (0.94, 1.0), (0.94, 1.0)][i] ) ) 




		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['mass_post_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_mass_post_SD'].append( MODHist(copy.deepcopy(track_mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 



		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", y_range=[(0., 4.5), (0., 4.5), (0., 4.5), (0., 4.5)][i], additional_text=additional_text, legend_location=[(0.96, 1.0), (0.94, 1.0), (0.94, 1.0), (0.94, 1.0)][i] ) )  

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", y_range=[(0., 4.5), (0., 4.5), (0., 4.5), (0., 4.5)][i], additional_text=additional_text, legend_location=[(0.96, 1.0), (0.94, 1.0), (0.94, 1.0), (0.94, 1.0)][i] ) ) 



		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 




		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 


		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 




		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} 2}$", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} 2}$", additional_text=additional_text ) ) 


		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} 2}$", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		all_hists['track_thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} 2}$", additional_text=additional_text ) ) 





		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['zg_10'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), y_range=(0, 9), additional_text=additional_text ) ) 
	
		additional_text = [ ( (0.50, 0.98), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['rg_10'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta}$", x_range=(0.0, 1.0), y_range=(0., 5.0), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['e1_10'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['e2_10'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['e05_10'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0, 0.6), additional_text=additional_text ) ) 



		additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_zg_10'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), y_range=(0, 9), additional_text=additional_text ) ) 
	
		additional_text = [ ( (0.50, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_rg_10'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 5.), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_e1_10'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_e2_10'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists['track_e05_10'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 







		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 
	
		additional_text = [ ( (0.50, 0.98), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$\\theta_g$", y_label="$\displaystyle \\fra\\theta{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.5), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 





		additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_zg_05'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 
	
		additional_text = [ ( (0.50, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_rg_05'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_e1_05'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_e2_05'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
		all_hists['track_e05_05'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 







		additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['zg_20'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0., 11.), x_range=(0.0, 0.6), additional_text=additional_text ) ) 
	
		additional_text = [ ( (1.0, 0.60), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['rg_20'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$\\theta_g$", y_label="$\displaystyle \\fra\\theta{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 6.), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['e1_20'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['e2_20'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['e05_20'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0, 0.6), additional_text=additional_text ) ) 



		additional_text = [ ( (-0.08, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_zg_20'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 
	
		additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_rg_20'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_e1_20'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_e2_20'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
		all_hists['track_e05_20'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 




	additional_text = [ ( (0.49, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['hardest_pT'].append( MODHist(Hist(200, 5, 2005, title="pT"), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$p_T$ [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T}$", y_scale='log', x_range=(0, 1000), y_range=(1e-9, 1e-1), additional_text=additional_text ) )

	additional_text = [ ( (0.49, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['hardest_pT'].append( MODHist(Hist(200, 0, 2000, title="pT"), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$p_T$ [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T}$", y_scale='log', x_range=(0, 1000), y_range=(1e-9, 1e-1), additional_text=additional_text ) )

	additional_text = [ ( (0.49, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['hardest_pT'].append( MODHist(Hist(200, 0, 2000, title="pT"), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$p_T$ [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T}$", y_scale='log', x_range=(0, 1000), y_range=(5e-7, 1e-1), additional_text=additional_text ) )




	additional_text = [ ( (0.47, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['uncor_hardest_pT'].append( MODHist(Hist(200, 5, 2005, title="pT"), conditions=[('hardest_pT', (85, None)), ('hardest_eta', (-2.4, 2.4))], use_prescale=True, x_label="$p_T$ [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T}$", x_scale='linear', y_scale='log', x_range=(0, 1000), y_range=(1e-9, 1e-1), additional_text=additional_text ) )

	additional_text = [ ( (0.47, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['uncor_hardest_pT'].append( MODHist(Hist(200, 0, 2000, title="pT"), conditions=[('hardest_pT', (150, None)), ('hardest_eta', (-2.4, 2.4))], use_prescale=True, x_label="$p_T$ [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T}$", x_scale='linear', y_scale='log', x_range=(0, 1000), y_range=(1e-9, 1e-1), additional_text=additional_text ) )

	additional_text = [ ( (0.47, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['uncor_hardest_pT'].append( MODHist(Hist(200, 0, 2000, title="pT"), conditions=[('hardest_pT', (250, None)), ('hardest_eta', (-2.4, 2.4))], use_prescale=True, x_label="$p_T$ [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T}$", x_scale='linear', y_scale='log', x_range=(0, 1000), y_range=(5e-7, 1e-1), additional_text=additional_text ) )




	additional_text = [ ( (0.49, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['frac_pT_loss'].append( MODHist(Hist(25, 0., 1.), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Fractional $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (0.49, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['frac_pT_loss'].append( MODHist(Hist(25, 0., 1.), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Fractional $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 


	additional_text = [ ( (0.49, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['frac_pT_loss'].append( MODHist(Hist(25, 0., 1.), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Fractional $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 




	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['jec'].append( MODHist(jec_hist, conditions=[('hardest_pT', (85, None)), ('hardest_eta', (-2.4, 2.4))], use_prescale=True, x_label="JEC", y_label="A.U.", x_range=(1.0, 1.15), x_scale='linear', additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['jec'].append( MODHist(jec_hist, conditions=[('hardest_pT', (150, None)), ('hardest_eta', (-2.4, 2.4))], use_prescale=True, x_label="JEC", y_label="A.U.", x_range=(1.0, 1.15), x_scale='linear', additional_text=additional_text ) ) 


	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['jec'].append( MODHist(jec_hist, conditions=[('hardest_pT', (250, None)), ('hardest_eta', (-2.4, 2.4))], use_prescale=True, x_label="JEC", y_label="A.U.", x_range=(1.0, 1.15), x_scale='linear', additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['hardest_area'].append( MODHist(hardest_area_hist, conditions=[('hardest_pT', (85, None)), ('hardest_eta', (-2.4, 2.4))], use_prescale=True, x_label="Jet Area", y_label="A.U.", x_range=(-0.5, 0.5), x_scale='linear', additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['hardest_area'].append( MODHist(hardest_area_hist, conditions=[('hardest_pT', (150, None)), ('hardest_eta', (-2.4, 2.4))], use_prescale=True, x_label="Jet Area", y_label="A.U.", x_range=(-0.5, 0.5), x_scale='linear', additional_text=additional_text ) ) 


	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['hardest_area'].append( MODHist(hardest_area_hist, conditions=[('hardest_pT', (250, None)), ('hardest_eta', (-2.4, 2.4))], use_prescale=True, x_label="Jet Area", y_label="A.U.", x_range=(-0.5, 0.5), x_scale='linear', additional_text=additional_text ) ) 




	additional_text = [ ( (0.49, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['pT_after_SD'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$p_T$", y_label="A.U.", y_scale='log', x_range=(0, 1000), y_range=(1e-7, 1e0), additional_text=additional_text ) )

	additional_text = [ ( (0.49, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['pT_after_SD'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$p_T$", y_label="A.U.", y_scale='log', x_range=(0, 1000), y_range=(1e-7, 1e0), additional_text=additional_text ) )

	additional_text = [ ( (0.49, 0.68), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['pT_after_SD'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$p_T$", y_label="A.U.", y_scale='log', x_range=(0, 1000), y_range=(1e-7, 1e0), additional_text=additional_text ) )




	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['hardest_eta'].append( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_pT', (85, None))], use_prescale=True, x_label="$\eta$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} \eta}$", x_range=(-5., 5.), y_range=(0., 0.35), mark_regions=[(-2.4, 0.24, 'right', 0.20), (2.4, 0.24, 'left', 0.20)], additional_text=additional_text ) ) 

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['hardest_eta'].append( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_pT', (150, None))], use_prescale=True, x_label="$\eta$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} \eta}$", x_range=(-5., 5.), y_range=(0., 0.40), mark_regions=[(-2.4, 0.24, 'right', 0.20), (2.4, 0.24, 'left', 0.20)], additional_text=additional_text ) ) 


	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['hardest_eta'].append( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_pT', (250, None))], use_prescale=True, x_label="$\eta$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} \eta}$", x_range=(-5., 5.), y_range=(0., 0.50), mark_regions=[(-2.4, 0.28, 'right', 0.23), (2.4, 0.28, 'left', 0.23)], additional_text=additional_text ) ) 






	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['hardest_phi'].append( MODHist(copy.deepcopy(hardest_phi_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$\phi$", x_range=(0, 2 * np.pi), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} \phi}$", y_range=(0.0, 0.5), additional_text=additional_text, axes_label_pi=True ) )

	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['hardest_phi'].append( MODHist(copy.deepcopy(hardest_phi_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$\phi$", x_range=(0, 2 * np.pi), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} \phi}$", y_range=(0.0, 0.3), additional_text=additional_text, axes_label_pi=True ) )


	additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['hardest_phi'].append( MODHist(copy.deepcopy(hardest_phi_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$\phi$", x_range=(0, 2 * np.pi), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} \phi}$", y_range=(0.0, 0.5), additional_text=additional_text, axes_label_pi=True ) )





	additional_text = [ ( (0.50, 0.68), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['softkill_pT_loss'].append( MODHist(Hist(25, 0., 0.2), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Softkill $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.68), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['softkill_pT_loss'].append( MODHist(Hist(25, 0., 0.2), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Softkill $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.68), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['softkill_pT_loss'].append( MODHist(Hist(25, 0., 0.2), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Softkill $p_T$ Loss", y_label="A.U.", additional_text=additional_text ) ) 




	


	additional_text = [ ( (0.995, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['mul_pre_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Constituent Multiplicity", x_range=(0, 50), y_range=(0., 0.10), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text, legend_location=(1.0, 1.0) ) ) 

	additional_text = [ ( (0.995, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['mul_pre_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Constituent Multiplicity", y_range=(0., 0.08), x_range=(0, 50), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text, legend_location=(1.0, 1.0) ) ) 

	additional_text = [ ( (0.995, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['mul_pre_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Constituent Multiplicity", x_range=(0, 50), y_range=(0., 0.10), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text, legend_location=(1.0, 1.0) ) ) 




	additional_text = [ ( (0.995, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_mul_pre_SD'].append( MODHist(copy.deepcopy(track_constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, y_range=(0, 0.10), x_label="Track Constituent Multiplicity", x_range=(0, 50), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text, legend_location=(1.0, 1.0) ) ) 

	additional_text = [ ( (0.995, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_mul_pre_SD'].append( MODHist(copy.deepcopy(track_constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, y_range=(0, 0.08), x_label="Track Constituent Multiplicity", x_range=(0, 50), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text, legend_location=(1.0, 1.0) ) ) 

	additional_text = [ ( (0.995, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_mul_pre_SD'].append( MODHist(copy.deepcopy(track_constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, y_range=(0, 0.10), x_label="Track Constituent Multiplicity", x_range=(0, 50), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text, legend_location=(1.0, 1.0) ) ) 






	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['mul_post_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Constituent Multiplicity", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['mul_post_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Constituent Multiplicity", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['mul_post_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Constituent Multiplicity", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_mul_post_SD'].append( MODHist(copy.deepcopy(track_constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track Constituent Multiplicity", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_mul_post_SD'].append( MODHist(copy.deepcopy(track_constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track Constituent Multiplicity", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_mul_post_SD'].append( MODHist(copy.deepcopy(track_constituent_mul_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track Constituent Multiplicity", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} N}$", additional_text=additional_text ) ) 






	additional_text = [ ( (0.99, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", y_range=(0, 10), additional_text=additional_text ) ) 

	additional_text = [ ( (0.99, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", y_range=(0, 10), additional_text=additional_text ) ) 


	additional_text = [ ( (0.99, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", y_range=(0, 10), additional_text=additional_text ) ) 






	additional_text = [ ( (0.99, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 

	additional_text = [ ( (0.99, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 

	additional_text = [ ( (0.99, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 


	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_pT_D_post_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $p_T^D$", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} p_T^D}$", additional_text=additional_text ) ) 






	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['mass_pre_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", x_range=(0.0, 40.0), y_range=(0, 0.11), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['mass_pre_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", x_range=(0.0, 60.0), y_range=(0, 0.08), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['mass_pre_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", x_range=(0.0, 80.0), y_range=(0, 0.09), additional_text=additional_text ) ) 




	
	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_mass_pre_SD'].append( MODHist(copy.deepcopy(track_mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track Mass [GeV]", x_range=(0.0, 40.0), y_range=(0, 0.11), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_mass_pre_SD'].append( MODHist(copy.deepcopy(track_mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track Mass [GeV]", x_range=(0.0, 35.0), y_range=(0, 0.08), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 


	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_mass_pre_SD'].append( MODHist(copy.deepcopy(track_mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track Mass [GeV]", x_range=(0.0, 40.0), y_range=(0, 0.09), y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['mass_post_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['mass_post_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['mass_post_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 




	
	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_mass_post_SD'].append( MODHist(copy.deepcopy(track_mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_mass_post_SD'].append( MODHist(copy.deepcopy(track_mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_mass_post_SD'].append( MODHist(copy.deepcopy(track_mass_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track Mass [GeV]", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} m}$", additional_text=additional_text ) ) 







	additional_text = [ ( (0.99, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", y_range=(0., 3.5), additional_text=additional_text ) ) 

	additional_text = [ ( (0.99, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", y_range=(0., 4.0), additional_text=additional_text ) ) 

	additional_text = [ ( (0.99, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", y_range=(0., 6.0), additional_text=additional_text ) ) 






	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", y_range=(0., 3.5), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", y_range=(0., 4.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", y_range=(0., 6.0), additional_text=additional_text ) ) 







	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 


	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_LHA_post_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track LHA", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1/2)}}$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_width_post_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track Width", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(1)}}$", additional_text=additional_text ) ) 







	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 







	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 




	additional_text = [ ( (1.0, 0.70), 'upper.right', "$\mathrm{PFC}~p_T > 0.5~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 








	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	all_hists['track_thrust_post_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track Thrust", y_label="$\\frac{1}{\sigma}\\frac{\mathrm{d} \sigma}{\mathrm{d} e^{(2)}}$", additional_text=additional_text ) ) 




	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['zg_10'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), y_range=(0, 9), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['rg_10'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta}$", x_range=(0.0, 1.0), y_range=(0., 5.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e1_10'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e2_10'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e05_10'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0, 0.6), additional_text=additional_text ) ) 



	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_zg_10'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), y_range=(0, 9), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_rg_10'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 5.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e1_10'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e2_10'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e05_10'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 







	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\fra\\theta{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.5), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_zg_05'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_rg_05'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e1_05'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e2_05'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e05_05'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 







	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['zg_20'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0., 11.), x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['rg_20'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\fra\\theta{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 6.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e1_20'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e2_20'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e05_20'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0, 0.6), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_zg_20'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_rg_20'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e1_20'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e2_20'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e05_20'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 
	





	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['zg_10'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), y_range=(0, 9), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['rg_10'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta}$", x_range=(0.0, 1.0), y_range=(0., 5.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e1_10'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e2_10'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e05_10'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0, 0.6), additional_text=additional_text ) ) 



	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_zg_10'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), y_range=(0, 9), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_rg_10'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 5.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e1_10'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e2_10'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e05_10'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 







	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\fra\\theta{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.5), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_zg_05'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_rg_05'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e1_05'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e2_05'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e05_05'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 







	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['zg_20'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0., 11.), x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['rg_20'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\fra\\theta{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 6.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e1_20'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e2_20'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e05_20'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0, 0.6), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_zg_20'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_rg_20'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e1_20'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e2_20'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e05_20'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 




	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['zg_10'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), y_range=(0, 9), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['rg_10'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta}$", x_range=(0.0, 1.0), y_range=(0., 5.0), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e1_10'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e2_10'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['e05_10'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0, 0.6), additional_text=additional_text ) ) 



	additional_text = [ ( (1.0, 0.60), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_zg_10'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), y_range=(0, 9), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_rg_10'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 5.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e1_10'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e2_10'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists['track_e05_10'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 







	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['zg_05'].append( MODHist(copy.deepcopy(zg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['rg_05'].append( MODHist(copy.deepcopy(rg_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\fra\\theta{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.5), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e1_05'].append( MODHist(copy.deepcopy(e1_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e2_05'].append( MODHist(copy.deepcopy(e2_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['e05_05'].append( MODHist(copy.deepcopy(e05_05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 





	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_zg_05'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_rg_05'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e1_05'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e2_05'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.05$" ) ]
	all_hists['track_e05_05'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 







	additional_text = [ ( (-0.08, 0.98), 'upper left', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['zg_20'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", y_range=(0., 11.), x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.60), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['rg_20'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$\\theta_g$", y_label="$\displaystyle \\fra\\theta{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 6.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e1_20'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e2_20'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['e05_20'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="$e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0, 0.6), additional_text=additional_text ) ) 



	additional_text = [ ( (-0.08, 0.98), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_zg_20'].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_rg_20'].append( MODHist(copy.deepcopy(rg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $\\theta_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.0, 1.0), y_range=(0., 4.), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e1_20'].append( MODHist(copy.deepcopy(e1_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $e_g^1$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e2_20'].append( MODHist(copy.deepcopy(e2_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $e_g^2$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "\\textbf{Track Only} \n $\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.2$" ) ]
	all_hists['track_e05_20'].append( MODHist(copy.deepcopy(e05_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], use_prescale=True, x_label="Track $e_g^{0.5}$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} e_g }$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 



	return all_hists










def trigger_hists():
	all_hists = {}

		
	pT_hist = Hist(50, 0, 300)

	# all_hists['corr_hardest_pT'] = []


	trigger_names = ["Jet15U_HcalNoiseFiltered", "Jet30U", "Jet50U", "Jet70U", "Jet100U", "Jet140U"]


	for trigger_name in trigger_names:
		additional_text = [ ( (-0.08, 0.98), 'upper left', "$ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
		all_hists[trigger_name] =  MODHist(copy.deepcopy(Hist(50, 0, 300)), conditions=[(['trigger_fired', 1], lambda x, y: y == x), (['jet_quality', 3], lambda x, y: y >= x), (['trig_jet_matched', 1], lambda x, y: y == x)], use_prescale=True, x_scale='log', x_label="Fractional $p_T$ Loss", y_label="A.U.", y_range=(0., 1.2), additional_text=additional_text )



	return all_hists


def get_pfc_hists():
	pfc_hists = {}

	pfc_hists['pfc_pT'] = []

	pT_boundaries = [85, 115, 150, 200, 250]



	# for i in range(len(pT_boundaries) - 1):
		
	# 	additional_text = [ ( (0.50, 0.70), 'upper left', "Charged PFCs \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
	# 	pfc_hists['pfc_pT'].append( MODHist(Hist(25, 0, 5), conditions=[('jet_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], mark_regions=[(1.0, 0.5, 'right', 0.45)], use_prescale=False, x_label="Charged PFC $p_T$ [GeV]", y_label="A.U.", x_range=(0, 5), additional_text=additional_text ) ) 


	# additional_text = [ ( (0.50, 0.70), 'upper left', "Charged PFCs \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	# pfc_hists['pfc_pT'].append( MODHist(Hist(25, 0, 5), conditions=[('jet_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], mark_regions=[(1.0, 0.5, 'right', 0.45)], use_prescale=True, x_label="Charged PFC $p_T$ [GeV]", y_label="A.U.", x_range=(0, 5), additional_text=additional_text ) ) 

	# additional_text = [ ( (0.50, 0.70), 'upper left', "Charged PFCs \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	# pfc_hists['pfc_pT'].append( MODHist(Hist(25, 0, 5), conditions=[('jet_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], mark_regions=[(1.0, 0.5, 'right', 0.45)], use_prescale=True, x_label="Charged PFC $p_T$ [GeV]", y_label="A.U.", x_range=(0, 5), additional_text=additional_text ) ) 

	# additional_text = [ ( (0.50, 0.70), 'upper left', "Charged PFCs \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	# pfc_hists['pfc_pT'].append( MODHist(Hist(25, 0, 5), conditions=[('jet_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], mark_regions=[(1.0, 0.5, 'right', 0.45)], use_prescale=True, x_label="Charged PFC $p_T$ [GeV]", y_label="A.U.", x_range=(0, 5), additional_text=additional_text ) ) 



	for i in range(len(pT_boundaries) - 1):
		
		additional_text = [ ( (0.50, 0.70), 'upper left', "Neutral PFCs \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$" ) ]
		pfc_hists['pfc_pT'].append( MODHist(Hist(25, 0, 5), conditions=[('jet_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], mark_regions=[(1.0, 0.8, 'right', 0.6)], use_prescale=False, x_label="Neutral PFC $p_T$ [GeV]", y_label="A.U.", x_range=(0, 5), additional_text=additional_text ) ) 


	additional_text = [ ( (0.50, 0.70), 'upper left', "Neutral PFCs \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 85~\mathrm{GeV}$" ) ]
	pfc_hists['pfc_pT'].append( MODHist(Hist(25, 0, 5), conditions=[('jet_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], mark_regions=[(1.0, 0.8, 'right', 0.6)], use_prescale=True, x_label="Neutral PFC $p_T$ [GeV]", y_label="A.U.", x_range=(0, 5), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.70), 'upper left', "Neutral PFCs \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 150~\mathrm{GeV}$" ) ]
	pfc_hists['pfc_pT'].append( MODHist(Hist(25, 0, 5), conditions=[('jet_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], mark_regions=[(1.0, 0.8, 'right', 0.6)], use_prescale=True, x_label="Neutral PFC $p_T$ [GeV]", y_label="A.U.", x_range=(0, 5), additional_text=additional_text ) ) 

	additional_text = [ ( (0.50, 0.70), 'upper left', "Neutral PFCs \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} > 250~\mathrm{GeV}$" ) ]
	pfc_hists['pfc_pT'].append( MODHist(Hist(25, 0, 5), conditions=[('jet_eta', (-2.4, 2.4)), ('hardest_pT', (250, None))], mark_regions=[(1.0, 0.8, 'right', 0.6)], use_prescale=True, x_label="Neutral PFC $p_T$ [GeV]", y_label="A.U.", x_range=(0, 5), additional_text=additional_text ) ) 



	return pfc_hists






def two_dim_hists():
	all_hists = {}

		
	zg_hist = Hist2D(50, 0.0, 0.5, 50, 0, 1.0)
	
	pT_boundaries = [85, 115, 150, 200, 250]

	all_hists[('zg_10', 'rg_10')] = []

	for i in range(len(pT_boundaries) - 1):

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists[('zg_10', 'rg_10')].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 




	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists[('zg_10', 'rg_10')].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists[('zg_10', 'rg_10')].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 



	return all_hists





def two_dim_log_hists():
	all_hists = {}

			
	# zg_hist = Hist2D(50, 0.0, 0.5, 50, 0, 1.0)
	# zg_10_hist = Hist(np.logspace(math.log(float(0.01), math.e), math.log(0.5, math.e), 50, base=np.e))
	zg_hist = Hist2D(np.logspace(math.log(float(0.1), math.e), math.log(0.5, math.e), 50, base=np.e), np.logspace(math.log(float(0.01), math.e), math.log(1.0, math.e), 50, base=np.e))

	pT_boundaries = [85, 115, 150, 200, 250]

	all_hists[('zg_10', 'rg_10')] = []

	for i in range(len(pT_boundaries) - 1):

		additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
		all_hists[('zg_10', 'rg_10')].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 


	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists[('zg_10', 'rg_10')].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (85, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.62), 'upper right', "$\mathrm{PFC}~p_T > 1.0~\mathrm{GeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5; \left| \eta \\right| < 2.4$ \n Hardest Jet $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n Soft Drop: $\\beta = 0; z_{\mathrm{cut}} = 0.1$" ) ]
	all_hists[('zg_10', 'rg_10')].append( MODHist(copy.deepcopy(zg_hist), conditions=[('hardest_eta', (-2.4, 2.4)), ('hardest_pT', (150, None))], use_prescale=True, x_label="$z_g$", y_label="$\displaystyle \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text ) ) 



	return all_hists






def get_hist_template(var):
	all_hists = all_hist_templates()
	return all_hists[var]


if __name__ == "__main__":
	

	hist = Hist2D(1, 0.1, 0.2, 1, 0, 0.2)

	hist.fill_array([[0.1, 0.3]])