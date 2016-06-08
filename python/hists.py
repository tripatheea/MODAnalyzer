from __future__ import division
import numpy as np
import copy


# RootPy
from rootpy.plotting import Hist


class MODHist:

	def __init__(self, hist, conditions=[], use_prescale=True, x_label="", y_label="", y_scale='linear', x_range=(0, -1), y_range=(0, -1), additional_text=[]):
		# hist should be a RootPy Hist object (TH1). 
		# conditions should be a list of lambda functions that need to be satisfied by 'value' for this hist to be filled. 

		self._hist = hist
		self._conditions = conditions
		self._use_prescale = use_prescale

		self._x_label = x_label
		self._y_label = y_label

		self._y_scale = y_scale

		self._x_range = x_range
		self._y_range = y_range

		self._additional_text = additional_text


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

	def y_scale(self):
		return self._y_scale

	def x_range(self):
		return self._x_range

	def y_range(self):
		return self._y_range

	def additional_text(self):
		return self._additional_text



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

	'''
	all_hists['hardest_pT'] = []
	all_hists['hardest_eta'] = []

	for i in range(len(pT_boundaries) - 1):
		all_hists['hardest_pT'].append( copy.deepcopy( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_pT', lambda x: x > pT_boundaries[i] and x < pT_boundaries[i + 1])], use_prescale=False, x_lims=(pT_boundaries[i], pT_boundaries[i + 1])) ) )
		all_hists['hardest_eta'].append( copy.deepcopy( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_pT', lambda x: x > pT_boundaries[i] and x < pT_boundaries[i + 1])], use_prescale=False) ) )


	# all_hists['hardest_pT'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_pT', lambda x: x > 85)], use_prescale=True) )
	all_hists['hardest_eta'].append( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_eta', lambda x: x > 85)], use_prescale=True) )
	
	# all_hists['hardest_pT'].append( MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_pT', lambda x: x > 150)], use_prescale=True) )
	all_hists['hardest_eta'].append( MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_eta', lambda x: x > 150)], use_prescale=True) )
	'''

	a = MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_pT', lambda x: x > 150)], use_prescale=True)
	b = MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_pT', lambda x: x > 150)], use_prescale=True)
	
	all_hists['hardest_pT'] = [a, b]

	return all_hists



def multi_page_plot_hist_templates():
	all_hists = {}

	hardest_pT_hist = Hist(100, 5, 1005, title="pT")
	hardest_phi_hist = Hist(25, 0, 2*np.pi, title="phi")
	hardest_eta_hist = Hist(25, -5, 5, title="eta")
	constituent_mul_hist = Hist(30, 0, 60, title="mul")
	pT_D_hist = Hist(25, 0, 1, title="pT_D")
	mass_hist = Hist(40, 0, 80, title="mass")

	lha_hist = Hist(25, 0, 1)
	width_hist = Hist(25, 0, 1)
	thrust_hist = Hist(25, 0, 1)


	pT_boundaries = [85, 115, 150, 200, 250, 325]

	all_hists['hardest_pT'] = []
	all_hists['hardest_phi'] = []
	all_hists['mul_pre_SD'] = []
	all_hists['pT_D_pre_SD'] = []
	all_hists['mass_pre_SD'] = []

	all_hists['LHA_pre_SD'] = []
	all_hists['width_pre_SD'] = []
	all_hists['thrust_pre_SD'] = []

	

	for i in range(len(pT_boundaries) - 1):
		
		additional_text = [ ( (-0.09, 0.80), 'upper left', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV};\eta<2.4$" ) ]
		all_hists['hardest_pT'].append( MODHist(Hist(25, pT_boundaries[i], pT_boundaries[i + 1]), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Jet $p_T$", y_label="A.U.", x_range=(pT_boundaries[i], pT_boundaries[i + 1]), y_range=(1e-3, 1e1), additional_text=additional_text ) ) 

		additional_text = [ ( (-0.09, 0.80), 'upper left', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV};\eta<2.4$" ) ]
		all_hists['hardest_phi'].append( MODHist(copy.deepcopy(hardest_phi_hist), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Jet $\phi$", y_label="A.U.", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV};\eta<2.4$" ) ]
		all_hists['mul_pre_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Constituent Multiplicity", y_label="A.U.", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV};\eta<2.4$" ) ]
		all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="$p_T^D$", y_label="A.U.", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV};\eta<2.4$" ) ]
		all_hists['mass_pre_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Jet Mass", y_label="A.U.", additional_text=additional_text ) ) 



		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV};\eta<2.4$" ) ]
		all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Jet LHA", y_label="A.U.", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV};\eta<2.4$" ) ]
		all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Jet Width", y_label="A.U.", additional_text=additional_text ) ) 

		additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} \in [" + str(pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV};\eta<2.4$" ) ]
		all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))], use_prescale=False, x_label="Jet Thrust", y_label="A.U.", additional_text=additional_text ) ) 




	additional_text = [ ( (0.49, 0.68), 'upper left', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 85~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['hardest_pT'].append( MODHist(Hist(100, 5, 1005, title="pT"), conditions=[('hardest_pT', (85, None))], use_prescale=True, x_label="Jet $p_T$", y_label="A.U.", y_range=(1e-8, 1e0), additional_text=additional_text ) )

	additional_text = [ ( (0.49, 0.68), 'upper left', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 150~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['hardest_pT'].append( MODHist(Hist(100, 0, 1000, title="pT"), conditions=[('hardest_pT', (150, None))], use_prescale=True, x_label="Jet $p_T$", y_label="A.U.", y_range=(1e-8, 1e0), additional_text=additional_text ) )


	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 85~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['mul_pre_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_pT', (85, None))], use_prescale=True, x_label="Constituent Multiplicity", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 150~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['mul_pre_SD'].append( MODHist(copy.deepcopy(constituent_mul_hist), conditions=[('hardest_pT', (150, None))], use_prescale=True, x_label="Constituent Multiplicity", y_label="A.U.", additional_text=additional_text ) ) 


	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 85~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_pT', (85, None))], use_prescale=True, x_label="$p_T^D$", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 150~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['pT_D_pre_SD'].append( MODHist(copy.deepcopy(pT_D_hist), conditions=[('hardest_pT', (150, None))], use_prescale=True, x_label="$p_T^D$", y_label="A.U.", additional_text=additional_text ) ) 



	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 85~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['mass_pre_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_pT', (85, None))], use_prescale=True, x_label="Jet Mass", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 150~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['mass_pre_SD'].append( MODHist(copy.deepcopy(mass_hist), conditions=[('hardest_pT', (150, None))], use_prescale=True, x_label="Jet Mass", y_label="A.U.", additional_text=additional_text ) ) 

	

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 85~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_pT', (85, None))], use_prescale=True, x_label="Jet LHA", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 150~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['LHA_pre_SD'].append( MODHist(copy.deepcopy(lha_hist), conditions=[('hardest_pT', (150, None))], use_prescale=True, x_label="Jet LHA", y_label="A.U.", additional_text=additional_text ) ) 


	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 85~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_pT', (85, None))], use_prescale=True, x_label="Jet Width", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 150~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['width_pre_SD'].append( MODHist(copy.deepcopy(width_hist), conditions=[('hardest_pT', (150, None))], use_prescale=True, x_label="Jet Width", y_label="A.U.", additional_text=additional_text ) ) 


	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 85~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_pT', (85, None))], use_prescale=True, x_label="Jet Thrust", y_label="A.U.", additional_text=additional_text ) ) 

	additional_text = [ ( (1.0, 0.70), 'upper right', "$\mathrm{PFC}~pT > 500~\mathrm{MeV}$ \n $ \mathrm{Anti-}k_{t}\mathrm{:}~R = 0.5$ \n $p_{T} > 150~\mathrm{GeV};\eta<2.4$" ) ]
	all_hists['thrust_pre_SD'].append( MODHist(copy.deepcopy(thrust_hist), conditions=[('hardest_pT', (150, None))], use_prescale=True, x_label="Jet Thrust", y_label="A.U.", additional_text=additional_text ) ) 

	



	return all_hists


def get_hist_template(var):
	all_hists = all_hist_templates()
	return all_hists[var]