from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


import parse








def parse_theory_file():
	
	input_dir = "/home/aashish/Dropbox (MIT)/Research/CMSOpenData/theory distributions/"

	pTs = [85, 115, 150, 200]

	
	input_files = []
	input_files.extend( ["zg/zg_" + str(pT) for pT in pTs] )
	input_files.extend( ["rg/rg_" + str(pT) for pT in pTs] )
	input_files.extend( ["e1/e1_" + str(pT) for pT in pTs] )
	input_files.extend( ["e2/e2_" + str(pT) for pT in pTs] )
	input_files.extend( ["e05/e05_" + str(pT) for pT in pTs] )

	input_files.extend( ["zg/zg_" + str(pT) + "log" for pT in pTs] )
	input_files.extend( ["rg/rg_" + str(pT) + "log" for pT in pTs] )
	input_files.extend( ["e1/e1_" + str(pT) + "log" for pT in pTs] )
	input_files.extend( ["e2/e2_" + str(pT) + "log" for pT in pTs] )
	input_files.extend( ["e05/e05_" + str(pT) + "log" for pT in pTs] )

	hists = {}


	for f in input_files:
	
		var_name = f.split("/")[1]
		hists[var_name] = []

		input_file = input_dir + f + ".dat"

		hists[var_name] = [ [], [], [], [] ]

		with open(input_file) as infile:
			
			line_number = 0


			for line in infile:
				line_components = line.split()
				
				hists[var_name][0].append(float(line_components[0]))
				
				hists[var_name][1].append(float(line_components[2]))
				hists[var_name][2].append(float(line_components[1]))
				hists[var_name][3].append(float(line_components[3]))

	return hists


def compile_sources(parsed_hists):

	data_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[0], parsed_hists[1], parsed_hists[2], parsed_hists[3]
	
	theory_hists = parse_theory_file()

	return [data_hists, theory_hists, pythia_hists, herwig_hists, sherpa_hists]






def compile_data_and_pythia(all_hists, variables):

	data_hists, pythia_hists = all_hists[0], all_hists[1]

	
	data, pythia = [], []
	for var in variables:
		data.append( data_hists[var] )
		pythia.append( pythia_hists[var] )


	compilation = []
	for i in range(len(data[0])):
		temp = [data[0][i], data[1][i], pythia[0][i], pythia[1][i]]
		compilation.append( temp )

	return compilation








def compile_hists(var, parsed_hists, x_scale='linear'):
	
	compilation = []


	data_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[0], parsed_hists[2], parsed_hists[3], parsed_hists[4]
	max_index = len(data_hists[var])
	
	for i in range(max_index):
		
		sub_list = [ data_hists[var][i], pythia_hists[var][i], herwig_hists[var][i], sherpa_hists[var][i] ]
		
		compilation.append( sub_list )

	return compilation


def compile_hists_with_theory(var, parsed_hists, x_scale='linear'):

	compilation = []


	data_hists, theory_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[0], parsed_hists[1], parsed_hists[2], parsed_hists[3], parsed_hists[4]
	max_index = len(data_hists[var])
	
	for i in range(max_index):

		# Get the correct variable name to use for theory.
		theory_var = var.split("_")[0]

		theory_var += "_" + str( data_hists[var][i].conditions()[1][1][0] )	# Don't hardcode position of pT condition. Find the pT condition using "" in.

		if x_scale == "log":
			theory_var += "log"

		sub_list = [ data_hists[var][i], theory_hists[theory_var], pythia_hists[var][i], herwig_hists[var][i], sherpa_hists[var][i] ]

		compilation.append( sub_list )

	return compilation


default_dir = "plots/Version 5_2/"


start = time.time()


parsed_linear = parse.load_root_files_to_hist(log=False)
parsed_log = parse.load_root_files_to_hist(log=True)

parsed_hists = compile_sources(parsed_linear)	
parsed_log_hists = compile_sources(parsed_log)




end = time.time()


print "Finished parsing all files in {} seconds. Now plotting them!".format(end - start)


start = time.time()


# create_multi_page_plot(filename=default_dir + "pT.pdf", hists=compile_hists('hardest_pT', parsed_hists))

# create_multi_page_plot(filename=default_dir + "phi.pdf", hists=compile_hists('hardest_phi', parsed_hists))

# create_multi_page_plot(filename=default_dir + "eta.pdf", hists=compile_hists('hardest_eta', parsed_hists))

# create_multi_page_plot(filename=default_dir + "frac_pT_loss.pdf", hists=compile_hists('frac_pT_loss', parsed_hists))

# create_multi_page_plot(filename=default_dir + "frac_pT_loss_log.pdf", hists=compile_hists('frac_pT_loss', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')

# create_multi_page_plot(filename=default_dir + "softkill_pT_loss.pdf", hists=compile_hists('softkill_pT_loss', parsed_hists))
# create_multi_page_plot(filename=default_dir + "softkill_pT_loss_log.pdf", hists=compile_hists('softkill_pT_loss', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')

# create_multi_page_plot(filename=default_dir + "constituent_multiplicity.pdf", hists=compile_hists('mul_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "track_constituent_multiplicity.pdf", hists=compile_hists('track_mul_pre_SD', parsed_hists))


# create_multi_page_plot(filename=default_dir + "pT_D.pdf", hists=compile_hists('pT_D_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "track_pT_D.pdf", hists=compile_hists('track_pT_D_pre_SD', parsed_hists))

# create_multi_page_plot(filename=default_dir + "mass.pdf", hists=compile_hists('mass_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "track_mass.pdf", hists=compile_hists('track_mass_pre_SD', parsed_hists))



# create_multi_page_plot(filename=default_dir + "lha.pdf", hists=compile_hists('LHA_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "track_lha.pdf", hists=compile_hists('track_LHA_pre_SD', parsed_hists))

# create_multi_page_plot(filename=default_dir + "width.pdf", hists=compile_hists('width_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "track_width.pdf", hists=compile_hists('track_width_pre_SD', parsed_hists))

# create_multi_page_plot(filename=default_dir + "thrust.pdf", hists=compile_hists('thrust_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "track_thrust.pdf", hists=compile_hists('track_thrust_pre_SD', parsed_hists))


# create_multi_page_plot(filename=default_dir + "lha_log.pdf", hists=compile_hists('LHA_pre_SD', parsed_log_hists), x_scale='log')
# create_multi_page_plot(filename=default_dir + "track_lha_log.pdf", hists=compile_hists('track_LHA_pre_SD', parsed_log_hists), x_scale='log')

# create_multi_page_plot(filename=default_dir + "width_log.pdf", hists=compile_hists('width_pre_SD', parsed_log_hists), x_scale='log')
# create_multi_page_plot(filename=default_dir + "track_width_log.pdf", hists=compile_hists('track_width_pre_SD', parsed_log_hists), x_scale='log')

# create_multi_page_plot(filename=default_dir + "thrust_log.pdf", hists=compile_hists('thrust_pre_SD', parsed_log_hists), x_scale='log')
# create_multi_page_plot(filename=default_dir + "track_thrust_log.pdf", hists=compile_hists('track_thrust_pre_SD', parsed_log_hists), x_scale='log')



# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/zg/zg_10.pdf", hists=compile_hists_with_theory('zg_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/rg/rg_10.pdf", hists=compile_hists_with_theory('rg_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e1/e1_10.pdf", hists=compile_hists_with_theory('e1_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e2/e2_10.pdf", hists=compile_hists_with_theory('e2_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e05/e05_10.pdf", hists=compile_hists_with_theory('e05_10', parsed_hists), theory=True)



create_multi_page_plot(filename=default_dir + "theta_g/log/all/zg/zg_10.pdf", hists=compile_hists_with_theory('zg_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/rg/rg_10.pdf", hists=compile_hists_with_theory('rg_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e1/e1_10.pdf", hists=compile_hists_with_theory('e1_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e2/e2_10.pdf", hists=compile_hists_with_theory('e2_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e05/e05_10.pdf", hists=compile_hists_with_theory('e05_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')






# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/zg/zg_10.pdf", hists=compile_hists('track_zg_10', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/rg/rg_10.pdf", hists=compile_hists('track_rg_10', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e1/e1_10.pdf", hists=compile_hists('track_e1_10', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e2/e2_10.pdf", hists=compile_hists('track_e2_10', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e05/e05_10.pdf", hists=compile_hists('track_e05_10', parsed_hists), theory=False)


# create_multi_page_plot(filename=default_dir + "theta_g/log/track/zg/zg_10.pdf", hists=compile_hists('track_zg_10', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/rg/rg_10.pdf", hists=compile_hists('track_rg_10', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e1/e1_10.pdf", hists=compile_hists('track_e1_10', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e2/e2_10.pdf", hists=compile_hists('track_e2_10', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e05/e05_10.pdf", hists=compile_hists('track_e05_10', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')




# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/zg/zg_05.pdf", hists=compile_hists('zg_05', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/rg/rg_05.pdf", hists=compile_hists('rg_05', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e1/e1_05.pdf", hists=compile_hists('e1_05', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e2/e2_05.pdf", hists=compile_hists('e2_05', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e05/e05_05.pdf", hists=compile_hists('e05_05', parsed_hists), theory=False)



# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/zg/zg_20.pdf", hists=compile_hists('zg_20', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/rg/rg_20.pdf", hists=compile_hists('rg_20', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e1/e1_20.pdf", hists=compile_hists('e1_20', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e2/e2_20.pdf", hists=compile_hists('e2_20', parsed_hists), theory=False)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e05/e05_20.pdf", hists=compile_hists('e05_20', parsed_hists), theory=False)




# create_multi_page_plot(filename=default_dir + "theta_g/log/all/zg/zg_05.pdf", hists=compile_hists('zg_05', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/rg/rg_05.pdf", hists=compile_hists('rg_05', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e1/e1_05.pdf", hists=compile_hists('e1_05', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e2/e2_05.pdf", hists=compile_hists('e2_05', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e05/e05_05.pdf", hists=compile_hists('e05_05', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')


# create_multi_page_plot(filename=default_dir + "theta_g/log/all/zg/zg_20.pdf", hists=compile_hists('zg_20', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/rg/rg_20.pdf", hists=compile_hists('rg_20', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e1/e1_20.pdf", hists=compile_hists('e1_20', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e2/e2_20.pdf", hists=compile_hists('e2_20', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e05/e05_20.pdf", hists=compile_hists('e05_20', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')




# create_data_only_plot(filename=default_dir + "all_track_zg_10.pdf", hists=compile_data_and_pythia([ parsed_log_hists[0], parsed_log_hists[2] ], variables=['zg_10', 'track_zg_10']), labels=["Everything", "Track", "Pythia 8.215", "Pythia 8.215 (Track)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [50, 30], [50, 30]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)


# create_data_only_plot(filename=default_dir + "mass_softdrop.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['mass_pre_SD', 'mass_post_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)
# create_data_only_plot(filename=default_dir + "pT_softdrop.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['hardest_pT', 'pT_after_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)

# create_data_only_plot(filename=default_dir + "mul_softdrop.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['mul_pre_SD', 'mul_post_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)
# create_data_only_plot(filename=default_dir + "pT_D_softdrop.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['pT_D_pre_SD', 'pT_D_post_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)
# create_data_only_plot(filename=default_dir + "lha_softdrop.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['LHA_pre_SD', 'LHA_post_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)
# create_data_only_plot(filename=default_dir + "width_softdrop.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['width_pre_SD', 'width_post_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)
# create_data_only_plot(filename=default_dir + "thrust_softdrop.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['thrust_pre_SD', 'thrust_post_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)


# create_data_only_plot(filename=default_dir + "lha_softdrop_log.pdf", hists=compile_data_and_pythia([ parsed_log_hists[0], parsed_log_hists[2] ], variables=['LHA_pre_SD', 'LHA_post_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)
# create_data_only_plot(filename=default_dir + "width_softdrop_log.pdf", hists=compile_data_and_pythia([ parsed_log_hists[0], parsed_log_hists[2] ], variables=['width_pre_SD', 'width_post_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)
# create_data_only_plot(filename=default_dir + "thrust_softdrop_log.pdf", hists=compile_data_and_pythia([ parsed_log_hists[0], parsed_log_hists[2] ], variables=['thrust_pre_SD', 'thrust_post_SD']), labels=["Before SoftDrop", "After SoftDrop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio\nto\nPythia", ratio_to_index=2)







##################################################################################################### DATA ONLY PLOTS ###################################################################################################



# create_data_only_plot(filename=default_dir + "data_pT.pdf", hists=[ [ parsed_linear[0]['hardest_pT'][i], parsed_linear[0]['uncor_hardest_pT'][i] ] for i in range(len(parsed_linear[0]['uncor_hardest_pT'])) ], labels=["Jet Energy Corrected", "Jet Energy Uncorrected"], types=["error", "error"], colors=["black", "orange"], line_styles=[1, 1], ratio_to_label="Ratio\nto\nCorrected", ratio_to_index=0)

# create_data_only_plot(filename=default_dir + "jec.pdf", hists=[ [x] for x in parsed_linear[0]['jec'] ], labels=["CMS Open Data 2010"], types=["error"], colors=["black"], line_styles=[[]], ratio_plot=False)


end = time.time()

print "Finished all plotting in {} seconds.".format(end - start)
