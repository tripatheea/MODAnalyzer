from __future__ import division

import math
import time
import sys
import hists

from MODPlot import *

from rootpy.io.pickler import dump, load

input_analysis_file = sys.argv[1]




def parse_file(input_file, all_hists):

	print "Parsing {}.".format(input_file)
	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	keywords = []
	keywords_set = False

	with open(input_file) as infile:
		
		line_number = 0

		for line in infile:


			# if line_number > 10000:	# Ideal length.
			# if line_number > 100000:	# Big enough.
			# if line_number > 100:		# Small tests.
			if False:
				break

			line_number += 1

			if line_number % 10000 == 0:
				print "At line number {}".format(line_number)

			try:
				numbers = line.split()

				if numbers[0] == "#" and (not keywords_set):
					keywords = numbers[2:]
					keywords_set = True

				elif numbers[0] == "Entry":

					prescale_index = keywords.index("prescale") + 1

					for i in range(len(keywords)):

						keyword = keywords[i]

						if keyword in all_hists.keys():
							
							for mod_hist in all_hists[keyword]:
								hist = mod_hist.hist()
								conditions = mod_hist.conditions()

								condition_satisfied = 1
								for condition_keyword, condition_boundaries in conditions:
									keyword_index = keywords.index(condition_keyword) + 1

									if condition_boundaries[0] == None and condition_boundaries[1] != None:
										condition_satisfied *= int( float(numbers[keyword_index]) < condition_boundaries[1] ) 
									elif condition_boundaries[0] != None and condition_boundaries[1] == None:
										condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] ) 
									elif condition_boundaries[0] == None and condition_boundaries[1] == None:
										condition_satisfied *= 1 
									elif condition_boundaries[0] != None and condition_boundaries[1] != None:
										condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] and float(numbers[keyword_index]) < condition_boundaries[1] )

								condition_satisfied = bool(condition_satisfied)


								if condition_satisfied:

									x = float(numbers[i + 1]) # + 1 because we ignore the first keyword "Entry".

									if (not mod_hist.use_prescale()) and input_file == input_analysis_file:	# For data file only.
										hist.fill_array( [x] )	 
									else:
										hist.fill_array( [x], [float(numbers[prescale_index])] )

			except:
				pass


	return all_hists


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


start = time.time()

def parse_general():

	hist_templates = hists.multi_page_plot_hist_templates()

	data_hists = parse_file(input_analysis_file, copy.deepcopy(hist_templates))
	pythia_hists = parse_file("/home/aashish/pythia_truth.dat", copy.deepcopy(hist_templates))
	herwig_hists = parse_file("/home/aashish/herwig_truth.dat", copy.deepcopy(hist_templates))
	sherpa_hists = parse_file("/home/aashish/sherpa_truth.dat", copy.deepcopy(hist_templates))
	
	theory_hists = parse_theory_file()

	return [data_hists, theory_hists, pythia_hists, herwig_hists, sherpa_hists]


def parse_log():

	hist_templates = hists.multi_page_log_plot_hist_templates()

	log_data_hists = parse_file(input_analysis_file, copy.deepcopy(hist_templates))
	log_pythia_hists = parse_file("/home/aashish/pythia_truth.dat", copy.deepcopy(hist_templates))
	log_herwig_hists = parse_file("/home/aashish/herwig_truth.dat", copy.deepcopy(hist_templates))
	log_sherpa_hists = parse_file("/home/aashish/sherpa_truth.dat", copy.deepcopy(hist_templates))

	theory_hists = parse_theory_file()

	return [log_data_hists, theory_hists, log_pythia_hists, log_herwig_hists, log_sherpa_hists]


def parse_data_only():

	hist_templates = hists.multi_page_data_only_plot_hist_templates()

	data_hists = parse_file("/home/aashish/data.dat", copy.deepcopy(hist_templates))

	corrected = data_hists['cor_hardest_pT']
	uncorrected = data_hists['uncor_hardest_pT']

	compiled = []
	for i in range(len(corrected)):
		temp = [corrected[i], uncorrected[i]]
		compiled.append(temp)


	return compiled



end = time.time()

print "Finished parsing all files in {} seconds. Now plotting them!".format(end - start)



def load_root_file_to_hists(root_filename):
	hlist = load(root_filename)

	print hlist

def save_hists_to_root_file(root_filename, hists_to_save):
	for var in hists_to_save:
		print hists_to_save[var]

	dump(hists_to_save, root_filename)


def compile_hists(var, hists, x_scale='linear'):
	
	compilation = []

	if x_scale == "log":
		log_data_hists, log_pythia_hists, log_herwig_hists, log_sherpa_hists = parsed_hists[0], parsed_hists[2], parsed_hists[3], parsed_hists[4]
		max_index = len(log_data_hists[var])
	else:
		data_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[0], parsed_hists[2], parsed_hists[3], parsed_hists[4]
		max_index = len(data_hists[var])
	
	for i in range(max_index):
		
		if x_scale == "linear":
			sub_list = [ data_hists[var][i], pythia_hists[var][i], herwig_hists[var][i], sherpa_hists[var][i] ]
		elif x_scale == "log":
			sub_list = [ log_data_hists[var][i], log_pythia_hists[var][i], log_herwig_hists[var][i], log_sherpa_hists[var][i] ]

		compilation.append( sub_list )

	return compilation


def compile_hists_with_theory(var, parsed_hists, x_scale='linear'):

	compilation = []

	if x_scale == "log":
		log_data_hists, theory_hists, log_pythia_hists, log_herwig_hists, log_sherpa_hists = parsed_hists[0], parsed_hists[1], parsed_hists[2], parsed_hists[3], parsed_hists[4]
		max_index = len(log_data_hists[var])
	else:
		data_hists, theory_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[0], parsed_hists[1], parsed_hists[2], parsed_hists[3], parsed_hists[4]
		max_index = len(data_hists[var])
	
	for i in range(max_index):

		# Get the correct variable name to use for theory.
		theory_var = var.split("_")[0]

		theory_var += "_" + str( data_hists[var][i].conditions()[0][1][0] )

		if x_scale == "log":
			theory_var += "log"

		if x_scale == "linear":
			sub_list = [ data_hists[var][i], theory_hists[theory_var], pythia_hists[var][i], herwig_hists[var][i], sherpa_hists[var][i] ]
		elif x_scale == "log":
			sub_list = [ log_data_hists[var][i], theory_hists[theory_var], log_pythia_hists[var][i], log_herwig_hists[var][i], log_sherpa_hists[var][i] ]

		compilation.append( sub_list )

	return compilation


default_dir = "plots/Version 5_2/"


start = time.time()



# parsed_hists = parse_general()
# parsed_log_hists = parse_log()

parsed_data_only_hists = parse_data_only()


create_data_only_plot(filename=default_dir + "data_pT.pdf", hists=parsed_data_only_hists)


# create_multi_page_plot(filename=default_dir + "pT.pdf", hists=compile_hists('hardest_pT', parsed_hists))

# create_multi_page_plot(filename=default_dir + "phi.pdf", hists=compile_hists('hardest_phi', parsed_hists))

# create_multi_page_plot(filename=default_dir + "eta.pdf", hists=compile_hists('hardest_eta', parsed_hists))


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



# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/zg/zg_10.pdf", hists=compile_hists_with_theory('zg_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/rg/rg_10.pdf", hists=compile_hists_with_theory('rg_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e1/e1_10.pdf", hists=compile_hists_with_theory('e1_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e2/e2_10.pdf", hists=compile_hists_with_theory('e2_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e05/e05_10.pdf", hists=compile_hists_with_theory('e05_10', parsed_hists), theory=True)



# create_multi_page_plot(filename=default_dir + "theta_g/log/all/zg/zg_10.pdf", hists=compile_hists_with_theory('zg_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
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


end = time.time()

print "Finished all plotting in {} seconds.".format(end - start)
