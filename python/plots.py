from __future__ import division

import time
import sys
import hists

from MODPlot import *

input_analysis_file = sys.argv[1]

def get_key(keyword):

	'''
	if 'eta' in keyword:
		return 'eta'
	elif 'pT' in keyword:
		return 'pT'
	elif 'phi' in keyword:
		return 'phi'

	return keyword.split('_')[0]
	'''
	return keyword




def parse_file(input_file):

	print "Parsing {}".format(input_file)
	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	# all_hists = hists.all_hist_templates()
	all_hists = hists.multi_page_plot_hist_templates()

	keywords = []
	keywords_set = False

	
	with open(input_file) as infile:
		
		line_number = 0

		for line in infile:

			if line_number > 10000:
			# if line_number > 100:
			# if False:
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

					pT_index = keywords.index("hardest_pT") + 1
					softdrop_pT_index = keywords.index("pT_after_SD") + 1
					eta_index = keywords.index("hardest_eta") + 1
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
									if (not mod_hist.use_prescale()) and input_file == input_analysis_file:	# For data file only.
										hist.fill_array( [float(numbers[i + 1])] )	 # + 1 because we ignore the first keyword "Entry".
									else:
										hist.fill_array( [float(numbers[i + 1])], [float(numbers[prescale_index])] )	 # + 1 because we ignore the first keyword "Entry".

						'''
						keyword = keywords[i]

						if key in all_hists.keys():

							print "tada"

							# Loop through each individual hist in all_hists[key] and fill it only if the current line satisfies all the accompanying conditions.
							for mod_hist in all_hists[key]:
								print mod_hist
								# all_hists[key].fill_array( [ float(numbers[i + 1]) ], [ float(numbers[prescale_index]) ] ) # + 1 because we ignore the first keyword "Entry".
						'''


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

data_hists = parse_file(input_analysis_file)
pythia_hists = parse_file("/home/aashish/pythia_truth.dat")
herwig_hists = parse_file("/home/aashish/herwig_truth.dat")
sherpa_hists = parse_file("/home/aashish/sherpa_truth.dat")

theory_hists = parse_theory_file()


end = time.time()

print "Finished parsing all files in {} seconds. Now plotting them!".format(end - start)


def get_hist_list(var):
	return [ data_hists[var], pythia_hists[var], herwig_hists[var], sherpa_hists[var] ]



def compile_hists(var):
	compilation = []
	for i in range(len(data_hists[var])):
		sub_list = [ data_hists[var][i], pythia_hists[var][i], herwig_hists[var][i], sherpa_hists[var][i] ]
		compilation.append( sub_list )

	return compilation


def compile_hists_with_theory(var):

	

	compilation = []
	for i in range(len(data_hists[var])):

		# Get the correct variable name to use for theory.
		theory_var = var.split("_")[0]

		theory_var += "_" + str( data_hists[var][i].conditions()[0][1][0] )

		sub_list = [ data_hists[var][i], theory_hists[theory_var], pythia_hists[var][i], herwig_hists[var][i], sherpa_hists[var][i] ]
		compilation.append( sub_list )

	return compilation


default_dir = "plots/Version 5.2/"


start = time.time()



# create_multi_page_plot(filename=default_dir + "pT.pdf", hists=compile_hists('hardest_pT'))

# create_multi_page_plot(filename=default_dir + "phi.pdf", hists=compile_hists('hardest_phi'))

# create_multi_page_plot(filename=default_dir + "eta.pdf", hists=compile_hists('hardest_eta'))


# create_multi_page_plot(filename=default_dir + "constituent_multiplicity.pdf", hists=compile_hists('mul_pre_SD'))
# create_multi_page_plot(filename=default_dir + "track_constituent_multiplicity.pdf", hists=compile_hists('track_mul_pre_SD'))


# create_multi_page_plot(filename=default_dir + "pT_D.pdf", hists=compile_hists('pT_D_pre_SD'))
# create_multi_page_plot(filename=default_dir + "track_pT_D.pdf", hists=compile_hists('track_pT_D_pre_SD'))

# create_multi_page_plot(filename=default_dir + "mass.pdf", hists=compile_hists('mass_pre_SD'))
# create_multi_page_plot(filename=default_dir + "track_mass.pdf", hists=compile_hists('track_mass_pre_SD'))

# create_multi_page_plot(filename=default_dir + "lha.pdf", hists=compile_hists('LHA_pre_SD'))
# create_multi_page_plot(filename=default_dir + "track_lha.pdf", hists=compile_hists('track_LHA_pre_SD'))

# create_multi_page_plot(filename=default_dir + "width.pdf", hists=compile_hists('width_pre_SD'))
# create_multi_page_plot(filename=default_dir + "track_width.pdf", hists=compile_hists('track_width_pre_SD'))

# create_multi_page_plot(filename=default_dir + "thrust.pdf", hists=compile_hists('thrust_pre_SD'))
# create_multi_page_plot(filename=default_dir + "track_thrust.pdf", hists=compile_hists('track_thrust_pre_SD'))



create_multi_page_plot(filename=default_dir + "theta_g/linear/all/zg/zg_10.pdf", hists=compile_hists_with_theory('zg_10'), theory=True)
create_multi_page_plot(filename=default_dir + "theta_g/linear/all/rg/rg_10.pdf", hists=compile_hists_with_theory('rg_10'), theory=True)
create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e1/e1_10.pdf", hists=compile_hists_with_theory('e1_10'), theory=True)
create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e2/e2_10.pdf", hists=compile_hists_with_theory('e2_10'), theory=True)
create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e05/e05_10.pdf", hists=compile_hists_with_theory('e05_10'), theory=True)


create_multi_page_plot(filename=default_dir + "theta_g/linear/track/zg/zg_10.pdf", hists=compile_hists('track_zg_10'), theory=False)
create_multi_page_plot(filename=default_dir + "theta_g/linear/track/rg/rg_10.pdf", hists=compile_hists('track_rg_10'), theory=False)
create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e1/e1_10.pdf", hists=compile_hists('track_e1_10'), theory=False)
create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e2/e2_10.pdf", hists=compile_hists('track_e2_10'), theory=False)
create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e05/e05_10.pdf", hists=compile_hists('track_e05_10'), theory=False)


'''
print "Plotting eta!"

eta_plot = MODPlot( get_hist_list('hardest_eta'), plot_types=plot_types, plot_colors=colors, plot_labels=labels, ratio_plot=True, ratio_to_index=1, ratio_label="Ratio\nto\nPythia", x_label="Jet $\eta$", y_label="A.U.", x_lims=(-5., 5.))
eta_plot.plot("hardest_eta.pdf")


print "Plotting phi!"

phi_plot = MODPlot( get_hist_list('hardest_phi'), plot_types=plot_types, plot_colors=colors, plot_labels=labels, ratio_plot=True, ratio_to_index=1, ratio_label="Ratio\nto\nPythia", x_label="Jet $\phi$", y_label="A.U.", x_lims=(0, 2*np.pi))
phi_plot.plot("hardest_phi.pdf")

'''

'''
print "Plotting pT!"


pT_plot = MODPlot( get_hist_list('hardest_pT'), plot_types=plot_types, plot_colors=colors, plot_labels=labels, y_scale='log', ratio_plot=True, ratio_to_index=1, ratio_label="Ratio\nto\nPythia", x_label="Jet $p_T$", y_label="A.U.")
pT_plot.plot()
pT_plot.save_plot("hardest_pT.pdf")
'''



'''
print "Plotting constituent multiplicity!"

constituent_multiplicity_plot = MODPlot( get_hist_list('mul_pre_SD'), plot_types=plot_types, plot_colors=colors, plot_labels=labels, ratio_plot=True, ratio_to_index=1, ratio_label="Ratio\nto\nPythia", x_label="Constituent Multiplicity", y_label="A.U.")
constituent_multiplicity_plot.plot("constituent_multiplicity.pdf")


print "Plotting zg!"

constituent_multiplicity_plot = MODPlot( get_hist_list('zg_10'), plot_types=plot_types, plot_colors=colors, plot_labels=labels, ratio_plot=True, ratio_to_index=1, ratio_label="Ratio\nto\nPythia", x_label="$z_g$", y_label="A.U.")
constituent_multiplicity_plot.plot("zg.pdf")
'''


end = time.time()

print "Finished all plotting in {} seconds.".format(end - start)
