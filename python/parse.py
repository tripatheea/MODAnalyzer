from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


from rootpy.io import File as TFile



import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt



# output_directory = "/home/aashish/root/macros/MODAnalyzer/parsed_root_files/"
# output_directory = "/media/aashish/My Files/Dropbox (MIT)/Research/data/June Generation (MC)/"
# output_directory = "/media/aashish/My Files/Dropbox (MIT)/"
# output_directory = "/home/aashish/Feb22/"
output_directory = "/media/aashish/opendata_mod/Feb5/histogrammed/pristine/"
# output_directory = "/home/aashish/Feb5/histogrammed/"





# data_file   = "/media/aashish/My Files/Dropbox (MIT)/Research/data/MC/analyzed/pristine.dat"
# data_file   = "/home/aashish/data.dat"
# data_file   = "/home/aashish/Feb5/analyzed/pristine.dat"
data_file   = "/media/aashish/opendata_mod/Feb5/analyzed/pristine.dat"
# data_file   = "/media/aashish/opendata_mod/Feb5/analyzed/pt/0004.dat"
# data_file   = "/home/aashish/root/macros/MODAnalyzer/abcde.dat"
# data_file   = "/media/aashish/opendata_mod/Feb5/analyzed/pt.dat"
# data_file   = "/home/aashish/root/macros/MODAnalyzer/test.dat"
pythia_file = "/media/aashish/opendata_mod/Feb5/analyzed/pythia.dat"
# herwig_file = "/media/aashish/My Files/Dropbox (MIT)/Research/data/MC/analyzed/herwig.dat"
# sherpa_file = "/media/aashish/My Files/Dropbox (MIT)/Research/data/MC/analyzed/sherpa.dat"


pfc_data_file = "/media/aashish/My Files/Dropbox (MIT)/Research/data/MC/analyzed/pfc.dat"
pfc_pythia_file = "/media/aashish/My Files/Dropbox (MIT)/Research/data/MC/analyzed/pythia_pfc.dat"
pfc_herwig_file = "/media/aashish/My Files/Dropbox (MIT)/Research/data/MC/analyzed/herwig_pfc.dat"
pfc_sherpa_file = "/media/aashish/My Files/Dropbox (MIT)/Research/data/MC/analyzed/sherpa_pfc.dat"

average_prescales = {}

average_prescales[(250, None)] = 1.	
average_prescales[(200, 250)] = 1.933420103
average_prescales[(150, 200)] = 5.361922609
average_prescales[(115, 150)] =  100.3122906
average_prescales[(85, 115)] =  851.3943491




my_prescales = defaultdict(float)
my_numbers = defaultdict(int)

def parse_file(input_file, all_hists, log_hists):

	f = open("./pt_output.dat", 'w')

	print "Parsing {}".format(input_file)
	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	keywords = []
	keywords_set = False

	with open(input_file) as infile:
		
		line_number = 0


		all_pT_s = []


		for line in infile:


			if line_number > 10000:	# Ideal length.
			# if line_number > 100000:	# Big enough.
			# if line_number > 100:		# Small tests.
			# if line_number > 30000:		# Small tests.
			# if line_number > 1000000:		# Small tests.
			# if line_number > 150000:		# Small tests.
			# if line_number > 100000:		# Small tests.
			# if False:
				break

			if len(line.strip()) == 0:
				continue

			line_number += 1

			if line_number % 1000 == 0:
				print "At line number {}".format(line_number)

			try:
				numbers = line.split()

				# print numbers

				# print numbers[2]

				if numbers[0] == "#" and (not keywords_set):
					keywords = numbers[2:]
					keywords_set = True

				elif numbers[0] == "Entry":

					prescale_index = keywords.index("prescale") + 1

					# print numbers

					for i in range(len(keywords)):

						keyword = keywords[i]
						
						if keyword == "hardest_pT":
							pT_of_this_event = float(numbers[i + 1]) # + 1 because we ignore the first keyword "Entry".

							all_pT_s.append(pT_of_this_event)

						if keyword in all_hists.keys():
							
							for mod_hist in all_hists[keyword]:
							# for mod_hist in [all_hists['hardest_pT'][0]]:
								hist = mod_hist.hist()
								conditions = mod_hist.conditions()

								# print conditions

								# print conditions

								condition_satisfied = True

								try:					
									for condition_keyword, condition_boundaries in conditions:
										keyword_index = keywords.index(condition_keyword) + 1

										
										
										if condition_boundaries[0] != None and float(numbers[keyword_index]) < float(condition_boundaries[0]):
											condition_satisfied = False
											

										if condition_boundaries[1] != None and float(numbers[keyword_index]) > float(condition_boundaries[1]):
											condition_satisfied = False
											


									
								except Exception as e:
									print "ASF", e
								
								# if keyword == "jec":
								# 	print condition_satisfied



								if condition_satisfied:

									x = float(numbers[i + 1]) # + 1 because we ignore the first keyword "Entry".


									'''			
									if i == len(keywords) - 1:
									


										# if pT_of_this_event > 130.4 and pT_of_this_event < 140.2:
										if pT_of_this_event > 148.6 and pT_of_this_event < 150.:
											# print numbers[2], numbers[4], numbers[8]
											# print "Test"
											f.write("{}\t{}\t{}\n".format(numbers[2], numbers[4], numbers[8]))
										# hardest_pt => 2; jec => 3


									continue
									'''

									if input_file == data_file:	# For data file only.

										if not mod_hist.use_prescale():
											hist.fill_array( [x] )
										else:
											
											#  Average prescale. 

											# To find which prescale to use, we need to find which trigger fired. 
											# To do that, we need to find the pT of the hardest jet.
											
											prescale_to_use = 0.0

											if pT_of_this_event > 250.:
												prescale_to_use = 1.0
											else:
												for pT_boundaries, prescale in average_prescales.items():
													
													lower, upper = pT_boundaries

													if upper != None:
														if pT_of_this_event > float(lower) and pT_of_this_event < float(upper):
															prescale_to_use = prescale
															break

											# print "Using the prescale", prescale_to_use
											hist.fill_array( [x], [prescale_to_use] )	 

									
									else:	# MC so always use prescales.
										# if type()
										hist.fill_array( [x], [float(numbers[prescale_index])] )

									# if keyword == "hardest_area":
										# print "Tada"
										# print [type(x)], [float(numbers[prescale_index])]


						if keyword in log_hists.keys():
							
							for mod_hist in log_hists[keyword]:
								hist = mod_hist.hist()
								conditions = mod_hist.conditions()

								condition_satisfied = True

								try:					
									for condition_keyword, condition_boundaries in conditions:
										keyword_index = keywords.index(condition_keyword) + 1

										
										
										if condition_boundaries[0] != None and float(numbers[keyword_index]) < float(condition_boundaries[0]):
											condition_satisfied = False
											

										if condition_boundaries[1] != None and float(numbers[keyword_index]) > float(condition_boundaries[1]):
											condition_satisfied = False
											


									
								except Exception as e:
									print "ASF", e
								
								# if keyword == "jec":
								# 	print condition_satisfied



								if condition_satisfied:

									x = float(numbers[i + 1]) # + 1 because we ignore the first keyword "Entry".

									if input_file == data_file:	# For data file only.

										if not mod_hist.use_prescale():
											hist.fill_array( [x] )
										else:
											
											#  Average prescale. 

											# To find which prescale to use, we need to find which trigger fired. 
											# To do that, we need to find the pT of the hardest jet.
											
											prescale_to_use = 0.0

											if pT_of_this_event > 250.:
												prescale_to_use = 1.0
											else:
												for pT_boundaries, prescale in average_prescales.items():
													
													lower, upper = pT_boundaries

													if upper != None:
														if pT_of_this_event > float(lower) and pT_of_this_event < float(upper):
															prescale_to_use = prescale
															break
											
											hist.fill_array( [x], [prescale_to_use] )	 
									
									else:	# MC so always use prescales.
										# if type()
										hist.fill_array( [x], [float(numbers[prescale_index])] )

									# if keyword == "hardest_area":
										# print "Tada"
										# print [type(x)], [float(numbers[prescale_index])]

			except Exception as e:
				print "Some exception occured!",
				print e

	# print "Largest pT was", max(all_pT_s)




	return all_hists, log_hists



def parse_to_root_file(input_filename, output_filename, hist_templates):

	print "Parsing {} to {}".format(input_filename, output_filename)
	
	
	parsed_hists, parsed_log_hists = parse_file( input_filename, copy.deepcopy( hist_templates[0] ), copy.deepcopy(hist_templates[1]) )


	f = TFile(output_filename[0], "RECREATE")

	for var in parsed_hists.keys():
		
		
	
		index = 0

		for mod_hist in parsed_hists[var]:
			hist = copy.deepcopy( mod_hist.hist() )
			hist.SetName("{}#{}".format(var, index))
			hist.Write()

			index += 1

	f.Close()




	f = TFile(output_filename[1], "RECREATE")

	for var in parsed_log_hists.keys():
		
		index = 0

		for mod_hist in parsed_log_hists[var]:
			hist = copy.deepcopy( mod_hist.hist() )
			hist.SetName("{}#{}".format(var, index))
			hist.Write()

			index += 1

	f.Close()


def root_file_to_hist(input_filename, hist_templates, is_this_data):

	hists = copy.deepcopy( hist_templates )

	root_file = TFile(input_filename, "read")

	for var in hists.keys():
		
		index = 0

		# if var in ['hardest_pT', 'uncor_hardest_pT', 'hardest_eta']:
		# if var not in ['uncor_hardest_pT']:
		if is_this_data:
			for mod_hist in hists[var]:
				hist_name = "{}#{}".format(var, index)

				# Get hist from ROOT file.
				hist = root_file.Get(hist_name)

				mod_hist.replace_hist(hist)

				index += 1

		else:
			if var != 'uncor_hardest_pT':
				for mod_hist in hists[var]:
					hist_name = "{}#{}".format(var, index)

					# Get hist from ROOT file.
					hist = root_file.Get(hist_name)

					mod_hist.replace_hist(hist)

					index += 1

	return hists




def parse_to_root_files():
	hist_templates = hists.multi_page_plot_hist_templates()
	log_hist_templates = hists.multi_page_log_plot_hist_templates()

	parse_to_root_file(input_filename=data_file, output_filename=(output_directory + "data.root", output_directory + "data_log.root"), hist_templates=(hist_templates, log_hist_templates))
	# parse_to_root_file(input_filename=data_file, output_filename=(output_directory + "pt.root", output_directory + "pt_log.root"), hist_templates=(hist_templates, log_hist_templates))
	# parse_to_root_file(input_filename=pythia_file, output_filename=(output_directory + "pythia.root", output_directory + "pythia_log.root"), hist_templates=(hist_templates, log_hist_templates))
	# parse_to_root_file(input_filename=herwig_file, output_filename=(output_directory + "herwig.root", output_directory + "herwig_log.root"), hist_templates=(hist_templates, log_hist_templates))
	# parse_to_root_file(input_filename=sherpa_file, output_filename=(output_directory + "sherpa.root", output_directory + "sherpa_log.root"), hist_templates=(hist_templates, log_hist_templates))


def load_root_files_to_hist(log=False):
	
	if not log:
		hist_templates = hists.multi_page_plot_hist_templates()
		# filenames = ["data.root", "pythia.root", "herwig.root", "sherpa.root"]
		# filenames = ["data.root", "pythia.root", "pythia.root", "pythia.root"]
		filenames = ["data.root", "data.root", "data.root", "data.root"]
		# filenames = ["pt.root", "pt.root", "pt.root", "pt.root"]
	else:
		hist_templates = hists.multi_page_log_plot_hist_templates()
		# filenames = ["data_log.root", "pythia_log.root", "herwig_log.root", "sherpa_log.root"]
		# filenames = ["data_log.root", "pythia_log.root", "pythia_log.root", "pythia_log.root"]
		filenames = ["data_log.root", "data_log.root", "data_log.root", "data_log.root"]
		# filenames = ["pt.root", "pt.root", "pt.root", "pt.root"]

	return  [ root_file_to_hist(output_directory + filename, hist_templates, is_this_data) for filename, is_this_data in zip(filenames, [True, False, False, False]) ] 




if __name__ == "__main__":

	parse_to_root_files()

	# load_root_files_to_hist()

	pass
