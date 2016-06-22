from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


from rootpy.io import File as TFile



import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt



output_directory = "/home/aashish/root/macros/MODAnalyzer/parsed_root_files/"

data_file = "/home/aashish/pristine.dat"
pythia_file = "/home/aashish/pythia_truth.dat"
herwig_file = "/home/aashish/herwig_truth.dat"
sherpa_file = "/home/aashish/sherpa_truth.dat"




def parse_file(input_file, all_hists):

	print "Parsing {}".format(input_file)
	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	keywords = []
	keywords_set = False

	with open(input_file) as infile:
		
		line_number = 0

		for line in infile:


			# if line_number > 10000:	# Ideal length.
			# if line_number > 100000:	# Big enough.
			# if line_number > 1000:		# Small tests.
			# if line_number > 30000:		# Small tests.
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

									if (not mod_hist.use_prescale()) and input_file == data_file:	# For data file only.
										hist.fill_array( [x] )	 
									else:
										hist.fill_array( [x], [float(numbers[prescale_index])] )

			except Exception as e:
				print "Some exception occured!"
				print e


	return all_hists



def parse_to_root_file(input_filename, output_filename, hist_templates):

	print "Parsing {} to {}".format(input_filename, output_filename)
	
	parsed_hists = parse_file( input_filename, copy.deepcopy( hist_templates ) )

	f = TFile(output_filename, "RECREATE")

	for var in parsed_hists.keys():
		
		index = 0

		for mod_hist in parsed_hists[var]:
			hist = copy.deepcopy( mod_hist.hist() )
			hist.SetName("{}#{}".format(var, index))
			hist.Write()

			index += 1

	f.Close()


def root_file_to_hist(input_filename, hist_templates):

	hists = copy.deepcopy( hist_templates )

	root_file = TFile(input_filename, "read")

	for var in hists.keys():
		
		index = 0

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

	# parse_to_root_file(input_filename=data_file, output_filename=output_directory + "data.root", hist_templates=hist_templates)
	# parse_to_root_file(input_filename=pythia_file, output_filename=output_directory + "pythia.root", hist_templates=hist_templates)
	# parse_to_root_file(input_filename=herwig_file, output_filename=output_directory + "herwig.root", hist_templates=hist_templates)
	parse_to_root_file(input_filename=sherpa_file, output_filename=output_directory + "sherpa.root", hist_templates=hist_templates)

	# parse_to_root_file(input_filename=data_file, output_filename=output_directory + "data_log.root", hist_templates=log_hist_templates)
	# parse_to_root_file(input_filename=pythia_file, output_filename=output_directory + "pythia_log.root", hist_templates=log_hist_templates)
	# parse_to_root_file(input_filename=herwig_file, output_filename=output_directory + "herwig_log.root", hist_templates=log_hist_templates)
	parse_to_root_file(input_filename=sherpa_file, output_filename=output_directory + "sherpa_log.root", hist_templates=log_hist_templates)



def load_root_files_to_hist(log=False):
	hist_templates = hists.multi_page_plot_hist_templates()

	if log:
		filenames = ["data.root", "pythia.root", "herwig.root", "sherpa.root"]
	else:
		filenames = ["data_log.root", "pythia_log.root", "herwig_log.root", "sherpa_log.root"]

	return  [ root_file_to_hist(output_directory + filename, hist_templates) for filename in filenames ] 




# parse_to_root_files()

load_root_files_to_hist()