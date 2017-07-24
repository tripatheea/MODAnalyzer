from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


from rootpy.io import File as TFile



import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt



# output_directory = "/media/aashish/opendata_mod/Feb5/histogrammed/Apr2/"
# output_directory = "/media/aashish/opendata_mod/Feb5/histogrammed/backup/wed_morning/"
output_directory = "/media/aashish/My Files/Dropbox (MIT)/histogrammed/"
# output_directory = "/home/aashish/Dropbox (MIT)/Research/data/June Generation (MC)/"
# output_directory = "/media/aashish/My Files/Dropbox (MIT)/Research/data/June Generation (MC)/"
# output_directory = "/home/aashish/"

# data_file = "/home/aashish/Dropbox (MIT)/Research/data/June Generation (MC)/analyzed/experiment/trig.dat"
data_file = "/media/aashish/opendata_mod/Feb5/analyzed/trigger/trigger.dat"



average_prescales = {}

average_prescales[(250, None)] = 1.	
average_prescales[(200, 250)] = 1.933420103
average_prescales[(150, 200)] = 5.361922609
average_prescales[(115, 150)] =  100.3122906
average_prescales[(85, 115)] =  851.3943491



def parse_file(input_file, all_hists):

	print "Parsing {}".format(input_file)
	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	keywords = []
	keywords_set = False

	with open(input_file) as infile:
		
		line_number = 0

		current_line_trig_name = ""

		for line in infile:


			# if line_number > 100000:	# Ideal length.
			# if line_number > 100000000:	# Very Huge.
			# if line_number > 10000000:	# Huge.
			# if line_number > 1000000:	# Big enough.
			# if line_number > 1000:		# Small tests.
			# if line_number > 30000:		# Small tests.
			if False:
				break

			line_number += 1

			if line_number % 1000000 == 0:
				print "At line number {}".format(line_number)

			try:
				numbers = line.split()

				if len(numbers) > 0 and numbers[0] == "#" and (not keywords_set):
					keywords = numbers[2:]
					keywords_set = True

				elif len(numbers) > 0 and numbers[0] == "Entry":



					prescale_index = keywords.index("prescale") + 1
					pT_index = keywords.index("corr_hardest_pT") + 1
					trig_name_index = keywords.index("trigger_name") + 1
					eta_index = keywords.index("hardest_eta") + 1

					pT_of_this_event = float(numbers[pT_index])
					eta_of_this_event = float(numbers[eta_index])
					prescale = float(numbers[prescale_index])
					current_line_trig_name = numbers[trig_name_index]

					if abs(eta_of_this_event) > 2.4:
						# print "oops, pseudorapidity too large.", eta_of_this_event
						continue

					for key in all_hists.keys():
						
						# See whether condition/s are satisfied or not.
						mod_hist = all_hists[key]
						hist = mod_hist.hist()

						if key in current_line_trig_name or ("prescale" in key and key.split("_prescale")[0] in current_line_trig_name):
													
							conditions = mod_hist.conditions()

							try:
								condition_satisfied = True

								for condition_keyword, condition_func in conditions:
									keyword_index = keyword_index = keywords.index(condition_keyword[0]) + 1
									condition_func_param = numbers[keyword_index]

									if (not condition_func(condition_keyword[1], int(condition_func_param))):
										condition_satisfied = False

							except Exception as excp:
								print "Some exception occured while processing conditions.", excp


							if condition_satisfied:
								x = pT_of_this_event
								y = prescale
								
								if "prescale" in key:
									#  Average prescale. 
									all_hists[key].hist().fill_array( [y] )
								else:
									# Trigger pT's.
									all_hists[key].hist().fill_array([x], [y])


			except Exception as exc:
				print "Some exception occured.", exc
				print line
				print "\n"
					

	return all_hists



def parse_to_root_file(input_filename, output_filename, hist_templates):

	print "Parsing {} to {}".format(input_filename, output_filename)
	
	parsed_hists = parse_file( input_filename, copy.deepcopy( hist_templates ) )

	f = TFile(output_filename, "RECREATE")

	for var in parsed_hists.keys():
		mod_hist = parsed_hists[var]
		hist = copy.deepcopy( mod_hist.hist() )
		hist.SetName("{}".format(var))
		hist.Write()

	f.Close()


def root_file_to_hist(input_filename, hist_templates):

	hists = copy.deepcopy( hist_templates )

	root_file = TFile(input_filename, "read")

	for var in hists.keys():
		mod_hist = hists[var]
		hist = root_file.Get(var)
		mod_hist.replace_hist(hist)

		
	return hists




def parse_to_root_files():
	
	hist_templates = hists.trigger_hists()

	parse_to_root_file(input_filename=data_file, output_filename=output_directory + "trig.root", hist_templates=hist_templates)

def load_root_files_to_hist(log=False):
	
	hist_templates = hists.trigger_hists()

	filenames = ["trig.root"]

	return  [ root_file_to_hist(output_directory + filename, hist_templates) for filename in filenames ] 




if __name__ == "__main__":

	
	parse_to_root_files()

	# load_root_files_to_hist()

	pass