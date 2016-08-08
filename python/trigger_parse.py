from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


from rootpy.io import File as TFile



import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt



output_directory = "./"
# output_directory = "/home/aashish/Dropbox (MIT)/Research/data/June Generation (MC)/"
# output_directory = "/media/aashish/My Files/Dropbox (MIT)/Research/data/June Generation (MC)/"

# data_file = "/home/aashish/Dropbox (MIT)/Research/data/June Generation (MC)/analyzed/experiment/trig.dat"
data_file = "/media/aashish/Transcend/experiment/trig.dat"





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


			# if line_number > 10000:	# Ideal length.
			# if line_number > 100000000:	# Huge.
			# if line_number > 1000000:	# Big enough.
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
					pT_index = keywords.index("corr_hardest_pT") + 1
					trig_name_index = keywords.index("trigger_name") + 1

					

						

					# print numbers[pT_index]

					pT_of_this_event = float(numbers[pT_index])
				
					prescale = float(numbers[prescale_index])
				
					current_line_trig_name = numbers[trig_name_index]

					# print pT_of_this_event, prescale, current_line_trig_name


					
					for key in all_hists.keys():
						if key in current_line_trig_name:
							
							# print key, current_line_trig_name
							# print all_hists
							# print all_hists[key]

							mod_hist = all_hists[key]
							hist = mod_hist.hist()

							conditions = mod_hist.conditions()

							try:
								condition_satisfied = 1
								for condition_keyword, condition_func in conditions:

									keyword_index = keywords.index(condition_keyword[0]) + 1
									condition_func_param = numbers[keyword_index]
									
									# print condition_keyword

									if condition_keyword[0] == "jet_quality":
										condition_satisfied *= int(condition_func(int(condition_keyword[1]), int(condition_func_param)))
									elif condition_keyword[0] == "trig_jet_matched":
										# print condition_keyword[1], condition_func_param, condition_func(int(condition_keyword[1]), int(condition_func_param))
										condition_satisfied *= int(condition_func(int(condition_keyword[1]), int(condition_func_param)))
									else:
										condition_satisfied *= int(condition_func(condition_keyword[1], condition_func_param))
								

									

								condition_satisfied = bool(condition_satisfied)
							except Exception as e:
								print "ASF", e
							
							if condition_satisfied:

								x = pT_of_this_event
								y = prescale

								# print "filling ", key, "with", x, y
								all_hists[key].hist().fill_array( [x], [y] )
								break


							# print condition_func_param, condition_func("Jet30U", condition_func_param)
						

						'''
						if keyword in all_hists.keys():
							
							for mod_hist in all_hists[keyword]:
								hist = mod_hist.hist()
								conditions = mod_hist.conditions()

								try:
									condition_satisfied = 1
									for condition_keyword, condition_func in conditions:

										
									
										keyword_index = keywords.index(condition_keyword[0]) + 1
										condition_func_param = numbers[keyword_index]
										
										# print condition_keyword

										if condition_keyword[0] == "jet_quality":
											condition_satisfied *= int(condition_func(int(condition_keyword[1]), int(condition_func_param)))
											# print condition_func_param, condition_keyword[1], condition_func(int(condition_keyword[1]), int(condition_func_param)), "; ",
										else:
											condition_satisfied *= int(condition_func(condition_keyword[1], condition_func_param))
									

										

									condition_satisfied = bool(condition_satisfied)
								except Exception as e:
									print "ASF", e
								
								if condition_satisfied:

									x = float(numbers[i + 1]) # + 1 because we ignore the first keyword "Entry".
									y = float(numbers[prescale_index])

									hist.fill_array( [x], [y] )

							# print condition_func_param, condition_func("Jet30U", condition_func_param)
						'''
			except Exception as e:
				print "Some exception occured!",
				print e
				print line
				continue


	return all_hists



def parse_to_root_file(input_filename, output_filename, hist_templates):

	print "Parsing {} to {}".format(input_filename, output_filename)
	
	parsed_hists = parse_file( input_filename, copy.deepcopy( hist_templates ) )

	f = TFile(output_filename, "RECREATE")

	for var in parsed_hists.keys():
		
		index = 0

		mod_hist = parsed_hists[var]
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

		mod_hist = hists[var]

	
		hist_name = "{}#{}".format(var, index)

		# Get hist from ROOT file.
		hist = root_file.Get(hist_name)

		mod_hist.replace_hist(hist)

		index += 1

	return hists




def parse_to_root_files():
	
	hist_templates = hists.trigger_hists()

	parse_to_root_file(input_filename=data_file, output_filename=output_directory + "trig_main.root", hist_templates=hist_templates)

def load_root_files_to_hist(log=False):
	
	hist_templates = hists.trigger_hists()

	filenames = ["trig.root"]

	return  [ root_file_to_hist("./"+ filename, hist_templates) for filename in filenames ] 




if __name__ == "__main__":

	
	parse_to_root_files()

	# load_root_files_to_hist()

	pass