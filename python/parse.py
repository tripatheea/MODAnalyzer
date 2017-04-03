from __future__ import division

import math
import time
import sys
import hists

import subprocess
import os.path

from MODPlot import *


from rootpy.io import File as TFile


import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt



# output_directory = "/home/aashish/root/macros/MODAnalyzer/parsed_root_files/"
# output_directory = "/media/aashish/My Files/Dropbox (MIT)/Research/data/June Generation (MC)/"
# output_directory = "/media/aashish/My Files/Dropbox (MIT)/"
# output_directory = "/home/aashish/Feb22/"
# output_directory = "/media/aashish/opendata_mod/Feb5/histogrammed/pristine/"
output_directory = "/media/aashish/opendata_mod/Feb5/histogrammed/"
# output_directory = "/media/aashish/7CA48778A48733A4/Mar-13-analysis/MC/sherpa/"





# data_file   = "/media/aashish/My Files/Dropbox (MIT)/Research/data/MC/analyzed/pristine.dat"
# data_file   = "/home/aashish/data.dat"
# data_file   = "/home/aashish/Feb5/analyzed/pristine.dat"
data_file   = "/media/aashish/opendata_mod/Feb5/analyzed/pristine.dat"
# data_file   = "/media/aashish/opendata_mod/Feb5/analyzed/pt/0004.dat"
# data_file   = "/home/aashish/root/macros/MODAnalyzer/abcde.dat"
# data_file   = "/media/aashish/opendata_mod/Feb5/analyzed/pt.dat"
# data_file   = "/home/aashish/root/macros/MODAnalyzer/test.dat"
pythia_file = "/media/aashish/My Files/Dropbox (MIT)/Research/data/MC/pythia/part_1_2_3_4.dat"
herwig_file = "/media/aashish/69B4875A78793A41/MC/herwig/herwig.dat"
# sherpa_file = "/media/aashish/7CA48778A48733A4/Mar-13-analysis/MC/sherpa/sherpa.dat"
sherpa_file = "/media/aashish/7CA48778A48733A4/Mar-13-analysis/MC/sherpa/sherpa_2_3_4_5_6_7.dat"


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



def parse_file(input_file, output_filename, all_hists, log_hists):

	print "Parsing {}".format(input_file)
	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	keywords = []
	keywords_set = False

	with open(input_file) as infile:
		
		line_number = 0


		largest_pT = 0.0


		for line in infile:

			start = time.time()

			# if line_number > 10000:		# Ideal length.
			# if line_number > 100000:		# Big enough.
			# if line_number > 1000:		# Small tests.
			# if line_number > 30000:		# Small tests.
			# if line_number > 1000000:		# Small tests.
			# if line_number > 150000:		# Small tests.
			# if line_number > 100000:		# Small tests.
			# if line_number > 2:				# Trivial test.
			if False:
				break

			if len(line.strip()) == 0:
				continue

			line_number += 1

			if line_number % 10000 == 0:
				print "At line number {}".format(line_number)
				print "Largest pT so far is", largest_pT
				write_to_root_files(all_hists, log_hists, output_filename)

			try:
			# if True:
				numbers = line.split()

				if numbers[0] == "#" and (not keywords_set):
					keywords = numbers[1:]
					keywords_set = True

					# Build a dictionary of keyword indices.
					keywords_index_dictionary = {}
					for keyword_index, keyword in enumerate(keywords):
						keywords_index_dictionary[keyword] = keyword_index



					# print keywords_index_dictionary['hardest_pT']

					prescale_index = keywords.index("prescale") 

					

				elif numbers[0] == "Entry":

					# Find what prescale to use.
					pT_of_this_event = float(numbers[keywords_index_dictionary['hardest_pT']]) # + 1 because we ignore the first keyword "Entry".


					# Find out what prescale to use.
					if input_file == data_file:	# For data file only.

						# print " i am data lol"
						#  Average prescale. 

						# To find which prescale to use, we need to find which trigger fired. 
						# To do that, we need to find the pT of the hardest jet.
						
						prescale_to_use = 0.0

						if pT_of_this_event > 250.:
							prescale_to_use = 1.0
						else:

							for pT_boundaries, prescale in average_prescales.items():
								
								lower, upper = pT_boundaries

								# print lower, upper, pT_of_this_event

								if upper != None:
									if pT_of_this_event > float(lower) and pT_of_this_event < float(upper):
										prescale_to_use = prescale
										break

					else:	# MC so always use prescales.
						# print "and i am a mc",
						prescale_to_use = float(numbers[prescale_index]) 
						# print prescale_to_use
					



					for some_hists in [all_hists, log_hists]:

						# for i in xrange(len(keywords)):
						for keyword in some_hists.keys():

							for mod_hist in some_hists[keyword]:
									
								
								hist = mod_hist.hist()

								# Now, see whether or not the conditions are met.
								condition_satisfied = True
								conditions = mod_hist.conditions()

								for condition_keyword, condition_boundaries in conditions:
									keyword_index = keywords_index_dictionary[condition_keyword]

									# print condition_keyword, keyword_index, numbers[keyword_index]

									
									if condition_boundaries[0] != None and float(numbers[keyword_index]) < float(condition_boundaries[0]):
										condition_satisfied = False
										
									if condition_boundaries[1] != None and float(numbers[keyword_index]) > float(condition_boundaries[1]):
										condition_satisfied = False
									
									# print condition_boundaries, float(numbers[keyword_index]), condition_satisfied	


								if condition_satisfied:			

									if pT_of_this_event > largest_pT:
										largest_pT = pT_of_this_event

									# if keyword == 'hardest_eta':
									# print conditions					
									
									# print "i am here"
									# print keyword, keywords_index_dictionary[keyword]
									x = float(numbers[keywords_index_dictionary[keyword]]) 
									
									hist.fill_array( [x], [prescale_to_use] )	 

													
										
									
			except Exception as e:
				print "Some exception occured!",
				print e

			end = time.time()

			# print "Took {} seconds for current line.".format(end - start)

	print "Largest pT was", largest_pT


	return all_hists, log_hists



def parse_from_conditional_files(output_directory, source, all_hists, log_hists):
	# print output_directory, source, all_hists


	for keyword, mod_hists in all_hists.items():
		for mod_hist in all_hists[keyword]:
			
			conditional_filename = get_conditional_filename(mod_hist.conditions())
			file_to_read = output_directory + "/" + source + "/conditions/" + conditional_filename + ".dat"

			# Now, read the file.

	
	pass



def get_conditional_filename(conditions):
	conditional_filename = ""
	for condition in conditions:
			condition_keyword = condition[0]
			condition_boundaries = ""
			for x in condition[1]:
				condition_boundaries += str(x) + "-"

			
			conditional_filename += condition_keyword + "-" + condition_boundaries[:-1] + "_"

	conditional_filename = conditional_filename[:-1]

	return conditional_filename    

def filter_events(input_filename, output_directory, source, all_mod_hists, all_mod_log_hists):

	# Get a list of unique conditions.

	all_conditions = set()

	for mod_hist_dicts in all_mod_hists:
		for keyword, mod_hists in mod_hist_dicts.items():
			for mod_hist in mod_hist_dicts[keyword]:
				
				conditions = mod_hist.conditions()

				all_conditions.add(tuple(conditions))


	# Build a keyword dictionary.
	
	keywords_index_dictionary = {}

	with open(input_filename) as infile:
		
		for line in infile:
			if len(line.strip()) == 0:
				continue

			numbers = line.split()

			
			if numbers[0] == "#":

				keyword_line = numbers[0]

				# Build a dictionary of keyword indices.
				keywords_index_dictionary = {}
				for keyword_index, keyword in enumerate(numbers):
					keywords_index_dictionary[keyword] = keyword_index + 1  # + 1 to offset for 0-counting.

				break

	# Create a new directory in the output_directory.
	subprocess.call(["mkdir", "-p", "{}/{}/conditions/".format(output_directory, source) ])
			
	for conditions in all_conditions:

		if_conditional = "if ( "

		for condition in conditions:
			
			# Build corresponding awk from condition.
			condition_keyword = condition[0]
			condition_boundaries = condition[1]

			keyword_index = keywords_index_dictionary[condition[0]]

			if condition_boundaries[0] != None and condition_boundaries[1] != None:
				if_conditional += "${} > {} && ${} < {}".format(keyword_index, float(condition_boundaries[0]), keyword_index, float(condition_boundaries[1])) 
			elif condition_boundaries[0] != None and condition_boundaries[1] == None:
				if_conditional += "${} > {}".format(keyword_index, float(condition_boundaries[0]))
			elif condition_boundaries[0] == None and condition_boundaries[1] != None:
				if_conditional += "${} < {}".format(keyword_index, float(condition_boundaries[0]))

			if_conditional += " && "

		if_conditional = if_conditional[:-3] + ")"

		# Now, build a filename.

		conditional_filename = get_conditional_filename(conditions)
		
		if_conditional = if_conditional + "{ print $0; }"
		
		proc = subprocess.Popen("echo '" + keyword_line + "' > " + "{}/{}/conditions/{}.dat".format(output_directory, source, conditional_filename), stdout=subprocess.PIPE, shell=True)

		command = "awk '{" + if_conditional + "}' " + input_filename + " >> " + "{}/{}/conditions/{}.dat".format(output_directory, source, conditional_filename)
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
		output = proc.stdout.read()
		
		print command
		print output
		

def write_to_root_files(parsed_hists, parsed_log_hists, output_filename):
	start = time.time()

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

	end = time.time()

	# print "\n" * 3
	# print "Took {} seconds to write down to files.".format(end - start)


def parse_to_root_file(input_filename, output_filename, hist_templates):

	# print "Parsing {} to {}".format(input_filename, output_filename)
	
	# First, get the already histogrammed hists.
	

	if os.path.exists(output_filename[0]) and os.path.exists(output_filename[1]):


		if "pythia" in input_filename or "herwig" in input_filename or "sherpa" in input_filename:
			is_this_data = False
		else:
			is_this_data = True

		hists = root_file_to_hist(output_filename[0], hist_templates[0], is_this_data)
		log_hists = root_file_to_hist(output_filename[1], hist_templates[1], is_this_data)
	else:
		hists, log_hists = copy.deepcopy(hist_templates[0]), copy.deepcopy(hist_templates[1])

	start = time.time()
	
	# print "Is this data?", is_this_data

	parsed_hists, parsed_log_hists = parse_file( input_filename, output_filename, hists, log_hists )
	
	# parsed_hists, parsed_log_hists = parse_from_conditional_files( output_directory, source, copy.deepcopy( hist_templates[0] ), copy.deepcopy(hist_templates[1]) )
	

	end = time.time()

	
	# print "Took {} seconds to parse everything- this is where I deepcopy.".format(end - start)
	write_to_root_files(parsed_hists, parsed_log_hists, output_filename)

	


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

	# parse_to_root_file(input_filename=data_file, output_filename=(output_directory + "data.root", output_directory + "data_log.root"), hist_templates=(hist_templates, log_hist_templates))
	# parse_to_root_file(input_filename=pythia_file, output_filename=(output_directory + "pythia.root", output_directory + "pythia_log.root"), hist_templates=(hist_templates, log_hist_templates))
	parse_to_root_file(input_filename=herwig_file, output_filename=(output_directory + "herwig.root", output_directory + "herwig_log.root"), hist_templates=(hist_templates, log_hist_templates))
	# parse_to_root_file(input_filename=sherpa_file, output_filename=(output_directory + "sherpa.root", output_directory + "sherpa_log.root"), hist_templates=(hist_templates, log_hist_templates))


def load_root_files_to_hist(log=False):
	
	if not log:
		hist_templates = hists.multi_page_plot_hist_templates()
		filenames = ["data.root", "pythia.root", "herwig.root", "sherpa.root"]
		# filenames = ["data.root", "pythia.root", "sherpa.root", "sherpa.root"]
		# filenames = ["data.root", "pythia.root", "herwig.root", "pythia.root"]
		# filenames = ["data.root", "data.root", "data.root", "data.root"]
		# filenames = ["pt.root", "pt.root", "pt.root", "pt.root"]
	else:
		hist_templates = hists.multi_page_log_plot_hist_templates()
		filenames = ["data_log.root", "pythia_log.root", "herwig_log.root", "sherpa_log.root"]
		# filenames = ["data_log.root", "pythia_log.root", "sherpa_log.root", "sherpa_log.root"]
		# filenames = ["data_log.root", "pythia_log.root", "pythia_log.root", "pythia_log.root"]
		# filenames = ["data_log.root", "data_log.root", "data_log.root", "data_log.root"]
		# filenames = ["pt.root", "pt.root", "pt.root", "pt.root"]

	return  [ root_file_to_hist(output_directory + filename, hist_templates, is_this_data) for filename, is_this_data in zip(filenames, [True, False, False, False]) ] 




if __name__ == "__main__":

	# filter_events(input_filename=pythia_file, output_directory=output_directory, source="pythia", all_mod_hists=(hists.multi_page_plot_hist_templates(), hists.multi_page_log_plot_hist_templates()))
	
	parse_to_root_files()

	# load_root_files_to_hist()


	pass
