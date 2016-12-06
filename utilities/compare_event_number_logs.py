from __future__ import division
import sys
import os
from shutil import copyfile
from subprocess import call

mod_log_path = sys.argv[1]
registry_log_path = sys.argv[2]

input_path = sys.argv[3]
output_path = sys.argv[4]

def compare(mod_log_path, registry_log_path):

	mod_numbers = {}
	with open(mod_log_path) as mod_log:
		for line in mod_log.readlines():
			elements = line.split()
			mod_numbers[elements[0]] = elements[1]

	registry_numbers = {}
	with open(registry_log_path) as registry_log:
		for line in registry_log.readlines():
			components = line.split()
			registry_numbers[components[0]] = components[1]

	files_not_found = []
	mismatch_files = []

	total_missing_event_numbers = 0
	total_missing_event_numbers_from_absent_mod_files = 0

	for filename, number in registry_numbers.items():
		if filename in mod_numbers.keys():
			if mod_numbers[filename] != registry_numbers[filename]:
				mismatch_files.append((filename, mod_numbers[filename], registry_numbers[filename]))
				total_missing_event_numbers += int(registry_numbers[filename]) - int(mod_numbers[filename])
		else:
			files_not_found.append(filename)
			total_missing_event_numbers_from_absent_mod_files += int(registry_numbers[filename])

	if len(files_not_found) != 0:
		print
		print "The following files were not found in the MOD directoy:"

		for f in files_not_found:
			# print f
			print "\"root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Jet/AOD/Apr21ReReco-v1/0004/" + str(f) + ".root\", ",

		print
		print "This accounts for missing {} events.".format(total_missing_event_numbers_from_absent_mod_files)
		print

	if len(mismatch_files) != 0:
		print
		print "There's a mismatch of event numbers for the following files:"

		print "\n[",
		for f in mismatch_files:
			# print "{} => {} for MOD vs. {} for registry".format(f[0], f[1], f[2])
			print "\"{}.root\", ".format(f[0]),

		print "]"
		print
		print "This accounts for a mismatch of {} events.".format(total_missing_event_numbers)
		print

	if (len(files_not_found) == 0 and len(mismatch_files) == 0):
		print "\nEverything matches perfectly! :D \n"

	return files_not_found, mismatch_files

def move_missing_files(input_path, output_path):
	a = compare(mod_log_path, registry_log_path)

	for f in os.listdir(input_path):
		# if f.endswith("mod"):
		# to_analyze.append(f)

		if f[:-5] in a[0]:
			print "Moving {} to {}".format(input_path + "/" + f, output_path + "/" + f)
			call("mv {} {}".format(input_path + "/" + f, output_path + "/" + f), shell=True)


compare(mod_log_path, registry_log_path)
# move_missing_files(input_path, output_path)