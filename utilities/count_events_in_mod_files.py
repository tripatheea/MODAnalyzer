from __future__ import division
import sys
import subprocess
import os

input_path = sys.argv[1]
output_file = sys.argv[2]


# proc = subprocess.Popen('ifconfig | grep inet', stdout=subprocess.PIPE, shell=True)
# output = proc.stdout.read().strip()

def get_all_number_of_events_in_files(files_path):
	for f in os.listdir(files_path):
		path_to_file = files_path + "/" + f

		proc = subprocess.Popen(["grep -o BeginEvent {} | wc -l".format(path_to_file)], shell=True, stdout=subprocess.PIPE)
		print proc.stdout()


get_all_number_of_events_in_files(input_path)