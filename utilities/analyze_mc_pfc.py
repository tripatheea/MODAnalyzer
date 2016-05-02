from subprocess import call
import os
from time import time
import sys
from collections import defaultdict


input_path = sys.argv[1]
output_file = sys.argv[2]

# def run_analyzer():
	
# 	location = "/Users/aashish/MODMonteCarlo/data/"
# 	all_mc = ["pythia", "herwig", "sherpa"]

# 	call(['make'])

# 	for mc in all_mc:

# 		# call(['rm', "/Users/aashish/" + mc + "_truth.dat"])
# 		# call(['rm', "/Users/aashish/" + mc + "_reco.dat"])

# 		call(['./bin/analyze', location + mc + "_truth.mod", "/Users/aashish/" + mc + "_truth_final.dat"])
# 		# call(['./bin/analyze', location + mc + "_reco.mod", "/Users/aashish/" + mc + "_reco.dat"])


def analyze_mc(input_path, output_file):
	
	files_to_process = []

	for f in os.listdir(input_path):
		if f.endswith('.mod'):
			files_to_process.append(input_path + '/' + f)

	for f in files_to_process:
		call(['./bin/analyze_pfc', f, output_file])		


start = time()

analyze_mc(input_path, output_file)

end = time()

print "Everything done in " + str(end - start) + " seconds!"