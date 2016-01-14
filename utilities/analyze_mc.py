from subprocess import call
import os
from time import time
import sys
from collections import defaultdict


def run_analyzer():
	
	location = "/Users/aashish/MODMonteCarlo/data/"
	all_mc = ["pythia", "herwig", "sherpa"]

	call(['make'])

	for mc in all_mc:

		# call(['rm', "/Users/aashish/" + mc + "_truth.dat"])
		# call(['rm', "/Users/aashish/" + mc + "_reco.dat"])

		call(['./bin/analyze', location + mc + "_truth.mod", "/Users/aashish/" + mc + "_truth_final.dat"])
		# call(['./bin/analyze', location + mc + "_reco.mod", "/Users/aashish/" + mc + "_reco.dat"])


start = time()

run_analyzer()

end = time()

print "Everything done in " + str(end - start) + " seconds!"