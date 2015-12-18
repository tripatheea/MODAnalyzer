from subprocess import call
import os
from time import time
import sys
from collections import defaultdict


def run_analyzer():
	
	location = "/home/aashish/MODMonteCarlo/data/"
	all_mc = ["pythia", "herwig", "sherpa"]

	call(['make'])

	for mc in all_mc:

		call(['rm', "/home/aashish/" + mc + "_truth.dat"])
		call(['rm', "/home/aashish/" + mc + "_reco.dat"])
		call(['rm', "/home/aashish/" + mc + "_truth_qcd.dat"])
		call(['rm', "/home/aashish/" + mc + "_reco_qcd.dat"])

		call(['./bin/analyze', location + mc + "_truth.mod", "/home/aashish/" + mc + "_truth.dat", "mc", "truth"])
		call(['./bin/analyze', location + mc + "_reco.mod", "/home/aashish/" + mc + "_reco.dat", "mc", "reco"])

		call(['./bin/analyze_beta', location + mc + "_truth.mod", "/home/aashish/" + mc + "_truth_qcd.dat", "mc", "truth"])
		call(['./bin/analyze_beta', location + mc + "_reco.mod", "/home/aashish/" + mc + "_reco_qcd.dat", "mc", "reco"])


start = time()

run_analyzer()

end = time()

print "Everything done in " + str(end - start) + " seconds!"