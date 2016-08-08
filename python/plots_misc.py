from __future__ import division




import sys
import math


# matplotlib
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt


def parse_weights(input_filename):
	
	weights = []

	with open(input_filename) as f:
		for line in f:
			
			components = line.split(";")
			
			try:	
				weights.append(float(components[3]))
			except ValueError as e:
				pass
				
	return weights


def plot_weights(input_filename):
	weights = parse_weights(input_filename)

	plt.hist(weights, bins=200, normed=1)

	print max(weights)


	print "Total number of events = ", len(weights)

	plt.yscale("log")

	plt.autoscale()

	plt.xlabel("Weight")
	plt.ylabel("# of Events")

	plt.savefig("plots/weights.pdf")
	plt.clf()

plot_weights(sys.argv[1])