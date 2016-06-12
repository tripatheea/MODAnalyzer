from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

rg_path = sys.argv[1]
zg_path = sys.argv[2]
output_path = sys.argv[3]


def compute_theta_g(rg_path, zg_path, output_path):
	
	

	for f in os.listdir(rg_path):
		zg_s = []
		rg_s = []
		with open(rg_path + f) as f2:
			for line in f2:
				components = line.split()
				c2 = [float(c) for c in components]
				rg_s.append( c2 )

		f3 = open(output_path + f.replace("rg", "theta_g"),'w')
		for line in rg_s:
			x = line[0] 
			output = "{} {} {} {}\n".format(x, line[1], line[2], line[3])
			f3.write(output)


start = time()

compute_theta_g(rg_path, zg_path, output_path)

end = time()

print "Everything done in " + str(end - start) + " seconds!"