from subprocess import call
import os
from time import time, sleep
import sys
from collections import defaultdict


input_path = sys.argv[1]
output_path = sys.argv[2]
registry_path = sys.argv[3]
completed_log = sys.argv[4]

def run_analyzer(input_path, output_path, registry_path, completed_log):
	to_analyze = []
	for f in os.listdir(input_path):
		if f.endswith("mod"):
			to_analyze.append(f)


	to_analyze.sort()

	for f in to_analyze:
		call(['./bin/move_events_to_correct_file', input_path, f, registry_path, output_path, completed_log])
		
		sleep(2)

		call(['rm', input_path + "/" + f])



start = time()

run_analyzer(input_path, output_path, registry_path, completed_log)

end = time()

print "Everything done in " + str(end - start) + " seconds!"