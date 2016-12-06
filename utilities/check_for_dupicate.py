from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

input_file = sys.argv[1]

def check_for_duplicate_events(input_file):
	
  duplicates = 0  
  total = 0
  unique = 0

  all_events = {}
  with open(input_file) as f:
	for line in f.readlines():

		if total % 10000 == 0:
			print "On line number {}".format(total)

		if len(line.strip()) != 0:
			run_number, event_number = line.split()
			
			if (run_number, event_number) in all_events.keys():
				print "Found a duplicate already."
				duplicates += 1
			else:
				all_events[(run_number, event_number)] = 1
				unique += 1

			total += 1 

  return total, unique, duplicates

start = time()

total, unique, duplicates = check_for_duplicate_events(input_file)

print "Parsed {} lines.".format(total)
print "Found {} unique events.".format(unique)
print "Found {} duplicate events.".format(duplicates)

end = time()

print "Everything done in " + str(end - start) + " seconds!"