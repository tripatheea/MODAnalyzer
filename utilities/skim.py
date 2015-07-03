from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

input_path = sys.argv[1]
output_file_path = sys.argv[2]
error_log_path = sys.argv[3]

def run_skimmer(input_path, output_file_path, error_log_path):
  to_analyze = []
  for f in os.listdir(input_path):
    if f.endswith("mod"):
      if not os.path.exists(output_file_path + "/" + f):
        to_analyze.append(f)

  to_analyze.sort()

  # Open the error log file to empty it.
  f = open(error_log_path + "skim_error_log.log", "w")
  f.write("\n")
  f.close()

  # Create a temporary file to count the total number of events. 
  # We'll delete it at the end.
  f = open(error_log_path + "skim_error_log.log" + ".num", "w")
  f.close()

  for f in to_analyze:
    call(['./bin/skim', input_path + f, output_file_path + f, error_log_path + "skim_error_log.log"])
  


start = time()

run_skimmer(input_path, output_file_path, error_log_path)

end = time()


# Find the total number of events processed.
number = sum([int(line.rstrip('\n')) for line in open(error_log_path + "skim_error_log.log" + ".num")])

print "\n\nSkimmed " + str(number) + " events in " + str(end - start) + " seconds!"

# Write the end result to the log file.
f = open(error_log_path + "skim_error_log.log", 'a')
f.write("\n\nSkimmed " + str(number) + " events in " + str(end - start) + " seconds!")
f.close()

# Delete the num file.
os.remove(error_log_path + "skim_error_log.log" + ".num")