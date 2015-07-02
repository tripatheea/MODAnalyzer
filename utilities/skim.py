from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

input_path = sys.argv[1]
output_file_path = sys.argv[2]

def run_skimmer(input_path, output_file_path):
  to_analyze = []
  for f in os.listdir(input_path):
    if f.endswith("mod"):
      if not os.path.exists(output_file_path + "/" + f):
        to_analyze.append(f)

  to_analyze.sort()

  for f in to_analyze:
    call(['./bin/skim', input_path + f, output_file_path + f])
  


start = time()

run_skimmer(input_path, output_file_path)

end = time()

print "Everything done in " + str(end - start) + " seconds!"