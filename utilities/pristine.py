from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

input_path = sys.argv[1]
output_path = sys.argv[2]

def run_analyzer(input_path, output_path):
  to_analyze = []
  for f in os.listdir(input_path):
    if f.endswith("mod"):
      to_analyze.append(f)

  to_analyze.sort()

  for f in to_analyze:
    call(['./bin/convert_to_pristine', input_path + f, output_path + f, output_path])


start = time()

run_analyzer(input_path, output_path)

end = time()

print "Everything done in " + str(end - start) + " seconds!"