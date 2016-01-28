from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

input_path = sys.argv[1]
output_file = sys.argv[2]

def run_lumi(input_path, output_file):
  to_analyze = []
  for f in os.listdir(input_path):
    if f.endswith("mod"):
      to_analyze.append(f)

  to_analyze.sort()

  for f in to_analyze:
    call(['./bin/luminosity', input_path + f, output_file])


start = time()

run_lumi(input_path, output_file)

end = time()

print "Everything done in " + str(end - start) + " seconds!"


# /media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/MOD/Apr21ReReco-v1/0000