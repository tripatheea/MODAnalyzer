from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

input_path = sys.argv[1]
output_path = sys.argv[2]

def run_analyzer(input_path, output_path):
  
    # Get files we already have in one jet form.
    already_analyzed = [f for f in os.listdir(output_path) if f.endswith("mod")]

    # print already_analyzed
    print len(already_analyzed), len(os.listdir(input_path))

    to_analyze = []
    for f in os.listdir(input_path):
        if f.endswith("mod") and f not in already_analyzed:
            to_analyze.append(f)

    to_analyze.sort()


    print to_analyze

    for f in to_analyze:
        # call(['./bin/convert_to_one_jet', input_path + f, output_path + f, output_path])
        pass


start = time()

run_analyzer(input_path, output_path)

end = time()

print "Everything done in " + str(end - start) + " seconds!"