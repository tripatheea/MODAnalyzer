from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

input_file = sys.argv[1]


def count_events(input_file, pt_cut):
  f = open(input_file, 'r')
  lines = f.read().split("\n")
  
  properties = defaultdict(list)

  number_of_events = 0

  for line in lines:
    try:
      numbers = line.split()
      if (float(numbers[3]) > pt_cut):
        number_of_events += 1;
    except:
      pass

  return number_of_events

print count_events(input_file, 150.0)
