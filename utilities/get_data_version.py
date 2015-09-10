import os

def get_version(input_file):

  f = open(input_file, 'r')
  lines = f.read().split("\n")

  properties = defaultdict(list)

  for line in lines:
    try:
      numbers = line.split()
      
      if not numbers[0] == "#":
        if (float(numbers[4]) > pT_lower_cut) and (float(numbers[17]) > pfc_pT_cut) and (float(numbers[4]) < pT_upper_cut):
          properties['event_number'].append( float( numbers[1] ) )
          properties['run_number'].append( float( numbers[2] ) )

          properties['uncorrected_hardest_pts'].append( float( numbers[3] ) )