from __future__ import division
import sys
import os

input_path = sys.argv[1]
compare_to_path = sys.argv[2]
output_path = sys.argv[3]


def move_done_files(input_path, compare_to_path, output_path):
	
	all_mod_files = [f for f in os.listdir(input_path) if (f.endswith("mod") and f in os.listdir(compare_to_path))  ]
	# all_mod_files = [f for f in os.listdir(input_path) if (f.endswith("mod"))  ]


	

	for f in all_mod_files:
		os.rename(input_path + "/" + f, output_path + "/" + f)

	print "Done moving {} files.".format(len(all_mod_files))


move_done_files(input_path, compare_to_path, output_path)