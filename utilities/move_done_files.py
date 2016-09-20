from __future__ import division
import sys
import os

roots_file_path = sys.argv[1]
done_root_files_path = sys.argv[2]
mod_file_path = sys.argv[3]


def move_done_files(roots_file_path, done_root_files_path, mod_file_path):
	
	all_mod_files = [f[:-4] for f in os.listdir(mod_file_path) if f.endswith("mod")  ]


	roots_file_path_to_move = [f for f in os.listdir(roots_file_path) if f.endswith("root") and f[:-5] in all_mod_files]

	for root_file in roots_file_path_to_move:
		os.rename(roots_file_path + "/" + root_file, done_root_files_path + "/" + root_file)

	print "Done moving {} files.".format(len(roots_file_path_to_move))


move_done_files(roots_file_path, done_root_files_path, mod_file_path)