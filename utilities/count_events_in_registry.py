from __future__ import division
import sys

registry_file = sys.argv[1]



all_event_numbers = set()
all_run_numbers = set()
all_filenames = set()

all_event_numbers_list = []
all_run_numbers_list = []
all_filenames_list = []

with open(registry_file) as f:
	for i, line in enumerate(f):
		line_contents = line.split(" ")
		event_number, run_number, file_name = line_contents[0], line_contents[1], line_contents[3]

		all_event_numbers.add(event_number)
		all_run_numbers.add(run_number)
		all_filenames.add(file_name)

		all_event_numbers_list.append(event_number)
		all_run_numbers_list.append(run_number)
		all_filenames_list.append(file_name)

		if i % 1000000 == 0:
			print "At line number {}".format(i)


print "=" * 50
print "According to the sets:"
print "Total number of lines : {}".format(i + 1)
print "Total number of events: {}".format(len(all_event_numbers_list))
print "Total number of runs  : {}".format(len(all_run_numbers))
print "Total number of files : {}".format(len(all_filenames))


