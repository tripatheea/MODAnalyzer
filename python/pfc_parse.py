from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


from rootpy.io import File as TFile


import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt

output_directory = "/media/aashish/opendata_mod/Feb5/histogrammed/jul8/"


pfc_data_file = "/media/aashish/opendata_mod/Feb5/analyzed/pfc/jul13/data_pfc.dat"
pfc_pythia_file = "/media/aashish/opendata_mod/Feb5/analyzed/pfc/jul13/pythia_pfc.dat"
pfc_herwig_file = "/media/aashish/opendata_mod/Feb5/analyzed/pfc/jul13/herwig_pfc.dat"
pfc_sherpa_file = "/media/aashish/opendata_mod/Feb5/analyzed/pfc/jul13/sherpa_pfc.dat"

average_prescales = {}

average_prescales[(250, None)] = 1.0
average_prescales[(200, 250)] = 1.933420103
average_prescales[(150, 200)] = 5.361922609
average_prescales[(115, 150)] = 100.3122906
average_prescales[(85, 115)] = 851.3943491


my_prescales = defaultdict(float)
my_numbers = defaultdict(int)


def parse_file(input_file, output_filename, all_hists):

    print "Parsing {}".format(input_file)

    # We read the file line by line, and for each line, we fill the corresponding histograms.
    # This is desirable to creating lists of values since this will not hold
    # anything in memory.

    keywords = []
    keywords_set = False

    with open(input_file) as infile:

        line_number = 0

        all_pT_s = []

        for line in infile:

            # if line_number > 10000:  # Ideal length.
            # if line_number > 100000:	# Big enough.
            # if line_number > 100:		# Small tests.
            # if line_number > 30000:		# Small tests.
            # if line_number > 1000000:		# Small tests.
            # if line_number > 150000:		# Small tests.
            # if line_number > 100000:		# Small tests.
            # if line_number > 100000:		# Small tests.
            if False:
                break

            if len(line.strip()) == 0:
                continue

            line_number += 1

            if line_number % 1000000 == 0 and line_number > 1:
                print "At line number {}".format(line_number)
                write_to_root_file(output_filename, all_hists)

            try:
                numbers = line.split()

                if numbers[0] == "#" and (not keywords_set):
                    keywords = numbers[2:]
                    keywords_set = True

                elif numbers[0] == "Entry":

                    prescale_index = keywords.index("prescale") + 1

                    for i in range(len(keywords)):

                        keyword = keywords[i]

                        if keyword == "hardest_pT":
                            # + 1 because we ignore the first keyword "Entry".
                            pT_of_this_event = float(numbers[i + 1])

                            all_pT_s.append(pT_of_this_event)

                        if keyword in all_hists.keys():

                            for mod_hist in all_hists[keyword]:
                                hist = mod_hist.hist()
                                conditions = mod_hist.conditions()

                                try:
                                    condition_satisfied = 1
                                    for condition_keyword, condition_boundaries in conditions:
                                        keyword_index = keywords.index(
                                            condition_keyword) + 1

                                        # print condition_boundaries[0],
                                        # condition_boundaries[1],
                                        # keyword_index

                                        # print numbers[keyword_index]

                                        if len(condition_boundaries) <= 2:
                                            if condition_boundaries[0] == None and condition_boundaries[1] != None:
                                                condition_satisfied *= int(
                                                    float(numbers[keyword_index]) < float(condition_boundaries[1]))
                                            elif condition_boundaries[0] != None and condition_boundaries[1] == None:
                                                condition_satisfied *= int(
                                                    float(numbers[keyword_index]) > float(condition_boundaries[0]))
                                            elif condition_boundaries[0] == None and condition_boundaries[1] == None:
                                                condition_satisfied *= 1
                                            elif condition_boundaries[0] != None and condition_boundaries[1] != None:
                                                condition_satisfied *= int(float(numbers[keyword_index]) > float(condition_boundaries[
                                                                           0]) and float(numbers[keyword_index]) < float(condition_boundaries[1]))
                                        else:

                                            if condition_boundaries[0] == "in":
                                                condition_satisfied *= int(
                                                    int(numbers[keyword_index]) in list(condition_boundaries[1:]))
                                            else:
                                                condition_satisfied *= int(
                                                    int(numbers[keyword_index]) not in list(condition_boundaries[1:]))

                                    condition_satisfied = bool(
                                        condition_satisfied)
                                except Exception as e:
                                    print "ASF", e

                                # if keyword == "jec":
                                # 	print condition_satisfied

                                if condition_satisfied:

                                    # + 1 because we ignore the first keyword "Entry".
                                    x = float(numbers[i + 1])

                                    if input_file == pfc_data_file:  # For data file only.

                                        if not mod_hist.use_prescale():
                                            hist.fill_array([x])
                                        else:

                                            #  Average prescale.

                                            # To find which prescale to use, we need to find which trigger fired.
                                            # To do that, we need to find the
                                            # pT of the hardest jet.

                                            prescale_to_use = 0.0

                                            if pT_of_this_event > 250.:
                                                prescale_to_use = 1.0
                                            else:
                                                for pT_boundaries, prescale in average_prescales.items():

                                                    lower, upper = pT_boundaries

                                                    if upper != None:
                                                        if pT_of_this_event > float(lower) and pT_of_this_event < float(upper):
                                                            prescale_to_use = prescale
                                                            break

                                            hist.fill_array(
                                                [x], [prescale_to_use])

                                    else:  # MC so always use prescales.
                                        # if type()
                                        hist.fill_array(
                                            [x], [float(numbers[prescale_index])])

                                    # if keyword == "hardest_area":
                                        # print "Tada"
                                        # print [type(x)],
                                        # [float(numbers[prescale_index])]

            except Exception as e:
                print "Some exception occured!",
                print e

    # print "Largest pT was", max(all_pT_s)

    return all_hists


def write_to_root_file(output_filename, parsed_hists):
    f = TFile(output_filename, "RECREATE")

    for var in parsed_hists.keys():

        index = 0

        for mod_hist in parsed_hists[var]:
            hist = copy.deepcopy(mod_hist.hist())
            hist.SetName("{}#{}".format(var, index))
            hist.Write()

            index += 1

    f.Close()


def parse_to_root_file(input_filename, output_filename, hist_templates):

    print "Parsing {} to {}".format(input_filename, output_filename)

    parsed_hists = parse_file(
        input_filename, output_filename, copy.deepcopy(hist_templates))

    write_to_root_file(output_filename, parsed_hists)


def root_file_to_hist(input_filename, hist_templates):

    hists = copy.deepcopy(hist_templates)

    root_file = TFile(input_filename, "read")

    for var in hists.keys():

        index = 0

        for mod_hist in hists[var]:
            hist_name = "{}#{}".format(var, index)

            # Get hist from ROOT file.
            hist = root_file.Get(hist_name)

            mod_hist.replace_hist(hist)

            index += 1

    return hists


def parse_pfc_to_root_files():
    hist_templates = hists.get_pfc_hists()

    # parse_to_root_file(input_filename=pfc_data_file, output_filename=output_directory +
    #                    "data_pfc.root", hist_templates=hist_templates)
    # parse_to_root_file(input_filename=pfc_pythia_file, output_filename=output_directory +
    #                    "pythia_pfc.root", hist_templates=hist_templates)
    # parse_to_root_file(input_filename=pfc_herwig_file,
    # output_filename=output_directory +
    # "herwig_pfc.root", hist_templates=hist_templates)
    parse_to_root_file(input_filename=pfc_sherpa_file, output_filename=output_directory +
                       "sherpa_pfc.root", hist_templates=hist_templates)


def load_pfc_root_files_to_hist():
    hist_templates = hists.get_pfc_hists()

    filenames = ['data_pfc.root', 'pythia_pfc.root',
                 'herwig_pfc.root', 'sherpa_pfc.root']
    # filenames = ['data_pfc.root', 'data_pfc.root', 'data_pfc.root', 'data_pfc.root']
    # filenames = ['data_pfc.root', 'pythia_pfc.root', 'pythia_pfc.root',
    # 'pythia_pfc.root']

    return [root_file_to_hist(output_directory + filename, hist_templates) for filename in filenames]


if __name__ == "__main__":

    # load_root_files_to_hist()

    parse_pfc_to_root_files()

    pass
