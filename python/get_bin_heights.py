from __future__ import division
import sys

import parse

average_prescales = {}

average_prescales[(250, None)] = 1. 
average_prescales[(200, 250)] = 1.933420103
average_prescales[(150, 200)] = 5.361922609
average_prescales[(115, 150)] =  100.3122906
average_prescales[(85, 115)] =  851.3943491


def get_bin_width(hist):
    if hist.x_scale() == "log":
        return (np.log(hist.hist().upperbound()) - np.log(hist.hist().lowerbound())) / hist.hist().nbins()
    
    return (hist.hist().upperbound() - hist.hist().lowerbound()) / hist.hist().nbins()

def normalize_hist(hist):
    bin_width = get_bin_width(hist)

    if hist.hist().GetSumOfWeights() != 0.0:
        hist.hist().Scale(1.0 / ( hist.hist().GetSumOfWeights() * bin_width ))
    else:
        hist.hist().Scale(0.0)

def get_number_of_events(height, bin_width, prescale_used):
    return height / (prescale_used)

def get_bin_heights(hist):
    root_hist = hist.hist()

    lower_bound, upper_bound = root_hist.bounds()
    n_bins = root_hist.nbins()
    bin_width = get_bin_width(hist)

    normalize_hist(hist)


    bin_heights = []
    try:
        # pT_bounds = filter(lambda x: x[0] == "hardest_pT", hist.conditions())[0][1]

        # prescale_used = average_prescales[pT_bounds]
    
        for bin_count in root_hist.bins_range():
            # bin_heights.append((root_hist.GetBinCenter(bin_count) - 0.5 * bin_width, get_number_of_events(root_hist.GetBinContent(bin_count), bin_width, prescale_used), root_hist.GetBinContent(bin_count)))
            bin_heights.append((root_hist.GetBinCenter(bin_count) - 0.5 * bin_width, root_hist.GetBinContent(bin_count)))
    except Exception as e:
        print "Some error occured.", e
        pass

    return bin_heights


def write_bin_heights(heights, f):
    
    for lower_edge, height in heights:
        f.write("{}\t{:.15f}\n".format(lower_edge, height))

if __name__ == '__main__':
    output_filename = sys.argv[1]
    
    f = open(output_filename, "w")
    f.close()

    all_filled_hists = parse.load_root_files_to_hist()

    

    f = open(output_filename, "a")

    # print filter(lambda x: x[0] == "hardest_pT", pT_hists[4].conditions())[0][1]
    pT_hists = all_filled_hists[0]['hardest_pT']
    f.write("Corrected pT:\n")
    for i in range(len(pT_hists)):
        pT_bounds = filter(lambda x: x[0] == "hardest_pT", pT_hists[i].conditions())[0][1]
        f.write(str(pT_bounds) + "\n")
        write_bin_heights(get_bin_heights(pT_hists[i]), f)
        f.write("\n")


    pT_hists = all_filled_hists[0]['uncor_hardest_pT']
    f.write("Uncorrected pT:\n")
    for i in range(len(pT_hists)):
        pT_bounds = filter(lambda x: x[0] == "hardest_pT", pT_hists[i].conditions())[0][1]
        f.write(str(pT_bounds) + "\n")
        write_bin_heights(get_bin_heights(pT_hists[i]), f)
        f.write("\n")