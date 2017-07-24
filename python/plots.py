from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


import parse
import pfc_parse


def parse_theory_file():

    input_dir = "/home/aashish/theory distributions/20/"

    pTs = [85, 115, 150, 200]

    input_files = []
    input_files.extend(["zg/zg_" + str(pT) for pT in pTs])
    input_files.extend(["rg/rg_" + str(pT) for pT in pTs])
    input_files.extend(["e1/e1_" + str(pT) for pT in pTs])
    input_files.extend(["e2/e2_" + str(pT) for pT in pTs])
    input_files.extend(["e05/e05_" + str(pT) for pT in pTs])

    input_files.extend(["zg/zg_log" + str(pT) for pT in pTs])
    input_files.extend(["rg/rg_log" + str(pT) for pT in pTs])
    input_files.extend(["e1/e1_log" + str(pT) for pT in pTs])
    input_files.extend(["e2/e2_log" + str(pT) for pT in pTs])
    input_files.extend(["e05/e05_log" + str(pT) for pT in pTs])

    input_files.extend(["zg/zg_85long", "zg/zg_150long", "zg/zg_250long",
                        "zg/zg_log85long", "zg/zg_log150long", "zg/zg_log250long"])
    input_files.extend(["rg/rg_85long", "rg/rg_150long", "rg/rg_250long",
                        "rg/rg_log85long", "rg/rg_log150long", "rg/rg_log250long"])
    input_files.extend(["e1/e1_85long", "e1/e1_150long", "e1/e1_250long",
                        "e1/e1_log85long", "e1/e1_log150long", "e1/e1_log250long"])
    input_files.extend(["e2/e2_85long", "e2/e2_150long", "e2/e2_250long",
                        "e2/e2_log85long", "e2/e2_log150long", "e2/e2_log250long"])
    input_files.extend(["e05/e05_85long", "e05/e05_150long", "e05/e05_250long",
                        "e05/e05_log85long", "e05/e05_log150long", "e05/e05_log250long"])

    hists = {}

    for f in input_files:

        # print f

        var_name = f.split("/")[1]

        hists[var_name] = []

        input_file = input_dir + f + ".dat"

        hists[var_name] = [[], [], [], []]

        with open(input_file) as infile:

            line_number = 0

            for line in infile:
                line_components = line.split()

                hists[var_name][0].append(float(line_components[0]))

                hists[var_name][1].append(float(line_components[2]))
                hists[var_name][2].append(float(line_components[1]))
                hists[var_name][3].append(float(line_components[3]))

    return hists


def compile_sources(parsed_hists):

    data_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[
        0], parsed_hists[1], parsed_hists[2], parsed_hists[3]

    theory_hists = parse_theory_file()

    return [data_hists, theory_hists, pythia_hists, herwig_hists, sherpa_hists]


def compile_data_and_pythia(all_hists, variables):

    data_hists, pythia_hists = all_hists[0], all_hists[1]

    data, pythia = [], []
    for var in variables:
        data.append(data_hists[var])
        pythia.append(pythia_hists[var])

    compilation = []
    for i in range(len(data[0])):
        temp = [data[0][i], data[1][i], pythia[0][i], pythia[1][i]]
        compilation.append(temp)

    return compilation


def compile_hists(var, parsed_hists, x_scale='linear'):

    compilation = []

    data_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[
        0], parsed_hists[2], parsed_hists[3], parsed_hists[4]

    # print data_hists, data_hists, var

    max_index = len(data_hists[var])

    for i in range(max_index):

        sub_list = [data_hists[var][i], pythia_hists[var][
            i], herwig_hists[var][i], sherpa_hists[var][i]]
        # sub_list = [ data_hists[var][i] ]

        compilation.append(sub_list)

    return compilation


def compile_hists_with_theory(var, parsed_hists, x_scale='linear'):

    compilation = []

    data_hists, theory_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[
        0], parsed_hists[1], parsed_hists[2], parsed_hists[3], parsed_hists[4]
    max_index = len(data_hists[var])

    for i in range(max_index):

        theory_var = var

        if "track" in theory_var:
            theory_var = theory_var.replace("track_", "")

        # Get the correct variable name to use for theory.
        theory_var = theory_var.split("_")[0]

        if data_hists[var][i].conditions()[1][1][1] == None:
            # Open pT boundaries. Use long.
            # Don't hardcode position of pT condition. Find the pT condition
            # using "" in.
            theory_var += "_" + \
                str(data_hists[var][i].conditions()[1][1][0]) + "long"
        else:
            # Closed pT boundaries.
            # Don't hardcode position of pT condition. Find the pT condition
            # using "" in.
            theory_var += "_" + str(data_hists[var][i].conditions()[1][1][0])

        if x_scale == "log":

            theory_var = theory_var.split(
                "_")[0] + "_log" + theory_var.split("_")[1]

            # theory_var += "log"

        sub_list = [data_hists[var][i], theory_hists[theory_var], pythia_hists[
            var][i], herwig_hists[var][i], sherpa_hists[var][i]]

        compilation.append(sub_list)

    return compilation


default_dir = "plots/Version 6/"


start = time.time()


parsed_linear = parse.load_root_files_to_hist(log=False)
parsed_hists = compile_sources(parsed_linear)

parsed_log = parse.load_root_files_to_hist(log=True)
parsed_log_hists = compile_sources(parsed_log)


# parsed_pfc = pfc_parse.load_pfc_root_files_to_hist()
# parsed_pfc_hists = compile_sources(parsed_pfc)


end = time.time()


print "Finished parsing all files in {} seconds. Now plotting them!".format(end - start)


start = time.time()


# create_multi_page_plot(filename=default_dir + "pfc_pT.pdf", hists=compile_hists('pfc_pT', parsed_pfc_hists))


# create_multi_page_plot(filename=default_dir + "pfc_eta.pdf",
#                        hists=compile_hists('pfc_eta', parsed_pfc_hists))


create_multi_page_plot(filename=default_dir + "hardest_jet_pT_all_linear.pdf",
                       hists=compile_hists('hardest_pT', parsed_hists))
# create_data_only_plot(filename=default_dir + "hardest_jet_pT_jec_all_linear.pdf", hists=[[parsed_linear[0]['hardest_pT'][i], parsed_linear[0]['uncor_hardest_pT'][i]] for i in range(len(parsed_linear[0]['uncor_hardest_pT']))], labels=[
#     "Jet Energy Corrected", "Jet Energy Uncorrected"], types=["error", "error"], colors=["black", "orange"], line_styles=[1, 1], ratio_to_label="Ratio to\nCorrected", ratio_to_index=0, text_outside_the_frame=True)


# create_multi_page_plot(filename=default_dir + "hardest_jet_phi_all_linear.pdf",
#                        hists=compile_hists('hardest_phi', parsed_hists))
# create_multi_page_plot(filename=default_dir + "hardest_jet_eta_all_linear.pdf",
#                        hists=compile_hists('hardest_eta', parsed_hists))


# create_multi_page_plot(filename=default_dir + "basic_sub_mul_all_linear.pdf",
#                        hists=compile_hists('mul_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "basic_sub_mul_track_linear.pdf",
#                        hists=compile_hists('track_mul_pre_SD', parsed_hists))


# create_data_only_plot(filename=default_dir + "basic_sub_mul_softdrop_all_linear.pdf", hists=compile_data_and_pythia([parsed_hists[0], parsed_hists[2]], variables=['mul_pre_SD', 'mul_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)
# create_data_only_plot(filename=default_dir + "basic_sub_mul_softdrop_track_linear.pdf", hists=compile_data_and_pythia([parsed_hists[0], parsed_hists[2]], variables=['track_mul_pre_SD', 'track_mul_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)


# create_multi_page_plot(filename=default_dir + "basic_sub_mass_all_linear.pdf",
#                        hists=compile_hists('mass_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "basic_sub_mass_track_linear.pdf",
#                        hists=compile_hists('track_mass_pre_SD', parsed_hists))

# create_data_only_plot(filename=default_dir + "basic_sub_mass_softdrop_all_linear.pdf", hists=compile_data_and_pythia([parsed_hists[0], parsed_hists[2]], variables=['mass_pre_SD', 'mass_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)
# create_data_only_plot(filename=default_dir + "basic_sub_mass_softdrop_track_linear.pdf", hists=compile_data_and_pythia([parsed_hists[0], parsed_hists[2]], variables=['track_mass_pre_SD', 'track_mass_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)


# create_multi_page_plot(filename=default_dir + "basic_sub_pTD_all_log.pdf",
#                        hists=compile_hists('pT_D_pre_SD', parsed_log_hists, x_scale='log'), x_scale='log')
# create_multi_page_plot(filename=default_dir + "basic_sub_pTD_track_log.pdf",
#                        hists=compile_hists('track_pT_D_pre_SD', parsed_log_hists,
#                                            x_scale='log'), x_scale='log')

# create_data_only_plot(filename=default_dir + "basic_sub_pTD_softdrop_all_log.pdf", hists=compile_data_and_pythia([parsed_log_hists[0], parsed_log_hists[2]], variables=['pT_D_pre_SD', 'pT_D_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)
# create_data_only_plot(filename=default_dir + "basic_sub_pTD_softdrop_track_log.pdf", hists=compile_data_and_pythia([parsed_log_hists[0], parsed_log_hists[2]], variables=['track_pT_D_pre_SD', 'track_pT_D_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)


# create_multi_page_plot(filename=default_dir + "basic_sub_lha_all_log.pdf",
#                        hists=compile_hists('LHA_pre_SD', parsed_log_hists), x_scale='log')
# create_multi_page_plot(filename=default_dir + "basic_sub_lha_track_log.pdf",
#                        hists=compile_hists(
#                            'track_LHA_pre_SD', parsed_log_hists),
#                        x_scale='log')

# create_data_only_plot(filename=default_dir + "basic_sub_lha_softdrop_all_log.pdf", hists=compile_data_and_pythia([parsed_log_hists[0], parsed_log_hists[2]], variables=['LHA_pre_SD', 'LHA_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)
# create_data_only_plot(filename=default_dir + "basic_sub_lha_softdrop_track_log.pdf", hists=compile_data_and_pythia([parsed_log_hists[0], parsed_log_hists[2]], variables=['track_LHA_pre_SD', 'track_LHA_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)


# create_multi_page_plot(filename=default_dir + "basic_sub_width_all_log.pdf",
#                        hists=compile_hists('width_pre_SD', parsed_log_hists), x_scale='log')
# create_multi_page_plot(filename=default_dir + "basic_sub_width_track_log.pdf",
#                        hists=compile_hists(
#                            'track_width_pre_SD', parsed_log_hists),
#                        x_scale='log')

# create_data_only_plot(filename=default_dir + "basic_sub_width_softdrop_all_log.pdf", hists=compile_data_and_pythia([parsed_log_hists[0], parsed_log_hists[2]], variables=['width_pre_SD', 'width_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)
# create_data_only_plot(filename=default_dir + "basic_sub_width_softdrop_track_log.pdf", hists=compile_data_and_pythia([parsed_log_hists[0], parsed_log_hists[2]], variables=['track_width_pre_SD', 'track_width_post_SD']), labels=[
#                       "Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)


# create_multi_page_plot(filename=default_dir + "basic_sub_thrust_all_log.pdf",
#                        hists=compile_hists('thrust_pre_SD', parsed_log_hists), x_scale='log')
# create_multi_page_plot(filename=default_dir + "basic_sub_thrust_track_log.pdf",
#                        hists=compile_hists(
#                            'track_thrust_pre_SD', parsed_log_hists),
#                        x_scale='log')

# create_data_only_plot(filename=default_dir + "basic_sub_thrust_softdrop_all_log.pdf", hists=compile_data_and_pythia([parsed_log_hists[0], parsed_log_hists[2]], variables=['thrust_pre_SD', 'thrust_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=[
#                       "error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)
# create_data_only_plot(filename=default_dir + "basic_sub_thrust_softdrop_track_log.pdf", hists=compile_data_and_pythia([parsed_log_hists[0], parsed_log_hists[2]], variables=['track_thrust_pre_SD', 'track_thrust_post_SD']), labels=[
#                       "Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index={0: 2, 1: 3, 2: 2, 3: 3}, text_outside_the_frame=True)


# create_multi_page_plot(filename=default_dir + "softdrop_frac_pT_loss_log.pdf", hists=compile_hists(
#     'frac_pT_loss', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')


# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/zg/zg_05.pdf",
#                        hists=compile_hists_with_theory('zg_05', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/zg/zg_05.pdf",
#                        hists=compile_hists_with_theory('track_zg_05', parsed_hists),
#                        theory=True)

# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/rg/rg_05.pdf",
#                        hists=compile_hists_with_theory('rg_05', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/rg/rg_05.pdf",
#                        hists=compile_hists_with_theory('track_rg_05', parsed_hists),
#                        theory=True)


# create_multi_page_plot(filename=default_dir + "theta_g/log/all/zg/zg_05.pdf",
#                        hists=compile_hists_with_theory('zg_05', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/zg/zg_05.pdf",
#                        hists=compile_hists_with_theory('track_zg_05', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/rg/rg_05.pdf",
#                        hists=compile_hists_with_theory('rg_05', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/rg/rg_05.pdf",
#                        hists=compile_hists_with_theory('track_rg_05', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e1/e1_05.pdf",
#                        hists=compile_hists_with_theory('e1_05', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e1/e1_05.pdf",
#                        hists=compile_hists_with_theory('track_e1_05', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e2/e2_05.pdf",
#                        hists=compile_hists_with_theory('e2_05', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e2/e2_05.pdf",
#                        hists=compile_hists_with_theory('track_e2_05', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e05/e05_05.pdf",
#                        hists=compile_hists_with_theory('e05_05', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e05/e05_05.pdf",
#                        hists=compile_hists_with_theory('track_e05_05', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')


# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/zg/zg_10.pdf",
#                        hists=compile_hists_with_theory('zg_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/zg/zg_10.pdf",
#                        hists=compile_hists_with_theory(
#                            'track_zg_10', parsed_hists),
#                        theory=True)

# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/rg/rg_10.pdf",
#                        hists=compile_hists_with_theory('rg_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/rg/rg_10.pdf",
#                        hists=compile_hists_with_theory(
#                            'track_rg_10', parsed_hists),
#                        theory=True)


# create_multi_page_plot(filename=default_dir + "theta_g/log/all/zg/zg_10.pdf",
#                        hists=compile_hists_with_theory('zg_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/zg/zg_10.pdf",
#                        hists=compile_hists_with_theory('track_zg_10', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/rg/rg_10.pdf",
#                        hists=compile_hists_with_theory('rg_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/rg/rg_10.pdf",
#                        hists=compile_hists_with_theory('track_rg_10', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e1/e1_10.pdf",
#                        hists=compile_hists_with_theory('e1_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e1/e1_10.pdf",
#                        hists=compile_hists_with_theory('track_e1_10', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e2/e2_10.pdf",
#                        hists=compile_hists_with_theory('e2_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e2/e2_10.pdf",
#                        hists=compile_hists_with_theory('track_e2_10', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e05/e05_10.pdf",
#                        hists=compile_hists_with_theory('e05_10', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e05/e05_10.pdf",
#                        hists=compile_hists_with_theory('track_e05_10', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')


# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/zg/zg_20.pdf",
#                        hists=compile_hists_with_theory('zg_20', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/zg/zg_20.pdf",
#                        hists=compile_hists_with_theory(
#                            'track_zg_20', parsed_hists),
#                        theory=True)


# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/rg/rg_20.pdf",
#                        hists=compile_hists_with_theory('rg_20', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/rg/rg_20.pdf",
#                        hists=compile_hists_with_theory(
#                            'track_rg_20', parsed_hists),
#                        theory=True)


# create_multi_page_plot(filename=default_dir + "theta_g/log/all/zg/zg_20.pdf",
#                        hists=compile_hists_with_theory('zg_20', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/zg/zg_20.pdf",
#                        hists=compile_hists_with_theory('track_zg_20', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/rg/rg_20.pdf",
#                        hists=compile_hists_with_theory('rg_20', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/rg/rg_20.pdf",
#                        hists=compile_hists_with_theory('track_rg_20', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e1/e1_20.pdf",
#                        hists=compile_hists_with_theory('e1_20', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e1/e1_20.pdf",
#                        hists=compile_hists_with_theory('track_e1_20', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e2/e2_20.pdf",
#                        hists=compile_hists_with_theory('e2_20', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e2/e2_20.pdf",
#                        hists=compile_hists_with_theory('track_e2_20', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')

# create_multi_page_plot(filename=default_dir + "theta_g/log/all/e05/e05_20.pdf",
#                        hists=compile_hists_with_theory('e05_20', parsed_log_hists, x_scale='log'), theory=True, x_scale='log')
# create_multi_page_plot(filename=default_dir + "theta_g/log/track/e05/e05_20.pdf",
#                        hists=compile_hists_with_theory('track_e05_20', parsed_log_hists,
#                                                        x_scale='log'), theory=True, x_scale='log')


##########################################################################


# create_data_only_plot(filename=default_dir + "jec.pdf", hists=[ [x] for x in parsed_linear[0]['jec'] ], labels=["CMS 2010 Open Data"], types=["error"], colors=["black"], line_styles=[[]], ratio_plot=False)
# create_data_only_plot(filename=default_dir + "area.pdf", hists=[ [x] for x in parsed_linear[0]['hardest_area'] ], labels=["CMS 2010 Open Data"], types=["error"], colors=["black"], line_styles=[[]], ratio_plot=False)


##########################################################################


# create_multi_page_plot(filename=default_dir + "hardest_jet_pT_all_log.pdf", hists=compile_hists('hardest_pT', parsed_log_hists, x_scale='log'), x_scale='log')
# create_data_only_plot(filename=default_dir + "hardest_jet_pT_jec_all_log.pdf", hists=[ [ parsed_log[0]['hardest_pT'][i], parsed_log[0]['uncor_hardest_pT'][i] ] for i in range(len(parsed_log[0]['uncor_hardest_pT'])) ], labels=["Jet Energy Corrected", "Jet Energy Uncorrected"], types=["error", "error"], colors=["black", "orange"], line_styles=[1, 1], ratio_to_label="Ratio to\nCorrected", ratio_to_index=0)


# create_multi_page_plot(filename=default_dir + "basic_sub_pTD_all_linear.pdf", hists=compile_hists('pT_D_pre_SD', parsed_hists))# create_multi_page_plot(filename=default_dir + "basic_sub_pTD_track_linear.pdf", hists=compile_hists('track_pT_D_pre_SD', parsed_hists))
#


#
# create_data_only_plot(filename=default_dir + "all_track_zg_10.pdf", hists=compile_data_and_pythia([ parsed_log_hists[0], parsed_log_hists[2] ], variables=['zg_10', 'track_zg_10']), labels=["Everything", "Track", "Pythia 8.215", "Pythia 8.215 (Track)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [50, 30], [50, 30]], ratio_to_label="Ratio to\nPythia", ratio_to_index=2)

# create_multi_page_plot(filename=default_dir + "softdrop_frac_pT_loss_linear.pdf", hists=compile_hists('frac_pT_loss', parsed_hists))
#
#
# create_multi_page_plot(filename=default_dir + "basic_sub_lha_all_linear.pdf", hists=compile_hists('LHA_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "basic_sub_lha_track_linear.pdf", hists=compile_hists('track_LHA_pre_SD', parsed_hists))

# create_multi_page_plot(filename=default_dir + "basic_sub_width_all_linear.pdf", hists=compile_hists('width_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "basic_sub_width_track_linear.pdf", hists=compile_hists('track_width_pre_SD', parsed_hists))

# create_multi_page_plot(filename=default_dir + "basic_sub_thrust_all_linear.pdf", hists=compile_hists('thrust_pre_SD', parsed_hists))
# create_multi_page_plot(filename=default_dir + "basic_sub_thrust_track_linear.pdf", hists=compile_hists('track_thrust_pre_SD', parsed_hists))

# create_multi_page_plot(filename=default_dir + "softkill_pT_loss_linear.pdf", hists=compile_hists('softkill_pT_loss', parsed_hists))
# create_multi_page_plot(filename=default_dir + "softkill_pT_loss_log.pdf", hists=compile_hists('softkill_pT_loss', parsed_log_hists, x_scale='log'), theory=False, x_scale='log')


#
# create_data_only_plot(filename=default_dir + "hardest_jet_pT_softdrop_all_linear.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['hardest_pT', 'pT_after_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index=2)
#
# create_data_only_plot(filename=default_dir + "basic_sub_lha_softdrop_track_linear.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['track_LHA_pre_SD', 'track_LHA_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index=2)
# create_data_only_plot(filename=default_dir + "basic_sub_width_softdrop_track_linear.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['track_width_pre_SD', 'track_width_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index=2)create_data_only_plot(filename=default_dir + "basic_sub_thrust_softdrop_track_linear.pdf", hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ], variables=['track_thrust_pre_SD', 'track_thrust_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error", "hist", "hist"], colors=["black", "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia", ratio_to_index=2)
# create_data_only_plot(filename=default_dir +
# "basic_sub_lha_softdrop_all_linear.pdf", hists=compile_data_and_pythia([
# parsed_hists[0], parsed_hists[2] ], variables=['LHA_pre_SD',
# 'LHA_post_SD']), labels=["Before Soft Drop", "After Soft Drop", "Pythia
# 8.215 (Before)", "Pythia 8.215 (After)"], types=["error", "error",
# "hist", "hist"], colors=["black", "red", "black", "red"],
# line_styles=[[], [], [7, 7], [7, 7]], ratio_to_label="Ratio to\nPythia",
# ratio_to_index=2)create_data_only_plot(filename=default_dir +
# "basic_sub_width_softdrop_all_linear.pdf",
# hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ],
# variables=['width_pre_SD', 'width_post_SD']), labels=["Before Soft
# Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215
# (After)"], types=["error", "error", "hist", "hist"], colors=["black",
# "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]],
# ratio_to_label="Ratio to\nPythia",
# ratio_to_index=2)create_data_only_plot(filename=default_dir +
# "basic_sub_thrust_softdrop_all_linear.pdf",
# hists=compile_data_and_pythia([ parsed_hists[0], parsed_hists[2] ],
# variables=['thrust_pre_SD', 'thrust_post_SD']), labels=["Before Soft
# Drop", "After Soft Drop", "Pythia 8.215 (Before)", "Pythia 8.215
# (After)"], types=["error", "error", "hist", "hist"], colors=["black",
# "red", "black", "red"], line_styles=[[], [], [7, 7], [7, 7]],
# ratio_to_label="Ratio to\nPythia", ratio_to_index=2)


# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e05/e05_10.pdf", hists=compile_hists_with_theory('e05_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e05/e05_10.pdf", hists=compile_hists_with_theory('track_e05_10', parsed_hists), theory=True)

# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e1/e1_10.pdf", hists=compile_hists_with_theory('e1_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e1/e1_10.pdf", hists=compile_hists_with_theory('track_e1_10', parsed_hists), theory=True)

# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e2/e2_10.pdf", hists=compile_hists_with_theory('e2_10', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e2/e2_10.pdf", hists=compile_hists_with_theory('track_e2_10', parsed_hists), theory=True)


# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e1/e1_05.pdf", hists=compile_hists_with_theory('e1_05', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e2/e2_05.pdf", hists=compile_hists_with_theory('e2_05', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e05/e05_05.pdf", hists=compile_hists_with_theory('e05_05', parsed_hists), theory=True)

# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e1/e1_05.pdf", hists=compile_hists_with_theory('track_e1_05', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e2/e2_05.pdf", hists=compile_hists_with_theory('track_e2_05', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e05/e05_05.pdf", hists=compile_hists_with_theory('track_e05_05', parsed_hists), theory=True)

# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e1/e1_20.pdf", hists=compile_hists_with_theory('e1_20', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e2/e2_20.pdf", hists=compile_hists_with_theory('e2_20', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/all/e05/e05_20.pdf", hists=compile_hists_with_theory('e05_20', parsed_hists), theory=True)

# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e1/e1_20.pdf", hists=compile_hists_with_theory('track_e1_20', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e2/e2_20.pdf", hists=compile_hists_with_theory('track_e2_20', parsed_hists), theory=True)
# create_multi_page_plot(filename=default_dir + "theta_g/linear/track/e05/e05_20.pdf", hists=compile_hists_with_theory('track_e05_20', parsed_hists), theory=True)


end = time.time()
print "Finished all plotting in {} seconds.".format(end - start)
