from __future__ import division
from pyPdf import PdfFileWriter, PdfFileReader
import os


def reshuffle_big_5(output_directory, input_file, filename):
    output_pdf = PdfFileWriter()

    with open(input_file, 'rb') as readfile:

        input_pdf = PdfFileReader(readfile)
        total_pages = input_pdf.getNumPages()

        if total_pages == 21:

            for j in [1, 0, 2]:

                output_pdf.addPage(input_pdf.getPage(5 + 7 * j))
                output_pdf.addPage(input_pdf.getPage(4 + 7 * j))
                output_pdf.addPage(input_pdf.getPage(6 + 7 * j))

                for i in range(0, 4):
                    output_pdf.addPage(input_pdf.getPage(i + 7 * j))

            with open(output_directory + filename, "wb") as writefile:
                output_pdf.write(writefile)


def reshuffle(output_directory, input_file, filename):
    output_pdf = PdfFileWriter()

    with open(input_file, 'rb') as readfile:

        input_pdf = PdfFileReader(readfile)
        total_pages = input_pdf.getNumPages()

        if total_pages == 7:

            print filename,

            if "hardest_jet_phi_all_linear" in filename or "hardest_jet_eta_all_linear" in filename or "hardest_jet_pT_all_linear" in filename or "hardest_jet_pT_jec_all_linear" in filename or "area" in filename or "jec" in filename or "pfc_neutral_0_100_pT" in filename or "pfc_charged_0_100_pT" in filename or "pfc_neutral_0_5_pT" in filename or "pfc_charged_0_5_pT" in filename:
                output_pdf.addPage(input_pdf.getPage(4))
                output_pdf.addPage(input_pdf.getPage(5))
                output_pdf.addPage(input_pdf.getPage(6))
                print "85"
            else:
                output_pdf.addPage(input_pdf.getPage(5))
                output_pdf.addPage(input_pdf.getPage(4))
                output_pdf.addPage(input_pdf.getPage(6))
                print "150"

            for i in range(0, 4):
                output_pdf.addPage(input_pdf.getPage(i))

            with open(output_directory + filename, "wb") as writefile:
                output_pdf.write(writefile)


if __name__ == "__main__":

    output_directory = "/home/aashish/root/macros/MODAnalyzer/plots/Version 6/reshuffled/"
    input_directory = "/home/aashish/root/macros/MODAnalyzer/plots/Version 6/"

    # output_directory = "/home/aashish/root/macros/MODAnalyzer/plots/Version 6/theta_g/reshuffled/"
    # input_directory = "/home/aashish/root/macros/MODAnalyzer/plots/Version 6/theta_g/"

    for file in os.listdir(input_directory):
        if file.endswith(".pdf"):
            # print(os.path.join("/mydir", file))

            reshuffle(output_directory, os.path.join(input_directory, file), file)
            # reshuffle_big_5(output_directory, os.path.join(input_directory, file), file)

    # reshuffle(output_directory, )
