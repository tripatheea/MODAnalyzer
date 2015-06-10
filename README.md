# MIT Open Data Analyzer

This repository has code to analyze CMS Open Data data produced using [https://github.com/tripatheea/MODProducer][1].

### How to Install/Use
 - Install [FastJet][2]. Note the path to FastJet config. 
 - Open `./Makefile` and paste the path to FastJet config to the **PATH\_TO\_FASTJET** variable.
 - Run [CMSOpenDataProducer][3]. It will output a `.dat` file. Copy it to the `./data/` directory.
 - Compile everything with `make`.
 - The analyzer takes two mandatory and one optional arguments- path to the input file, path to the output file and the number of events to process, respectively. If the third argument is absent, all events found in the input file will be processed.

   ``./bin/validate data/CMS_JetSample.dat data/CMS_JetSample_ validated.dat (57)``
   ``./bin/filter data/CMS_JetSample.dat data/CMS_JetSample_filtered.dat (29)``
   ``./bin/analyze data/CMS_JetSample.dat data/CMS_JetSample_analyzed.dat (33)``

 - The code will output a `DAT` file `./data/output.dat`.
 - To make the plots, move to the directory `./root`.
 - Run root and compile/execute the file `plots.cc`

  ``root ``

  ``.x plots.cc+`` 

## TODO
- Add jet energy correction factors to AK5/AK7
- Write a simple analysis that double checks that PFCandidates yield the same AK5/7 jets as CMS finds.
- Get Pythia 8 up and running and generate simulated data (with simulated triggers)
- Figure out if there is any way to get luminosity information.
- Find out the full list of triggers associated with the Jet primary dataset.  The list might be here:  https://fwyzard.web.cern.ch/fwyzard/hlt/2010/dataset 

## TODO (Plots)
- Create an AK5 hardest pT spectrum, color coded by which trigger is being used.
- Create an AK7 hardest pT spectrum, again color coded by trigger.
- Create a corrected AK5 spectrum, and compared to the published CMS data.
- Create an Ntilde distribution.
- Create a quasi-corrected Ntilde distribution (I need to teach you how to do this.)


[1]: https://github.com/tripatheea/MODProducer
[2]: http://www.fastjet.fr/
[3]: https://github.com/tripatheea/MODProducer
