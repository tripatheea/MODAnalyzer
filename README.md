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
- Look into gzip thingy for size issues.
- Talk to Sal about the jet areas being 0.

## TODO
- Get Pythia 8 up and running and generate simulated data (with simulated triggers)



[1]: https://github.com/tripatheea/MODProducer
[2]: http://www.fastjet.fr/
[3]: https://github.com/tripatheea/MODProducer
