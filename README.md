# MIT Open Data Analyzer

This repository has code to analyze MOD (MIT Open Data) files produced using [MODProducer](https://github.com/tripatheea/MODProducer "MODProducer").

## Workflow

We adopt the following workflow for extracting MOD files out of the given AOD files.

1.  Download all the ROOT files and arrange them in the same directory structure as they live on the CMS server.

2. Create a registry that maps each event and run number to a certain ROOT file. This is done so that things can be done one file at a time as we're dealing with multiple TeraBytes of data here and it's a pain to have to do everything at once. 

3. Run the Producer on those AOD files. This reads the download directory and processes only the files in there. This produces N MOD files. 

4. Filter those N MOD files to get only those files for which the correct trigger fired. This process is called Skimming in this workflow. This will produce M <= N MOD files. For a certain AOD file, if none of the events in there have the correct trigger fired, a corresponding skimmed MOD file will not be written. That's why M might be less than N.

5. Read in those "skimmed" M <= N output files one by one and calculate stuff to produce a single DAT file. 

6. Produce plots using the DAT file produced in step (5).

This repository is concerned with steps (4) to (6) only. Steps (1) to (3) are carried out by the [MODProducer](https://github.com/tripatheea/MODProducer/ "MODProducer") package.

## Usage Instruction

### Preparation

 - Install [FastJet](http://www.fastjet.fr/ "FastJet"). Note the path to FastJet config. 
 
 - Open `./Makefile` and paste the path to FastJet config to the **PATH\_TO\_FASTJET** variable on line 5.
 
 -  Note the path to the directory that contains the MOD files you've produced from the [MODProducer](https://github.com/tripatheea/MODProducer/ "MODProducer") package.
 
 - Compile everything with `make`.

### Workflow Instructions
 
 - First, we run the skimmer to filter those N MOD files you produced to get only those files for which the correct trigger fired. This is accomplished by the Python script `utilities/skim.py`. This script takes a single argument,  a path to the directory that contains all the MOD files.
   
   ```
   python ./utilities/skim.py /media/aashish/opendata/MIT_CMS/eos/opendata/cms/Run2010B/Jet/MOD/Apr21ReReco-v1/0000/
   ```
 
    This step maintains the same directory structure as the input directory except MOD replaced with SKIM. That's why you do not need to enter an output directory. It will also output an error log in the same directory.


 - Next, we run the analyzer. We use the Python script `utilities/analyze.py` for data analysis. This script will run the executable `bin/analyze` M times for M "skimmed" MOD files. This script takes two arguments:
 
   1. path to the directory that holds the skimmed files.
   2. path to a filename to write the analyzed data into. 

     ```
     python ./utilities/analyze.py /media/aashish/opendata/MIT_CMS/eos/opendata/cms/Run2010B/Jet/SKIM/Apr21ReReco-v1/0000/ ~/analyzed_data.dat
     ```

  - Finally, we use the output file to produce plots using the Python script `python/plots.py`. This script uses [matplotlib](http://matplotlib.org/ "matplotlib") to produce the plots but feel free to explore the output file and use any plotting routine you want. If you run `python/plots.py`, it will produce plots as PDFs in a directory `plots`. It takes one argument- path to the output file produced in the previous step.
  
    ```
    python ./python/plots.py ~/analyzed_data.dat
    ```

## TODO
- Talk to Sal about what the rapidity range is for AK5.
- Figure out migration in charged to neutral.


### Low Priority TODOs
- Look into gzip thingy for size issues.


