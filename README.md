# CMSOpenDataAnalyzer

This repository has code to analyze CMS Open Data data produced using [https://github.com/tripatheea/CMSOpenDataProducer][1].

### How to Install/Use
 - Install [FastJet][2]. Note the path to FastJet config. 
 - Open ./Makefile and paste the path to FastJet config to the **PATH\_TO\_FASTJET** variable.
 - Run [CMSOpenDataProducer][3]. It will output a `.dat` file. Copy it to the `./data/' directory.
 - Compile everything with `make`.
 - Move to the 'bin' directory and run the analyzer.

```cd bin```

```./analysis```

 - The code will output a `DAT` file `antikt_multiplicities.dat`.
 - To make the plots, move to the directory `./root'.
 - Run root and compile/execute the file `plots.cc`

```root```

```.x plots.cc+``` 

  

[1]: https://github.com/tripatheea/CMSOpenDataProducer
  [2]: http://www.fastjet.fr/
  [3]: https://github.com/tripatheea/CMSOpenDataProducer