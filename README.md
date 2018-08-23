# PlotProcessedData
Repository to save code for making analysis plots on the new data

Current workflow is:

Processed files are downloaded in some format to /data/snoplus/OfficialProcessing

If ntuple files are available, use these for fastest plotting.

Otherwise the code does support plotting from ratds files and will also allow plotting of PMT related variables.

If desired, conversion of ratds files into ntuples can be done using scripts/root_converter.py.
    
If only zdabs are available, processing of these into both ratds and root files can be run on the batch using process.py script.
    
The plotting takes either a list of ratds or ntuple files, one file per line (as e.g. result of `ls /data/snoplus/OfficialProcessing/processed_6.2.3/*_ntuple.root &> filelist_ntuple.dat`. Code assumes default is ratds files and ntuples will have string "ntuple" in the name of the filelist.

Filelists for plotting can also be pruned according to some good run list using prune_filelist.py

The code is built with `make` and you must also run source scripts/setup_libs.sh to link the created library correctly. The code is then run with e.g. `./bin/makeRunPlots --filelist=filelists/filelist_rat625_waterFit_ntuple.dat --directory=plots/`

Plots can be uploaded to a webpage after creation. The most useful format requires the use of a script called scripts/generate_jsons.py which creates some json files linking together plots from each run on the webpage. The index.html file is then generated running `python gallery/gallery.py plots/`.
