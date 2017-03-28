# PlotProcessedData
Repository to save code for making analysis plots on the new data

Current workflow is:

Processed files are downloaded in some format to /data/snoplus/OfficialProcessing
    
    - if ntuple files are available, use these for fastest plotting.
    
    - Otherwise the code does support plotting from ratds files.
    
    - Conversion of ratds files into ntuples can be done using root_converter.py.
    
    - If only zdabs are available, processing of these into both ratds and root files can be run on the batch using process.py script.
    
The plotting takes either a list of ratds or ntuple files, one file per line (as e.g. result of `ls /data/snoplus/OfficialProcessing/processed_6.2.3/*_ntuple.root &> filelist_ntuple.dat`. Code assumes default is ratds files and ntuples will have string "ntuple" in the name of the filelist.

Filelists for plotting can also be pruned according to some good run list using prune_filelist.py

Plotting is then run with e.g. `root -l -q -b 'makeRunPlots.C+("filelists/filelist_goodruns_ntuple.dat")`

Plots can be uploaded to a webpage after creation. The most useful format requires the use of a script called generate_jsons.py which creates some json files linking together plots from each run on the webpage. The index.html file is then generated running `python gallery/gallery.py plots/`.
