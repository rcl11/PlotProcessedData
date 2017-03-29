#!/usr/bin/python

import json
import os
import glob
import re
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-d", "--directory", dest="dirname",
                  help="Directory where the plots are kept")

(options, args) = parser.parse_args()


#Extract known information from filename
for filename in glob.glob(options.dirname+"/"+"*.png"):
    data = {}
    json_name = filename.replace(".png",".json")
    json_file = open(json_name, "w")
    p = re.compile("r([0-9]*)_s")
    if not p.findall(filename): 
        data['plot type'] = "vs Run"
        json.dump(data,json_file)
        json_file.close()
        continue
        
    #run number    
    run_number = p.findall(filename)[0]
    j = re.compile("_s([0-9]*)")
    subrun_number = j.findall(filename)[0]

    q = re.compile("/(.*)_r")
    #plot type    
    plot_type = q.findall(filename)[0]

    data['run number'] = run_number+"_"+subrun_number
    data['plot type'] = plot_type

    json.dump(data,json_file)
    json_file.close()



#Also can add special info - certain type of run, whether or not marked as "good" etc
