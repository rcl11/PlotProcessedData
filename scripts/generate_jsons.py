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


nice_labels = {
    'nhits' : 'Number of hits',
    'nhits_log' : 'Number of hits log scale',
    'posR' : 'R position',
    'posR3' : 'R cubed position',
    'posRz' : 'R vs z position',
    'posxy' : 'x vs y position',
    'posx' : 'x position',
    'posy' : 'y position',
    'posz' : 'z position',
    'rpmt' : 'PMT R position',
    'totalQ' : 'Total charge',
    'tpmt' : 'PMT time',
    'xpmt' : 'PMT x position',
    'ypmt' : 'PMT y position',
    'zpmt' : 'PMT z position',
    'fitValid' : 'Fit validity',
    'itr' : 'ITR',
}

maintenance_runs = ['14713','14714','14715','15091']
experimental_runs = ['15143','15144','15145','15153','15154','15155','15156','15157','15158','15159']

#Extract known information from filename
for filename in glob.glob(options.dirname+"/"+"*.png"):
    data = {}
    json_name = filename.replace(".png",".json")
    json_file = open(json_name, "w")
    p = re.compile("r([0-9]*)_s")
    if not p.findall(filename): 
        data['plot type'] = "vs Run"
        data['trigger type'] = "Any trigger"
        json.dump(data,json_file)
        json_file.close()
        continue
        
    #run number    
    run_number = p.findall(filename)[0]
    if run_number in maintenance_runs :
        data['run type'] = 'maintenance'
    elif run_number in experimental_runs:
        data['run type'] = 'experimental'
    j = re.compile("_s([0-9]*)")
    subrun_number = j.findall(filename)[0]

    q = re.compile("/(.*)_r")
    #plot type    
    plot_type = q.findall(filename)[0]
    data['run number'] = run_number+"_"+subrun_number

    #trigger type    
    trig_type = "Any trigger"
    if "nhits" in plot_type:
        if not "nhits_" in plot_type or "nhits_log" in plot_type:
            trig_type = "Any trigger"
        else:    
            trig_type = plot_type.split("nhits_",1)[1]
        if "log" in plot_type:
            plot_type = 'nhits_log'
            trig_type = trig_type.replace("_log","")
        else:
            plot_type = 'nhits'
    
    if "totalQ" in plot_type:
        if not "totalQ_" in plot_type:
            trig_type = "Any trigger"
        else:
            trig_type = plot_type.split("totalQ_",1)[1]
        plot_type = 'totalQ'
    data['trigger type'] = trig_type

    data['plot type'] = nice_labels[plot_type]

    json.dump(data,json_file)
    json_file.close()


#Also can add special info - certain type of run, whether or not marked as "good" etc
