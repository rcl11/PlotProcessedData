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
    'nhitsz' : 'Number of hits vs z',
    'nhitstemp' : 'Number of hits vs temperature',
    'nhits_log' : 'Number of hits log scale',
    'posR' : 'R position',
    'posR3' : 'R cubed position',
    'posRz' : 'R vs z position',
    'posrhoz' : 'rho vs z position',
    'posxy' : 'x vs y position',
    'posx' : 'x position',
    'posy' : 'y position',
    'posz' : 'z position',
    'totalQ' : 'Total charge',
    'duration' : 'Run duration',
    'trigger' : 'Trigger word',
    'dataclean' : 'Data clean',
    'temp': 'cavity temp',
    'fitValid' : 'Fit validity',
    'rpmt' : 'PMT R position',
    'tpmt' : 'PMT time',
    'xpmt' : 'PMT x position',
    'ypmt' : 'PMT y position',
    'zpmt' : 'PMT z position',
    'itr' : 'ITR',
    'udotr' : 'U dot R',
    'udotrz' : 'U dot R vs z',
    'beta14': 'beta14',
    'energy': 'energy',
    'errposx': 'x err vs x',
    'errposy': 'y err vs y',
    'errposz': 'z err vs z',
    'errposxnhits': 'x err vs nhits',
    'errposxitr': 'x err vs ITR',
    'errposynhits': 'y err vs nhits',
    'errposyitr': 'y err vs ITR',
    'errposznhits': 'z err vs nhits',
    'errposzitr': 'z err vs ITR',
    'errtimex': 'time err vs x',
    'errtimey': 'time err vs y',
    'errtimez': 'time err vs z',
    'time': 'fitted time',
    'timeposx': 'fitted time vs x',
    'timeposy': 'fitted time vs y',
    'timeposz': 'fitted time vs z',
    'errenergy': 'energy err vs energy',
}

maintenance_runs = ['14713','14714','14715','15091']
experimental_runs = ['15143','15144','15145','15153','15154','15155','15156','15157','15158','15159']
physics_runs = ['15060','15061','15459','16071','16072','16073','15664','15665','15666']
runs1hr = ['15060','15061','15459']

#Extract known information from filename
for filename in glob.glob(options.dirname+"/"+"*.png"):
    print filename
    data = {}
    json_name = filename.replace(".png",".json")
    json_file = open(json_name, "w")
    p = re.compile("r([0-9]*)_")
    if "_vs_" in filename: 
        data['plot type'] = "vs Run"
        data['trigger type'] = "Any trigger"
        json.dump(data,json_file)
        json_file.close()
        continue
        
    #run number    
    if not "_to_" in filename and p.findall(filename):
        run_number = p.findall(filename)[0]
    if run_number in maintenance_runs :
        data['run type'] = 'maintenance'
    elif run_number in experimental_runs:
        data['run type'] = 'experimental'
    elif run_number in physics_runs:
        data['run type'] = 'physics'
    
    if run_number in runs1hr:
        data['run duration'] = '1 hour'
        
    
 #   j = re.compile("_s([0-9]*)")
 #   if j.findall(filename):
 #     subrun_number = j.findall(filename)[0]

    t = re.compile("[0-9*]/(.*)_r[0-9]*_to_")
    q = re.compile("[0-9*]/(.*)_r")
    #plot type    
    if "_to_" in filename:
        data['run number'] = "sum of runs"
        plot_type = t.findall(filename)[0]
    elif not "_vs_" in filename:    
      if q.findall(filename):
        plot_type = q.findall(filename)[0]
        data['run number'] = run_number

    #trigger type    
    if "goodfit" in plot_type:
        data['fit'] = 'good'
        plot_type = plot_type.replace("goodfit","")
    else: data['fit'] = 'all'    
    
    trig_type = "Any trigger"
    if "nhitsz" in plot_type:
        plot_type = 'nhitsz'
    if "nhitstemp" in plot_type:
        plot_type = 'nhitstemp'
    elif "errposxnhits" in plot_type:
        plot_type = 'errposxnhits'
    elif "errposynhits" in plot_type:
        plot_type = 'errposynhits'
    elif "errposznhits" in plot_type:
        plot_type = 'errposznhits'
    elif "nhits" in plot_type:
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
