#!/usr/bin/python

import os
import glob

with open('/data/snoplus/lane/ProcessedData/PlotProcessedData/good_run_list.txt') as f:
    interesting_runs = f.read().splitlines()

for run_num in interesting_runs:
    for subrunfile in glob.glob("/data/snoplus/OfficialProcessing/raw/good_runs_all/SNOP*"+str(run_num)+"*.l2.zdab"):
        
        print subrunfile
        subrun_num = subrunfile.split(str(run_num)+"_",1)[1] 
        subrun_num = subrun_num.split(".l2",1)[0] 

        script = open("process_"+str(run_num)+"_"+str(subrun_num)+".sh", "w")

        
        script.write("#!/bin/sh\n")
       
        script.write("cd /data/lane/RAT/rat\n")
        script.write('source /data/lane/RAT/env_rat-6.2.0.sh\n')
        script.write('rat -b postgres://snoplus:dontestopmenow@pgsql.snopl.us:5400/ratdb -i /data/snoplus/OfficialProcessing/raw/good_runs_all/SNOP_00000{0}_{1}.l2.zdab mac/processing/partial_water/first_pass_data_cleaning.mac\n'.format(run_num, subrun_num)) 
        script.write('rat -b postgres://snoplus:dontestopmenow@pgsql.snopl.us:5400/ratdb -i /data/snoplus/OfficialProcessing/raw/good_runs_all/SNOP_00000{0}_{1}.l2.zdab -o /data/snoplus/OfficialProcessing/processed_pre6.2.5_localprocess/Processing_r{0}_s{1}_p000.root mac/processing/partial_water/processing.mac\n'.format(run_num, subrun_num))
        script.close()

        os.system("qsub -cwd -l h_rss=4G,h_vmem=4G -q *SL6 process_"+str(run_num)+"_"+str(subrun_num)+".sh")
