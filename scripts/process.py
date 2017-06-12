#!/usr/bin/python

import os
import glob

from optparse import OptionParser


parser = OptionParser()
parser.add_option("-g", "--good_runs", dest="good_runs",
                  help="List of interesting runs for which we want to run processing")
parser.add_option("-i", "--input_dir", dest="input_dir",
                  help="Directory where the zdabs files are located")
parser.add_option("-o", "--output_dir", dest="output_dir",
                  help="Directory where the root files are to be stored")
parser.add_option("-r", "--rat_dir", dest="rat_dir",
                  default = "/data/lane/RAT/rat", help="Path of rat release")
parser.add_option("-e", "--rat_env", dest="rat_env",
                  default = "/data/lane/RAT/env_rat-6.2.0.sh", help="Environment file for rat")
parser.add_option("-s", "--submit", default = False, action = "store_true", dest="submit",
                  help="Actually submit the batch jobs! Please do a dry run first without this!")

(options, args) = parser.parse_args()

good_run = str(options.good_runs)
input_dir = str(options.input_dir) + "/"
output_dir = str(options.output_dir) + "/"
rat_dir = str(options.rat_dir)
rat_env = str(options.rat_env)
submit = options.submit



with open(good_run) as f:
    interesting_runs = f.read().splitlines()

count=0

for run_num in interesting_runs:
    for subrunfile in glob.glob(input_dir+"SNOP*"+str(run_num)+"*.l2.zdab"):
        
        subrun_num = subrunfile.split(str(run_num)+"_",1)[1] 
        subrun_num = subrun_num.split(".l2",1)[0] 

        script = open("process_"+str(run_num)+"_"+str(subrun_num)+".sh", "w")

        print 'Generating script for run number:', run_num, 'subrun number', subrun_num
        
        script.write("#!/bin/sh\n")
       
        script.write("cd "+rat_dir+"\n")
        script.write("source "+rat_env+"\n")
        script.write("mkdir /tmp/{0}_{1}\n".format(run_num,subrun_num))
        script.write("cd /tmp/{0}_{1}\n".format(run_num,subrun_num))
        script.write('rat -i {0}/SNOP_0000{1}_{2}.l2.zdab mac/processing/water/first_pass_data_cleaning.mac\n'.format(input_dir, run_num, subrun_num))
        script.write('rat -i {0}/SNOP_0000{1}_{2}.l2.zdab -o {3}/output_r{1}_s{2} mac/processing/water/second_pass_processing.mac\n'.format(input_dir, run_num, subrun_num, output_dir))
        script.write('rat -i {2}/output_r{0}_s{1}.root -o {2}/Processing_r{0}_s{1}_p000 mac/processing/water/third_pass_analysis_processing.mac\n'.format(run_num, subrun_num, output_dir))
        script.write("cp /tmp/{0}_{1}/Processing_r{0}_s{1}_p000* {2}/\n".format(run_num,subrun_num,output_dir))
        script.close()

        if submit:
            count+=1
            #attempt not to overload the servers
            if(count%10==0): os.system("sleep 1800")
            os.system("qsub -cwd -l h_rss=4G,h_vmem=4G -q SL6 process_"+str(run_num)+"_"+str(subrun_num)+".sh")
