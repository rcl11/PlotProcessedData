#!/usr/bin/python

import os
import glob
import re

from optparse import OptionParser


parser = OptionParser()
parser.add_option("-o", "--output_dir", dest="output_dir",
                  help="Directory where the root files are to be stored")
parser.add_option("-r", "--rat_dir", dest="rat_dir",
                  default = "/data/lane/RAT/rat", help="Path of rat release")
parser.add_option("-e", "--rat_env", dest="rat_env",
                  default = "/data/lane/RAT/env_rat-6.2.0.sh", help="Environment file for rat")
parser.add_option("-f", "--filelist", dest="filelist",
                  default = "filelists/filelist_test.dat", help="Filelist containing runs to be simulated")
parser.add_option("-s", "--submit", default = False, action = "store_true", dest="submit",
                  help="Actually submit the batch jobs! Please do a dry run first without this!")

(options, args) = parser.parse_args()

output_dir = str(options.output_dir) + "/"
rat_dir = str(options.rat_dir)
rat_env = str(options.rat_env)
submit = options.submit
filelist = str(options.filelist)



with open(filelist) as f:
    interesting_runs = f.read().splitlines()

q = re.compile("_r(.*)_s")

#simulations = ['Tl208','Tl208_av','Bi214','Bi214_av']
simulations = ['Tl208','Tl208_exwater','Bi214','Bi214_exwater']
#simulations = ['Tl208']
count = 0
for simulation in simulations:
    for run in interesting_runs:
        run_num = q.findall(run)[0]

        script = open("produce_"+str(simulation)+"_"+str(run_num)+".sh", "w")

        print 'Generating script for run number:', run_num
        
        script.write("#!/bin/sh\n")
       
        script.write("cd "+rat_dir+"\n")
        script.write("source "+rat_env+"\n")
        script.write("mkdir /tmp/{0}_{1}\n".format(run_num,simulation))
        script.write("cd /tmp/{0}_{1}\n".format(run_num,simulation))
        script.write('rat -r {0} {2}/mac/run-by-run-production/water/{1}.mac -o {1}_r{0}_s000_p000\n'.format(run_num,simulation,rat_dir))
        script.write("cp /tmp/{0}_{1}/{1}_r{0}_s000_p000* {2}/\n".format(run_num,simulation,output_dir))
        script.close()

        if submit:
            count+=1
            #attempt not to overload the servers
            if(count%10==0): os.system("sleep 600")
            os.system("qsub -cwd -l h_rss=4G,h_vmem=4G -q SL6 produce_"+str(simulation)+"_"+str(run_num)+".sh")
