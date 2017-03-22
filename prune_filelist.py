#!/usr/bin/python

with open("good_run_list.txt") as f1:
    good_runs = [str(int(line)) for line in f1]


    
newfile = open("filelist_goodruns_ntuple.dat", "w")

with open("filelist_ntuple.dat") as oldfile:
    for line in oldfile:
        if any(good_run in line for good_run in good_runs):
            newfile.write(line)
    
