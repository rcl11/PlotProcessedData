#!/usr/bin/python

import os
import glob
import re

for subrunfile in glob.glob("/data/snoplus/OfficialProcessing/processed_prewater_freija/nino.lbl.gov/~freija/processed-latest/*.root"):
    print subrunfile
    j = re.compile("/([0-9]+)_")
    run = j.findall(subrunfile)[0]
    k = re.compile("_([0-9]+)-processed.root")
    subrun = k.findall(subrunfile)[0]
    newname = "Processing_r"+run+"_s"+subrun+"_p000.root"
    newfullpath = subrunfile.replace(run+"_"+subrun+"-processed.root",newname)
    print newfullpath
    os.system("mv " + subrunfile + " " + newfullpath) 

