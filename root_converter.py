#!/usr/bin/python

import os
import glob

for filename in glob.glob("*[0-9].root"):
    print filename
    filename = filename.replace(".root","")
    print filename
    print 'rat -i /data/snoplus/OfficialProcessing/processed_6.2.3/%(filename)s.root -o /data/snoplus/OfficialProcessing/processed_6.2.3/%(filename)s_ntuple.root root_converter.mac' % vars()
    os.system('rat -i /data/snoplus/OfficialProcessing/processed_6.2.3/%(filename)s.root -o /data/snoplus/OfficialProcessing/processed_6.2.3/%(filename)s_ntuple.root root_converter.mac' % vars()) 
