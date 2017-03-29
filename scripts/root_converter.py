#!/usr/bin/python

import os
import glob
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-i", "--input", dest="dirname",
                  help="Directory where the files to convert are kept")

(options, args) = parser.parse_args()

directory = str(options.dirname) + "/"

for filename in glob.glob(directory+"/"+"*[0-9].root"):
    print filename
    filename = filename.replace(".root","")
    if not os.path.exists(directory+filename+"_ntuple.root"):
        print 'rat -i %(directory)s%(filename)s.root -o %(directory)s%(filename)s_ntuple.root scripts/root_converter.mac' % vars()
        os.system('rat -i %(directory)s%(filename)s.root -o %(directory)s%(filename)s_ntuple.root scripts/root_converter.mac' % vars()) 
