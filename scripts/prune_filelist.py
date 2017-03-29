#!/usr/bin/python
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-g", "--goodruns", dest="good_runs",
                  help="Good run list")
parser.add_option("-i", "--input", dest="input_file",
                  help="Input file to be pruned")
parser.add_option("-o", "--output", dest="output_file",
                  help="destination file to store pruned list")

(options, args) = parser.parse_args()

good_run = str(options.good_runs)
input_file = str(options.input_file)
output_file = str(options.output_file)


with open(good_run) as f1:
    good_runs = [str(int(line)) for line in f1]
    
newfile = open(output_file, "w")

with open(input_file) as oldfile:
    for line in oldfile:
        if any(good_run in line for good_run in good_runs):
            newfile.write(line)
    
