import commands

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input files")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
parser.add_argument("-t", "--ttree", dest="ttree", default="Ana/passedEvents", help="TTree Name")
args = parser.parse_args()

file = commands.getstatusoutput('ls'+str(args.inputfiles))
for filename in file:
    print "filename="+str(filename)
