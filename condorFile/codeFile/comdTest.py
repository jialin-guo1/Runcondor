import commands

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input files")
args = parser.parse_args()

file = commands.getstatusoutput('xrdfs root://cmsio5.rc.ufl.edu/ ls  '+str(args.inputfiles))
for filename in file:
    print "filename="+str(filename)
