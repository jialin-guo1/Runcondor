import commands
import os

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input files")
args = parser.parse_args()

outputfile = os.popen('xrdfs root://cmsio5.rc.ufl.edu/ ls  '+str(args.inputfiles))
for line in outputfile:
    print "filename="+str(line.strip('\n'))
