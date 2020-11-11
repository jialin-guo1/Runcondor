
import os

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input files")
args = parser.parse_args()

def ifROOT(line):
    line=line.strip('\n')
#    print line[-3:]
    if line[-4:] != "root":
        return False

outputfile = os.popen('xrdfs root://cmsio5.rc.ufl.edu/ ls  '+str(args.inputfiles))
for line in outputfile:
    if(ifROOT(line)==False):
        continue
    print "root://cms-xrd-global.cern.ch/"+args.inputfiles+"/*.root"
    print "root://cms-xrd-global.cern.ch/"+line
#    print "filename="+str(line.strip('\n'))
