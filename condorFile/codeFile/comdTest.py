
import os

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input files")
parser.add_argument("-t", "--ttree", dest="ttree", default="Ana/passedEvents", help="TTree Name")
args = parser.parse_args()

import ROOT
chain = ROOT.TChain(args.ttree)
dataset = "root://cms-xrd-global.cern.ch/"+args.inputfiles+"/*.root"
print dataset
chain.Add(dataset)
print 'Total number of events: ' + str(chain.GetEntries())

def ifROOT(line):
    line=line.strip('\n')
    if line[-4:] != "root":
        return False

outputfile = os.popen('xrdfs root://cmsio5.rc.ufl.edu/ ls  '+str(args.inputfiles))
for line in outputfile:
    if(ifROOT(line)==False):
        continue
    line=line.strip('\n')
    filename = "root://cms-xrd-global.cern.ch/"+str(line)
    files = ROOT.TFile(filename)
    print "root://cms-xrd-global.cern.ch/"+line
#    print "filename="+str(line.strip('\n'))
