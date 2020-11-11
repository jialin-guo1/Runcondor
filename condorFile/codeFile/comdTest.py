
import os

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input files")
parser.add_argument("-t", "--ttree", dest="ttree", default="Ana/passedEvents", help="TTree Name")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
args = parser.parse_args()


import ROOT
chain = ROOT.TChain(args.ttree)
#dataset = "root://cms-xrd-global.cern.ch/"+args.inputfiles+"/*.root"
#print dataset
#chain.Add(dataset)
#print 'Total number of events: ' + str(chain.GetEntries())

def ifROOT(line):
    line=line.strip('\n')
    if line[-4:] != "root":
        return False

Nevent = 0

outputfile = os.popen('xrdfs root://cmsio5.rc.ufl.edu/ ls  '+str(args.inputfiles))
for line in outputfile:
    if(ifROOT(line)==False):
        continue
    line=line.strip('\n')
    filename = "root://cms-xrd-global.cern.ch/"+str(line)
    chain.Add(filename)
    print "chain "+filename+" file!"
    files = ROOT.TFile.Open(filename)
    Nevent_h = files.Ana.Get('nEvents')
    Nevent += Nevent_h.GetBinContent(1)
    print "root://cms-xrd-global.cern.ch/"+line

print 'Total number of events: ' + str(chain.GetEntries())
print "Total Nevent = "+str(Nevent)

from array import array
weight = array('f',[0.])
file_out = ROOT.TFile(args.outputfile, 'recreate')
passedEvents = ROOT.TTree("passedEvents","passedEvents")
passedEvents.Branch("weight",weight,"weight/F")
for ievent,event in enumerate(chain):
    if(not event.passedTrig): continue
    if(not event.passedFullSelection): continue
    weight[0] = event.eventWeight/Nevent
    passedEvents.Fill()

file_out.Write()
file_out.Close()
