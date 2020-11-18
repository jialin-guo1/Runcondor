#!/usr/bin/env python
import ROOT
import os
import numpy as np
from array import array

import argparse
parser = argparse.ArgumentParser(description = "A simple ttree plotter")
parser.add_argument("-t", "--ttree", dest="ttree", default="Ana/passedEvents", help="TTree Name")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input root files")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
parser.add_argument("-s","--substring", dest="substring",default="",help='only submit datasets with this string in the name')
args = parser.parse_args()

print("get file " + str(args.inputfiles))
#input Ntuple
chain = ROOT.TChain(args.ttree)
chain.Add(args.inputfiles)

passedZXCRSelection = array('i',[0])
nZXCRFailedLeptons = array('i',[0])

file_out = ROOT.TFile(args.outputfile, 'recreate')
passedEvents = ROOT.TTree("passedEvents","passedEvents")
passedEvents.Branch("passedZXCRSelection",passedZXCRSelection,"passedZXCRSelection/I")
passedEvents.Branch("nZXCRFailedLeptons",nZXCRFailedLeptons,"nZXCRFailedLeptons/I")

for ievent,event in enumerate(chain):
    if ievent == 50000: break
    if(not event.passedTrig): continue
    if(not event.passedZXCRSelection): continue
    passedZXCRSelection[0] = event.passedZXCRSelection
    nZXCRFailedLeptons[0] = event.nZXCRFailedLeptons
    print "passedZXCRSelection = " + str(event.passedZXCRSelection)
    print "nZXCRFailedLeptons = " + str(event.nZXCRFailedLeptons)

    passedEvents.Fill()

file_out.Write()
file_out.Close()
