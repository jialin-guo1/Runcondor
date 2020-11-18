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

passedFullSelection = array('i',[-90])
passedTrig = array('i',[-90])
passedZXCRSelection = array('i',[-90])
nZXCRFailedLeptons = array('i',[-90])
H_FSR = array('f',[0])

file_out = ROOT.TFile(args.outputfile, 'recreate')
passedEvents = ROOT.TTree("passedEvents","passedEvents")
passedEvents.Branch("passedZXCRSelection",passedZXCRSelection,"passedZXCRSelection/I")
passedEvents.Branch("nZXCRFailedLeptons",nZXCRFailedLeptons,"nZXCRFailedLeptons/I")
passedEvents.Branch("passedFullSelection",passedFullSelection,"passedFullSelection/I")
passedEvents.Branch("passedTrig",passedTrig,"passedTrig/I")
passedEvents.Branch("H_FSR",H_FSR,"H_FSR/F")

for ievent,event in enumerate(chain):
    passedZXCRSelection[0] = event.passedZXCRSelection
    nZXCRFailedLeptons[0] = event.nZXCRFailedLeptons
    passedFullSelection[0] = event.passedFullSelection
    passedTrig[0] = event.passedTrig
    l1FSR = ROOT.TLorentzVector()
    l2FSR = ROOT.TLorentzVector()
    l3FSR = ROOT.TLorentzVector()
    l4FSR = ROOT.TLorentzVector()
    l1FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[0]],event.lepFSR_eta[event.lep_Hindex[0]],event.lepFSR_phi[event.lep_Hindex[0]],event.lepFSR_mass[event.lep_Hindex[0]])
    l2FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[1]],event.lepFSR_eta[event.lep_Hindex[1]],event.lepFSR_phi[event.lep_Hindex[1]],event.lepFSR_mass[event.lep_Hindex[1]])
    l3FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[2]],event.lepFSR_eta[event.lep_Hindex[2]],event.lepFSR_phi[event.lep_Hindex[2]],event.lepFSR_mass[event.lep_Hindex[2]])
    l4FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[3]],event.lepFSR_eta[event.lep_Hindex[3]],event.lepFSR_phi[event.lep_Hindex[3]],event.lepFSR_mass[event.lep_Hindex[3]])
    H4massFSR = ROOT.TLorentzVector()
    H4massFSR = l1FSR+l2FSR+l3FSR+l4FSR
    H_FSR[0] = H4massFSR.M()
    passedEvents.Fill()

file_out.Write()
file_out.Close()
