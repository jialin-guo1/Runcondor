#!/user/bin/env python

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input files")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
parser.add_argument("-t", "--ttree", dest="ttree", default="Ana/passedEvents", help="TTree Name")
args = parser.parse_args()

import numpy as np
import ROOT
from array import array

chain = ROOT.TChain(args.ttree)
#chain.Add("root://cms-xrd-global.cern.ch/"+args.inputfiles+"/*.root")
#print 'Total number of events: ' + str(chain.GetEntries())

#get nevent
SumW = 0

import os
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
    chain.Add(filename)
    print "chain "+filename+" file!"
    files = ROOT.TFile.Open(filename)

#    Nevent_h = files.Ana.Get('nEvents')
#    Nevent += Nevent_h.GetBinContent(1)

    SumW_h = files.Ana.Get('sumWeights')
    SumW += SumW_h.GetBinContent(1)

print 'Total number of events: ' + str(chain.GetEntries())
print "Total Nevent = "+str(SumW)
#    SumWPU_h = files.Ana.Get('sumWeightsPU')
#    SumWPU += SumWPU_h.GetBinContent(1)

#    NV_h = files.Ana.Get('nVtx')
#    NV +=  NV_h.GetBinContent(1)

#    NV_ReW_h = files.Ana.Get('nVtx_ReWeighted')
#    NV_ReW += NV_ReW_h.GetBinContent(1)

#    NInter_h = files.Ana.Get('nInteractions')
#    NInter += NInter_h.GetBinContent(1)

#    NInter_ReW_h = files.Ana.Get('nInteraction_ReWeighted')
#    NInter_ReW += NInter_ReW_h.GetBinContent(1)

#variables
lep1_pt = array('f',[0.])
lep1_eta = array('f',[0.])
lep1_phi = array('f',[0.])
lep1_mass = array('f',[0.])
lep2_pt = array('f',[0.])
lep2_eta = array('f',[0.])
lep2_phi = array('f',[0.])
lep2_mass = array('f',[0.])
lep3_pt = array('f',[0.])
lep3_eta = array('f',[0.])
lep3_phi = array('f',[0.])
lep3_mass = array('f',[0.])
lep4_pt = array('f',[0.])
lep4_eta = array('f',[0.])
lep4_phi = array('f',[0.])
lep4_mass = array('f',[0.])
lep1FSR_pt = array('f',[0.])
lep1FSR_eta = array('f',[0.])
lep1FSR_phi = array('f',[0.])
lep1FSR_mass = array('f',[0.])
lep2FSR_pt = array('f',[0.])
lep2FSR_eta = array('f',[0.])
lep2FSR_phi = array('f',[0.])
lep2FSR_mass = array('f',[0.])
lep3FSR_pt = array('f',[0.])
lep3FSR_eta = array('f',[0.])
lep3FSR_phi = array('f',[0.])
lep3FSR_mass = array('f',[0.])
lep4FSR_pt = array('f',[0.])
lep4FSR_eta = array('f',[0.])
lep4FSR_phi = array('f',[0.])
lep4FSR_mass = array('f',[0.])
lep_4mass = array('f',[0.])
lepnoFSR_4mass = array('f',[0])
ledZ_mass = array('f',[0])
subledZ_mass  = array('f',[0])
h_pt = array('f',[0.])
h_eta = array('f',[0.])
h_phi = array('f',[0.])
h_mass = array('f',[0.])
H = array('f',[0.])
H_FSR = array('f',[0.])
EMCweight = array('f',[0.0])
weight = array('f',[0.])


#Output file and any Branch we want
file_out = ROOT.TFile(args.outputfile, 'recreate')
passedEvents = ROOT.TTree("passedEvents","passedEvents")
passedEvents.Branch("lep1_pt",lep1_pt,"lep1_pt/F")
passedEvents.Branch("lep1_eta",lep1_eta,"lep1_eta/F")
passedEvents.Branch("lep1_phi",lep1_phi,"lep1_phi/F")
passedEvents.Branch("lep1_mass",lep1_mass,"lep1_mass/F")
passedEvents.Branch("lep2_pt",lep2_pt,"lep2_pt/F")
passedEvents.Branch("lep2_eta",lep2_eta,"lep2_eta/F")
passedEvents.Branch("lep2_phi",lep2_phi,"lep2_phi/F")
passedEvents.Branch("lep2_mass",lep2_mass,"lep2_mass/F")
passedEvents.Branch("lep3_pt",lep3_pt,"lep3_pt/F")
passedEvents.Branch("lep3_eta",lep3_eta,"lep3_eta/F")
passedEvents.Branch("lep3_phi",lep3_phi,"lep3_phi/F")
passedEvents.Branch("lep3_mass",lep3_mass,"lep3_mass/F")
passedEvents.Branch("lep4_pt",lep4_pt,"lep4_pt/F")
passedEvents.Branch("lep4_eta",lep4_eta,"lep4_eta/F")
passedEvents.Branch("lep4_phi",lep4_phi,"lep4_phi/F")
passedEvents.Branch("lep4_mass",lep4_mass,"lep4_mass/F")
passedEvents.Branch("lep1FSR_pt",lep1FSR_pt,"lep1FSR_pt/F")
passedEvents.Branch("lep1FSR_eta",lep1FSR_eta,"lep1FSR_eta/F")
passedEvents.Branch("lep1FSR_phi",lep1FSR_phi,"lep1FSR_phi/F")
passedEvents.Branch("lep1FSR_mass",lep1FSR_mass,"lep1FSR_mass/F")
passedEvents.Branch("lep2FSR_pt",lep2FSR_pt,"lep2FSR_pt/F")
passedEvents.Branch("lep2FSR_eta",lep2FSR_eta,"lep2FSR_eta/F")
passedEvents.Branch("lep2FSR_phi",lep2FSR_phi,"lep2_phi/F")
passedEvents.Branch("lep2FSR_mass",lep2FSR_mass,"lep2FSR_mass/F")
passedEvents.Branch("lep3FSR_pt",lep3FSR_pt,"lep3FSR_pt/F")
passedEvents.Branch("lep3FSR_eta",lep3FSR_eta,"lep3FSR_eta/F")
passedEvents.Branch("lep3FSR_phi",lep3FSR_phi,"lep3FSR_phi/F")
passedEvents.Branch("lep3FSR_mass",lep3FSR_mass,"lep3FSR_mass/F")
passedEvents.Branch("lep4FSR_pt",lep4FSR_pt,"lep4FSR_pt/F")
passedEvents.Branch("lep4FSR_eta",lep4FSR_eta,"lep4FSR_eta/F")
passedEvents.Branch("lep4FSR_phi",lep4FSR_phi,"lep4FSR_phi/F")
passedEvents.Branch("lep4FSR_mass",lep4FSR_mass,"lep4FSR_mass/F")
passedEvents.Branch("h_pt",h_pt,"h_pt/F")
passedEvents.Branch("h_eta",h_eta,"h_eta/F")
passedEvents.Branch("h_phi",h_phi,"h_phi/F")
passedEvents.Branch("h_mass",h_mass,"h_mass/F")
passedEvents.Branch("lep_4mass",lep_4mass,"lep_4mass/F")
passedEvents.Branch("lepnoFSR_4mass",lepnoFSR_4mass,"lepnoFSR_4mass/F")
passedEvents.Branch("ledZ_mass",ledZ_mass,"ledZ_mass/F")
passedEvents.Branch("subledZ_mass",subledZ_mass,"subledZ_mass/F")
passedEvents.Branch("H",H,"H/F")
passedEvents.Branch("H_FSR",H_FSR,"H_FSR/F")
passedEvents.Branch("weight",weight,"weight/F")
passedEvents.Branch("EMCweight",EMCweight,"EMCweight/F")

#Loop over all the events in the input ntuple
for ievent,event in enumerate(chain):
    if(not event.passedTrig): continue
    if(not event.passedFullSelection): continue

    lep_4mass[0] = event.mass4l
    lepnoFSR_4mass[0] = event.mass4l_noFSR
    ledZ_mass[0] = event.massZ1
    subledZ_mass[0] = event.massZ2
    weight[0] = event.eventWeight/SumW
    EMCweight[0] = event.dataMCWeight/SumW


    Nlep = event.lep_pt.size()
    for i in range(Nlep):

    #fill tree

       lep1_pt[0] = event.lep_pt[event.lep_Hindex[0]]
       lep1_eta[0] = event.lep_eta[event.lep_Hindex[0]]
       lep1_phi[0] = event.lep_phi[event.lep_Hindex[0]]
       lep1_mass[0] = event.lep_mass[event.lep_Hindex[0]]


       lep2_pt[0] = event.lep_pt[event.lep_Hindex[1]]
       lep2_eta[0] = event.lep_eta[event.lep_Hindex[1]]
       lep2_phi[0] = event.lep_phi[event.lep_Hindex[1]]
       lep2_mass[0] = event.lep_mass[event.lep_Hindex[1]]

       lep3_pt[0] = event.lep_pt[event.lep_Hindex[2]]
       lep3_eta[0] = event.lep_eta[event.lep_Hindex[2]]
       lep3_phi[0] = event.lep_phi[event.lep_Hindex[2]]
       lep3_mass[0] = event.lep_mass[event.lep_Hindex[2]]

       lep4_pt[0] = event.lep_pt[event.lep_Hindex[3]]
       lep4_eta[0] = event.lep_eta[event.lep_Hindex[3]]
       lep4_phi[0] = event.lep_phi[event.lep_Hindex[3]]
       lep4_mass[0] = event.lep_mass[event.lep_Hindex[3]]

       lep1FSR_pt[0] = event.lepFSR_pt[event.lep_Hindex[0]]
       lep1FSR_eta[0] = event.lepFSR_eta[event.lep_Hindex[0]]
       lep1FSR_phi[0] = event.lepFSR_phi[event.lep_Hindex[0]]
       lep1FSR_mass[0] = event.lepFSR_mass[event.lep_Hindex[0]]

       lep2FSR_pt[0] = event.lepFSR_pt[event.lep_Hindex[1]]
       lep2FSR_eta[0] = event.lepFSR_eta[event.lep_Hindex[1]]
       lep2FSR_phi[0] = event.lepFSR_phi[event.lep_Hindex[1]]
       lep2FSR_mass[0] = event.lepFSR_mass[event.lep_Hindex[1]]

       lep3FSR_pt[0] = event.lepFSR_pt[event.lep_Hindex[2]]
       lep3FSR_eta[0] = event.lepFSR_eta[event.lep_Hindex[2]]
       lep3FSR_phi[0] = event.lepFSR_phi[event.lep_Hindex[2]]
       lep3FSR_mass[0] = event.lepFSR_mass[event.lep_Hindex[2]]

       lep4FSR_pt[0] = event.lepFSR_pt[event.lep_Hindex[3]]
       lep4FSR_eta[0] = event.lepFSR_eta[event.lep_Hindex[3]]
       lep4FSR_phi[0] = event.lepFSR_phi[event.lep_Hindex[3]]
       lep4FSR_mass[0] = event.lepFSR_mass[event.lep_Hindex[3]]

       l1 = ROOT.TLorentzVector()
       l2 = ROOT.TLorentzVector()
       l3 = ROOT.TLorentzVector()
       l4 = ROOT.TLorentzVector()
       l1.SetPtEtaPhiM(event.lep_pt[event.lep_Hindex[0]],event.lep_eta[event.lep_Hindex[0]],event.lep_phi[event.lep_Hindex[0]],event.lep_mass[event.lep_Hindex[0]])
       l2.SetPtEtaPhiM(event.lep_pt[event.lep_Hindex[1]],event.lep_eta[event.lep_Hindex[1]],event.lep_phi[event.lep_Hindex[1]],event.lep_mass[event.lep_Hindex[1]])
       l3.SetPtEtaPhiM(event.lep_pt[event.lep_Hindex[2]],event.lep_eta[event.lep_Hindex[2]],event.lep_phi[event.lep_Hindex[2]],event.lep_mass[event.lep_Hindex[2]])
       l4.SetPtEtaPhiM(event.lep_pt[event.lep_Hindex[3]],event.lep_eta[event.lep_Hindex[3]],event.lep_phi[event.lep_Hindex[3]],event.lep_mass[event.lep_Hindex[3]])
       H4mass = ROOT.TLorentzVector()
       H4mass = l1+l2+l3+l4
       H[0] = H4mass.M()

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



    Hhiggs = event.H_pt.size()
    for i in range(Hhiggs):
        h_pt[0] = event.H_pt[i]
        h_eta[0] =event.H_eta[i]
        h_phi[0] = event.H_phi[i]
        h_mass[0] = event.H_mass[i]
#      lep_Hindex[i] = event.lep_Hindex[i]
    passedEvents.Fill()

file_out.Write()
file_out.Close()
