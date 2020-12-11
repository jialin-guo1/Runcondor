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
#outputfile0 = os.popen('xrdfs root://cmsio5.rc.ufl.edu/ ls /store/user/ferrico/2018data/UFHZZAnalysisRun2/myTask_MC/ZZTo4L_13TeV_powheg_pythia8_ext1/crab_ZZTo4L_13TeV_powheg_pythia8_ext1_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/191119_094706/0000')

for line in outputfile:
    if(ifROOT(line)==False):
        continue
    line=line.strip('\n')
    filename = "root://cms-xrd-global.cern.ch/"+str(line)
    chain.Add(filename)
    print "chain "+filename+" file!"
    files = ROOT.TFile.Open(filename)
    SumW_h = files.Ana.Get('sumWeights')
    SumW += SumW_h.GetBinContent(1)

#print 'Total number of events before: ' + str(chain.GetEntries())
#print "Total Nevent before = "+str(SumW)

#for line in outputfile0:
#    if(ifROOT(line)==False):
#        continue
#    line=line.strip('\n')
#    filename = "root://cms-xrd-global.cern.ch/"+str(line)
#    chain.Add(filename)
#    files = ROOT.TFile.Open(filename)
#    SumW_h = files.Ana.Get('sumWeights')
#    SumW += SumW_h.GetBinContent(1)

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
lep1_id = array("i",[0])
lep1_pt = array('f',[0.])
lep1_eta = array('f',[0.])
lep1_phi = array('f',[0.])
lep1_mass = array('f',[0.])
lep2_id = array("i",[0])
lep2_pt = array('f',[0.])
lep2_eta = array('f',[0.])
lep2_phi = array('f',[0.])
lep2_mass = array('f',[0.])
lep3_id = array("i",[0])
lep3_pt = array('f',[0.])
lep3_eta = array('f',[0.])
lep3_phi = array('f',[0.])
lep3_mass = array('f',[0.])
lep4_id = array("i",[0])
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
EMCweight = array('f',[0.0])
weight = array('f',[0.])
k_gg = array('f',[0.])
k_qq_qcd_dPhi = array('f',[0.])
k_qq_qcd_M = array('f',[0.])
k_qq_ewk = array('f',[0.])
k_qq_qcd_pt = array('f',[0.])
cross = array('f',[0.])
Cat = array('i',[0])
#ZX
SS4e = array('f',[0.])
SS4mu = array('f',[0.])
SS2e2mu = array('f',[0.])
SS2mu2e = array('f',[0.])
lep_RelIsoNoFSR1 = array('f',[0.])
lep_RelIsoNoFSR2 = array('f',[0.])
lep_RelIsoNoFSR3 = array('f',[0.])
lep_RelIsoNoFSR4 = array('f',[0.])
passedZXCRSelection = array('i',[0])
nZXCRFailedLeptons = array('i',[0])
passedFullSelection = array('i',[0])
passedSSCRSelection = array('i',[0])
passedTrig = array('i',[0])
#massWindow test
Run = array('l',[0])
Event = array('l',[0])
LumiSect = array('l',[0])




#Output file and any Branch we want
file_out = ROOT.TFile(args.outputfile, 'recreate')
passedEvents = ROOT.TTree("passedEvents","passedEvents")
passedEvents.Branch("lep1_id",lep1_id,"lep1_id/I")
passedEvents.Branch("lep1_pt",lep1_pt,"lep1_pt/F")
passedEvents.Branch("lep1_eta",lep1_eta,"lep1_eta/F")
passedEvents.Branch("lep1_phi",lep1_phi,"lep1_phi/F")
passedEvents.Branch("lep1_mass",lep1_mass,"lep1_mass/F")
passedEvents.Branch("lep2_id",lep2_id,"lep2_id/I")
passedEvents.Branch("lep2_pt",lep2_pt,"lep2_pt/F")
passedEvents.Branch("lep2_eta",lep2_eta,"lep2_eta/F")
passedEvents.Branch("lep2_phi",lep2_phi,"lep2_phi/F")
passedEvents.Branch("lep2_mass",lep2_mass,"lep2_mass/F")
passedEvents.Branch("lep3_id",lep3_id,"lep3_id/I")
passedEvents.Branch("lep3_pt",lep3_pt,"lep3_pt/F")
passedEvents.Branch("lep3_eta",lep3_eta,"lep3_eta/F")
passedEvents.Branch("lep3_phi",lep3_phi,"lep3_phi/F")
passedEvents.Branch("lep3_mass",lep3_mass,"lep3_mass/F")
passedEvents.Branch("lep4_id",lep4_id,"lep4_id/I")
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
passedEvents.Branch("lep_4mass",lep_4mass,"lep_4mass/F")
passedEvents.Branch("lepnoFSR_4mass",lepnoFSR_4mass,"lepnoFSR_4mass/F")
passedEvents.Branch("ledZ_mass",ledZ_mass,"ledZ_mass/F")
passedEvents.Branch("subledZ_mass",subledZ_mass,"subledZ_mass/F")
passedEvents.Branch("weight",weight,"weight/F")
passedEvents.Branch("EMCweight",EMCweight,"EMCweight/F")
passedEvents.Branch("k_gg",k_gg,"k_gg/F")
passedEvents.Branch("k_qq_qcd_dPhi",k_qq_qcd_dPhi,"k_qq_qcd_dPhi/F")
passedEvents.Branch("k_qq_qcd_M",k_qq_qcd_M,"k_qq_qcd_M/F")
passedEvents.Branch("k_qq_ewk",k_qq_ewk,"k_qq_ewk/F")
passedEvents.Branch("k_qq_qcd_pt",k_qq_qcd_pt,"k_qq_qcd_pt/F")
passedEvents.Branch("cross",cross,"cross/F")
passedEvents.Branch("Cat",Cat,"Cat/I")
passedEvents.Branch("lep_RelIsoNoFSR1",lep_RelIsoNoFSR1,"lep_RelIsoNoFSR1/F")
passedEvents.Branch("lep_RelIsoNoFSR2",lep_RelIsoNoFSR2,"lep_RelIsoNoFSR1/F")
passedEvents.Branch("lep_RelIsoNoFSR3",lep_RelIsoNoFSR3,"lep_RelIsoNoFSR1/F")
passedEvents.Branch("lep_RelIsoNoFSR4",lep_RelIsoNoFSR4,"lep_RelIsoNoFSR1/F")
passedEvents.Branch("passedZXCRSelection",passedZXCRSelection,"passedZXCRSelection/I")
passedEvents.Branch("nZXCRFailedLeptons",nZXCRFailedLeptons,"nZXCRFailedLeptons/I")
passedEvents.Branch("passedFullSelection",passedFullSelection,"passedFullSelection/I")
passedEvents.Branch("passedSSCRSelection",passedSSCRSelection,"passedSSCRSelection/I")
passedEvents.Branch("passedTrig",passedTrig,"passedTrig/I")
passedEvents.Branch("Run",Run,"Run/l")
passedEvents.Branch("Event",Event,"Event/l")
passedEvents.Branch("LumiSect",LumiSect,"LumiSect/l")
passedEvents.Branch("SS4e",SS4e,"SS4e/F")
passedEvents.Branch("SS4mu",SS4mu,"SS4mu/F")
passedEvents.Branch("SS2e2mu",SS2e2mue,"SS2e2mue/F")
passedEvents.Branch("SS2mu2e",SS2mu2e,"SS2mu2e/F")


#Loop over all the events in the input ntuple
for ievent,event in enumerate(chain):


    n_Zs = 0
    Z_pt = []
    Z_eta = []
    Z_phi = []
    Z_mass = []
    Z_lepindex1 = []
    Z_lepindex2 = []
    lep_index = [0,0,0,0]

    Zmass = 91.1876
    minZ1DeltaM = 99999.99
#    if(ievent==1000): break
    passedZXCRSelection[0] = event.passedZXCRSelection
    nZXCRFailedLeptons[0] = event.nZXCRFailedLeptons
    passedFullSelection[0] = event.passedFullSelection
    passedTrig[0] = event.passedTrig
    passedSSCRSelection[0] = False


    lep_RelIsoNoFSR1[0] = event.lep_RelIsoNoFSR[event.lep_Hindex[0]]
    lep_RelIsoNoFSR2[0] = event.lep_RelIsoNoFSR[event.lep_Hindex[1]]
    lep_RelIsoNoFSR3[0] = event.lep_RelIsoNoFSR[event.lep_Hindex[2]]
    lep_RelIsoNoFSR4[0] = event.lep_RelIsoNoFSR[event.lep_Hindex[3]]

    lep_4mass[0] = event.mass4l
    lepnoFSR_4mass[0] = event.mass4l_noFSR
    ledZ_mass[0] = event.massZ1
    subledZ_mass[0] = event.massZ2
    weight[0] = event.eventWeight
    EMCweight[0] = event.dataMCWeight
    k_gg[0] = event.k_ggZZ
    k_qq_qcd_dPhi[0] = event.k_qqZZ_qcd_dPhi
    k_qq_qcd_M[0] = event.k_qqZZ_qcd_M
    k_qq_ewk[0] = event.k_qqZZ_ewk
    k_qq_qcd_pt[0] = event.k_qqZZ_qcd_Pt
    cross[0] = event.crossSection
    Cat[0] = event.EventCat
    Run[0] = event.Run
    Event[0] = event.Event
    LumiSect[0] = event.LumiSect

    lep1_pt[0] = event.lep_pt[event.lep_Hindex[0]]
    lep1_eta[0] = event.lep_eta[event.lep_Hindex[0]]
    lep1_phi[0] = event.lep_phi[event.lep_Hindex[0]]
    lep1_mass[0] = event.lep_mass[event.lep_Hindex[0]]
    lep1_id[0] = event.lep_id[event.lep_Hindex[0]]


    lep2_pt[0] = event.lep_pt[event.lep_Hindex[1]]
    lep2_eta[0] = event.lep_eta[event.lep_Hindex[1]]
    lep2_phi[0] = event.lep_phi[event.lep_Hindex[1]]
    lep2_mass[0] = event.lep_mass[event.lep_Hindex[1]]
    lep2_id[0] = event.lep_id[event.lep_Hindex[1]]

    lep3_pt[0] = event.lep_pt[event.lep_Hindex[2]]
    lep3_eta[0] = event.lep_eta[event.lep_Hindex[2]]
    lep3_phi[0] = event.lep_phi[event.lep_Hindex[2]]
    lep3_mass[0] = event.lep_mass[event.lep_Hindex[2]]
    lep3_id[0] = event.lep_id[event.lep_Hindex[2]]

    lep4_pt[0] = event.lep_pt[event.lep_Hindex[3]]
    lep4_eta[0] = event.lep_eta[event.lep_Hindex[3]]
    lep4_phi[0] = event.lep_phi[event.lep_Hindex[3]]
    lep4_mass[0] = event.lep_mass[event.lep_Hindex[3]]
    lep4_id[0] = event.lep_id[event.lep_Hindex[3]]

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

    # method call for Same sign control regions
    Nlep = event.lep_pt.size()

    for i in range(Nlep):
        for j in range(i+1,Nlep):
            # make all Z candidates including any FSR photons,
            lifsr = ROOT.TLorentzVector()
            ljfsr = ROOT.TLorentzVector()
            lifsr.SetPtEtaPhiM(event.lepFSR_pt[i], event.lepFSR_eta[i], event.lepFSR_phi[i], event.lepFSR_mass[i])
            ljfsr.SetPtEtaPhiM(event.lepFSR_pt[j], event.lepFSR_eta[j], event.lepFSR_phi[j], event.lepFSR_mass[j])
            Z = ROOT.TLorentzVector()
            Z = (lifsr + ljfsr)
            if(Z.M()>0):
                n_Zs = n_Zs + 1
                Z_pt.append(Z.Pt())
                Z_eta.append(Z.Eta())
                Z_phi.append(Z.Phi())
                Z_mass.append(Z.M())
                Z_lepindex1.append(i)
                Z_lepindex2.append(j)

    for i in range(n_Zs):
        for j in range(i+1,n_Zs):
            i1 = Z_lepindex1[i]
            i2 = Z_lepindex2[i]
            j1 = Z_lepindex1[j]
            j2 = Z_lepindex2[j]
            if(i1 == j1 or i1==j2 or i2==j1 or i2==j2): continue
            # opposite sign and same flavour for Z ; same sign same flavour for fack leptons
            if(event.lep_id[i1]+event.lep_id[i2]!=0 or event.lep_id[j1]*event.lep_id[j2]<0 or abs(event.lep_id[j1])!=abs(event.lep_id[j2])): continue
            lep1 = ROOT.TLorentzVector()
            lep2 = ROOT.TLorentzVector()
            lep3 = ROOT.TLorentzVector()
            lep4 = ROOT.TLorentzVector()
            lep1.SetPtEtaPhiM(event.lepFSR_pt[i1],event.lepFSR_eta[i1],event.lepFSR_phi[i1],event.lepFSR_mass[i1])
            lep2.SetPtEtaPhiM(event.lepFSR_pt[i2],event.lepFSR_eta[i2],event.lepFSR_phi[i2],event.lepFSR_mass[i2])
            lep3.SetPtEtaPhiM(event.lepFSR_pt[j1],event.lepFSR_eta[j1],event.lepFSR_phi[j1],event.lepFSR_mass[j1])
            lep4.SetPtEtaPhiM(event.lepFSR_pt[j2],event.lepFSR_eta[j2],event.lepFSR_phi[j2],event.lepFSR_mass[j2])

            Z1 = ROOT.TLorentzVector()
            SSleps = ROOT.TLorentzVector()
            SSleps = lep3+lep4
            Z1 = lep1 + lep2
            Z1_lepindex = [0,0]
            Z1DeltaM = abs(Z1.M()-Zmass)

            # Check Leading and Subleading pt Cut
            if(lep1.Pt()>lep2.Pt()):
                Z1_lepindex[0] = i1
                Z1_lepindex[1] = i2
                if(lep1.Pt()<20 or lep2.Pt()<10): continue
            else:
                Z1_lepindex[0] = i2
                Z1_lepindex[1] = i1
                if(lep2.Pt()<20 or lep1.Pt()<10): continue
            npassedpt += 1
            # Check dR(li,lj)>0.02 for any i,j
            if (deltaR(lep1.Eta(),lep1.Phi(),lep2.Eta(),lep2.Phi())<0.02): continue
            if (deltaR(lep1.Eta(),lep1.Phi(),lep3.Eta(),lep3.Phi())<0.02): continue
            if (deltaR(lep1.Eta(),lep1.Phi(),lep4.Eta(),lep4.Phi())<0.02): continue
            if (deltaR(lep2.Eta(),lep2.Phi(),lep3.Eta(),lep3.Phi())<0.02): continue
            if (deltaR(lep2.Eta(),lep2.Phi(),lep4.Eta(),lep4.Phi())<0.02): continue
            if (deltaR(lep3.Eta(),lep3.Phi(),lep4.Eta(),lep4.Phi())<0.02): continue
            npassedDR += 1
            # Check M(l+,l-)>4.0 GeV for any OS pair , Do not include FSR photons
            allM = []
            lep1_noFSR = ROOT.TLorentzVector()
            lep2_noFSR = ROOT.TLorentzVector()
            lep3_noFSR = ROOT.TLorentzVector()
            lep4_noFSR = ROOT.TLorentzVector()
            lep1_noFSR.SetPtEtaPhiM(event.lep_pt[i1],event.lep_eta[i1],event.lep_phi[i1],event.lep_mass[i1])
            lep2_noFSR.SetPtEtaPhiM(event.lep_pt[i2],event.lep_eta[i2],event.lep_phi[i2],event.lep_mass[i2])
            lep3_noFSR.SetPtEtaPhiM(event.lep_pt[j1],event.lep_eta[j1],event.lep_phi[j1],event.lep_mass[j1])
            lep4_noFSR.SetPtEtaPhiM(event.lep_pt[j2],event.lep_eta[j2],event.lep_phi[j2],event.lep_mass[j2])
            i1i2 = ROOT.TLorentzVector()
            i1i2 = lep1_noFSR+lep2_noFSR
            allM.append(i1i2.M())
            if(event.lep_id[i1]*event.lep_id[j1]<0):
                i1j1 = ROOT.TLorentzVector()
                i1j2 = ROOT.TLorentzVector()
                i1j2 = lep1_noFSR+lep4_noFSR
                i1j1 = lep1_noFSR+lep3_noFSR
                allM.append(i1j1.M())
                allM.append(i1j2.M())
            else:
                i2j1 = ROOT.TLorentzVector()
                i2j2 = ROOT.TLorentzVector()
                i2j1 = lep2_noFSR+lep3_noFSR
                i2j2 = lep2_noFSR+lep4_noFSR
                allM.append(i2j1.M())
                allM.append(i2j2.M())
            allM.sort()
            # print " list of OS pair mass after sort : " + str(allM)
            if(allM[0]<4.0): continue

            # Check isolation cut (without FSR ) for Z1 leptons
            #if(abs(event.lep_id[Z1_lepindex[0]])==13):
            if( event.lep_RelIsoNoFSR[Z1_lepindex[0]]>0.35): continue
            #if(abs(event.lep_id[Z1_lepindex[1]])==13):
            if( event.lep_RelIsoNoFSR[Z1_lepindex[1]]>0.35): continue
            npassedios +=1
            # Check tight ID cut for Z1 leptons
            if (not event.lep_tightId[Z1_lepindex[0]]): continue
            if (not event.lep_tightId[Z1_lepindex[1]]): continue
            # check the masswindow
            if(Z1.M()<40 or Z1.M()>120 or SSleps.M()<12 or SSleps.M()>120): continue
            # Check if this candidate has the best Z1 and highest scalar sum of Z2 lepton pt
            if(Z1DeltaM<minZ1DeltaM):
                minZ1DeltaM = Z1DeltaM
                lep_index[0] = Z1_lepindex[0]
                lep_index[1] = Z1_lepindex[1]
                lep_index[2] = j1
                lep_index[3] = j2
                passedSSCRselection[0] = True


    if(passedSSCRSelection[0]):
        l1 = ROOT.TLorentzVector()
        l2 = ROOT.TLorentzVector()
        l3 = ROOT.TLorentzVector()
        l4 = ROOT.TLorentzVector()
        l1.SetPtEtaPhiM(event.lepFSR_pt[lep_index[0]],event.lepFSR_eta[lep_index[0]],event.lepFSR_phi[lep_index[0]],event.lepFSR_mass[lep_index[0]])
        l2.SetPtEtaPhiM(event.lepFSR_pt[lep_index[1]],event.lepFSR_eta[lep_index[1]],event.lepFSR_phi[lep_index[1]],event.lepFSR_mass[lep_index[1]])
        l3.SetPtEtaPhiM(event.lepFSR_pt[lep_index[2]],event.lepFSR_eta[lep_index[2]],event.lepFSR_phi[lep_index[2]],event.lepFSR_mass[lep_index[2]])
        l4.SetPtEtaPhiM(event.lepFSR_pt[lep_index[3]],event.lepFSR_eta[lep_index[3]],event.lepFSR_phi[lep_index[3]],event.lepFSR_mass[lep_index[3]])
        mass4l = (l1+l2+l3+l4).M()

        if(abs(event.lep_id[lep_index[0]])==abs(event.lep_id[lep_index[1]])==abs(event.lep_id[lep_index[2]])==abs(event.lep_id[lep_index[3]])==11):
            SS4e[0] = mass4l
        if(abs(event.lep_id[lep_index[0]])==abs(event.lep_id[lep_index[1]])==abs(event.lep_id[lep_index[2]])==abs(event.lep_id[lep_index[3]])==13):
            SS4mu[0] = mass4l
        if(abs(event.lep_id[lep_index[0]])==abs(event.lep_id[lep_index[1]])==11 and abs(event.lep_id[lep_index[2]])==abs(event.lep_id[lep_index[3]])==13):
            SS2e2mu[0] = mass4l
        if(abs(event.lep_id[lep_index[0]])==abs(event.lep_id[lep_index[1]])==13 and abs(event.lep_id[lep_index[2]])==abs(event.lep_id[lep_index[3]])==11):
            SS2mu2e[0] = mass4l


    passedEvents.Fill()


file_out.Write()
file_out.Close()
