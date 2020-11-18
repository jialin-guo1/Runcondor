#!/usr/bin/env python
import ROOT
f = ROOT.TFile("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/2016_DATA_new.root")
t = f.Get("passedEvents")
h = ROOT.TH1D("h","h",50,70,170)
for ievent,event in enumerate(t):
    if(not event.passedFullSelection):
        h.Fill(event.H_FSR)

h.Draw("E1")
