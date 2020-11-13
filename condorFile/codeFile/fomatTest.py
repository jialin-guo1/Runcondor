import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input files")
parser.add_argument("-t", "--ttree", dest="ttree", default="Ana/passedEvents", help="TTree Name")
#parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
args = parser.parse_args()

import ROOT
chain = ROOT.TChain(args.ttree)

for i in range(1,10):
    filename = "root://cms-xrd-global.cern.ch//store/user/ferrico/2018data/UFHZZAnalysisRun2/myTask_MC/ZZTo4L_13TeV_powheg_pythia8_ext1/crab_ZZTo4L_13TeV_powheg_pythia8_ext1_RunIISummer16MiniAODv2/191202_132458/0000/ZZTo4L_13TeV_powheg_pythia8_ext1_RunIISummer16MiniAODv2_{0:s}.root".format(str(i))
    print filename
    chain.Add(filename)

print 'Total number of events: ' + str(chain.GetEntries())
