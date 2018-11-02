#!/usr/bin/env python

import itertools
from math import *
from collections import OrderedDict
import multiprocessing 
from multiprocessing import Pool
from string import *
from array import array
import ROOT
from ROOT import * 
import optparse 
import random

#global Preselection 
#Preselection = 1

def SortBySrootB(v):
    return v[0]

def quickplot(File, tree, plot, var, Cut, Weight):
    temp = plot.Clone("temp")
    chain = TChain(tree)
    chain.Add(File)
    chain.Draw(var + ">>" + "temp", "(" + Weight + ")*(" + Cut + ")", "goff") 
    plot.Add(temp)

window = {'twentyfive':"Evt_aM>0&Evt_aM<40", 'fifty':"Evt_aM>20&Evt_aM<70", "seventyfive":"Evt_aM>60&Evt_aM<80"}

def optimization(L):
    A = L[0]
    B = L[1]
    C = L[2]
    D = L[3]
   
    sig = TH1F("SIG" + str(L), "", 1000, 0, 500)
    qcd = TH1F("QCD" + str(L), "", 1000, 0, 500)
    ttbar = TH1F("TTBAR" + str(L), "", 1000, 0, 500)
#    Cut = "&" + Preselection
    Cut = " "
    Cut += window["fifty"] 
    Cut += "&(J1sdm - J2sdm)/(J1sdm + J2sdm)<" + str(B)
    Cut += "&abs(J1eta - J2eta)<" + str(A)
    Cut += "&J1tau21<" + str(C) + "&J2tau21<" + str(C)
    Cut += "&J1bb>" + str(D) + "&J2bb>" + str(D)
    quickplot("/home/rek81/userArea/CMSSW_8_0_20/src/nano_to_tree/signal/X1000a50_gxsd330p1.root", "tree", sig, "Evt_aM", Cut, "hlt_JJ300200m30*weight")
    quickplot("/home/rek81/userArea/CMSSW_8_0_20/src/nano_to_tree/QCD_newWeight/2016_QCD_HT.root", "tree", qcd, "Evt_aM", Cut, "hlt_JJ300200m30*weight")
    quickplot("/home/rek81/userArea/CMSSW_8_0_20/src/nano_to_tree/TT/TTandJets.root", "tree", ttbar, "Evt_aM", Cut, "hlt_JJ300200m30*weight")
    sigI = sig.Integral()
    qcdI = qcd.Integral()
    ttbarI = ttbar.Integral()
    if sigI > 0 and qcdI > 0 and ttbarI > 0:
        srootb = sigI/TMath.Sqrt(qcdI+ttbarI)
        print srootb, L
        return [srootb, L]
    else:
        return [0, L]

if __name__ == '__main__':
    from optparse import OptionParser    
    parser = OptionParser()
    parser.add_option('--p', '--preselection', metavar='PRESELECTION', type='string', dest='preselections', help="Preselection cuts for optimization")
    (option, args) = parser.parse_args()
    preselection = option.preselections
    
    A = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0)
    B = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    C = (0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)
    D = (0.6, 0.8)
    arrays = [A, B, C, D]
    L = list(itertools.product(*arrays))
    
    pool = Pool(24)

    output = pool.map(optimization, L)
    output_sorted = sorted(output, key=SortBySrootB)

    print output_sorted
    G = TGraph()
    G.GetXaxis().SetTitle("X Mass ordered by optimization")
    G.GetYaxis().SetTitle("s/sqrt{b}")
    for i in range(len(output_sorted)):
        G.SetPoint(G.GetN(), i, float(output_sorted[i][0]))
    C = TCanvas()
    C.cd()
    G.Draw("AP")
    C.Print("optimization_plots/Optimiziation_X1000a50_gxsd330p1_hlt_JJ300200m30.root")


