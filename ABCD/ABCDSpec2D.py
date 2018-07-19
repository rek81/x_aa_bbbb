#2D Spectator
import ROOT
import os
from ROOT import *
import sys
import math
from array import array
import functools
import scipy
import random
import string
import ABCD
from ABCD import *

FitFile = TFile("2D_Test_v1.root")
Fit = FitFile.Get("pol22/FitWithErrors")

def Fill2DHist(F, N):


	MSR = TH1F("msr_"+N, ";Average a mass (GeV);events", 60, 0, 300)
	EST = TH1F("est_"+N, ";Average a mass (GeV);events", 60, 0, 300)
	EST.GetYaxis().SetTitleOffset(1.05)
	EST_S = TH1F("stat_"+N, ";Average a mass (GeV);events", 60, 0, 300)
	EST_E = TH1F("err"+N, ";Average a mass (GeV);events", 60, 0, 300)
	MSR.SetStats(0)
	EST.SetStats(0)


	T = F.Get("tree")
	n = T.GetEntries()
	for j in range(0, n):
		T.GetEntry(j)
		if T.doublecsv_1 > 0.6 and math.fabs(T.SD_0-T.SD_1)/(T.SD_0+T.SD_1) < 0.15:
			if T.weight<15 and T.SD_0>10. and T.SD_1>10. and T.Tau21_1<0.55 and T.Tau21_0<0.55 and T.Tau32_1>0.57 and T.Tau32_0>0.57 and T.pT_0>475 and math.fabs(T.eta_0-T.eta_1)<2.0:
				AVG = (T.SD_0+T.SD_1)/2.
				X = T.PhiM/1000.
				W = T.weight*T.puW
				if T.doublecsv_0 > 0.6:
					MSR.Fill(AVG, W)
				else:
					b = Fit.FindBin(AVG, X)
					fitW = Fit.GetBinContent(b)
					fitWE = Fit.GetBinError(b)
					EST.Fill(AVG, W*fitW)
					EST_S.Fill(AVG, W)
					EST_E.Fill(AVG, W*(fitW+fitWE))
	FindAndSetMax(MSR, EST_E, EST)
	return (MSR, EST, EST_S, EST_E)


F_qcd = TFile("QCD.root")
QCD = Fill2DHist(F_qcd, "QCD")
Pull = QCD[0].Clone("pull")

TBoxes = []
PBoxes = []
for a in range(QCD[0].GetXaxis().GetNbins()):
	#errors
	BinLow = QCD[0].GetXaxis().GetBinLowEdge(a+1)
	BinHigh = QCD[0].GetXaxis().GetBinUpEdge(a+1)
	E = QCD[3].GetBinContent(a+1) - QCD[1].GetBinContent(a+1)
	if QCD[2].GetBinContent(a+1) > 0:
		E = math.sqrt(E*2 + (QCD[2].GetBinError(a+1)/QCD[2].GetBinContent(a+1))**2)
	teBox = TBox(BinLow, max(0,QCD[1].GetBinContent(a+1)-E) , BinHigh, QCD[1].GetBinContent(a+1)+E)
	teBox.SetFillColor(kAzure-4)
	teBox.SetFillStyle(3324)
	TBoxes.append(teBox)
	#pulls
	A = QCD[0].GetBinContent(a+1)
	AE = QCD[0].GetBinError(a+1)
	if AE < 1.: AE = 1.4
	p = (A-QCD[1].GetBinContent(a+1))/AE
	Pull.SetBinContent(a+1, p)
	Pull.SetBinError(a+1, 1.)
	teBox = TBox(BinLow, -E/AE , BinHigh, E/AE)
	teBox.SetFillColor(kAzure-4)
	teBox.SetFillStyle(3324)
	PBoxes.append(teBox)


GoodPlotFormat(QCD[0], "e0markers", kBlack, 20)
GoodPlotFormat(QCD[1], "thickline", kAzure, 20)


Pull.GetXaxis().SetNdivisions(0)
Pull.GetYaxis().SetNdivisions(6)
Pull.GetYaxis().SetTitle("#frac{Msr - Est}{#sigma}")
Pull.GetYaxis().CenterTitle(True)
Pull.GetYaxis().SetLabelSize(85/22*Pull.GetYaxis().GetLabelSize())
Pull.GetYaxis().SetTitleSize(4.4*Pull.GetYaxis().GetTitleSize())
Pull.GetYaxis().SetTitleOffset(0.175)
Pull.GetYaxis().SetRangeUser(-3., 3.)
Pull.SetLineColor(kGray+3)


GoodPlotFormat(Pull, "e0markers", kBlack, 20)
CenterLine = TLine(0, 0., 300, 0.)
CenterLine.SetLineColor(kRed)
CenterLine.SetLineStyle(1)
CenterLine.SetLineWidth(1)
Linem1 = TLine(0, -1., 300, -1.)
Linep1 = TLine(0, 1., 300, 1.)
Linem2 = TLine(0, -2., 300, -2.)
Linep2 = TLine(0, 2., 300, 2.)
for l in [Linep2,Linem2,Linep1,Linem1]:
	l.SetLineColor(kPink-4)
	l.SetLineWidth(1)



C = TCanvas()
plot = TPad("plot", "The pad 80% of the height",0,0.15,1,1)
pull = TPad("pull", "The pad 20% of the height",0,0,1.0,0.15)
plot.Draw()
pull.Draw()
plot.cd()
QCD[1].Draw("hist")
for b in TBoxes: b.Draw("same")
QCD[0].Draw("e0same")
pull.cd()
Pull.Draw("e")
for l in [Linep2,Linem2,Linep1,Linem1,CenterLine]:
	l.Draw("same")
for b in PBoxes:
	b.Draw("same")