# Two-Dimensional version of the ABCD code:
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


class ABCD2D:
	def __init__(self, printname, BKG, SUB, presel, lumi):
		print "\n\n\n\033[0;37;44m-- == -- == -- == -- == -- == -- == -- == -- == -- == -- == --"+ "\033[0m"
		print "\033[0;37;44m-- == -- == -- =      ABCD (2D) METHOD        = -- == -- == --"+ "\033[0m"
		print "\033[0;37;44m-- == -- == -- == -- == -- == -- == -- == -- == -- == -- == --"+ "\033[0m\n"
		self.printname = printname
		self.BKG = BKG
		self.SUB = SUB
		self.presel = presel
		self.lumi = lumi
		print "Backgrounds (main): "
		for b in BKG:
			print "   \033[0;32;40m" + b.name+ "\033[0m"
		print "Backgrounds (to be Subtracted): "
		for b in SUB:
			print "   \033[0;32;40m" + b.name+ "\033[0m"

	def SetFitVariables(self, xaxis, yaxis): # What two variables is our fit evaluated on:
		self.fitX = xaxis
		self.fitY = yaxis
		self.presel += "&("+xaxis.var+">"+str(xaxis.bins[0])+"&"+yaxis.var+">"+str(yaxis.bins[0])+"&"+xaxis.var+"<"+str(xaxis.bins[-1])+"&"+yaxis.var+"<"+str(yaxis.bins[-1])+")"

	def SetABCDVariables(self, xaxis, yaxis, xV, yV): # How to chop up our regions:
		####  my method looks like this:
		####    B C
		####    A D
		####
		self.ABvCD = xaxis
		self.ACvBD = yaxis
		print "\n\n\n\033[0;37;44m-- == -- == -- == -- == -- == -- == -- == -- == -- == -- == --"+ "\033[0m\n"
		self.CUT1 = ""
		for c in xV:
			self.CUT1 += self.ABvCD.var+c
			if c != xV[-1]: self.CUT1 += " & "
		self.CUT2 = ""
		for c in yV:
			self.CUT2 += self.ACvBD.var+c
			if c != yV[-1]: self.CUT2 += " & "
		print "abscissa cut   =   \033[0;31;40m" + self.CUT1 + "\033[0m\n"
		print "ordinate cut   =   \033[0;31;40m" + self.CUT2 + "\033[0m\n"

		self.cutA = self.presel+"&("+self.CUT1+")"+"&("+self.CUT2+")"
		self.cutB = self.presel+"&("+self.CUT1+")"+"&!("+self.CUT2+")"
		self.cutC = self.presel+"&!("+self.CUT1+")"+"&("+self.CUT2+")"
		self.cutD = self.presel+"&!("+self.CUT1+")"+"&!("+self.CUT2+")"

	def Estimate(self, fit):
		print "Preselection   =   \033[0;31;40m" + self.presel + "\033[0m\n"
		self.F = fit
		self.OUTPUT = TFile(self.printname + ".root", "recreate")
		self.OUTPUT.cd()
		print "  \n\n\n\033[0;37;44m-- == -- == -- == -- == -- == -- == -- == -- == -- == -- == --"+ "\033[0m\n"
		self.A= TH2F("A", ";"+self.fitX.name+";"+self.fitY.name, self.fitX.nbins, scipy.array(self.fitX.bins), self.fitY.nbins, scipy.array(self.fitY.bins))
		self.B= TH2F("B", ";"+self.fitX.name+";"+self.fitY.name, self.fitX.nbins, scipy.array(self.fitX.bins), self.fitY.nbins, scipy.array(self.fitY.bins))
		self.C= TH2F("C", ";"+self.fitX.name+";"+self.fitY.name, self.fitX.nbins, scipy.array(self.fitX.bins), self.fitY.nbins, scipy.array(self.fitY.bins))
		self.D= TH2F("D", ";"+self.fitX.name+";"+self.fitY.name, self.fitX.nbins, scipy.array(self.fitX.bins), self.fitY.nbins, scipy.array(self.fitY.bins))
		self.A.SetStats(0)
		self.B.SetStats(0)
		self.C.SetStats(0)
		self.D.SetStats(0)
		for bkg in self.BKG:
			print "Adding " + bkg.printname + " to regions A,B,C and D."
			self.chn = TChain(bkg.treename)
			self.chn.Add(bkg.name)
			if bkg.MC: 
				W = str(self.lumi)+"*"+bkg.weight
			else: W = "1.0"
			self.ADD_HIST(self.A, "("+W+")*("+self.cutA+")")
			self.ADD_HIST(self.B, "("+W+")*("+self.cutB+")")
			self.ADD_HIST(self.C, "("+W+")*("+self.cutC+")")
			self.ADD_HIST(self.D, "("+W+")*("+self.cutD+")")
		for bkg in self.SUB:
			print "Removing " + bkg.printname + " contribution from regions B,C and D."
			chn = TChain(bkg.treename)
			chn.Add(bkg.name)
			if bkg.MC: 
				W = str(self.lumi)+"*"+bkg.weight
			else: W = "1.0"
			self.ADD_HIST(self.B, "-1*("+W+")*("+self.cutB+")")
			self.ADD_HIST(self.C, "-1*("+W+")*("+self.cutC+")")
			self.ADD_HIST(self.D, "-1*("+W+")*("+self.cutD+")")
		self.CoverD = self.C.Clone("CoverD")
		self.CoverD.Divide(self.D)
		EstNoFit = self.B.Clone("EST_fit_NONE")
		EstNoFit.Multiply(self.CoverD)
		Ex = EstNoFit.ProjectionX("Eproj_NOFIT"+self.fitX.name)
		Ey = EstNoFit.ProjectionY("Eproj_NOFIT"+self.fitY.name)
		GoodPlotFormat(Ex, "thickline", self.BKG[0].color)
		GoodPlotFormat(Ey, "thickline", self.BKG[0].color)
		self.A_x = self.A.ProjectionX("Aproj_"+self.fitX.name)
		GoodPlotFormat(self.A_x, "e0markers", kBlack, 20)
		self.A_y = self.A.ProjectionY("Aproj_"+self.fitY.name)
		GoodPlotFormat(self.A_y, "e0markers", kBlack, 20)

		Cx = TCanvas("Check_NOFIT"+self.fitX.name, "", 800, 800)
		Cx.cd()
		self.A_x.Draw("e")
		Ex.Draw("samehist")
		Cx.Write()
		Cy = TCanvas("Check_NOFIT"+self.fitY.name, "", 800, 800)
		Cy.cd()
		self.A_y.Draw("e")
		Ey.Draw("samehist")
		Cy.Write()


		for fit in self.F:
			FOLDER = self.OUTPUT.mkdir(fit[1])
			FOLDER.cd()
			print "Using function:	   \033[0;31;40m"  + fit[1] + "\033[0m\n"
			self.FitResults = self.CoverD.Fit(fit[0])
			F = fit[0].Clone("CoverD_fit")
			F.Write()
			Ndof = fit[0].GetNDF()
			Chi2 = fit[0].GetChisquare()
			print "chi2/Ndof = "+ '{:4.2f}'.format(Chi2)+ "/" + str(Ndof) + " = " + str(Chi2/Ndof)
			# GET THE ERRORS:
			FWE = TH2F("FitWithErrors", "", 250, self.fitX.bins[0], self.fitX.bins[-1], 250, self.fitY.bins[0], self.fitY.bins[-1])
			EG = TGraph2DErrors(1)
			for i in range(250):
				for j in range(250):
					x = self.fitX.bins[0] + i*(self.fitX.bins[-1]- self.fitX.bins[0])/250.
					y = self.fitY.bins[0] + j*(self.fitY.bins[-1]- self.fitY.bins[0])/250.
					EG.SetPoint(0, x, y, 0)
					TVirtualFitter.GetFitter().GetConfidenceIntervals(EG)
					ze = EG.GetErrorZ(0)
					z = F.Eval(x,y)
					b = FWE.FindBin(x,y)
					FWE.SetBinContent(b, z)
					FWE.SetBinError(b, ze)
					print str(x) + ", " + str(y) + " ---> " + str(z) + "+/-" + str(ze)
			FWE.SetStats(0)
			ESTFit = self.B.Clone("EST_fit")
			FWE.Write()
			ESTFit.Multiply(fit[0])
			ESTFit.Write()
			Ex = ESTFit.ProjectionX("Eproj_"+self.fitX.name)
			Ey = ESTFit.ProjectionY("Eproj_"+self.fitY.name)
			GoodPlotFormat(Ex, "thickline", self.BKG[0].color)
			GoodPlotFormat(Ey, "thickline", self.BKG[0].color)
			Cx = TCanvas("Check_"+self.fitX.name, "", 800, 800)
			Cx.cd()
			self.A_x.Draw("e")
			Ex.Draw("samehist")
			Cx.Write()
			Cy = TCanvas("Check_"+self.fitY.name, "", 800, 800)
			Cy.cd()
			self.A_y.Draw("e")
			Ey.Draw("samehist")
			Cy.Write()



		self.OUTPUT.Write()
		self.OUTPUT.Save()
		self.OUTPUT.Close()

	def ADD_HIST(self, H, WC):
		tmp = H.Clone("tmp")
		self.chn.Draw(self.fitY.var+":"+self.fitX.var+">>"+H.GetName(), WC, "goff")
		H.Add(tmp)
		del tmp

if __name__ == '__main__':


	XBins = []
	for i in range(9):
		XBins.append(0.5 + 0.5*i)
	aBins = []
	for i in range(31):
		aBins.append(10.*i)
	DetaBins = [0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5]
	MasymBins = [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
	BBBins = []
	for i in range(11):
		BBBins.append(i*0.1)


	deta = VARIABLE("abs(eta_0-eta_1)", DetaBins, "#Delta#eta")
	masym = VARIABLE("abs(SD_0-SD_1)/(SD_0+SD_1)", MasymBins, "Mass Asymmetry")
	XM = VARIABLE("PhiM/1000.", XBins, "X Mass (TeV)")
	aM = VARIABLE("(SD_0+SD_1)/2", aBins, "Average Mass (GeV)")
	bb1 = VARIABLE("doublecsv_0", BBBins, "Jet 1 double b tag")
	bb2 = VARIABLE("doublecsv_1", BBBins, "Jet 2 double b tag")

	presel = "weight<15&SD_0>10.&SD_1>10.&Tau21_1<0.55&Tau21_0<0.55&Tau32_1>0.57&Tau32_0>0.57&pT_0>475&abs(eta_0-eta_1)<2.0&doublecsv_1 > 0.6"


	qcd = DISTRIBUTION("QCD.root", "tree", "weight*puW", kAzure, "QCD", False)
	ttbar = DISTRIBUTION("ttbar.root", "tree", "weight*puW", kViolet, "t#bar{t}", False)

	Poly33 = "[0] + [1]*x + [2]*y + [3]*x*y + [4]*x*x + [5]*y*y + [6]*x*x*y + [7]*x*y*y + [8]*x*x*x + [9]*y*y*y"
	fit33 = TF2("PolyFit33", Poly33, aBins[0], aBins[-1], XBins[0], XBins[-1])

	Poly22 = "[0] + [1]*x + [2]*y + [3]*x*y + [4]*x*x + [5]*y*y"
	fit22 = TF2("PolyFit22", Poly22, aBins[0], aBins[-1], XBins[0], XBins[-1])

	ABCD = ABCD2D("2D_Test_v1", [qcd], [], presel, 10.0)
	ABCD.SetFitVariables(aM, XM)
	ABCD.SetABCDVariables(masym, bb1, ["< 0.15"], ["> 0.6"])
	ABCD.Estimate([[fit33, "pol33"], [fit22, "pol22"]])