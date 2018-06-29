import ROOT
from ROOT import *
import os
import sys
import math
from array import array
import functools
import scipy
import random
import string

def Normalize(H):
	Int = H.Integral()
	print "Integra: ", Int
	H.Scale(1./Int)
	return Int

def quickplot(File, tree, plot, var, Cut, Weight):
        temp = plot.Clone("temp")
        chain = ROOT.TChain(tree)
        chain.Add(File)
        chain.Draw(var+">>"+"temp", "("+Weight+")*("+Cut+")", "goff")
        plot.Add(temp)

def quick2dplot(File, tree, plot, var, var2, Cut, Weight):
        temp = plot.Clone("temp")
        chain = ROOT.TChain(tree)
        chain.Add(File)
        chain.Draw(var2+":"+var+">>"+"temp", "("+Weight+")*("+Cut+")", "goff")
	plot.GetYaxis().SetTitleOffset(1.45)
        plot.Add(temp)

def FindAndSetMax(*args):
	if len(args) == 1: args = args[0]
        maximum = 0.0
        for i in args:
                i.SetStats(0)
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in args:
                j.GetYaxis().SetRangeUser(0,maximum*1.35)
	return maximum*1.35


def HistSortBySize(H):
	return H.Integral()

def GoodPlotFormat(H, *args):
	try: H.SetStats(0)
	except: print " ------------ [  No stats box found!  ]"
	if args[0] == 'thickline':
		H.SetLineColor(args[1])
		H.SetLineWidth(2)
	if args[0] == 'thinline':
		H.SetLineColor(args[1])
		H.SetLineWidth(1)
	if args[0] == 'fill':
		H.SetLineColor(args[1])
		H.SetFillColor(args[1])
		H.SetFillStyle(args[2])
	if args[0] == 'e0markers':
		H.SetLineColor(args[1])
		H.SetMarkerColor(args[1])
		H.SetMarkerStyle(args[2])
	H.GetXaxis().SetTitleSize(0.04)

def AddCMSLumi(pad, fb, extra):

	cmsText     = "CMS";
	cmsTextFont   = 61  

	extraText   = extra
	extraTextFont = 52 

	lumiTextSize     = 0.4
	lumiTextOffset   = 0.15

	cmsTextSize      = 0.5
	cmsTextOffset    = 0.15


	H = pad.GetWh()
	W = pad.GetWw()
	l = pad.GetLeftMargin()
	t = pad.GetTopMargin()
	r = pad.GetRightMargin()
	b = pad.GetBottomMargin()
	e = 0.025

	pad.cd()

	lumiText = str(fb)+" fb^{-1} (13 TeV)"

	latex = TLatex()
	latex.SetNDC()
	latex.SetTextAngle(0)
	latex.SetTextColor(kBlack)	
	
	extraTextSize = 0.76*cmsTextSize
	
	latex.SetTextFont(42)
	latex.SetTextAlign(31) 
	latex.SetTextSize(lumiTextSize*t)	

	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)

	pad.cd()

	latex.SetTextFont(cmsTextFont)
	latex.SetTextSize(cmsTextSize*t)
	latex.SetTextAlign(11)
	latex.DrawLatex(l, 1-t+cmsTextOffset*t, cmsText)
	latex.SetTextFont(extraTextFont)
	latex.SetTextAlign(11)
	latex.SetTextSize(extraTextSize*t)
	latex.DrawLatex(l+0.11, 1-t+cmsTextOffset*t, extraText)
 

	pad.Update()

def GetQuantileProfiles(Th2f, cut, name):
		q1 = []
		nxbins = Th2f.GetXaxis().GetNbins();
		xlo = Th2f.GetXaxis().GetBinLowEdge(1);
		xhi = Th2f.GetXaxis().GetBinUpEdge(Th2f.GetXaxis().GetNbins() );
		for i in range(nxbins):
				H = Th2f.ProjectionY("ProjY"+str(i),i+1,i+1)
				probSum = array('d', [cut])
				q = array('d', [0.0]*len(probSum))
				H.GetQuantiles(len(probSum), q, probSum)
				q1.append(q[0])
		H1 = TH1F(name, "", nxbins,xlo,xhi)
		for i in range(nxbins):
				H1.SetBinContent(i+1,q1[i])
		return H1


class VARIABLE:
	def __init__(self, var, bins, name):
		self.name = name
		self.var = var
		self.bins = bins
		self.nbins = len(bins)-1

class DISTRIBUTION:
	def __init__(self, name, treename, weight, color, printname, isData):
		self.name = name
		self.treename = treename
		self.weight = weight
		self.color = color
		self.printname = printname
		self.MC = not isData

class ABCDEST:
	def __init__(self, printname, BKG, SUBBKG, SIG, var1, var2, varf, presel, lumi, extra):
		gROOT.ProcessLine( "gErrorIgnoreLevel = 3001;")
		print "\n\n\n\033[0;37;44m-- == -- == -- == -- == -- == -- == -- == -- == -- == -- == --"+ "\033[0m"
		print "\033[0;37;44m-- == -- == -- =         ABCD METHOD          = -- == -- == --"+ "\033[0m"
		print "\033[0;37;44m-- == -- == -- == -- == -- == -- == -- == -- == -- == -- == --"+ "\033[0m\n"
		self.printname = printname
		self.BKG = BKG
		self.SUBBKG = SUBBKG
		self.SIG = SIG
		self.varf = varf
		self.var1 = var1
		self.var2 = var2
		self.presel = presel
		self.lumi = lumi
		self.extra = extra
		print "Creating framework with: "
		print "	  " + var1.name + "	\033[0;31;40m("+var1.var+")"+ "\033[0m"
		print "		   & "
		print "	  " + var2.name + "	\033[0;31;40m("+var2.var+")"+ "\033[0m\n"
		print 
		print "Preselection   =   \033[0;31;40m" + presel + "\033[0m\n"

		try:
			os.stat(printname)
		except:
			os.mkdir(printname)  

		print "Backgrounds (Main): "
		for b in BKG:
			print "	  \033[0;32;40m" + b.name+ "\033[0m"
		print "Backgrounds (Subtracted): "
		for b in SUBBKG:
			print "	  \033[0;32;40m" + b.name+ "\033[0m"
		print "Signals: "
		for s in SIG:
			print "	  \033[0;32;40m" + s.name+ "\033[0m\n"

		print "Saving output to: \033[0;37;40m" + os.path.dirname(os.path.realpath(__file__)) + "/"+printname + "\033[0m\n"

		plot = self.MakeSimple2DPlots(BKG, var1, var2, "bkg")
		sigplot = self.MakeSimple2DPlots(SIG, var1, var2, "sig")
		subplot = self.MakeSimple2DPlots(SUBBKG, var1, var2, "subbkg")
		plot.Add(subplot, -1.)
		C1 = TCanvas("C_Plane", "", 500, 500)
		C1.cd()
		plot.Draw("col")
		AddCMSLumi(gPad, self.lumi, self.extra)
		C1.Print(self.printname+"/"+"ABCD_regions_bkg.png")
		C1.Print(self.printname+"/"+"ABCD_regions_bkg.root")
		ProfsM5 = []		
		ProfsM10 = []
		for i in [9,8,7,6,5,4,3,2,1]:
			ProfsM10.append(GetQuantileProfiles(plot, 0.1*i, "PlaneProf_"+str(i)))
		for i in [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]:
			ProfsM5.append(GetQuantileProfiles(plot, 0.1*i, "PlaneProf_"+str(i)))
		C2 = TCanvas("C2_Plane", "", 500, 500)
		C2.cd()
		plot.SetMarkerColor(kWhite)
		plot.Draw("")
		for i in ProfsM10:
			i.SetLineColor(kBlue)
			i.SetLineStyle(2)
			i.Draw("same")
		for i in ProfsM5:
			i.SetLineColor(kBlue)
			i.SetLineStyle(3)
			i.Draw("same")
		AddCMSLumi(gPad, self.lumi, self.extra)
		C2.Print(self.printname+"/"+"ABCD_quantiles_bkg.png")
		C2.Print(self.printname+"/"+"ABCD_quantiles_bkg.root")
		C3 = TCanvas("C3_Plane", "", 500, 500)
		C3.cd()
		sigplot.Draw("col")
		AddCMSLumi(gPad, self.lumi, self.extra)
		C3.Print(self.printname+"/"+"ABCD_regions_sig.png")
		C3.Print(self.printname+"/"+"ABCD_regions_sig.root")
		del plot, sigplot, subplot
		del C1, C2, C3

	def MakeSimple2DPlots(self, DIST, var1, var2, name):
		Plane = TH2F("Plane"+name, ";"+var1.name+";"+var2.name, var1.nbins, scipy.array(var1.bins), var2.nbins, scipy.array(var2.bins))
		Plane.SetStats(0)
		for d in DIST:
			tmp = TH2F("Plane"+name+d.name, ";"+var1.name+";"+var2.name, var1.nbins, scipy.array(var1.bins), var2.nbins, scipy.array(var2.bins))
			chn = TChain(d.treename)
			chn.Add(d.name)
			chn.Draw(var2.var+":"+var1.var+">>Plane"+name+d.name, "("+str(self.lumi)+"*"+d.weight+")*("+self.presel+")", "goff")
			Plane.Add(tmp)
		return Plane

	def Estimate(self, name, C1, C2, FitFunc):
		self.C2 = C2
		self.C1 = C1
		self.FitFunc = FitFunc
		print "\n\n\n\033[0;37;44m-- == -- == -- == -- == -- == -- == -- == -- == -- == -- == --"+ "\033[0m\n"
		self.CUT1 = ""
		for c in C1:
			self.CUT1 += self.var1.var+c
		self.CUT2 = ""
		for c in C2:
			self.CUT2 += self.var2.var+c
		print "abscissa cut   =   \033[0;31;40m" + self.CUT1 + "\033[0m\n"
		print "ordinate cut   =   \033[0;31;40m" + self.CUT2 + "\033[0m\n"


		self.cutA = self.presel+"&("+self.CUT1+")"+"&("+self.CUT2+")"
		self.cutB = self.presel+"&("+self.CUT1+")"+"&!("+self.CUT2+")"
		self.cutC = self.presel+"&!("+self.CUT1+")"+"&("+self.CUT2+")"
		self.cutD = self.presel+"&!("+self.CUT1+")"+"&!("+self.CUT2+")"


		L = TLegend(0.53,0.62,0.89,0.89)
		L.SetLineColor(0)
		L.SetFillColor(0)


		S = TH1F("Signal", ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
		for s in self.SIG:
			tS = TH1F("TSignal"+s.name, ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
			Schn = TChain(s.treename)
			Schn.Add(s.name)
			if s.MC: 
				W = str(self.lumi)+"*"+s.weight
			else: W = "1."
			Schn.Draw(self.varf.var+">>TSignal"+s.name, "("+W+")*("+self.cutA+")", "goff")
			S.Add(tS)
		GoodPlotFormat(S, "fill", self.SIG[0].color, 3002)
		L.AddEntry(S, self.SIG[0].printname, "F")

		AKeep = []

		print "Filling A regions:"
		A = TH1F("BkgA", ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
		Pull = TH1F("Pull", ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
		GoodPlotFormat(Pull, "e0markers", kBlack, 20)
		self.B = TH1F("BkgB", ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))

		GoodPlotFormat(A, "e0markers", kBlack, 20)
		L.AddEntry(A, self.BKG[0].printname, "PL")
		GoodPlotFormat(self.B, "thickline", self.BKG[0].color)
		L.AddEntry(self.B, "Data-Driven Estimate (B #times C/D)", "L")

		self.R = TH1F("BkgR", ";"+self.varf.name+";C/D rate", self.varf.nbins, scipy.array(self.varf.bins))

		for b in self.BKG:
			randomthing = random.choice(string.letters)
			tA = TH1F("TBkgA"+b.name+randomthing, ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
			Bchn = TChain(b.treename)
			Bchn.Add(b.name)
			if b.MC: 
				W = str(self.lumi)+"*"+b.weight
			else: W = "1."
			Bchn.Draw(self.varf.var+">>TBkgA"+b.name+randomthing, "("+W+")*("+self.cutA+")", "goff")
			A.Add(tA)
		for b in self.SUBBKG:
			randomthing = random.choice(string.letters)
			tA = TH1F("sTBkgA"+b.name+randomthing, ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
			Bchn = TChain(b.treename)
			Bchn.Add(b.name)
			if b.MC: 
				W = str(self.lumi)+"*"+b.weight
			else: W = "1."
			Bchn.Draw(self.varf.var+">>sTBkgA"+b.name+randomthing, "("+W+")*("+self.cutA+")", "goff")
			AKeep.append(tA)
			GoodPlotFormat(tA, "thickline", b.color)
			L.AddEntry(tA, b.printname, "L")

		self.Max = 0.14
		print "Filling B, C and D regions:"
		for b in range(len(self.varf.bins)-1):
			self.FillMassBin(b)

		if self.FitFunc != "none":
			print "Using function:	   \033[0;31;40m"  + self.FitFunc + "\033[0m\n"
			self.F = TF1("fit", self.FitFunc, self.varf.bins[0], self.varf.bins[-1])
			self.F.SetLineColor(kRed)
			self.FitResults = self.R.Fit(self.F)
			self.B.Multiply(self.F)
		else:
			self.B.Multiply(self.R)

		for b in AKeep:
			self.B.Add(b)
		PlotKeep = sorted(AKeep, key=HistSortBySize)
		#PlotKeep.reverse()

		EBox = []
		PEBox = []
		SigBox = []
		for b in range(len(self.varf.bins)-1):
			blow = self.varf.bins[b]
			bhigh = self.varf.bins[b+1]
			c = self.B.GetBinContent(b+1)
			ce = self.B.GetBinError(b+1)
			tBox = TBox(blow, max(0,c - ce) , bhigh, c + ce)
			tBox.SetFillColor(kAzure-4)
			tBox.SetFillStyle(3324)
			EBox.append(tBox)
			a = A.GetBinContent(b+1)
			ae = A.GetBinError(b+1)
			if ae < 1.: ae = 1.4
			p = (a-c)/ae
			Pull.SetBinContent(b+1, p)
			Pull.SetBinError(b+1, 1.)
			teBox = TBox(blow, -ce/ae , bhigh, ce/ae)
			teBox.SetFillColor(kAzure-4)
			teBox.SetFillStyle(3324)
			PEBox.append(teBox)
			if S.GetBinContent(b+1)/ae > 0.25:
				sBox = TBox(blow, 0 , bhigh, S.GetBinContent(b+1)/ae)
				sBox.SetFillColor(self.SIG[0].color)
				sBox.SetFillStyle(3002)
				SigBox.append(sBox)


		CenterLine = TLine(self.varf.bins[0], 0., self.varf.bins[-1], 0.)
		CenterLine.SetLineColor(kRed)
		CenterLine.SetLineStyle(1)
		CenterLine.SetLineWidth(1)
		Linem1 = TLine(self.varf.bins[0], -1., self.varf.bins[-1], -1.)
		Linep1 = TLine(self.varf.bins[0], 1., self.varf.bins[-1], 1.)
		Linem2 = TLine(self.varf.bins[0], -2., self.varf.bins[-1], -2.)
		Linep2 = TLine(self.varf.bins[0], 2., self.varf.bins[-1], 2.)
		for l in [Linep2,Linem2,Linep1,Linem1]:
			l.SetLineColor(kPink-4)
			l.SetLineWidth(1)

		Pull.GetXaxis().SetNdivisions(0)
		Pull.GetYaxis().SetNdivisions(6)
		Pull.GetYaxis().SetTitle("#frac{Msr - Est}{#sigma}")
		Pull.GetYaxis().CenterTitle(True)
		Pull.GetYaxis().SetLabelSize(85/22*Pull.GetYaxis().GetLabelSize())
		Pull.GetYaxis().SetTitleSize(4.4*Pull.GetYaxis().GetTitleSize())
		Pull.GetYaxis().SetTitleOffset(0.175)
		Pull.GetYaxis().SetRangeUser(-3., 3.)
		Pull.SetLineColor(kGray+3)

		L.AddEntry(EBox[0], "Estimate Uncertainty", "F")
		FindAndSetMax(S,A,self.B)

		C = TCanvas("C_"+name, "", 600, 600)
		plot = TPad("plot", "The pad 80% of the height",0,0.15,1,1)
		pull = TPad("pull", "The pad 20% of the height",0,0,1.0,0.15)
		plot.Draw()
		pull.Draw()
		plot.cd()
		self.B.Draw("hist")
		A.Draw("e0same")
		for i in range(len(PlotKeep)):
			PlotKeep[i].Draw("histsame")
			if i == len(PlotKeep)-1: break
			PlotKeep[i+1].Add(PlotKeep[i])

		S.Draw("histsame")
		for b in EBox:
			b.Draw("same")
		L.Draw("same")
		AddCMSLumi(gPad, self.lumi, self.extra)
		pull.cd()
		Pull.Draw("e0")
		for b in SigBox:
			b.Draw("same")
		for l in [Linep2,Linem2,Linep1,Linem1,CenterLine]:
			l.Draw("same")
		for b in PEBox:
			b.Draw("same")
		C.Print(self.printname+"/"+name+"_EST.png")
		C.Print(self.printname+"/"+name+"_EST.root")

		if self.FitFunc != "none":

			self.EG = TGraphErrors(1000)
			for i in range(1000):
				self.EG.SetPoint(i, self.varf.bins[0] + i*(self.varf.bins[-1]- self.varf.bins[0])/1000., 0)
			TVirtualFitter.GetFitter().GetConfidenceIntervals(self.EG)
			Ndof = self.F.GetNDF()
			Chi2 = self.F.GetChisquare()
			GoodPlotFormat(self.EG, "e0markers", kRed-10, 1)
		GoodPlotFormat(self.R, "e0markers", kBlack, 23)

		if self.FitFunc != "none":
			L2 = TLegend(0.13,0.65,0.45,0.8975)
		else:
			L2 = TLegend(0.13,0.75,0.45,0.85)
		L2.SetLineColor(0)
		L2.SetFillColor(0)
		L2.AddEntry(self.R, "C/D", "PL")
		self.R.Draw("e")
		if self.FitFunc != "none":
			L2.SetHeader("#chi^{2}/Ndof = "+ '{:4.2f}'.format(Chi2)+ "/" + str(Ndof))
			L2.AddEntry(self.F, "fit (pol3)", "L")
		self.R.GetYaxis().SetRangeUser(0.,self.Max*1.3)
		C2 = TCanvas("C2_"+name, "", 500, 500)
		C2.cd()
		self.R.Draw("e")
		if self.FitFunc != "none":
			self.EG.Draw("sameP")
			self.F.Draw("sameL")
			self.R.Draw("esame")
		L2.Draw("same")
		AddCMSLumi(gPad, self.lumi, self.extra)
		C2.Print(self.printname+"/"+name+"_FIT.png")
		C2.Print(self.printname+"/"+name+"_FIT.root")

		if self.FitFunc != "none":
			FitF = TFile(self.printname+"/ABCDFit.root", "recreate")
			FitF.cd()
			F1 = self.F.Clone("FIT")
			F1.Write()
			E1 = self.EG.Clone("FIT_ERROR")
			E1.Write()
			R1 = self.R.Clone("BINS")
			R1.Write()
			FitF.Write()
			FitF.Close()
			return [self.F, self.EG]
		else:
			FitF = TFile(self.printname+"/ABCDBins.root", "recreate")
			FitF.cd()
			R1 = self.R.Clone("BINS")
			R1.Write()
			FitF.Write()
			FitF.Close()
			return [self.R]

	def FillMassBin(self, b):
		blow = self.varf.bins[b]
		bhigh = self.varf.bins[b+1]
		massbincut = self.varf.var+"<"+str(bhigh)+"&"+self.varf.var+">"+str(blow)
		tB = TH1F("BkgtB"+str(b), ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
		tC = TH1F("BkgtC"+str(b), ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
		tD = TH1F("BkgtD"+str(b), ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
		for bb in self.BKG:
			randomthing = random.choice(string.letters)
			Bchn = TChain(bb.treename)
			Bchn.Add(bb.name)
			ttB = TH1F("BkgtB"+str(b)+bb.name+randomthing, ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
			ttC = TH1F("BkgtC"+str(b)+bb.name+randomthing, ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
			ttD = TH1F("BkgtD"+str(b)+bb.name+randomthing, ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
			if bb.MC: 
				W = str(self.lumi)+"*"+bb.weight
			else: W = "1."
			Bchn.Draw(self.varf.var+">>BkgtB"+str(b)+bb.name+randomthing, "("+W+")*("+massbincut+"&"+self.cutB+")", "goff")
			Bchn.Draw(self.varf.var+">>BkgtC"+str(b)+bb.name+randomthing, "("+W+")*("+massbincut+"&"+self.cutC+")", "goff")
			Bchn.Draw(self.varf.var+">>BkgtD"+str(b)+bb.name+randomthing, "("+W+")*("+massbincut+"&"+self.cutD+")", "goff")
			tB.Add(ttB)
			tC.Add(ttC)
			tD.Add(ttD)
		for bb in self.SUBBKG:
			randomthing = random.choice(string.letters)
			Bchn = TChain(bb.treename)
			Bchn.Add(bb.name)
			ttB = TH1F("sBkgtB"+str(b)+bb.name+randomthing, ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
			ttC = TH1F("sBkgtC"+str(b)+bb.name+randomthing, ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
			ttD = TH1F("sBkgtD"+str(b)+bb.name+randomthing, ";"+self.varf.name, self.varf.nbins, scipy.array(self.varf.bins))
			if bb.MC: 
				W = str(self.lumi)+"*"+bb.weight
			else: W = "1."
			Bchn.Draw(self.varf.var+">>sBkgtB"+str(b)+bb.name+randomthing, "("+W+")*("+massbincut+"&"+self.cutB+")", "goff")
			Bchn.Draw(self.varf.var+">>sBkgtC"+str(b)+bb.name+randomthing, "("+W+")*("+massbincut+"&"+self.cutC+")", "goff")
			Bchn.Draw(self.varf.var+">>sBkgtD"+str(b)+bb.name+randomthing, "("+W+")*("+massbincut+"&"+self.cutD+")", "goff")
			tB.Add(ttB,-1)
			tC.Add(ttC,-1)
			tD.Add(ttD,-1)
			tC.Sumw2()
			tD.Sumw2()
		if tD.Integral() > 0:
			tC.Divide(tD)
			self.B.Add(tB)
			if tC.Integral() < 2.5:
				if tC.Integral() > self.Max:
					self.Max = tC.Integral()
				self.R.Add(tC)

