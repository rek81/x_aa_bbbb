import os
import ROOT
from ROOT import *
from array import array
import math
from math import *
import sys

def getPUPPIweight( puppipt, puppieta):
	genCorr  = 1.
	recoCorr = 1.
	totalWeight = 1.
	
	corrGEN = ROOT.TF1("corrGEN","[0]+[1]*pow(x*[2],-[3])",200,3500)
	corrGEN.SetParameter(0,1.00626)
	corrGEN.SetParameter(1, -1.06161)
	corrGEN.SetParameter(2,0.0799900)
	corrGEN.SetParameter(3,1.20454)
	
	corrRECO_cen = ROOT.TF1("corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",200,3500)
	corrRECO_cen.SetParameter(0,1.09302)
	corrRECO_cen.SetParameter(1,-0.000150068)
	corrRECO_cen.SetParameter(2,3.44866e-07)
	corrRECO_cen.SetParameter(3,-2.68100e-10)
	corrRECO_cen.SetParameter(4,8.67440e-14)
	corrRECO_cen.SetParameter(5,-1.00114e-17)
	
	corrRECO_for = ROOT.TF1("corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",200,3500)
	corrRECO_for.SetParameter(0,1.27212)
	corrRECO_for.SetParameter(1,-0.000571640)
	corrRECO_for.SetParameter(2,8.37289e-07)
	corrRECO_for.SetParameter(3,-5.20433e-10)
	corrRECO_for.SetParameter(4,1.45375e-13)
	corrRECO_for.SetParameter(5,-1.50389e-17)
	genCorr =  corrGEN.Eval( puppipt )
	if( abs(puppieta)  < 1.3 ):
		recoCorr = corrRECO_cen.Eval( puppipt )
	else:
		recoCorr = corrRECO_for.Eval( puppipt );
	totalWeight = genCorr*recoCorr
	#print totalWeight
	return totalWeight

class MicroTree:
	def __init__(self, name, File,TREE):
		self.name = name
		self.File = File
		self.__book__()
		self.FillMicroTree(TREE)
	def __book__(self):
		self.f = ROOT.TFile( self.name + ".root", "recreate" )
		self.f.cd()
		self.tree = ROOT.TTree("tree", "tree")
		self.weight = array('f', [0.0])
		self.addBranch('weight', self.weight, self.tree)
		self.puW = array('f', [0.0])
		self.addBranch('puW', self.puW, self.tree)
		self.kf = array('f', [0.0])
		self.addBranch('kf', self.kf, self.tree)
		self.kfNLO = array('f', [0.0])
		self.addBranch('kfNLO', self.kfNLO, self.tree)
		self.pT_0 = array('f', [-1.0])
		self.addBranch('pT_0', self.pT_0, self.tree)
		self.pT_1 = array('f', [-1.0])
		self.addBranch('pT_1', self.pT_1, self.tree)
		self.Tau21_0 = array('f', [-1.0])
		self.addBranch('Tau21_0', self.Tau21_0, self.tree)
		self.Tau21_1 = array('f', [-1.0])
		self.addBranch('Tau21_1', self.Tau21_1, self.tree)
		self.Tau32_0 = array('f', [-1.0])
		self.addBranch('Tau32_0', self.Tau32_0, self.tree)
		self.Tau32_1 = array('f', [-1.0])
		self.addBranch('Tau32_1', self.Tau32_1, self.tree)
		self.SD_0 = array('f', [-1.0])
		self.addBranch('SD_0', self.SD_0, self.tree)
		self.SD_1 = array('f', [-1.0])
		self.addBranch('SD_1', self.SD_1, self.tree)
		self.eta_0 = array('f', [100.0])
		self.addBranch('eta_0', self.eta_0, self.tree)
		self.eta_1 = array('f', [100.0])
		self.addBranch('eta_1', self.eta_1, self.tree)
		self.phi_0 = array('f', [100.0])
		self.addBranch('phi_0', self.phi_0, self.tree)
		self.phi_1 = array('f', [100.0])
		self.addBranch('phi_1', self.phi_1, self.tree)
		self.doublecsv_0 = array('f', [-5.0])
		self.addBranch("doublecsv_0", self.doublecsv_0, self.tree)
		self.doublecsv_1 = array('f', [-5.0])
		self.addBranch("doublecsv_1", self.doublecsv_1, self.tree)
		self.minsubcsv_0 = array('f', [-5.0])
		self.addBranch("minsubcsv_0", self.minsubcsv_0, self.tree)
		self.minsubcsv_1 = array('f', [-5.0])
		self.addBranch("minsubcsv_1", self.minsubcsv_1, self.tree)
		self.puppet = array('f', [-1.0])
		self.addBranch('puppet', self.puppet, self.tree)
		self.pfmet = array('f', [-1.0])
		self.addBranch('pfmet', self.pfmet, self.tree)
		self.csv_0 = array('f', [-5.0])
		self.addBranch("csv_0", self.csv_0, self.tree)
		self.csv_1 = array('f', [-5.0])
		self.addBranch("csv_1", self.csv_1, self.tree)
		self.aM = array('f', [-5.0])
		self.addBranch("aM", self.aM, self.tree)
		self.Masym = array('f', [-5.0])
		self.addBranch("Masym", self.Masym, self.tree)
		self.XM = array('f', [-5.0])
		self.addBranch("XM", self.XM, self.tree)
		self.XBoost = array('f', [-5.0])
		self.addBranch("XBoost", self.XBoost, self.tree)
		self.dR = array('f', [-5.0])
		self.addBranch("dR", self.dR, self.tree)
		self.dEta = array('f', [-5.0])
		self.addBranch("dEta", self.dEta, self.tree)
		self.dPhi = array('f', [-5.0])
		self.addBranch("dPhi", self.dPhi, self.tree)
		# VARS:
	def FillMicroTree(self, TREE):
		File = TFile(self.File)
		print " "
		print " "
		print "---- Starting:" + self.name
		print "Path = " + self.File
				#print "Nevents at production: " + str(File.Get("NEvents").GetEntries())
		Tree = File.Get(TREE)
		n = Tree.GetEntries()
		for j in range(0, n): # Here is where we loop over all events.
			if j % 5000 == 0 or j == 1:
				percentDone = float(j) / float(n) * 100.0
				print 'Processing '+self.name+' {0:10.0f}/{1:10.0f} : {2:5.2f} %'.format(j, n, percentDone )
			Tree.GetEntry(j)
			if Tree.neleLoose == 0 and Tree.nmuLoose == 0 and Tree.ntau==0 and Tree.pfmet < 100.:
				TheaW0 = getPUPPIweight(Tree.AK8Puppijet0_pt, Tree.AK8Puppijet0_eta)
				TheaW1 = getPUPPIweight(Tree.AK8Puppijet1_pt, Tree.AK8Puppijet1_eta)
				self.pT_0[0] = Tree.AK8Puppijet0_pt
				self.pT_1[0] = Tree.AK8Puppijet1_pt
				self.Tau21_0[0] = Tree.AK8Puppijet0_tau21
				self.Tau21_1[0] = Tree.AK8Puppijet1_tau21
				self.Tau32_0[0] = Tree.AK8Puppijet0_tau32
				self.Tau32_1[0] = Tree.AK8Puppijet1_tau32
				self.SD_0[0] =Tree.AK8Puppijet0_msd * TheaW0
				self.SD_1[0] =Tree.AK8Puppijet1_msd * TheaW1
				self.doublecsv_0[0] = Tree.AK8Puppijet0_doublecsv
				self.doublecsv_1[0] = Tree.AK8Puppijet1_doublecsv
				self.minsubcsv_0[0] = Tree.AK8Puppijet0_minsubcsv
				self.minsubcsv_1[0] = Tree.AK8Puppijet1_minsubcsv
				self.eta_0[0] = Tree.AK8Puppijet0_eta
				self.eta_1[0] = Tree.AK8Puppijet1_eta
				self.phi_0[0] = Tree.AK8Puppijet0_phi
				self.phi_1[0] = Tree.AK8Puppijet1_phi
				self.weight[0] = Tree.scale1fb
				self.puW[0] = Tree.puWeight
				self.kf[0] = Tree.kfactor
				self.kfNLO[0] = Tree.kfactorNLO
				self.puppet[0] = Tree.puppet
				self.pfmet[0] = Tree.pfmet
				self.csv_0[0] = Tree.AK8Puppijet0_csv
				self.csv_1[0] = Tree.AK8Puppijet1_csv
				J1 = TLorentzVector()
				J2 = TLorentzVector()
				J1.SetPtEtaPhiM(Tree.AK8Puppijet0_pt, Tree.AK8Puppijet0_eta, Tree.AK8Puppijet0_phi, Tree.AK8Puppijet0_mass*TheaW0)
				J2.SetPtEtaPhiM(Tree.AK8Puppijet1_pt, Tree.AK8Puppijet1_eta, Tree.AK8Puppijet1_phi, Tree.AK8Puppijet1_mass*TheaW1)
				self.XM[0] = (J1 + J2).M()
				self.aM[0] = (self.SD_0[0] + self.SD_1[0])/2
				self.Masym[0] = math.fabs(self.SD_0[0] - self.SD_1[0])/ (self.SD_0[0] + self.SD_1[0])
				self.XBoost[0] = (J1 + J2).Pt()
				self.dR[0] = J1.DeltaR(J2)
				self.dEta[0] = math.fabs(J1.Eta() - J2.Eta())
				self.dPhi[0] = J1.DeltaPhi(J2)
				self.tree.Fill()
		self.f.cd()
		self.f.Write()
		self.f.Close()
	def addBranch(self, name, var, T):
		T.Branch(name, var, name+'/F')
	def __del__(self):
				print "done!"
for y in ["650", "800", "1200"]:
	for x in ["25", "50", "125"]:
		if y == "1200" and x == "125": continue
		Sig = MicroTree("Xaa_"+y+"_"+x, "/uscms/home/rkowalsk/nobackup/CMSSW_8_0_20/src/BB_BB_optimization/signal_files/y"+y+"x"+x+".root", "Events")
#TT = MicroTree("ttbar", "TTBAR.root", "otree")
#QCD = MicroTree("qcd", "QCD.root", "otree")
#JetHT = MicroTree("NewJetHT", "FullJETHT.root", "Events")
