#!/usr/bin/env python

# This little tool is an easy to use check tool to see if new values are correct.
# Same usage as SusyAna just put my. in front of the value
# have fun Klaas Padeken

import ROOT as ro
from math import pi, fabs, cos, sqrt,fmod
from array import array
from ROOT import gDirectory,TFile,TGraph, TH1I , gROOT , TCanvas, TPad, TGraph2D, TH1D, TF1, gStyle, TClonesArray, TParticle, gSystem,ProcInfo_t,TLorentzVector
#from histos import *


allHists=[]
h1_truthpdgid=TH1D("h1_truthpdgid","h1_truthpdgid",100,-50,50);allHists.append(h1_truthpdgid)
h1_jetn=TH1D("h1_jetn","h1_jetn",20,0,20);allHists.append(h1_jetn)
h1_muon=TH1D("h1_muon","h1_muon",10,0,10);allHists.append(h1_muon)
h1_elen=TH1D("h1_elen","h1_elen",10,0,10);allHists.append(h1_elen)
h1_muo_pt=TH1D("h1_muo_pt","h1_muo_pt",100,0,100);allHists.append(h1_muo_pt)
h1_muo_minv=TH1D("h1_muo_minv","h1_muo_minv",100,0,100);allHists.append(h1_muo_minv)
h1_muo_RelTrkIso=TH1D("h1_muo_RelTrkIso","h1_muo_RelTrkIso",100,0,10);allHists.append(h1_muo_RelTrkIso)
h1_muo_TrkIso=TH1D("h1_muo_TrkIso","h1_muo_TrkIso",100,0,10);allHists.append(h1_muo_TrkIso)
h1_muo_ECalIso=TH1D("h1_muo_ECalIso","h1_muo_ECalIso",100,0,10);allHists.append(h1_muo_ECalIso)
h1_muo_HCalIso=TH1D("h1_muo_HCalIso","h1_muo_HCalIso",100,0,10);allHists.append(h1_muo_HCalIso)
h1_muo_TrkIsoDep=TH1D("h1_muo_TrkIsoDep","h1_muo_TrkIsoDep",100,0,10);allHists.append(h1_muo_TrkIsoDep)
h1_muo_ECalIsoDep=TH1D("h1_muo_ECalIsoDep","h1_muo_ECalIsoDep",100,0,10);allHists.append(h1_muo_ECalIsoDep)
h1_muo_HCalIsoDep=TH1D("h1_muo_HCalIsoDep","h1_muo_HCalIsoDep",100,0,10);allHists.append(h1_muo_HCalIsoDep)
h1_muo_AllIso=TH1D("h1_muo_AllIso","h1_muo_AllIso",100,0,10);allHists.append(h1_muo_AllIso)
h1_met_et3=TH1D("h1_met_et3","h1_met_et3",20,0,20);allHists.append(h1_met_et3)
h1_met_et5=TH1D("h1_met_et5","h1_met_et5",20,0,20);allHists.append(h1_met_et5)
h1_met_et8=TH1D("h1_met_et8","h1_met_et8",20,0,20);allHists.append(h1_met_et8)

myfile = TFile('out.root')
myfile.cd("ACSkimAnalysis")
my = gDirectory.Get('allData')
entries = my.GetEntriesFast()
info= ProcInfo_t()
print( "ACAna: Running over %i events"%(entries))

for jentry in range(entries):
    ientry = my.LoadTree(jentry)
    if ientry < 0:
        break
    nb = my.GetEntry(jentry)
    if (jentry%10000)==0:
        print( "Event %i" %(jentry))

    if (jentry%50000)==0:
        gSystem.GetProcInfo(info);
        print(r" -> Memory in MB : %.2f (resident)  %.2f (virtual)"%( info.fMemResident/1000.,info.fMemVirtual/1000.))
    if jentry==50000:
        break
    
    muoN=0
    eleN=0
    for i in range(my.truthl_n):
        if fabs(my.truthl_pdgid[i])==13:
            muoN+=1
        if fabs(my.truthl_pdgid[i])==11:
            eleN+=1
        if fabs(my.truthl_pdgid[i])!=21:
            h1_truthpdgid.Fill(my.truthl_pdgid[i])
    jetN=0
    for i in range(my.pfjet_n):
        if(fabs(my.pfjet_eta[i])<2.5):
            jetN+=1
    h1_jetn.Fill(my.pfjet_n)
    h1_muon.Fill(my.muo_n)
    h1_elen.Fill(my.ele_n)
    muon=[]
    for i in range(my.muo_n):
        h1_muo_pt.Fill(my.muo_pt[i])
        h1_muo_pt.Fill(my.muo_pt[i])
        h1_muo_RelTrkIso.Fill(my.muo_RelTrkIso[i])
        h1_muo_TrkIso.Fill(my.muo_TrkIso[i])
        h1_muo_ECalIso.Fill(my.muo_ECalIso[i])
        h1_muo_HCalIso.Fill(my.muo_HCalIso[i])
        h1_muo_TrkIsoDep.Fill(my.muo_TrkIsoDep[i])
        h1_muo_ECalIsoDep.Fill(my.muo_ECalIsoDep[i])
        h1_muo_HCalIsoDep.Fill(my.muo_HCalIsoDep[i])
        h1_muo_AllIso.Fill(my.muo_AllIso[i])
        muon.append(TLorentzVector())
        muon[-1].SetPtEtaPhiM(my.muo_pt[i],my.muo_eta[i],my.muo_phi[i],0.105658)
    if my.muo_n>=2:
        h1_muo_minv.Fill((muon[0]+muon[1]).M())
    h1_met_et3.Fill(my.met_et[3])
    h1_met_et5.Fill(my.met_et[5])

    
    
    

myfileout = TFile('outTree.root',"RECREATE")
myfileout.cd()
for i in allHists:
    i.Write()

