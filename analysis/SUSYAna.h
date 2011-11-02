#ifndef SUSYAna_h
#define SUSYAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <set>

#include "Combinations.h"

#include "TreeContent.h"

using std::cout;
using std::endl;
using namespace std;

class SUSYAna : public TreeContent {
  
public :
    
  SUSYAna(TTree *tree=0);
  void Loop(TString fout, bool debug, TString type);
 
  bool DEBUG;
  
  void init(int stages);
  void write(TString fout, int stages);

  void BasicDump(int i);
  void TriggerDump(TString sel);
  void MuonDump(bool full=1);
  void CaloJetDump();
  void PFJetDump();
  void TruthJetDump();
  void TruthDump();
  void VertexDump();
  void METDump();
  void SCDump();
  void EleDump(bool full=1);
  void PFEleDump(bool full=1);

  void doJESandRecalculateMET(TString corr);

  bool find_duplicate(int run, int evt, double x1, double x2);
  
  Double_t DeltaPhi(double a, double b);
  Double_t mT(double et1, double phi1, double et2, double phi2);

  Double_t AlphaT(const std::vector<TLorentzVector> & objects);
  Double_t MinDEt(const std::vector<TLorentzVector> & objects, 
		  std::vector<UInt_t> * lista = NULL, 
		  std::vector<UInt_t> * listb = NULL);
  Double_t SumET(const std::vector<TLorentzVector> & objects);
  Double_t MT(const std::vector<TLorentzVector> & objects);
  double M2(double E1, double px1, double py1, double pz1, double E2, double px2, double py2, double pz2);
  double M2v2(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2);
  double AngleMuons(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2);

  typedef std::pair< pair<int,int> , pair<double,double> > Key;
  typedef std::set<Key> KeySet;
  typedef KeySet::const_iterator KeyIter;
  
  KeySet _keys;

  std::vector<TH1F*> h1_mu_pt;
  std::vector<TH1F*> h1_mu_Minv;
  
};
#endif // #ifdef SUSYAna_cxx
