#ifndef SUSYAna_h
#define SUSYAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <set>

#include "TreeContent.h"

using std::cout;
using std::endl;
using namespace std;

class SUSYAna : public TreeContent {
  
public :
    
  SUSYAna(TTree *tree=0);
  void Loop();
 
  char _fname[100];
  void setOutFile(char fn[100]){strcpy(_fname,fn);}

  bool DEBUG;
  void setMode(bool bset){DEBUG = bset;}
  
  void init(int stages);
  void write(int stages);

  bool find_duplicate(int run, int evt, double x1, double x2);
  
  typedef std::pair< pair<int,int> , pair<double,double> > Key;
  typedef std::set<Key> KeySet;
  typedef KeySet::const_iterator KeyIter;
  
  KeySet _keys;

  std::vector<TH1F*> h1_mu_pt;
  std::vector<TH1F*> h1_mu_TrkIso;
  
};
#endif // #ifdef SUSYAna_cxx
