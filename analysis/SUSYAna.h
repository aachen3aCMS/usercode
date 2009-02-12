#ifndef SUSYAna_h
#define SUSYAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <algorithm>

#include <iostream>

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
  
  std::vector<TH1F*> h1_mu_pt;
  std::vector<TH1F*> h1_mu_TrkIso;
  
};
#endif // #ifdef SUSYAna_cxx
