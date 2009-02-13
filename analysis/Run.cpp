#include "SUSYAna.h"

#include <iostream>
#include <TString.h>
#include <TSystem.h>
#include <TStopwatch.h>

using std::cout;
using std::endl;

char inFile[1000];
char outFile[100];
bool useInFile;
bool useList;
bool deb;

void parseCommandLine(int argc, char *argv[]);

int main(int argc, char *argv[]) {

  ProcInfo_t info;

  TStopwatch timer;
  timer.Start();

  useInFile = false;
  useList = false;

  deb = false;

  cout << endl;
  cout << "  ------------------------------" << endl;
  cout << "  |                            |" << endl;
  cout << "  |  Start executing RunSUSY   |" << endl;
  cout << "  |                            |" << endl;
  cout << "  ------------------------------" << endl;
  cout << endl;
  gSystem->Exec("TEMPDATE=`date`; echo 'RunSUSY:' $TEMPDATE");
  cout << endl;

  parseCommandLine(argc,argv);

  SUSYAna *fth;

  TString s(inFile);
 
  if (useList) {

    cout << "RunSUSY: Collecting input files" << endl;

    int nfiles;
    char tex[100];
    TString ttex;
    TChain *fChain;
    FILE *filelist;
    
    filelist = fopen(inFile, "r");
    nfiles = 0;
    fChain = new TChain("ACSkimAnalysis/allData");
    while( !feof(filelist) ) {
      
      fscanf(filelist, "%s\n", tex);
      ttex = tex;
      if (ttex.Contains("pnfs"))
	ttex = "dcache:" + ttex;
      if( strcmp(tex,"//")==0 ) break;
      cout << "          " << ttex << endl;
      fChain->Add(ttex);
      nfiles++;
      //  if( nfiles==nmax ) break;
    }
    cout << "            " << "(" << nfiles << " files)" << endl;
    fth = new SUSYAna(fChain);
  }
  
  else if (useInFile) {
    TFile *f;
    if (s.Contains("pnfs")) {
      s = "dcache:" + s;
      f = TFile::Open(s);
    }
    else
      f = new TFile(inFile);

    TTree *tree = (TTree*)f->Get("ACSkimAnalysis/allData");

    fth = new SUSYAna(tree);
  }
  else {
    cout << endl;
    cout << " Usage: RunSUSY [-in <input> | -inlist <inputlist>] -out <outputfile>  [ -debug ] " << endl;
    cout << endl;
    exit(-1);
  }
  fth->setOutFile(outFile);
  fth->setMode(deb);
  cout << "RunSUSY: Starting Event Loop" << endl;
  gSystem->GetProcInfo(&info);
  cout << "RunSUSY: resident mem (MB) : " << info.fMemResident/1000. << endl;
  cout << "RunSUSY: virtual  mem (MB) : " << info.fMemVirtual/1000. << endl;
  cout << endl;

  fth->Loop();

  timer.Stop();

  gSystem->GetProcInfo(&info);

  cout << endl;
  cout << "RunSUSY: DONE " << endl;
  cout << "RunSUSY: REAL time (s)     : " << timer.RealTime() << endl;
  cout << "RunSUSY: CPU  time (s)     : " << timer.CpuTime() << endl; 
  cout << "RunSUSY: resident mem (MB) : " << info.fMemResident/1000. << endl;
  cout << "RunSUSY: virtual  mem (MB) : " << info.fMemVirtual/1000. << endl;
  cout << endl;
  gSystem->Exec("TEMPDATE=`date`; echo 'RunSUSY:' $TEMPDATE");
  cout << endl;
  cout << endl;

}


// stolen from O. Peter's DataFilter

void parseCommandLine(int argc, char *argv[]) {
  for (int arg=1; arg<argc; arg++) {
    if (strcmp(argv[arg], "-out") == 0) {
      if (argc < arg+1) {
	cout << "RunSUSY: Missing output filename" << endl;
	exit(-1);
      }
      strcpy(outFile,argv[arg+1]);
      cout << "RunSUSY: Outputfile " << outFile << endl;
      arg++;
      continue;
    } 
    else if (strcmp(argv[arg], "-in") == 0) {
      if (argc < arg+1) {
	cout << "RunSUSY: Missing filename" << endl;
	exit(-1);
      }
      strcpy(inFile,argv[arg+1]);
      cout << "RunSUSY: Inputfile " << inFile << endl;
      arg++;
      useInFile = true;
      useList = false;
      continue; 
    } 
    else if (strcmp(argv[arg], "-inlist") == 0) {
      if (argc < arg+1) {
	cout << "RunSUSY: Missing filename" << endl;
	exit(-1);
      }
      strcpy(inFile,argv[arg+1]);
      cout << "RunSUSY: Inputfilelist " << inFile << endl;
      arg++;
      useInFile = false;
      useList = true;
      continue;
    }
    else if (strcmp(argv[arg], "-debug") == 0) {
      deb = true;
      cout << "RunSUSY: Debug mode " << endl;
      arg++;
      continue;
    }
  }
}



