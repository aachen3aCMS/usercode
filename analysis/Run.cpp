#include "SUSYAna.h"

#include <iostream>
#include <TString.h>
#include <TSystem.h>
#include <TStopwatch.h>

using std::cout;
using std::endl;

char inFile[1000];
char outFile[100];
char ag[100];
bool useInFile;
bool useList;
bool deb;

void parseCommandLine(int argc, char *argv[]);
void usage();

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

  SUSYAna *susy;

  TString s(inFile);
  TString type(ag);
  TString outf(outFile);

  if (type == "\0") {
    cout << "RunSUSY: You MUST provide the input type !" << endl;
    usage();
  }
  if (outf == "\0") {
    cout << "RunSUSY: You MUST provide an output file !" << endl;
    usage();
  }

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
    susy = new SUSYAna(fChain);
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

    susy = new SUSYAna(tree);
  }
  else {
    cout << "RunSUSY: You MUST provide the input files or file lists !" << endl;
    usage();
  }

  cout << "RunSUSY: Starting Event Loop" << endl;

  gSystem->GetProcInfo(&info);

  cout << "RunSUSY: resident mem (MB) : " << info.fMemResident/1000. << endl;
  cout << "RunSUSY: virtual  mem (MB) : " << info.fMemVirtual/1000. << endl;
  cout << endl;

  susy->Loop(outf, deb, type);

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

void usage() {
  cout << endl;
  cout << " Usage: RunSUSY -in     <inputfile> -out <outputfile> -type <data|mc|signal> [ -debug ] " << endl;
  cout << " Usage: RunSUSY -inlist <inputlist> -out <outputfile> -type <data|mc|signal> [ -debug ] " << endl;
  cout << endl;
  exit(-1);
}

// stolen from O. Peter's DataFilter

void parseCommandLine(int argc, char *argv[]) {
  if (argc==1)
    usage();

  for (int arg=1; arg<argc; arg++) {
    if (strcmp(argv[arg], "-out") == 0) {
      strcpy(outFile,argv[arg+1]);
      cout << "RunSUSY: Outputfile " << outFile << endl;
      arg++;
      continue;
    } 
    else if (strcmp(argv[arg], "-in") == 0) {
      strcpy(inFile,argv[arg+1]);
      cout << "RunSUSY: Inputfile " << inFile << endl;
      arg++;
      useInFile = true;
      useList = false;
      continue; 
    } 
    else if (strcmp(argv[arg], "-inlist") == 0) {
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
      continue;
    }
    else if (strcmp(argv[arg], "-type") == 0) {
      strcpy(ag,argv[arg+1]);
      cout << "RunSUSY: Input type " << ag << endl;
      arg++;
      continue; 
    } 
  }
}



