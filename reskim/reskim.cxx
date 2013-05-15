// local includes
#include "Utilities.h"
#include "Analysis.h"

// ROOT includes
#include "TStopwatch.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"

// C/C++ includes
#include <stdlib.h>
#include <iostream>

using namespace std;
using namespace ROOT;

enum ErrorCodes { E_WRONG_PARAMS, E_OUTPUTFILE, E_NO_TREE, E_NO_ANALYSIS };

void deactivateBranches(TTree* inputTree){
  INFO("Setting output branch addresses");
  TObjArray * branchlist = inputTree->GetListOfBranches();
  ///######################################################################################################
  ///   Change the status of Branches not needed for your analysis
  ///
  inputTree->SetBranchStatus("noise_*",1);
  inputTree->SetBranchStatus("truthjet_*",1);
  inputTree->SetBranchStatus("SystJet_*",1);
  inputTree->SetBranchStatus("calojet_*",1);
  inputTree->SetBranchStatus("pfjet_*",1);
  inputTree->SetBranchStatus("fatjet_*",1);
  inputTree->SetBranchStatus("SC_*",1);
  inputTree->SetBranchStatus("ele_*",1);
  inputTree->SetBranchStatus("pfele_*",1);
  inputTree->SetBranchStatus("tau_*",1);
  inputTree->SetBranchStatus("pho_*",1);
  inputTree->SetBranchStatus("susy*",1);
}

int main(int argc, char *argv[])
{
  // process and timing information
  ProcInfo_t info;
  TStopwatch timer;
  timer.Start();

  // get date and time
  time_t rawtime;
  struct tm * timeinfo;
  time(&rawtime);
  timeinfo = localtime ( &rawtime );

  // user greeting
  INFO("");
  INFO("######################################################################");
  INFO("#                                                                    #");
  INFO("# Start executing analyzer at " << asctime(timeinfo)		       );
  INFO("#                                                                    #");
  INFO("######################################################################");
  INFO("");

  // check command line arguments
  if (argc < 3) {
    ERROR("Usage: reskim pdfset outputfile.root inputfile [inputfile ... ]");
    return E_WRONG_PARAMS;
  }

  // read configuration file, configure values from command line
  gLogLevel = 3;
  const char * outputFileName = argv[2];
  TString pdfset = argv[1];

  TChain * chain = new TChain("ACSkimAnalysis/allData");
  for (int n = 3; n < argc; n++) {
    // try to find out if inputfile is a directory.
    Text_t * basepath  = gSystem->ExpandPathName(argv[n]);
    void   * dirhandle = gSystem->OpenDirectory(basepath);
    if (dirhandle != 0) {
      const Text_t * basename;
      while ((basename = gSystem->GetDirEntry(dirhandle))) {
	// Skip non-ROOT files
	if (!strstr(basename, ".root"))
	  continue;
	string fullname = basepath + string("/") + basename;
	INFO("Adding file " << fullname);
	chain->Add(fullname.c_str());
      }
    }
    else {
      // OK, this is not a directory, add the given file(s)
      INFO("Adding file " << argv[n]);
      chain->Add(argv[n]);
    }
    delete basepath;
  }

  // create output file
  INFO("Creating output file and cloning tree");
  TFile * outFile = TFile::Open(outputFileName, "RECREATE");
  if (outFile == 0 || !outFile->IsOpen()) {
    ERROR("Could not open output file " << outputFileName);
    return E_OUTPUTFILE;
  }
  outFile->SetCompressionLevel(6);
  outFile->mkdir("ACSkimAnalysis");
  outFile->cd("ACSkimAnalysis");
  //chain->SetBranchStatus("*",0);
  // clone tree structure, do not copy events...
  deactivateBranches(chain);
  TTree * outTree = chain->CloneTree(0);
  if (outTree == 0) {
    ERROR("Could not clone tree from input file - maybe no tree in input file?");
    return E_NO_TREE;
  }
  // histograms belong to ROOT directory...
  outFile->cd();
  
  // create analysis object
  Analysis * analysis = 0;
  try {
    analysis = new Analysis(*chain, *outTree);
  }
  CATCH;
  if (analysis == 0) {
    ERROR("Could not create analysis class - Either exception catched or out of memory");
    return E_NO_ANALYSIS;
  }

  // mem info
  gSystem->GetProcInfo(&info);
// //   INFO("resident mem (MB) : " << info.fMemResident/1000.);
  INFO("virtual  mem (MB) : " << info.fMemVirtual/1000.);

  // start loop
//   INFO("starting event loop");
  try {
    analysis->Loop(pdfset);
  }
  CATCH;
  INFO("end of event loop");

  // save data in file
//   INFO("Saving data to file and closing...");
  // get currrent file - needed because of output tree spanning different
  // files, via setting TTree::SetMaxTreeSize
  outFile = outTree->GetCurrentFile();
  outFile->Write();
  outFile->Close();
  delete outFile;
  // no, we do not need to delete the objects in the file, this is done on the file delete... 
  
  // delete analysis object
//   INFO("Deleting analysis");
  delete analysis;

  // time and memory info
  timer.Stop();
  gSystem->GetProcInfo(&info);
  INFO("");
  INFO("real time (s)     : " << timer.RealTime());
  INFO("CPU time (s)      : " << timer.CpuTime()); 
  INFO("resident mem (MB) : " << info.fMemResident/1000.);
  INFO("virtual  mem (MB) : " << info.fMemVirtual/1000.);
  INFO("");

  // info at end of program
  time(&rawtime);
  timeinfo = localtime ( &rawtime );
  INFO("");
  INFO("######################################################################");;
  INFO("#                                                                    #");;
  INFO("# End executing analyzer at " << asctime(timeinfo)		       );
  INFO("#                                                                    #");;
  INFO("######################################################################");;
  INFO("");
  
  return 0;
}
