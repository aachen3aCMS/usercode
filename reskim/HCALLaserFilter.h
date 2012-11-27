#ifndef HCALLaserFilter_h
#define HCALLaserFilter_h

// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <TString.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifndef __CINT__
#include <zlib.h>
#endif

using namespace std;


class HCALLaserFilter{
  
    public :
        explicit HCALLaserFilter(TString eventFileName);
        ~HCALLaserFilter();
        bool filter(int run, int lumisection, int event);


    private:

        void readEventListFile(const string & eventFileName);
        void addEventString(const string & eventString);
        
        typedef std::vector< std::string >::iterator strVecI;
        //typedef struct gzFile_s *gzFile;    /* semi-opaque gzip file descriptor */
        
        std::vector< std::string > EventList_;  // vector of strings representing bad events, with each string in "run:LS:event" format
        bool verbose_;
        // Set run range of events in the BAD LASER LIST.  
        // The purpose of these values is to shorten the length of the EventList_ vector when running on only a subset of data
        int minrun_;
        int maxrun_;  // if specified (i.e., values > -1), then only events in the given range will be filtered
        int minRunInFile, maxRunInFile;
};
#endif 
