//
// Package:    UserCode/aachen3a/ACSusyAnalysis
// Class:      SusyACSkimAnalysis
// 
// Description: Trigger helper tools (packing/unpacking)
//
// Original Author:  Carsten Magass
//         Created:  May 2009
//

#include <iostream>
#include <string.h>
#include <algorithm>

using std::cout;
using std::endl;
using namespace std;

#define INTSIZE 4

namespace ACSusyAnalysis {
  union u64 {
    char c[INTSIZE];
    int i;
  };
  
  int get_size(const int* a) {
    int size=0;
    u64 tmp;
    for(int i=0; ; ++i ) {
      tmp.i = a[i];
      size++;
      for(int j=0; j < INTSIZE ; ++j )
	if(tmp.c[j] == '\x0')
	  return size;
    }
    return -1;
  }
  
  std::string unpack(const int* a) {
    int size = get_size(a);
    u64 tmp;
    std::string ret;
    char* c = new char[size*INTSIZE];
    for(int i=0; i < size ; ++i) {
      tmp.i = a[i];
      for(int j=0; j < INTSIZE; ++j )
	c[i*INTSIZE+j] = tmp.c[j];
    }
    ret = c;
    delete c;
    return ret;
  }

  int* pack(const char* c) {
    u64 tmp;
    int j=0,count=0;
    int size = strlen(c)/INTSIZE+1;
    int* ii=new int[size];
    for(int i=0 ; ; i++) {
      tmp.c[j++]=c[i];
      ii[count]=tmp.i;
      if(j==INTSIZE) {
	j=0;
	count++;
      }
      if(c[i] == '\x0')
	break;
    }
    return ii;
  }
}
