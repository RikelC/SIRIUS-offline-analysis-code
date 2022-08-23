#ifndef ReadDataFiles_h
#define ReadDataFiles_h 1
#include <iostream>
using namespace std;
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include "stdio.h"
#include "stdlib.h"
#include <iostream> 
#include "MFMAllFrames.h"
#include "myHistograms.hh"
#include "global.hh"
#include "signalTreatment.hh"
class ReadDataFiles
{
 private:
  
  MFMCommonFrame *  fFrame;
  MFMSiriusFrame *  fSiriusframe;
  myHistograms *hist;
  const char * filename = "/data/siriusX/test/acquisition/run/run_0052.dat.02-02-22_12h19m54s";
  int fVerbose=0;
  int fDumpsize =32;
  int frameCounter[8];
  string frameName[8]={"Cobo","Eby","Hello","xml header",
                       "rea trace","rea gen","sirius","else"};


  signalTreatment *pData;
  int metatype ;
  int blocksize ;
  bool endiness ;
  bool blobness ;

  int framesize ;
  int datasource;
  int type;
  int rivision;

  int headersize;
  uint32_t eventnumber;
  uint64_t timestamp;
  int  dumpsize;
  int maxdump =128;
  int channel=0, board=0,iboard=0;
  int NbItems,i;
  uint16_t value;
 

  void ReadUserFrame(MFMCommonFrame* commonframe);
  void ReadDefaultFrame(MFMCommonFrame* commonframe) ;
  void ReadSiriusFrame(MFMCommonFrame* commonframe) ;

 public:
  ReadDataFiles();
  ~ReadDataFiles();

  void perform();
};

#endif
