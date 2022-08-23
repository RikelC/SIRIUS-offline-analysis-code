#ifndef signalTreatment_h
#define signalTreatment_h 1
#include "global.hh"
#include "digitalFilters.hh"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TAttLine.h"
#include "TF1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include <stdio.h>
#include <algorithm>    // std::sort
#include <exception> 
#include <fstream>
#include <iomanip> 
#include <Math/Interpolator.h>
#include "global.hh"
#include "signalInfo.hh"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TAttFill.h"
#include <sys/types.h>
#include <dirent.h>
#include  "TRandom3.h"
#include "TChain.h"
#include<unistd.h>
#include "TList.h"
#include "TDirectory.h"
#include "../source/digitalFilters.cc"
//class digitalFilters;
class signalTreatment
{
 private:

  digitalFilters<UShort_t> *myFilter;
  ULong64_t ftimestamp;
  ULong64_t ftime;
  UInt_t feventNo;
  UShort_t ftrace[TRACE_SIZE];
  UInt_t fgain;
  UShort_t fboardID;
  UShort_t fchannelID;
  UShort_t boardIdx;
  //
  UShort_t kStart, kStop, kStep;
  UShort_t mStart, mStop, mStep;
  UShort_t k, l, m;
  UShort_t kTimes;
  UShort_t mTimes;
  Double_t trapAmplitude;
  TTree* trapTree; // store trapezoidal amplitudes
  TRandom3 * rand;
  signalInfo info;
    //private methods
  //Double_t signal(Int_t , Int_t , Int_t , UShort_t* , Double_t );
  template<class T>
  void invert_signal(T* v);
 public:
  signalTreatment();
  ~signalTreatment();
  void drawEntry(int, std::string);
  void drawRun(int, std::string);
  void perform();//get amp, calibrate, get time

  void set_eventnumber(uint32_t i){feventNo =i;}
  void set_timestamp(uint64_t i){ftimestamp =i;}
  void set_boardID(UShort_t i){fboardID =i;}
  void set_channelID(UShort_t i){fchannelID = i;}
  void set_gain(UInt_t i){fgain = i;}
void set_trace_value(int i, uint16_t v){ftrace[i] = v;}
  //
  uint32_t get_eventnumber()const {return feventNo;}
  uint64_t get_timestamp()const {return ftimestamp;}
  UShort_t get_boardID()const {return fboardID;}
  UShort_t get_channelID()const{return fchannelID;}
  UInt_t get_gain()const{return fgain;}
  uint16_t get_trace_value(int i)const {return ftrace[i];}
  uint16_t* get_trace(){return ftrace;}

};
#endif
