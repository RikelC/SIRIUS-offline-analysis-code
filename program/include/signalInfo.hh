#ifndef signalInfo_h
#define signalInfo_h 1
#include "TMath.h"
#include "global.hh"

class signalInfo{
 public:
  signalInfo();
 ~signalInfo();
  Double_t Baseline;
  Double_t Noise;//sigma
  Double_t Amplitude;
  UShort_t RiseTime;
  UShort_t Max_pos;
  UShort_t Trigger;
  void GetSignalInfo(UShort_t* v);
};

#endif
