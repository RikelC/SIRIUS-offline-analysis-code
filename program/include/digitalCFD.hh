#ifndef digitalCFD_h
#define digitalCFD_h 1

#include "global.hh"
#include "signalInfo.hh"
#include <Math/Interpolator.h>
class digitalCFD
{
private:
  
public:
digitalCFD();
~digitalCFD();
Double_t getCFD_Time1(UShort_t* , signalInfo* ,  Double_t );
Double_t getCFD_Time2(UShort_t*, signalInfo*,  Int_t, Double_t);
Double_t getCFD_Time3(UShort_t*, signalInfo*,  Double_t);
};

#endif
