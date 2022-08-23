#ifndef myHistograms_h
#define myHistograms_h 1

#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "global.hh"
#include "TFile.h"
#include "TDirectory.h"
class myHistograms
{

 private:
  TDirectory *dFile;
  TDirectory *dBoard[NBOARDS];
   
  
 public:
  myHistograms();
  ~ myHistograms();
  TH1I *** Gain ; //!
  TH1I *** FeedBack; //!
  TH1I *** Trace  ; //!
  TH2I *** Trace_Sum ;  //!
  TH1D *** hRaw;
  // TH1F *** hTrap;
  void rawHistInitialization(TFile* f);
  void deleteRawHistograms();
  
};

#endif
