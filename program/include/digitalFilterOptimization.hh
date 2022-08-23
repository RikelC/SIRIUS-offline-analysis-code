#ifndef digitalFilterOptimization_h
#define digitalFilterOptimization_h 1

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TList.h"
#include "TDirectory.h"

#include <stdio.h>//C library to perform Input/Output operations
#include <algorithm>    // std::sort
#include <exception> 
#include <fstream>
#include <iomanip> //setw,setprecision
#include <Math/Interpolator.h>



#include "global.hh"
#include "signalInfo.hh"
#include "digitalFilters.hh"
#include "../source/digitalFilters.cc"
#include "trap_parameter.hh"
#include "resolutionStrip.hh"
class digitalFilterOptimization
{
private:
digitalFilters<UShort_t> *filter;

// MyGlobal *gObj = MyGlobal::GetInstance();
//variables for the branches of the input tree
ULong64_t ftime;
UInt_t feventNo;
UShort_t ftrace[TRACE_SIZE];
UInt_t fgain;
UShort_t fboardID;
UShort_t fchannelID;
//variables for the branches of the trap tree
Double_t trapAmplitude; 
UShort_t fboardIdx;
UShort_t k, m,l, kIdx, mIdx;
//
UShort_t kStart, kStop, kStep;
UShort_t mStart, mStop, mStep;
UShort_t ktimes;
UShort_t mtimes;

TFile* file_tree;
TFile* file_spectra;
TTree* trapTree; // store trapezoidal amplitudes
TRandom3 * rand;
signalInfo signal;
  
//create histograms
int nPeaks;
Double_t *fPositionX;
Double_t *fPositionY;
//Fit functions
TF1 *fLinear;
TF1 *G0;
TF1 *G1;
TF1 *G2;
// graph for calibration
TGraph* gr_cal;
Double_t Energy[2];
TSpectrum *spectrum;
TCanvas* canvas;
TCanvas* canvas2;


TList* list_res;
  
TH1F ***** hRaw;
TH1F ***** hCalib;
TH2F *** hResolution;

TList* list_calib;
TH2F *** hCalibStrip;
TH2F* h_m;
TH2F* h_c;

//calibration parameters
Double_t ****gain;
Double_t ****offset;
void initialize_TObjects();
void delete_TObjects();
void setKnM_multiplicity(UShort_t k, UShort_t m);
public:
digitalFilterOptimization();
~digitalFilterOptimization();

void optimize_trapezoidal_filter_parameters(int, vector<string>);
void sort_data_with_min_trapParametes(int, std::string);
/*   void sort_data_with_min_trapParametes2(int, std::string);
     void sort_data_with_min_trapParametesNCalibParameters(int, std::string);*/

    
void set_Trapezoidal_optimization_limits(UShort_t* par){
kStart= par[0];
kStop = par[1];
kStep = par[2];
mStart = par[3];
mStop = par[4];
mStep = par[5];

ktimes = (kStop - kStart)/kStep + 1;
mtimes = (mStop - mStart)/mStep + 1;
};



};

#endif
