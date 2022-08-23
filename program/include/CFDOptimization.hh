#ifndef CFDOptimization_h
#define CFDOptimization_h 1
#include "global.hh"
#include "signalInfo.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TAttLine.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TSpectrum.h"
#include <iostream>
#include <string>
#include <stdio.h>
#include <algorithm>    // std::sort
#include <exception> 
#include <fstream>
#include <iomanip> 
#include <Math/Interpolator.h>
#include "CFDparameters.hh"
#include "TRandom3.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "digitalCFD.hh"

using namespace std;
class DataInfo{
public:
    unsigned long long int TimeStamp;
    double Baseline;
    double Noise;//sigma
    double Amplitude;
    double RiseTime;
    unsigned short Max_pos;
    unsigned short Trigger;
    double **CFD_Time;
    int col_num;
    int row_num;
  void clear_cfd_array(){
   for(int i = 0; i < row_num;i++) if(CFD_Time[i]) delete[] CFD_Time[i];

    if(CFD_Time)delete []CFD_Time;}
  void initialize_cfd_array(int row, int col){
    CFD_Time = new double*[row];
    for(int i = 0; i < row;i++)
        CFD_Time[i] = new double[col];

    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            CFD_Time[i][j] = 0.;
        }
    }
  }

};
class CFDOptimization
{
 private:
   //variables for the branches of the input tree
    ULong64_t ftime;
    UInt_t feventNo;
    UShort_t ftrace[TRACE_SIZE];
    UInt_t fgain;
    UShort_t fboardID;
    UShort_t fchannelID;
    //

  Bool_t SiFlag, SedFlag, PmtFlag;
    Bool_t draw_cfd_components;
    Bool_t eventsProcessFlag;

    Int_t si_cfd_delay_low, si_cfd_delay_high;
    Int_t sed_cfd_delay_low, sed_cfd_delay_high;
  Int_t pmt_cfd_delay_low, pmt_cfd_delay_high;    
  Int_t n_events_to_process;
    
  Int_t Si_CFDmethod;
    Int_t Sed_CFDmethod;

  Int_t Pmt_CFDmethod;

    Int_t sed_delay_step;
Int_t pmt_delay_step;    
Int_t si_delay_step;
    Int_t si_n_delay;
    Int_t sed_n_delay;
  Int_t pmt_n_delay;
    Double_t si_cfd_fract_step;
    Double_t sed_cfd_fract_step;
  Double_t pmt_cfd_fract_step;

    Double_t si_cfd_fract_start;
    Double_t si_cfd_fract_stop;
   
 Double_t sed_cfd_fract_start;
    Double_t sed_cfd_fract_stop;

 Double_t pmt_cfd_fract_start;
    Double_t pmt_cfd_fract_stop;

    Int_t nCFD_si;
    Int_t nCFD_sed;
    Int_t nCFD_pmt;
    
    Double_t sed_noise_max;
    Double_t sed_baseline_max;

    Double_t pmt_noise_max;
    Double_t pmt_baseline_max;

    Double_t si_noise_max;
    Double_t si_baseline_max;

    //To calculate the max value from a polynomial fit
    TGraph* gr;
    TF1* ampFitFunc;
    //for the Draw method to visualize the CFD components
    TGraph* hcfd_fracvin;
    TGraph*  hcfd_vin;
    TGraph* hcfd_vindelayed;
    TGraph* hcfd;

    //Methods
    void fillSignal(TH2F *, UShort_t* );//to fill the traces in a 2D histogram

  TRandom3 *rand;
   TRandom3 *rand2;
   TRandom3 *rand3;
    double **CFD_Time;
    int col_num;
    int row_num;
    void initialize_cfd_array(int row, int col);
 public:
  CFDOptimization();
  ~CFDOptimization();

   //setters
    void set_draw_cfd_components(Bool_t flag){draw_cfd_components = flag;}
    void set_Si_cfd_delay_parameters(Int_t n1, Int_t n2, Int_t n3){
        si_cfd_delay_low = n1;
        si_cfd_delay_high = n2;
        si_delay_step = n3;
        si_n_delay = (n2-n1)/n3 +1;
    }
    void set_Sed_cfd_delay_parameters(Int_t n1, Int_t n2, Int_t n3 ){
        sed_cfd_delay_low = n1;
        sed_cfd_delay_high = n2;
        sed_delay_step = n3;
        sed_n_delay = (Int_t)((n2-n1)/n3 +1);
    }
 void set_Pmt_cfd_delay_parameters(Int_t n1, Int_t n2, Int_t n3 ){
        pmt_cfd_delay_low = n1;
        pmt_cfd_delay_high = n2;
        pmt_delay_step = n3;
        pmt_n_delay = (Int_t)((n2-n1)/n3 +1);
    }

    void set_Si_cfd_attenuation_parameters(Double_t n1, Double_t n2, Double_t n3 ){
        si_cfd_fract_start = n1;
        si_cfd_fract_stop = n2;
        si_cfd_fract_step = n3;
        nCFD_si = (Int_t)((n2-n1)/n3 +1);
    }

    void set_Sed_cfd_attenuation_parameters(Double_t n1, Double_t n2, Double_t n3 ){
        sed_cfd_fract_start = n1;
        sed_cfd_fract_stop = n2;
        sed_cfd_fract_step = n3;
        nCFD_sed = (n2-n1)/n3 +1;
    }


    void set_Pmt_cfd_attenuation_parameters(Double_t n1, Double_t n2, Double_t n3 ){
        pmt_cfd_fract_start = n1;
        pmt_cfd_fract_stop = n2;
        pmt_cfd_fract_step = n3;
        nCFD_pmt = (n2-n1)/n3 +1;
    }

    void set_Sed_Time_calcul_method(Int_t n){Sed_CFDmethod =n;}
    void set_Si_Time_calcul_method(Int_t n){Si_CFDmethod =n;}
 void set_Pmt_Time_calcul_method(Int_t n){Pmt_CFDmethod =n;}
        
void set_n_events_to_process(Int_t n){eventsProcessFlag = 1; n_events_to_process =n;}
  void set_baseline_limit(Double_t pmt,Double_t sed, Double_t si){pmt_baseline_max = pmt;sed_baseline_max = sed; si_baseline_max = si;}
  void set_noise_limit(Double_t pmt, Double_t sed,  Double_t si){pmt_noise_max = pmt; sed_noise_max = sed;  si_noise_max = si;}
    //
    void perform(std::string );
    void do_optimization(std::string );
    void drawTrace(Int_t, std::string);
    void test(Int_t, std::string);
  void monteCarlo();
 void monteCarlo2();
};

#endif
