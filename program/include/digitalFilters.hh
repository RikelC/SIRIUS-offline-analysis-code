#ifndef digitalFilters_h
#define digitalFilters_h 1
#include <math.h>
#include <iostream>
#include <chrono>
#include <string>
#include "TH1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "global.hh"


using namespace std;
template<class T>
class digitalFilters
{
 private:
  Double_t Capacitance;//in pF
  Double_t RC_constant;
  double  trapezoidal_shaper_M ;
  double  trapezoidal_shaper_m2;// Gain of the Shaper Amplifier
  double  trapezoidal_shaper_m1; 
  //arrays for trapezoidal filters
  Double_t* array_u ;
  Double_t* array_v ;
  Double_t* array_d ;
  Double_t* array_p;
  Double_t* Rectangular;
  Double_t* Trapezoidal;
 Double_t signal(Int_t , Int_t , Int_t , T* , Double_t );
  Double_t get_max_val_trapezoidal(Double_t, Double_t* , UShort_t, Double_t);
 public:
  digitalFilters();
  ~digitalFilters();

  
  Double_t* get_array_u(){return array_u;}
  Double_t* get_array_v(){return array_v;}
  Double_t* get_array_d(){return array_d;}
  Double_t* get_array_p(){return array_p;}
  Double_t* get_array_Rectangular(){return Rectangular;}
  Double_t* get_array_Trapezoidal(){return Trapezoidal;}
  
  void generate_simulted_signal(int, int, double,  double, T*);
  Double_t CR_RC4_filter(UShort_t , T *, Double_t );
  Double_t trapezoidal_filter_algorithm1( UShort_t , UShort_t, T* , Double_t);
  Double_t trapezoidal_filter_algorithm2( UShort_t , UShort_t , T* , Double_t);
  Double_t trapezoidal_filter_algorithm3( UShort_t , UShort_t , T* , Double_t);
  Double_t trapezoidal_filter_algorithm4( UShort_t , UShort_t , T* , Double_t);
  Double_t trapezoidal_filter_algorithm5( UShort_t , UShort_t , T* , Double_t);
  void pole_zero_correction(T*);
  void compare_trap_algorithms(int, std::string, T*);


  
 
  //set gain
  void set_DSSD_gain(Double_t);
  void set_RC_constant(Double_t);

};


#endif
