#include "signalInfo.hh"

signalInfo::signalInfo()
{
  Baseline = 0.;
  Noise =  0.;
  Amplitude = 0.;
  RiseTime = 0;
  Max_pos = 0;
  Trigger = 0;
}

signalInfo::~signalInfo()
{

}


void signalInfo::GetSignalInfo(UShort_t* v){
  UShort_t step = 3;
  //calculate the trigger, baseline and noise
  Double_t accu =0.,accu2 =0., mv_avg =0., mv_sigma =0.;
  for (UShort_t i = 1; i < TRACE_SIZE; i++) {
    accu += (Double_t)v[i];
    mv_avg = accu/(Double_t)i;

    accu2 += ((Double_t)v[i] - mv_avg)*((Double_t)v[i] - mv_avg);
    mv_sigma = TMath::Sqrt(accu2/(Double_t)(i-1));

    if(v[i] > (mv_avg + 5.*mv_sigma)){
      Baseline  = mv_avg;
      Noise     = mv_sigma;
      Trigger   = i;
      break;
    }

  }

  //determine the polarity of the signal
  Double_t signal_is = 0.;
  for(UShort_t n = Trigger; n < 900; n++)
   if( (v[n] > Baseline + 10.*Noise) || (v[n] < Baseline- 10.*Noise)){
      signal_is+=v[n]- Baseline;
    }

  
  //Get max value and max pos
  Double_t max_val = (Double_t)v[Trigger];
  UShort_t max_pos = 0;
  for (UShort_t i = Trigger; i < Trigger+100; i++) {
    if(signal_is > 0.){

    if((Double_t)v[i] > max_val){
      max_val = (Double_t)v[i];
      max_pos = i;
    }
    }
    else{

    if((Double_t)v[i] < max_val){
      max_val = (Double_t)v[i];
      max_pos = i;
    }
    }
  }
  
  Amplitude = max_val - Baseline;
  Max_pos = max_pos;

  //calculate risetime
  //using namespace std;
  UShort_t  interval =0.;
  Double_t temp = 0.;
  for (UShort_t i = max_pos - 1; i > 0; i--) {
    temp = v[i] - Baseline;
    if(signal_is > 0.){
if(temp < 0.1*max_val) break;
    if(temp <= .9* max_val && temp >= 0.1*max_val){
      interval++;
    }
    }
    else{
if(temp > 0.1*max_val) break;
    if(temp >= .9* max_val && temp <= 0.1*max_val){
      interval++;
    }
    }
    
  }
  RiseTime = interval * sampling_period;
  Amplitude = std::abs(Amplitude);
}
