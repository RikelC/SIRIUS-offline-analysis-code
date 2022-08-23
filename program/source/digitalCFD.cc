#include "digitalCFD.hh"

digitalCFD::digitalCFD()
{

}

digitalCFD::~digitalCFD()
{

}


Double_t digitalCFD::getCFD_Time1(UShort_t* v, signalInfo* info,  Double_t fract){
  Double_t amplitude = info->Amplitude;
  Int_t max_pos = (Int_t)info->Max_pos;
  Double_t baseline = info->Baseline;
  Double_t temp = 0.;

  Double_t y2= 0., y1= 0., x2= 0., x1 = 0., y = 0., x = 0., F = 0.;

  temp = (Double_t)v[max_pos] - baseline;
  y = fract*amplitude;
  while ( temp > y && (max_pos>0)){
    max_pos--;
    temp = v[max_pos] - baseline;
  }


  //cfd[i] = (maxbin -1) * sampling_period;//-1 because the first bin is empty
  // a linear interpolation for more precise value
  // y = y1 + (x - x1)* (y2 - y1)/(x2 -x1)
  // F = (y2 - y1)/(x2 -x1)
  // x = x1 + (y -y1)/F
  x2 = (Double_t)(max_pos + 1)*sampling_period;
  y2 = v[max_pos + 1];

  x1 = (Double_t)(max_pos)*sampling_period;
  y1 = v[max_pos];

  F = (y2 - y1)/(x2 - x1);

  x = x1 + (y -y1)/F;

  return x;
}

Double_t digitalCFD::getCFD_Time2(UShort_t* v, signalInfo* info, Int_t delay,  Double_t fract){

  Int_t max_pos = info->Max_pos;
  Double_t amplitude = info->Amplitude;
  Double_t baseline = info->Baseline;
  
  Double_t* vout = new Double_t[TRACE_SIZE];
  vout[0] = 0.;
  Double_t* vin = new Double_t[TRACE_SIZE];
  vin[0] = 0.;
  Double_t* vindelayed = new Double_t[TRACE_SIZE];
  vindelayed[0] = 0.;

  Double_t y2= 0., y1= 0.;
  Int_t x2= 0, x1 = 0;
  Double_t y = 0., x = 0., M = 0., C = 0.;
  Double_t val = 0.;

  for (Int_t i = 1; i < TRACE_SIZE; i++) {
    vin[i] = v[i] - baseline;
  }


  for (Int_t j = 1; j < TRACE_SIZE; j++) {
    //forward delay
    if(j <= delay)
      vindelayed[j] = 0.;
    else
      vindelayed[j] = vin[j - delay];
  }

  for (Int_t t =1; t < TRACE_SIZE; t++) {
    vout[t] = vindelayed[t] - fract*vin[t];
    /* if(draw_cfd_components){
       hcfd_vin->SetPoint(t-1,(Double_t)(t-1), vin[t]);
       hcfd_vindelayed->SetPoint(t-1,(Double_t)(t-1), vindelayed[t]);
       hcfd_fracvin->SetPoint(t-1,(Double_t)(t-1), fract*vin[t]);
       hcfd->SetPoint(t-1,(Double_t)(t-1), vout[t]);
       }*/
  }

  //find max pos in the interval (max_pos - 100, max_pos+100)
  x2 = max_pos - 101;
  y2 = vout[x2];
  for(Int_t n = max_pos - 100; n < (max_pos + 100);n++){
    if(vout[n] > y2){
      y2 = vout[n];
      x2 = n;
    }
  }
  //find zero crossing point via linear interpolation
  while ( y2 > 0.){
    x2--;
    y2 = vout[x2];
  }
  if(y2 < 0.){x2++;y2 = vout[x2];}
  x1 = x2 - 1;
  y1 = vout[x1];
  M = (y2 - y1)/(Double_t)(x2 - x1);
  C = y1 - M*x1;
  y = 0.;
  x = (y -C)/M;

  Double_t cfdTime = x* sampling_period;

    
  
  delete [] vout;
  delete [] vin;
  delete [] vindelayed;

  return cfdTime;
}

Double_t digitalCFD::getCFD_Time3(UShort_t* v, signalInfo* info,  Double_t fract){

  Double_t amplitude = info->Amplitude;
  Int_t max_pos = info->Max_pos;
  Double_t baseline = info->Baseline;
  Double_t temp = 0.;
  unsigned int pts = 4;

  Double_t *x = new Double_t[pts];
  Double_t *y = new Double_t[pts];

  Double_t y_interp = 0., x_interp =0.;
  y_interp = fract*amplitude;
  temp = v[max_pos] - baseline;

  while ( temp > y_interp && (max_pos>0)){
    max_pos--;
    temp = v[max_pos] - baseline;
  }

  max_pos--;//
  for(Int_t p = 0; p < pts;p++){
    x[p] = (double) max_pos;
    y[p] = v[max_pos] -baseline;
    max_pos++;
  }

  ROOT::Math::Interpolator inter (pts, ROOT::Math::Interpolation::kCSPLINE);
  //ROOT::Math::Interpolation::Type type  : kLINEAR, kPOLYNOMIAL, kCSPLINE, kAKIMA, kCSPLINE_PERIODIC, kAKIMA_PERIODIC

  inter.SetData (pts, y, x);
  x_interp = inter.Eval(y_interp);


  delete [] x;
  delete [] y;
  return x_interp;
}
