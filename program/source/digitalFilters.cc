#include "digitalFilters.hh"
template <class T>
digitalFilters<T>::digitalFilters(){
  Capacitance = 1.;//in pF
  RC_constant = 550.;
  trapezoidal_shaper_M = 0. ;
  trapezoidal_shaper_m2 = 0.;
  trapezoidal_shaper_m1 = 0.;
  array_u = new Double_t[TRACE_SIZE];
  array_v = new Double_t[TRACE_SIZE];
  array_d = new Double_t[TRACE_SIZE];
  array_p = new Double_t[TRACE_SIZE];
  Rectangular = new Double_t[TRACE_SIZE];
  Trapezoidal = new Double_t[TRACE_SIZE];

}
template <class T>
digitalFilters<T>::~digitalFilters(){
  delete [] array_u;
  delete [] array_v;
  delete [] array_d;
  delete [] array_p;
  delete [] Rectangular;
  delete [] Trapezoidal;
}

template <class T>
void digitalFilters<T>::set_DSSD_gain(Double_t gain){
  Capacitance = gain; RC_constant =  Resistance*Capacitance;
  trapezoidal_shaper_M = pow((exp(sampling_period/RC_constant) - 1.), -1.);
  trapezoidal_shaper_m2 = 1.;// Gain of the Shaper Amplifier
  trapezoidal_shaper_m1 =  trapezoidal_shaper_m2 *  trapezoidal_shaper_M;
}
template <class T>
void digitalFilters<T>::set_RC_constant(Double_t  r){
  RC_constant = r;
  trapezoidal_shaper_M = pow((exp(sampling_period/RC_constant) - 1.), -1.);
  trapezoidal_shaper_m2 = 1.;// Gain of the Shaper Amplifier
  trapezoidal_shaper_m1 =  trapezoidal_shaper_m2 *  trapezoidal_shaper_M;
}


template<class T>
Double_t digitalFilters<T>::signal(Int_t n, Int_t nStart, Int_t m, T* v, Double_t baseline){
 if((n - m) < nStart) return baseline;
  else return (Double_t)v[n-m];
}

template <class T>
Double_t digitalFilters<T>::get_max_val_trapezoidal( Double_t signal_is, Double_t* v , UShort_t pos, Double_t max){
  //average on both sides of the max position
  Double_t maxR =0.; Double_t maxL = 0.;Double_t rR =0.; Double_t rL =0.;
  UShort_t j2 =pos + 1;
  if(signal_is > 0.){
    while(v[j2] > 0.999 * max){
      maxR += v[j2];
      j2++;
      rR++;
    }
    UShort_t j1 = pos -1;
    while(v[j1] > 0.999 * max){
      maxL += v[j1];
      j1--;
      rL++;
    }

  }else{
    while(v[j2] < 0.999 * max){
      maxR += v[j2];
      j2++;
      rR++;
    }
    UShort_t j1 = pos -1;
    while(v[j1] < 0.9999 * max){
      maxL += v[j1];
      j1--;
      rL++;
    }

  }
  
  //compute average
Double_t max_val = max + maxL + maxR;
  max_val = max_val /(rR+rL+1.);
  return max_val;
}
template <class T>
void digitalFilters<T>::generate_simulted_signal(int dt, int t0, double offset,  double Amplitude, T* ftrace){

  //double Amplitude = 500.;

  // int dt = 5;//sampling period
  // int t0 = 1000;
  // int t1 = 550;
  // double offset = 1000.;
  //Signal Generation with pile up
  for(int n = 0; n <TRACE_SIZE; n++){
    int t = n * dt;
    // ftrace[n] = Amplitude * exp(- t/ RC_constant); // Without delay
    double noise =  0;//rand->Uniform(0,50);
    if(t < t0){
      ftrace[n] = offset+ noise ; // With delay
    }

    else {
      ftrace[n]  = offset + Amplitude * exp(- (t - t0)/ RC_constant) + noise;
    }
    //if (t >= t1) { ftrace[n] += 200 * exp(- (t - t1)/ RC_constant);}//Pile Up Signal
  }
}

template <class T>
Double_t digitalFilters<T>::CR_RC4_filter(UShort_t S, T *x, Double_t G ){

  Double_t *y = Rectangular;
  Double_t *y_p0 = Trapezoidal;

  for(int i = 0; i < 10;i++){
    y[i] = 0.;
    y_p0[i] = 0.;
  }
  double period = sampling_period;
  double tau = RC_constant;
  Double_t alpha = (1.- exp(-period/tau));
  Double_t beta = (1+exp(-period/tau));
  //Ref. Ksenia thesis, page 86
  for(int n = 10; n < TRACE_SIZE; n++){
    y[n] = G*pow(alpha,4.)*beta*(Double_t)(x[n] - x[n-S]) + 5.*exp(-period/tau)*y[n-1] - 10.*exp(-2.*period/tau)*y[n-2] + 10.*exp(-3.*period/tau)*y[n-3] - 5.*exp(-4.*period/tau)*y[n-4] + 1.*exp(-5.*period/tau)*y[n-5];
    y_p0[n] = y_p0[n-1] + y[n] - exp(-period/tau)*y[n-1];
    
  }

  return 0;
}
template <class T>
Double_t digitalFilters<T>::trapezoidal_filter_algorithm1( UShort_t k , UShort_t m, T* trace , Double_t baseline){
  Double_t max_val =0.;
  UShort_t max_pos = 0;
  Double_t signal_is = 0.;
  UShort_t l = k +m;
  for(UShort_t n = 0; n <TRACE_SIZE; n++){
    if((n > 10 && n < 900) && (trace[n] > baseline + 100) || (trace[n] < baseline -100)){
      signal_is+=trace[n]-baseline;
    }

    if(n < k){
      array_u[n] =(Double_t) trace[n] - baseline;
      array_v[n] = (Double_t)trace[n] - baseline;
    }

    else{
      array_u[n]= (Double_t)(trace[n]- trace[n -k]);
    }

    if(n >= l+k){
      array_v[n] = (Double_t)(trace[n-l] - trace[n-l-k]);
    }

    array_d[n] = array_u[n] - array_v[n];

    Rectangular[n] = Rectangular[n-1] + array_d[n] - (exp(- 5 /RC_constant) *array_d[n-1]);

    Trapezoidal[n] = Trapezoidal[n-1] + Rectangular[n];


    //get max position here 

    if(signal_is > 0.){
      if(Trapezoidal[n] > max_val){
        max_val = Trapezoidal[n];
        max_pos = n;
      }
    }else{
      if(Trapezoidal[n] < max_val){
        max_val = Trapezoidal[n];
        max_pos = n;
      }
    }
  }

  max_val = get_max_val_trapezoidal(signal_is,Trapezoidal, max_pos, max_val);
  // cout<<"range "<<rR+rL+1<<endl;
  //return signal height
  Double_t signalAmplitude = max_val /(Double_t)(std::min(k,l));
  //cout<<"amplitude "<<signalAmplitude<<endl;
  return std::abs(signalAmplitude);
}

template<class T>
Double_t digitalFilters<T>::trapezoidal_filter_algorithm2( UShort_t k , UShort_t m, T* trace , Double_t baseline){
  Double_t max_val =0.;
  int max_pos = 0;
  Double_t signal_is = 0.;
  UShort_t l = k +m;
  for(UShort_t n = 0; n <TRACE_SIZE; n++){
    if((n > 10 && n < 900) && (trace[n] > baseline + 100) || (trace[n] < baseline -100)){
      signal_is+=trace[n]-baseline;
    }

    if( n < k){
      array_d[n] =(Double_t) trace[n] - baseline;
    }
    if(n >= k){
      array_d[n] = (Double_t)trace[n] - baseline;
      array_d[n] -= (Double_t)trace[n-k] - baseline;
    }
    if(n >= l +k){
      array_d[n]-= (Double_t)trace[n- l] - baseline;
      array_d[n] +=  (Double_t)trace[n - k-l] - baseline;
    }

    Rectangular[n] = Rectangular[n-1] + array_d[n] - (exp(- sampling_period /RC_constant) *array_d[n-1]);

    Trapezoidal[n] = Trapezoidal[n-1] + Rectangular[n];
    //get max position here
     
    if(signal_is > 0.){
      if(Trapezoidal[n] > max_val){
        max_val = Trapezoidal[n];
        max_pos = n;
      }
    }else{
      if(Trapezoidal[n] < max_val){
        max_val = Trapezoidal[n];
        max_pos = n;
      }
    }
  }
   max_val = get_max_val_trapezoidal(signal_is,Trapezoidal, max_pos,max_val);
  //return signal height
  Double_t signalAmplitude = max_val/std::min(k,l);
  return std::abs(signalAmplitude);
}

template<class T>
Double_t digitalFilters<T>::trapezoidal_filter_algorithm3( UShort_t k , UShort_t m, T* trace , Double_t baseline){
  Double_t max_val =0.;
  UShort_t max_pos = 0;
  Double_t signal_is = 0.;
  UShort_t l = k +m;
  for(UShort_t n = 0; n <TRACE_SIZE; n++){
    if((n > 10 && n < 900) && (trace[n] > baseline + 100) || (trace[n] < baseline -100)){
      signal_is+=trace[n]-baseline;
    }
    
    if (n < k)
      array_u[n] =  (Double_t)trace[n]- baseline;
    else
      array_u[n] = (Double_t)trace[n]- trace[n -k] ;

    array_d[n] = array_u[n];

    if(n >= k+l)
      array_d[n] -= array_u[n -k-l];
        
    array_p[n] = array_p[n-1] +trapezoidal_shaper_m2 * array_d[n];
        
    Rectangular[n] = array_p[n] +  array_d[n] * trapezoidal_shaper_m1;
        
    Trapezoidal[n] =Trapezoidal[n-1] + Rectangular[n];
    //get max position here 
   
    if(signal_is > 0.){
      if(Trapezoidal[n] > max_val){
        max_val = Trapezoidal[n];
        max_pos = n;
      }
    }else{
      if(Trapezoidal[n] < max_val){
        max_val = Trapezoidal[n];
        max_pos = n;
      }
    }
  }
  
  max_val = get_max_val_trapezoidal(signal_is,Trapezoidal, max_pos,max_val);  //return signal height
  Double_t signalAmplitude = max_val*sampling_period /(RC_constant *k);
  return std::abs(signalAmplitude);
}
template<class T>
Double_t digitalFilters<T>::trapezoidal_filter_algorithm4( UShort_t k , UShort_t m, T* trace , Double_t baseline){
  Double_t max_val =0.;
  UShort_t max_pos = 0;
  Double_t signal_is = 0.;
  UShort_t l = k +m;
  for(UShort_t n = 0; n <TRACE_SIZE; n++){
    if((n > 10 && n < 900) && (trace[n] > baseline + 100) || (trace[n] < baseline -100)){
      signal_is+=trace[n]-baseline;
    } 
    //equation 1: dkl[n] = v[n] - v[n-k] - v[n-l] + v[n-k-l];// signal(Int_t n, Int_t m, UShort_t* v, Double_t baseline) v[n] is the signal

    array_d[n] = (Double_t)trace[n] - signal(n,0, k, trace, baseline)
                 - signal(n,0, l, trace, baseline)
                 + signal(n,0,k+l, trace, baseline);

    //equation 2: p[n] = p[n-1] + dkl[n], n>= 0
    if(n == 0) array_p[n] = array_d[n];
    else array_p[n] = array_p[n-1] + array_d[n];

    //equation 3: r[n] = p[n] + M*dkl[n]
    Rectangular[n] = array_p[n] + (1./(TMath::Exp(sampling_period/RC_constant) - 1.))*array_d[n];

    //equation 4: s[n] = s[n-1] + r[n], n>= 0
    if(n == 0) Trapezoidal[n] = Rectangular[n];
    else Trapezoidal[n] = Trapezoidal[n-1] + Rectangular[n];

    //get max position here 
   
    if(signal_is > 0.){
      if(Trapezoidal[n] > max_val){
        max_val = Trapezoidal[n];
        max_pos = n;
      }
    }else{
      if(Trapezoidal[n] < max_val){
        max_val = Trapezoidal[n];
        max_pos = n;
      }
    }
  }

  max_val = get_max_val_trapezoidal(signal_is,Trapezoidal, max_pos,max_val);
  //return signal height
  Double_t signalAmplitude = max_val*sampling_period /(RC_constant * std::min(k,l));
  return std::abs(signalAmplitude);
}
template<class T>
Double_t digitalFilters<T>::trapezoidal_filter_algorithm5( UShort_t k , UShort_t m, T* trace , Double_t baseline){
  Double_t max_val =0.;
  UShort_t max_pos = 0;
  Double_t signal_is = 0.;
  UShort_t l = k+m;
  for(UShort_t n = 0; n <TRACE_SIZE; n++){
    if((n > 10 && n < 900) && (trace[n] > baseline + 100) || (trace[n] < baseline -100)){
      signal_is+=trace[n]-baseline;
    }
  }
  
  //  int m = abs(l - k);
  for(int i=m+k+m;i<TRACE_SIZE;i++){
    Trapezoidal[i] = Trapezoidal[i-1] + (trace[i]-trace[i-m])-(trace[i-m-k]-trace[i-m-k-m]);
    //get max position here 
    if(Trapezoidal[i] > max_val){
      max_val = Trapezoidal[i];
      max_pos = i;
    }
  }
    max_val = get_max_val_trapezoidal(signal_is,Trapezoidal, max_pos,max_val);
  //return signal height
  Double_t signalAmplitude = max_val/( std::min(k,l));
  return signalAmplitude;
}

//Check the code
template<class T>
void digitalFilters<T>::pole_zero_correction(T *v){
  double p1 = 0.;
  for(int n = 0; n < 100; n++){
    p1 +=(Double_t) v[n];
  }
  double offset1 = p1/100.;

  
  Double_t d = exp(sampling_period/RC_constant);
  double g = d/(1-d);
  double m = 2.*(1-d)/(1+d);
  for(int n = 0; n < TRACE_SIZE; n++)
    array_u[n] = (double)v[n]- offset1;

  array_v[0] = array_u[0]; 
 
  for(int n = 1; n < TRACE_SIZE; n++){
    double acc =0.;
    for(int i = 0; i <= n; i++){
      acc += array_u[i];
    }
    array_v[n] = (g*array_u[n] + acc)*m;  
  }
  
}
template<class T>
void digitalFilters<T>::compare_trap_algorithms(int entry, std::string input_file_name, T *ftrace){
  UShort_t k = 0, kStart = 20, kStop = 250, kStep = 10; 
  UShort_t m = 0, mStart = 0, mStop =0, mStep = 10;
 
  Double_t ampDiff = 0.;
  UShort_t ktimes = (kStop - kStart)/kStep + 1;
  UShort_t mtimes = (mStop - mStart)/mStep + 1;
  TH1I* h = new TH1I("hsignal","signal",TRACE_SIZE,0,TRACE_SIZE);
  TH1F** h_amp = new TH1F*[5];
  TH1F** h_compT = new TH1F*[5];
  int color[5] = { 1, 3, 4, 6, 2};// black, green, blue, pink
  int style[5] = { 3002, 3003, 3004, 3005, 3006};
  for(int i = 0; i < 5; i++){
    h_amp[i] = new TH1F(Form("h_amp%d", i+1),Form("TrapAlg_%d;k; abs (real - computed)",i+1),ktimes,kStart,kStop);
    h_amp[i]->SetLineColor(color[i]);
    h_amp[i]->SetFillColor(color[i]);
    h_amp[i]->SetFillStyle(style[i]);
    h_compT[i] = new TH1F(Form("h_compT%d", 1+ i),Form("TrapAlg_%d ;k; computational time [#mu sec]", i+1),ktimes,kStart,kStop);
    h_compT[i]->SetLineColor(color[i]);
    h_compT[i]->SetFillColor(color[i]);
    h_compT[i]->SetFillStyle(style[i]);
  }
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::microseconds;
  
  UShort_t l = 0;
  double signalHeight = 500.; double trapHeight = 0.;
  /* TFile* inFile = new TFile(input_file_name.data(), "READ");
     if(!inFile->IsOpen()){
     std::cerr << "ERROR: could not open input file\n";
     std::terminate();}
     TTree* inTree = (TTree*)inFile->Get("rawDataTree");
     inTree->SetBranchAddress("Time", &ftime);
     inTree->SetBranchAddress("EventNo", &feventNo);
     inTree->SetBranchAddress("Trace", ftrace);
     inTree->SetBranchAddress("Gain", &fgain);
     inTree->SetBranchAddress("BoardID", &fboardID);
     inTree->SetBranchAddress("ChannelID", &fchannelID);
     inTree->GetEntry(entry);
     // for(int i = 0; i < 991;i++)ftrace[i] = ftrace[i+1];
     */

  generate_simulted_signal(5,2500,1000.,signalHeight, ftrace);


  //Fill signal
  for(int n = 1; n < TRACE_SIZE; n++){
    h->SetBinContent(n+1, ftrace[n]);
  }
  //compute baseline
  int M = 50;
  double p1 = 0.;
  for(int n = 0; n < M; n++){
    p1 += ftrace[n];
  }
  double offset1 = p1/M;

  //Get amplitude
  signalHeight = h->GetMaximum() - offset1;
  cout<<"signal height "<<signalHeight<<endl;
  for(UShort_t k = kStart; k <= kStop;k = k + kStep ){
    for(UShort_t m = mStart; m <= mStop; m = m + mStep){
      //      l = k +m;
      //Algorithm 1
      auto t1_1 = high_resolution_clock::now();
      trapHeight = trapezoidal_filter_algorithm1( k , m , ftrace,offset1);
      auto t2_1 = high_resolution_clock::now();
      duration<double, std::micro> us_double_1 = t2_1 - t1_1;
      ampDiff = std::abs(signalHeight - trapHeight);
      h_amp[0]->Fill(k, ampDiff);
      h_compT[0]->Fill(k,us_double_1.count());
      cout<<"Computational time of Trap 1 "<<us_double_1.count()<<endl;
      //Algorithm 2
      auto t1_2 = high_resolution_clock::now();
      trapHeight = trapezoidal_filter_algorithm2( k , m , ftrace,offset1);
      auto t2_2 = high_resolution_clock::now();
      duration<double, std::micro> us_double_2 = t2_2 - t1_2;
      ampDiff = std::abs(signalHeight - trapHeight);
      h_amp[1]->Fill(k, ampDiff);
      h_compT[1]->Fill(k,us_double_2.count());
      cout<<"Computational time of Trap 2 "<<us_double_2.count()<<endl;
      //Algorithm 3
      auto t1_3 = high_resolution_clock::now();
      trapHeight = trapezoidal_filter_algorithm3( k , m , ftrace,offset1);
      auto t2_3 = high_resolution_clock::now();
      duration<double, std::micro> us_double_3 = t2_3 - t1_3;
      ampDiff = std::abs(signalHeight - trapHeight);
      h_amp[2]->Fill(k, ampDiff);
      h_compT[2]->Fill(k,us_double_3.count());
      cout<<"Computational time of Trap 3 "<<us_double_3.count()<<endl;
      //Algorithm 4
      auto t1_4 = high_resolution_clock::now();
      trapHeight = trapezoidal_filter_algorithm4( k , m , ftrace,offset1);
      auto t2_4 = high_resolution_clock::now();
      duration<double, std::micro> us_double_4 = t2_4 - t1_4;
      ampDiff = std::abs(signalHeight - trapHeight);
      h_amp[3]->Fill(k, ampDiff);
      h_compT[3]->Fill(k,us_double_4.count());
      cout<<"Computational time of Trap 4 "<<us_double_4.count()<<endl;
      //CR-RC4
      auto t1_5 = high_resolution_clock::now();
      trapHeight = CR_RC4_filter( 5 ,  ftrace,k);
      auto t2_5 = high_resolution_clock::now();
      duration<double, std::micro> us_double_5 = t2_5 - t1_5;
      ampDiff = std::abs(signalHeight - trapHeight);
      h_amp[4]->Fill(k, ampDiff);
      h_compT[4]->Fill(k,us_double_5.count());
      cout<<"Computational time of CR-RC4 "<<us_double_5.count()<<endl;
    }
  }
  TCanvas* c = new TCanvas;
  h_amp[0]->Draw("");
  h_amp[1]->Draw("same");
  h_amp[2]->Draw("same");
  h_amp[3]->Draw("same");
  c->BuildLegend();
  TCanvas *c1 = new TCanvas;
  h_compT[0]->Draw("");
  h_compT[1]->Draw("same");
  h_compT[2]->Draw("same");
  h_compT[3]->Draw("same");
  h_compT[4]->Draw("same");
  c1->BuildLegend();
} 
