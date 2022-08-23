#include "CFDOptimization.hh"
bool sortFunc_heavy (CFDparameters i, CFDparameters j) { return (i.sigma_HF<j.sigma_HF); }

bool sortFunc_light (CFDparameters i, CFDparameters j) { return (i.sigma_LF<j.sigma_LF); }

bool sortFunc_tof (CFDparameters i, CFDparameters j) { return (i.delta_tof_mt<j.delta_tof_mt); }

CFDOptimization::CFDOptimization()
{
  draw_cfd_components = 0;
  gr = new TGraph();
  ampFitFunc = new TF1("pol3", "pol3",0,992);
  eventsProcessFlag = 0;
  n_events_to_process = 0;
      col_num =0;
    row_num =0;
  si_cfd_delay_low = 20;
  sed_cfd_delay_low = 10;
  pmt_cfd_delay_low = 10;

  si_cfd_delay_high = 20;
  sed_cfd_delay_high = 10;
  pmt_cfd_delay_high = 10;

  pmt_baseline_max = 1E6;
  sed_baseline_max = 1E6;
  si_baseline_max = 1E6;

  pmt_noise_max = 1E6;
  sed_noise_max = 1E6;
  si_noise_max = 1E6;
  Sed_CFDmethod = 2;
  Si_CFDmethod = 2;
  Pmt_CFDmethod = 2;
  si_delay_step =sampling_period;//ns
  sed_delay_step =sampling_period;//ns
  pmt_delay_step =sampling_period;//ns
  si_n_delay = (si_cfd_delay_high - si_cfd_delay_low)/si_delay_step + 1;
  sed_n_delay = (sed_cfd_delay_high - sed_cfd_delay_low)/sed_delay_step + 1;
  pmt_n_delay = (pmt_cfd_delay_high - pmt_cfd_delay_low)/pmt_delay_step + 1;
    
  si_cfd_fract_step = 0.5;
  sed_cfd_fract_step = 0.5;
  pmt_cfd_fract_step = 0.5;


  si_cfd_fract_start = 0.;
  si_cfd_fract_stop = 0.2;
    
  sed_cfd_fract_start = 0.;
  sed_cfd_fract_stop = 0.2;

  pmt_cfd_fract_start = 0.;
  pmt_cfd_fract_stop = 0.2;
    
  
  nCFD_si = (Int_t)((si_cfd_fract_stop - si_cfd_fract_start)/si_cfd_fract_step + 1);
  nCFD_sed = (Int_t)((sed_cfd_fract_stop - sed_cfd_fract_start)/sed_cfd_fract_step + 1);

  nCFD_pmt = (Int_t)((pmt_cfd_fract_stop - pmt_cfd_fract_start)/pmt_cfd_fract_step + 1);
  rand = new TRandom3;
  rand2 = new TRandom3;
  rand3 = new TRandom3;
}

CFDOptimization::~CFDOptimization()
{
  delete gr;
  delete ampFitFunc;
  delete rand;
  delete rand2;
  delete rand3;
    for(int i = 0; i < row_num;i++) if(CFD_Time[i]) delete[] CFD_Time[i];

    if(CFD_Time)delete []CFD_Time;

}


void CFDOptimization::initialize_cfd_array(int row, int col){
    col_num = col;
    row_num = row;


    CFD_Time = new double*[row];
    for(int i = 0; i < row;i++)
        CFD_Time[i] = new double[col];

    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            CFD_Time[i][j] = 0.;
        }
    }

}

void CFDOptimization::fillSignal(TH2F *htrace_sum, UShort_t* trace){
  for (UShort_t i = 0; i < TRACE_SIZE; i++) {
    htrace_sum->Fill(i+1, trace[i]);
  }
}

void CFDOptimization::perform(std::string input_file_name){

  //Extract run numbers
  std::size_t pos = input_file_name.find_first_of("_");
  string runNo = input_file_name;
  runNo.erase(0,pos+1);
  pos = runNo.length();
  runNo.erase(4,pos);

  string runNo2 = input_file_name;
  pos = runNo2.find_first_of("_");
  runNo2.erase(0,pos+2);
  pos = runNo2.length();
  runNo2.erase(1,pos);
  //output file
  TString oFileName;
  oFileName.Form("../../Analysis_plots/run_%s_r%s-Numexo2_withminPar.root", runNo.data(), runNo2.data());
  TFile* fout = new TFile(oFileName,"RECREATE");
  if(Sed_CFDmethod == 2 && Si_CFDmethod ==2)
  {TTree* hTree = new TTree("header","header");
    //attenuation info
    hTree->Branch("sed_CFD_frac_start", &sed_cfd_fract_start, "sed_CFD_frac_start/D");
    hTree->Branch("sed_CFD_frac_stop", &sed_cfd_fract_stop,"sed_CFD_frac_stop/D");
    hTree->Branch("sed_CFD_frac_step", &sed_cfd_fract_step,"sed_CFD_frac_step/D");
    hTree->Branch("si_CFD_frac_start", &si_cfd_fract_start,"si_CFD_frac_start/D");
    hTree->Branch("si_CFD_frac_stop", &si_cfd_fract_stop,"si_CFD_frac_stop/D");
    hTree->Branch("si_CFD_frac_step", &si_cfd_fract_step,"si_CFD_frac_step/D");
    //delay info
    hTree->Branch("sed_CFD_delay_start", &sed_cfd_delay_low,"sed_CFD_delay_start/I");
    hTree->Branch("sed_CFD_delay_stop", &sed_cfd_delay_high,"sed_CFD_delay_stop/I");
    hTree->Branch("sed_CFD_delay_step", &sed_delay_step,"sed_CFD_delay_step/I");
    hTree->Branch("si_CFD_delay_start", &si_cfd_delay_low,"si_CFD_delay_start/I");
    hTree->Branch("si_CFD_delay_stop", &si_cfd_delay_high,"si_CFD_delay_stop/I");
    hTree->Branch("si_CFD_delay_step", &si_delay_step,"si_CFD_delay_step/I");

    hTree->Fill();

    delete hTree;
  }
  //Graphs
  //Sed graphs
  TGraph* gr_sed_amp = new TGraph();
  gr_sed_amp->SetNameTitle("gr_sed_amp","Sed signal amplitude as a function of event number");
  TGraph* gr_sed_maxpos = new TGraph();
  gr_sed_maxpos->SetNameTitle("gr_sed_maxpos","Sed signal maxpos as a function of event number");
  TGraph* gr_sed_baseline = new TGraph;
  gr_sed_baseline->SetNameTitle("gr_sed_baseline","baseline(evt number) in Sed");
  TGraph* gr_sed_noise = new TGraph;
  gr_sed_noise->SetNameTitle("gr_sed_noise","noise(evt number) in Sed");
  TGraph* gr_sed_risetime = new TGraph;
  gr_sed_risetime->SetNameTitle("gr_sed_risetime","risetime(evt number) in Sed");

  TGraph* gr_sed_trigger = new TGraph;
  gr_sed_trigger->SetNameTitle("gr_sed_trigger","trigger(evt number) in Si");

  TGraph **gr_cfd_sed = new TGraph*[nCFD_sed];
  TGraph **gr_cfd_sed_corr = new TGraph*[nCFD_sed];
  for(int l=0;l<nCFD_sed;l++){
    gr_cfd_sed[l] = new TGraph();
    gr_cfd_sed[l]->SetNameTitle(Form("gr_cfd_sed_%d", l),Form("Sed CFD value(event number) for fraction = 0.5x%d",(l+1)));
    gr_cfd_sed[l]->GetXaxis()->SetTitle("event number");
    gr_cfd_sed[l]->GetYaxis()->SetTitle("T_{cfd} (ns)");
    gr_cfd_sed_corr[l] = new TGraph();
    gr_cfd_sed_corr[l]->SetNameTitle(Form("gr_cfd_sed_corr_%d", l),Form("Sed CFD value(event number) for fraction = 0.5x%d",(l+1)));
    gr_cfd_sed_corr[l]->GetXaxis()->SetTitle("event number");
    gr_cfd_sed_corr[l]->GetYaxis()->SetTitle("T_{cfd} (ns)");
  }

  //PMT graphs
  TGraph* gr_pmt_amp = new TGraph();
  gr_pmt_amp->SetNameTitle("gr_pmt_amp","Pmt signal amplitude as a function of event number");
  TGraph* gr_pmt_maxpos = new TGraph();
  gr_pmt_maxpos->SetNameTitle("gr_pmt_maxpos","Pmt signal maxpos as a function of event number");
  TGraph* gr_pmt_baseline = new TGraph;
  gr_pmt_baseline->SetNameTitle("gr_pmt_baseline","baseline(evt number) in Pmt");
  TGraph* gr_pmt_noise = new TGraph;
  gr_pmt_noise->SetNameTitle("gr_pmt_noise","noise(evt number) in Pmt");
  TGraph* gr_pmt_risetime = new TGraph;
  gr_pmt_risetime->SetNameTitle("gr_pmt_risetime","risetime(evt number) in Pmt");

  TGraph* gr_pmt_trigger = new TGraph;
  gr_pmt_trigger->SetNameTitle("gr_pmt_trigger","trigger(evt number) in Si");

  TGraph **gr_cfd_pmt = new TGraph*[nCFD_pmt];
  TGraph **gr_cfd_pmt_corr = new TGraph*[nCFD_pmt];
  for(int l=0;l<nCFD_pmt;l++){
    gr_cfd_pmt[l] = new TGraph();
    gr_cfd_pmt[l]->SetNameTitle(Form("gr_cfd_pmt_%d", l),Form("Pmt CFD value(event number) for fraction = 0.5x%d",(l+1)));
    gr_cfd_pmt[l]->GetXaxis()->SetTitle("event number");
    gr_cfd_pmt[l]->GetYaxis()->SetTitle("T_{cfd} (ns)");
    gr_cfd_pmt_corr[l] = new TGraph();
    gr_cfd_pmt_corr[l]->SetNameTitle(Form("gr_cfd_pmt_corr_%d", l),Form("Pmt CFD value(event number) for fraction = 0.5x%d",(l+1)));
    gr_cfd_pmt_corr[l]->GetXaxis()->SetTitle("event number");
    gr_cfd_pmt_corr[l]->GetYaxis()->SetTitle("T_{cfd} (ns)");
  }

  //Si graphs
  TGraph* gr_si_amp = new TGraph();
  gr_si_amp->SetNameTitle("gr_si_amp","Si signal amplitude as a function of event number");
  TGraph* gr_si_maxpos = new TGraph();
  gr_si_maxpos->SetNameTitle("gr_si_maxpos","Si signal maxpos as a function of event number");
  TGraph* gr_si_baseline = new TGraph;
  gr_si_baseline->SetNameTitle("gr_si_baseline","baseline(evt number) in Si");
  TGraph* gr_si_noise = new TGraph;
  gr_si_noise->SetNameTitle("gr_si_noise","noise(evt number) in Si");
  TGraph* gr_si_risetime = new TGraph;
  gr_si_risetime->SetNameTitle("gr_si_risetime","risetime(evt number) in Si");

  TGraph* gr_si_trigger = new TGraph;
  gr_si_trigger->SetNameTitle("gr_si_trigger","trigger(evt number) in Si");

  TGraph **gr_cfd_si = new TGraph*[nCFD_si];

  for(int l=0;l<nCFD_si;l++){
    gr_cfd_si[l] = new TGraph();
    gr_cfd_si[l]->SetNameTitle(Form("gr_cfd_si_%d",l),Form("Si CFD value(event number) for fraction = 0.5x%d",(l+1)));
    gr_cfd_si[l]->GetXaxis()->SetTitle("event number");
    gr_cfd_si[l]->GetYaxis()->SetTitle("T_{cfd} (ns)");
  }
  //Histograms
  //Sed histograms
  TH1F *h_ampSed= new TH1F("ampSed","Sed signal amplitude; ADC ch",512,0,16384);
  TH1F *h_baselineSed= new TH1F("baselineSed","Sed signal baseline; ADC ch",200,0,2000);
  TH1F *h_noiseSed= new TH1F("noiseSed","Sed signal noise; ADC ch",100,0,100);
  TH2F* htrace_sumSed = new TH2F("htrace_sumSed","Sed traces",TRACE_SIZE, 0, TRACE_SIZE, 2148,-800.,16384);
  TH2F* hCFD_sed = new TH2F("hCFD_sed","Sed CFD value 0.5(i+1);i; CFD value;",nCFD_sed,0,nCFD_sed, 3000, 0,3000);

  TH2F** hsed_cfdTvsE = new TH2F*[nCFD_sed];
  for(int l=0;l<nCFD_sed;l++){
    hsed_cfdTvsE[l] = new TH2F(Form("sed_cfdTvsE%d",l),Form("Sed T vs E fraction = 0.5x%d; E_{Sed}; T_{cfd}(ns);",l), 512,0,16384, 1000, 2000,3000);
  }
  TH2F* hCFD_sed_corr = new TH2F("hCFD_sed_corr","Sed CFD value 0.5(i+1); in coincidence with Si;i; CFD value;",nCFD_sed,0,nCFD_sed, 3000, 0,3000);
  TH2F** hsed_cfdTvsE_corr = new TH2F*[nCFD_sed];
  for(int l=0;l<nCFD_sed;l++){
    hsed_cfdTvsE_corr[l] = new TH2F(Form("sed_cfdTvsE_corr%d",l),Form("Sed T vs E fraction = 0.5x%d; E_{Sed}; T_{cfd}(ns);",l), 250,0,6000, 1000, 2000,3000);
  }

  //Pmt histograms
  TH1F *h_ampPmt= new TH1F("ampPmt","Pmt signal amplitude; ADC ch",512,0,16384);
  TH1F *h_baselinePmt= new TH1F("baselinePmt","Pmt signal baseline; ADC ch",200,0,2000);
  TH1F *h_noisePmt= new TH1F("noisePmt","Pmt signal noise; ADC ch",100,0,100);
  TH2F* htrace_sumPmt = new TH2F("htrace_sumPmt","Pmt traces",TRACE_SIZE, 0, TRACE_SIZE, 2148,-800.,16384);
  TH2F* hCFD_pmt = new TH2F("hCFD_pmt","Pmt CFD value 0.5(i+1);i; CFD value;",nCFD_pmt,0,nCFD_pmt, 3000, 0,3000);

  TH2F** hpmt_cfdTvsE = new TH2F*[nCFD_pmt];
  for(int l=0;l<nCFD_pmt;l++){
    hpmt_cfdTvsE[l] = new TH2F(Form("pmt_cfdTvsE%d",l),Form("Pmt T vs E fraction = 0.5x%d; E_{Pmt}; T_{cfd}(ns);",l), 512,0,16384, 1000, 2000,3000);
  }
  TH2F* hCFD_pmt_corr = new TH2F("hCFD_pmt_corr","Pmt CFD value 0.5(i+1); in coincidence with Si;i; CFD value;",nCFD_pmt,0,nCFD_pmt, 3000, 0,3000);
  TH2F** hpmt_cfdTvsE_corr = new TH2F*[nCFD_pmt];
  for(int l=0;l<nCFD_pmt;l++){
    hpmt_cfdTvsE_corr[l] = new TH2F(Form("pmt_cfdTvsE_corr%d",l),Form("Pmt T vs E fraction = 0.5x%d; E_{Pmt}; T_{cfd}(ns);",l), 250,0,6000, 1000, 2000,3000);
  }


  //Si histograms
  TH1F *h_ampSi= new TH1F("ampSi","Si signal amplitude; ADC ch",512,0,16384);
  TH1F *h_baselineSi= new TH1F("baselineSi","Si signal baseline; ADC ch",200,0,2000);
  TH1F *h_noiseSi= new TH1F("noiseSi","Si signal noise; ADC ch",100,0,100);
  TH2F* htrace_sumSi = new TH2F("htrace_sumSi","Si traces ",TRACE_SIZE, 0, TRACE_SIZE, 2148,-800,16384);
  TH2F* hCFD_si = new TH2F("hCFD_si","Si CFD value 0.5(i+1);i; CFD value;",nCFD_si,0,nCFD_si, 3000, 0,3000);
  TH2F** hsi_cfdTvsE = new TH2F*[nCFD_si];
  for(int l=0;l<nCFD_si;l++){
    hsi_cfdTvsE[l] = new TH2F(Form("si_cfdTvsE%d",l),Form("Si T vs E fraction = 0.5x%d; E_{Si}; T_{cfd}(ns);",l), 512,0,16384, 1000, 2000,3000);
  }
  //time stamp diff
  TH1F *h_delT_si_sed= new TH1F("h_delT_si_sed","#Delta Timestamp(Si-SED); ns",100,0,1000);
  //time stamp diff
  TH1F *h_delT_si_pmt= new TH1F("h_delT_si_pmt","#Delta Timestamp(Si-PMT); ns",100,0,1000);
  //time stamp diff
  TH1F *h_delT_sed_pmt= new TH1F("h_delT_sed_pmt","#Delta Timestamp(SED-PMT); ns",1000,-10000,10000);
  



  //ToF cfd sed-si
  Int_t i1,i2,i3,i4;
  TH1F *****h_tof_si_sed = new TH1F****[sed_n_delay];
  for(Int_t d1 = 0; d1 < sed_n_delay; d1++){
    h_tof_si_sed[d1] = new TH1F***[nCFD_sed];
    for(Int_t n1 = 0; n1 < nCFD_sed; n1++){
      h_tof_si_sed[d1][n1] = new TH1F**[si_n_delay];
      for(Int_t d2 = 0; d2 < si_n_delay; d2++){
        h_tof_si_sed[d1][n1][d2] = new TH1F*[nCFD_si];
        for(Int_t n2=0; n2< nCFD_si; n2++){
          i1 = d1*sed_delay_step +sed_cfd_delay_low;
          i2 = n1* sed_cfd_fract_step*100 + sed_cfd_fract_start*100;
          i3 = d2*si_delay_step +si_cfd_delay_low;
          i4 = n2*(si_cfd_fract_step*100) + si_cfd_fract_start*100;
                    
          h_tof_si_sed[d1][n1][d2][n2] =new TH1F(Form("tof_sed_d%df%d_si_d%df%d",i1,i2,i3,i4),Form("#Delta T(Si-SED) Sed:d=%dns, f= %d%, Si:d=%dns, f=%d%; ToF(ns);",i1,i2,i3,i4),400,0,200);
        }
      }
    }
  }
  
  //ToF cfd pmt-si

  TH1F *****h_tof_si_pmt = new TH1F****[pmt_n_delay];
  for(Int_t d1 = 0; d1 < pmt_n_delay; d1++){
    h_tof_si_pmt[d1] = new TH1F***[nCFD_pmt];
    for(Int_t n1 = 0; n1 < nCFD_pmt; n1++){
      h_tof_si_pmt[d1][n1] = new TH1F**[si_n_delay];
      for(Int_t d2 = 0; d2 < si_n_delay; d2++){
        h_tof_si_pmt[d1][n1][d2] = new TH1F*[nCFD_si];
        for(Int_t n2=0; n2< nCFD_si; n2++){
          i1 = d1*pmt_delay_step +pmt_cfd_delay_low;
          i2 = n1* pmt_cfd_fract_step*100 + pmt_cfd_fract_start*100;
          i3 = d2*si_delay_step +si_cfd_delay_low;
          i4 = n2*si_cfd_fract_step*100 + si_cfd_fract_start*100;
                    
          h_tof_si_pmt[d1][n1][d2][n2] =new TH1F(Form("tof_pmt_d%df%d_si_d%df%d",i1,i2,i3,i4),Form("#Delta T(Si-PMT) Pmt:d=%dns, f= %d%, Si:d=%dns, f=%d%; ToF(ns);",i1,i2,i3,i4),800,0,200);
        }
      }
    }
  }

  //ToF cfd pmt-sed

  TH1F *****h_tof_pmt_sed = new TH1F****[pmt_n_delay];
  for(Int_t d1 = 0; d1 < pmt_n_delay; d1++){
    h_tof_pmt_sed[d1] = new TH1F***[nCFD_pmt];
    for(Int_t n1 = 0; n1 < nCFD_pmt; n1++){
      h_tof_pmt_sed[d1][n1] = new TH1F**[sed_n_delay];
      for(Int_t d2 = 0; d2 < sed_n_delay; d2++){
        h_tof_pmt_sed[d1][n1][d2] = new TH1F*[nCFD_sed];
        for(Int_t n2=0; n2< nCFD_sed; n2++){
          i1 = d1*pmt_delay_step +pmt_cfd_delay_low;
          i2 =n1* pmt_cfd_fract_step*100 + pmt_cfd_fract_start*100;
          i3 = d2*sed_delay_step +sed_cfd_delay_low;
          i4 = n2*sed_cfd_fract_step*100 + sed_cfd_fract_start*100;
                    
          h_tof_pmt_sed[d1][n1][d2][n2] =new TH1F(Form("tof_pmt_d%df%d_sed_d%df%d",i1,i2,i3,i4),Form("#Delta T(PMT-SED) PMT:d=%dns, f= %d%, Sed:d=%dns, f=%d%; ToF(ns);",i1,i2,i3,i4),1600,0,100);
        }
      }
    }
  }

  //Correlation 2D histograms
  TH2F *h_EsedEsi = new TH2F("h_EsedEsi","E_{Si} vs E_{Sed};E_{Si}(ADC ch); E_{Sed}(ADC ch);", 512,0,16384,512,0,16384);

  TH2F *h_Esitof_SiSed= new TH2F("h_Esitof_SiSed","E_{Si} vs ToF (Si-Sed); E_{Si}(ADC ch); ToF(ns);",512,0,16384, 800,0,200);

  TH2F *h_Esedtof= new TH2F("h_Esedtof","E_{Sed} vs ToF; E_{Sed}(ADC ch); ToF(ns);",512,0,16384, 800,0,200);


  //si-pmt
  TH2F *h_EpmtEsi = new TH2F("h_EpmtEsi","E_{Si} vs E_{Pmt};E_{Si}(ADC ch); E_{Pmt}(ADC ch);", 512,0,16384,512,0,16384);

  TH2F *h_Esitof_SiPmt= new TH2F("h_Esitof_SiPmt","E_{Si} vs ToF(Si-PMT); E_{Si}(ADC ch); ToF(ns);",512,0,16384, 800,0,200);

  TH2F *h_Epmttof_SiPmt= new TH2F("h_Epmttof_SiPmt","E_{Pmt} vs ToF(Si-PMT); E_{Pmt}(ADC ch); ToF(ns);",512,0,16384, 800,0,200);

  //sed-pmt
  TH2F *h_EsedEpmt = new TH2F("h_EsedEpmt","E_{Pmt} vs E_{Sed};E_{Pmt}(ADC ch); E_{Sed}(ADC ch);", 512,0,16384,512,0,16384);

  TH2F *h_Epmttof_SedPmt= new TH2F("h_Epmttof_SedPmt","E_{Pmt} vs ToF(Sed-PMT); E_{Pmt}(ADC ch); ToF(ns);",512,0,16384, 800,0,200);
  TH2F *h_Esedtof_SedPmt= new TH2F("h_Esedtof_SedPmt","E_{Sed} vs ToF(Sed-PMT); E_{Sed}(ADC ch); ToF(ns);",512,0,16384, 800,0,200);
  

  //
  TGraph* gr_tof_si_sed = new TGraph();
  gr_tof_si_sed->SetNameTitle("gr_tof_si_sed","ToF(si_sed) as a function of event number");


  TGraph* gr_tof_si_pmt = new TGraph();
  gr_tof_si_pmt->SetNameTitle("gr_tof_si_pmt","ToF(si-pmt) as a function of event number");

  TGraph* gr_tof_sed_pmt = new TGraph();
  gr_tof_sed_pmt->SetNameTitle("gr_tof_sed_pmt","ToF(sed-pmt) as a function of event number");

  //input file
  TFile* inFile = new TFile(input_file_name.data(), "READ");
  if(!inFile->IsOpen()){
    fout->Delete("T;*");
    fout->Close();
    remove(oFileName.Data());
    std::cerr << "ERROR: could not open input file\n";
    std::terminate();}
  TTree* inTree = (TTree*)inFile->Get("rawDataTree");
  inTree->SetBranchAddress("Time", &ftime);
  inTree->SetBranchAddress("EventNo", &feventNo);
  inTree->SetBranchAddress("Trace", ftrace);
  inTree->SetBranchAddress("Gain", &fgain);
  inTree->SetBranchAddress("BoardID", &fboardID);
  inTree->SetBranchAddress("ChannelID", &fchannelID);
  // number of entries in the tree
  Int_t nEntries = inTree->GetEntries();
  //variables
  Int_t siCounter =0;
  Int_t sedCounter =0;
  Int_t pmtCounter =0;

  Int_t siCounter2 =0;
  Int_t siCounter3 =0;
  Int_t sedCounter2 = 0;

  DataInfo sed_signal;
  sed_signal.initialize_cfd_array(sed_n_delay,nCFD_sed);
  DataInfo pmt_signal;
  pmt_signal.initialize_cfd_array(pmt_n_delay,nCFD_pmt);

  DataInfo si_signal;
  si_signal.initialize_cfd_array(si_n_delay,nCFD_si);
  UShort_t detector;
  Double_t Si_Time, Sed_Time, Pmt_Time;
  Int_t deltaTimestamp_sed_pmt =0;
  Int_t deltaTimestamp_sed_si =0;
  Int_t deltaTimestamp_pmt_si =0; 

  Double_t deltaT_si_sed =0;
  Double_t deltaT_pmt_sed =0;
  Double_t deltaT_si_pmt =0;
  signalInfo signal_info;
  if(eventsProcessFlag && n_events_to_process < nEntries){nEntries = n_events_to_process;}
  //Start event loops
  for(Int_t evt = 0; evt < nEntries; evt++){
    SiFlag = 0; SedFlag =0;PmtFlag =0;
    inTree->GetEntry(evt);
    //inTree->Show(evt);
    printf ("\rno. of events : %d treated out of: %d events", evt, nEntries);
    fflush (stdout);

    detector = fchannelID;
    ftime = ftime*timestamp_unit;

    switch (detector) {
      case PMT:

         
        PmtFlag = 1;
        signal_info.GetSignalInfo(ftrace);
       
        pmt_signal.Amplitude = signal_info.Amplitude;
        pmt_signal.Max_pos = (UShort_t)signal_info.Max_pos;
        pmt_signal.Baseline= signal_info.Baseline;
        pmt_signal.Noise = signal_info.Noise;
        pmt_signal.RiseTime = signal_info.RiseTime;
        pmt_signal.Trigger = (UShort_t)signal_info.Trigger;
        pmt_signal.TimeStamp = ftime;
        
        if(pmt_signal.RiseTime > 10 && pmt_signal.RiseTime < 20 ) {
          pmtCounter++;
          fillSignal(htrace_sumPmt, ftrace);
          //histo
          h_ampPmt->Fill(pmt_signal.Amplitude);
          h_baselinePmt->Fill(pmt_signal.Baseline);
          h_noisePmt->Fill(pmt_signal.Noise);

          //graphs
          gr_pmt_amp->SetPoint(pmtCounter -1, pmtCounter-1, pmt_signal.Amplitude);
          gr_pmt_maxpos->SetPoint(pmtCounter -1, pmtCounter-1, pmt_signal.Max_pos);
          gr_pmt_noise->SetPoint(pmtCounter -1, pmtCounter-1, pmt_signal.Noise);
          gr_pmt_baseline->SetPoint(pmtCounter -1,pmtCounter-1, pmt_signal.Baseline);
          gr_pmt_risetime->SetPoint(pmtCounter -1,pmtCounter-1, pmt_signal.RiseTime);
          gr_pmt_trigger->SetPoint(pmtCounter -1,pmtCounter-1, pmt_signal.Trigger);

          /*
          if(Pmt_CFDmethod == 1) getCFD_Time1(ftrace, signal_info, pmt_signal.CFD_Time);
          else if(Pmt_CFDmethod == 2) getCFD_Time2(ftrace, signal_info, pmt_signal.CFD_Time);
          else if(Pmt_CFDmethod == 3) getCFD_Time3(ftrace, signal_info, pmt_signal.CFD_Time );*/
   
          for (Int_t i  = 0; i<nCFD_pmt; i++) {
            hCFD_pmt->Fill(i,pmt_signal.CFD_Time[0][i]);
            hpmt_cfdTvsE[i]->Fill(pmt_signal.Amplitude,pmt_signal.CFD_Time[0][i] );
            gr_cfd_pmt[i]->SetPoint(pmtCounter -1, pmtCounter-1, pmt_signal.CFD_Time[0][i]);
          }
        }
        break;
      case SED:
        SedFlag = 1;
        deltaTimestamp_sed_pmt = (Int_t)(ftime - pmt_signal.TimeStamp);

        signal_info.GetSignalInfo(ftrace);
        //  if(signal_info[5] > 470){
        fillSignal(htrace_sumSed, ftrace);
        sed_signal.Amplitude = signal_info.Amplitude;
        sed_signal.Max_pos = (UShort_t)signal_info.Max_pos;
        sed_signal.Baseline= signal_info.Baseline;
        sed_signal.Noise = signal_info.Noise;
        sed_signal.RiseTime = signal_info.RiseTime;
        sed_signal.Trigger = (UShort_t)signal_info.Trigger;
        sed_signal.TimeStamp = ftime;
        //	if(sed_signal.Baseline > sed_baseline_max ||sed_signal.Noise > sed_noise_max) break;
        sedCounter++;

        //histo
        h_delT_sed_pmt ->Fill(deltaTimestamp_sed_pmt);
        h_ampSed->Fill(sed_signal.Amplitude);
        h_baselineSed->Fill(sed_signal.Baseline);
        h_noiseSed->Fill(sed_signal.Noise);

        //graphs
        gr_sed_amp->SetPoint(sedCounter -1, sedCounter-1, sed_signal.Amplitude);
        gr_sed_maxpos->SetPoint(sedCounter -1, sedCounter-1, sed_signal.Max_pos);
        gr_sed_noise->SetPoint(sedCounter -1, sedCounter-1, sed_signal.Noise);
        gr_sed_baseline->SetPoint(sedCounter -1,sedCounter-1, sed_signal.Baseline);
        gr_sed_risetime->SetPoint(sedCounter -1,sedCounter-1, sed_signal.RiseTime);
        gr_sed_trigger->SetPoint(sedCounter -1,sedCounter-1, sed_signal.Trigger);
        /*
        if(Sed_CFDmethod == 1) getCFD_Time1(ftrace, signal_info, sed_signal.CFD_Time);
        else if(Sed_CFDmethod == 2) getCFD_Time2(ftrace, signal_info, sed_signal.CFD_Time);
        else if(Sed_CFDmethod == 3) getCFD_Time3(ftrace, signal_info, sed_signal.CFD_Time );*/
        for (Int_t i  = 0; i<nCFD_sed; i++) {
          hCFD_sed->Fill(i,sed_signal.CFD_Time[0][i]);
          hsed_cfdTvsE[i]->Fill(sed_signal.Amplitude,sed_signal.CFD_Time[0][i] );
          gr_cfd_sed[i]->SetPoint(sedCounter -1, sedCounter-1, sed_signal.CFD_Time[0][i]);
        }

        if(deltaTimestamp_sed_pmt >= 0 && deltaTimestamp_sed_pmt < 200){
          sedCounter2++;

          //correlations plots
          h_EsedEpmt->Fill(pmt_signal.Amplitude, sed_signal.Amplitude);
          for (Int_t d1 = 0; d1< pmt_n_delay; d1++) {
            for (Int_t n1 = 0; n1<nCFD_pmt; n1++) {
              for (Int_t d2 = 0; d2< sed_n_delay; d2++) {
                for (Int_t n2 = 0; n2<nCFD_sed; n2++) {
                  Pmt_Time = pmt_signal.CFD_Time[d1][n1] + pmt_signal.TimeStamp - pre_trig_buffer*sampling_period;
                  Sed_Time  =  sed_signal.CFD_Time[d2][n2]  + sed_signal.TimeStamp - pre_trig_buffer*sampling_period;

                  deltaT_pmt_sed =  Sed_Time - Pmt_Time;
                  h_tof_pmt_sed[d1][n1][d2][n2]->Fill(deltaT_pmt_sed );
                }
              }
            }
          } 





          Pmt_Time = pmt_signal.CFD_Time[0][0] + pmt_signal.TimeStamp - pre_trig_buffer*sampling_period;
          Sed_Time =  sed_signal.CFD_Time[0][0] +sed_signal.TimeStamp - pre_trig_buffer*sampling_period;
          deltaT_pmt_sed = Sed_Time - Pmt_Time;

          h_Epmttof_SedPmt  ->Fill(pmt_signal.Amplitude, deltaT_pmt_sed);
          h_Esedtof_SedPmt ->Fill(sed_signal.Amplitude, deltaT_pmt_sed);
          gr_tof_sed_pmt ->SetPoint(sedCounter2 -1, sedCounter2 -1, deltaT_pmt_sed);
        }
        // }
      
        break;
   
      
      case Si:
        SiFlag = 1;
        deltaTimestamp_sed_si = (Int_t)(ftime - sed_signal.TimeStamp);
        deltaTimestamp_pmt_si = (Int_t)(ftime - pmt_signal.TimeStamp);
        signal_info.GetSignalInfo(ftrace);
        if(signal_info.Trigger > 470){

          fillSignal(htrace_sumSi, ftrace);
          si_signal.Amplitude = signal_info.Amplitude;
          si_signal.Max_pos = (UShort_t)signal_info.Max_pos;
          si_signal.Baseline= signal_info.Baseline;
          si_signal.Noise = signal_info.Noise;
          si_signal.RiseTime = signal_info.RiseTime;
          si_signal.TimeStamp = ftime;
          si_signal.Trigger = (UShort_t)signal_info.Trigger;
          //	if(si_signal.Baseline > si_baseline_max ||si_signal.Noise > si_noise_max) break;
          siCounter++;
     
          //histograms
          h_delT_si_sed ->Fill(deltaTimestamp_sed_si);
          h_ampSi->Fill(si_signal.Amplitude);
          h_baselineSi->Fill(si_signal.Baseline);
          h_noiseSi->Fill(si_signal.Noise);
          //graphs
          gr_si_amp->SetPoint(siCounter -1, siCounter-1, si_signal.Amplitude);
          gr_si_maxpos->SetPoint(siCounter -1, siCounter-1, si_signal.Max_pos);
          gr_si_noise->SetPoint(siCounter -1, siCounter-1, si_signal.Noise);
          gr_si_baseline->SetPoint(siCounter-1,siCounter-1, si_signal.Baseline);
          gr_si_risetime->SetPoint(siCounter-1,siCounter-1, si_signal.RiseTime);
          gr_si_trigger->SetPoint(siCounter -1,siCounter-1, si_signal.Trigger);
          /*
          if(Si_CFDmethod ==1)getCFD_Time1(ftrace, signal_info, si_signal.CFD_Time);
          else if(Si_CFDmethod ==2)getCFD_Time2(ftrace, signal_info, si_signal.CFD_Time);
          else if(Si_CFDmethod ==3)getCFD_Time3(ftrace, signal_info, si_signal.CFD_Time);
*/
          for (Int_t i = 0; i<nCFD_si; i++) {
            //Si
            hCFD_si->Fill(i,si_signal.CFD_Time[0][i]);
            hsi_cfdTvsE[i]->Fill(si_signal.Amplitude, si_signal.CFD_Time[0][i] );
            gr_cfd_si[i]->SetPoint(siCounter2 -1, siCounter2 -1, si_signal.CFD_Time[0][i]);
          }
                
          if(deltaTimestamp_sed_si >= 0 && deltaTimestamp_sed_si < 200){
        
            siCounter2++;
        
            //correlation plots
            h_EsedEsi->Fill(si_signal.Amplitude, sed_signal.Amplitude);
            for (Int_t d1 = 0; d1< sed_n_delay; d1++) {
              for (Int_t n1 = 0; n1<nCFD_sed; n1++) {
                for (Int_t d2 = 0; d2< si_n_delay; d2++) {
                  for (Int_t n2 = 0; n2<nCFD_si; n2++) {
                    
                    Sed_Time = sed_signal.CFD_Time[d1][n1] + sed_signal.TimeStamp - pre_trig_buffer*sampling_period;
                    Si_Time =  si_signal.CFD_Time[d2][n2]  + si_signal.TimeStamp - pre_trig_buffer*sampling_period;
                  
                    deltaT_si_sed = Si_Time - Sed_Time;
                    h_tof_si_sed[d1][n1][d2][n2]->Fill(deltaT_si_sed );
                   
                  }
                }
              }
            }
           
            //choose for a given delay
          

            for (Int_t i = 0; i<nCFD_sed; i++) {
              //Sed
              hCFD_sed_corr->Fill(i,sed_signal.CFD_Time[0][i]);
              hsed_cfdTvsE_corr[i]->Fill(sed_signal.Amplitude,sed_signal.CFD_Time[0][i] );
              gr_cfd_sed_corr[i]->SetPoint(siCounter2 -1, siCounter2 -1, sed_signal.CFD_Time[0][i]);
            }
           
            Si_Time = si_signal.CFD_Time[0][0] + si_signal.TimeStamp - pre_trig_buffer*sampling_period;
            Sed_Time =  sed_signal.CFD_Time[0][0] +sed_signal.TimeStamp - pre_trig_buffer*sampling_period;
            deltaT_si_sed = Si_Time - Sed_Time;

            h_Esitof_SiSed  ->Fill(si_signal.Amplitude, deltaT_si_sed);
            h_Esedtof->Fill(sed_signal.Amplitude, deltaT_si_sed);
            gr_tof_si_sed->SetPoint(siCounter2 -1, siCounter2 -1, deltaT_si_sed);
            // if(deltaT_si_sed > 45) fillSignal(htrace_sumSi, ftrace);

          }
        }
        // pmt-si
      	if(deltaTimestamp_pmt_si >= 0 && deltaTimestamp_pmt_si < 200){
        
          siCounter3++;
	
          //correlation plots
          h_EpmtEsi->Fill(si_signal.Amplitude, pmt_signal.Amplitude);
          for (Int_t d1 = 0; d1< pmt_n_delay; d1++) {
            for (Int_t n1 = 0; n1<nCFD_pmt; n1++) {
              for (Int_t d2 = 0; d2< si_n_delay; d2++) {
                for (Int_t n2 = 0; n2<nCFD_si; n2++) {
                  Pmt_Time = pmt_signal.CFD_Time[d1][n1] + pmt_signal.TimeStamp - pre_trig_buffer*sampling_period;
                  Si_Time =  si_signal.CFD_Time[d2][n2]  + si_signal.TimeStamp - pre_trig_buffer*sampling_period;
 
                  deltaT_si_pmt = Si_Time - Pmt_Time;
                  h_tof_si_pmt[d1][n1][d2][n2]->Fill(deltaT_si_pmt );
                }
              }
            }
          }
          //choose for a given delay
	

          for (Int_t i = 0; i<nCFD_pmt; i++) {
            //PMT
            hCFD_pmt_corr->Fill(i,pmt_signal.CFD_Time[0][i]);
            hpmt_cfdTvsE_corr[i]->Fill(pmt_signal.Amplitude,pmt_signal.CFD_Time[0][i] );
            gr_cfd_pmt_corr[i]->SetPoint(siCounter3 -1, siCounter3 -1, pmt_signal.CFD_Time[0][i]);
          }

          Si_Time = si_signal.CFD_Time[0][0] + si_signal.TimeStamp - pre_trig_buffer*sampling_period;
          Pmt_Time =  pmt_signal.CFD_Time[0][0] +pmt_signal.TimeStamp - pre_trig_buffer*sampling_period;
          deltaT_si_pmt = Si_Time - Pmt_Time;

          h_Esitof_SiPmt  ->Fill(si_signal.Amplitude, deltaT_si_pmt);
          h_Epmttof_SiPmt->Fill(pmt_signal.Amplitude, deltaT_si_pmt);
          gr_tof_si_pmt->SetPoint(siCounter2 -3, siCounter3 -1, deltaT_si_pmt);

        }
    
        break;
      default:
        break;
    }
  }
  //cout <<"end of event loop "<<endl;
  inFile->cd();
  inFile->Close();
  cout<<"input file closed"<<endl;


  cout<<"Saving plots......"<<endl;
  fout->cd();
  fout->Write();
  //graphs pmt
  gr_pmt_amp->Write();
  gr_pmt_maxpos->Write();
  gr_pmt_baseline->Write();
  gr_pmt_noise->Write();
  gr_pmt_risetime->Write();
  gr_pmt_trigger->Write();
  for(int i = 0; i< nCFD_sed;i++){
    gr_cfd_pmt[i]->Write();
  }
  //correlated with respect to si
  for(int i = 0; i< nCFD_sed;i++){
    gr_cfd_pmt_corr[i]->Write();
  }
  //graphs sed
  gr_sed_amp->Write();
  gr_sed_maxpos->Write();
  gr_sed_baseline->Write();
  gr_sed_noise->Write();
  gr_sed_risetime->Write();
  gr_sed_trigger->Write();
  for(int i = 0; i< nCFD_sed;i++){
    gr_cfd_sed[i]->Write();
  }
  //correlated with respect to si
  for(int i = 0; i< nCFD_sed;i++){
    gr_cfd_sed_corr[i]->Write();
  }
  //graphs si
  gr_si_amp->Write();
  gr_si_maxpos->Write();
  gr_si_baseline->Write();
  gr_si_noise->Write();
  gr_si_risetime->Write();
  gr_si_trigger->Write();
  for(int i = 0; i< nCFD_si;i++){
    gr_cfd_si[i]->Write();
  }
  gr_tof_si_sed ->Write();
  gr_tof_sed_pmt->Write();
  gr_tof_si_pmt->Write();
  fout->Close();
  cout<<"output file closed."<<endl;

  //clear memory
  cout<<"Clearing memory.................";
  delete fout; 
  delete gr_sed_amp; 
  delete gr_sed_maxpos;
  delete gr_sed_baseline;
  delete gr_sed_noise; 
  delete gr_sed_risetime;
  delete gr_sed_trigger;
  delete [] gr_cfd_sed; 
  delete [] gr_cfd_sed_corr; 
  delete gr_pmt_amp;
  delete gr_pmt_maxpos;
  delete gr_pmt_baseline;
  delete gr_pmt_noise; 
  delete gr_pmt_risetime;
  delete gr_pmt_trigger;
  

  delete [] gr_cfd_pmt;
  delete [] gr_cfd_pmt_corr;
  delete gr_si_amp;
  delete gr_si_maxpos;
  delete gr_si_baseline;
  delete gr_si_noise;
  delete gr_si_risetime;
  delete gr_si_trigger;
  delete[] gr_cfd_si;

  delete[] hsed_cfdTvsE ;

  delete [] hsed_cfdTvsE_corr;
  delete [] hpmt_cfdTvsE;
  delete [] hpmt_cfdTvsE_corr;
  delete[] hsi_cfdTvsE;
   
  for(Int_t d1 = 0; d1 < sed_n_delay; d1++){
    for(Int_t n1 = 0; n1 < nCFD_sed; n1++){
      for(Int_t d2 = 0; d2 < si_n_delay; d2++){
        delete [] h_tof_si_sed[d1][n1][d2];
      }
    }
  }
  for(Int_t d1 = 0; d1 < sed_n_delay; d1++){
    for(Int_t n1 = 0; n1 < nCFD_sed; n1++){
      delete [] h_tof_si_sed[d1][n1];
    }
  }  
  for(Int_t d1 = 0; d1 < sed_n_delay; d1++){
    delete [] h_tof_si_sed[d1];
  } 
  delete [] h_tof_si_sed;


  for(Int_t d1 = 0; d1 < pmt_n_delay; d1++){
    for(Int_t n1 = 0; n1 < nCFD_pmt; n1++){
      for(Int_t d2 = 0; d2 < si_n_delay; d2++){
        delete [] h_tof_si_pmt[d1][n1][d2];
      }
    }
  }
  
  for(Int_t d1 = 0; d1 < pmt_n_delay; d1++){
    for(Int_t n1 = 0; n1 < nCFD_pmt; n1++){
      delete [] h_tof_si_pmt[d1][n1];
    }
  }
  
  for(Int_t d1 = 0; d1 < pmt_n_delay; d1++){
    delete [] h_tof_si_pmt[d1];
  }
  delete [] h_tof_si_pmt;
 

  for(Int_t d1 = 0; d1 < pmt_n_delay; d1++){
    for(Int_t n1 = 0; n1 < nCFD_pmt; n1++){
      for(Int_t d2 = 0; d2 < sed_n_delay; d2++){
        delete [] h_tof_pmt_sed[d1][n1][d2];
      }
    }
  }
  for(Int_t d1 = 0; d1 < pmt_n_delay; d1++){
    for(Int_t n1 = 0; n1 < nCFD_pmt; n1++){
      delete [] h_tof_pmt_sed[d1][n1];
    }
  }
  for(Int_t d1 = 0; d1 < pmt_n_delay; d1++){
    delete [] h_tof_pmt_sed[d1];
  }
  delete [] h_tof_pmt_sed;

  delete inFile;
  /*
  //root deletes these histograms 
  delete h_ampSed;
  delete h_baselineSed;
  delete h_noiseSed;
  delete htrace_sumSed;
  delete hCFD_sed;
  delete hCFD_sed_corr;
  delete h_ampPmt;
  delete h_baselinePmt;
  delete h_noisePmt;
  delete htrace_sumPmt;
  delete hCFD_pmt;
  delete hCFD_pmt_corr;
  delete h_ampSi;
  delete h_baselineSi;
  delete h_noiseSi;
  delete htrace_sumSi;
  delete hCFD_si;
  delete h_delT_si_sed;
  delete h_delT_si_pmt;
  delete h_delT_sed_pmt;
  delete h_EsedEsi;
  delete h_Esitof_SiSed;
  delete h_Esedtof;
  delete h_EpmtEsi;
  delete h_Esitof_SiPmt;
  delete h_Epmttof_SiPmt;
  delete h_EsedEpmt;
  delete h_Epmttof_SedPmt;
  delete h_Esedtof_SedPmt;
  delete gr_tof_si_sed;
  delete gr_tof_si_pmt;
  delete gr_tof_sed_pmt;

  */

  cout<<"done\n";
  
  
}

void CFDOptimization::do_optimization(std::string input_file_name){
  //W. L. Stanley, Phys Rev 113(3), 1232, (1963)
  //velocity(H) = 0.9152 cm/ns          Mass(H) = 143.61 amu
  //velocity(L) = 1.232 cm/ns          Mass(L) = 108.39 amu

  //sigma(VH) = 0.0795 cm/ns           sigma(MH) = 6.77 amu
  //sigma(VL) = 0.0665 cm/ns           sigma(LH) = 6.77 amu
  //
  // v'= v + sigma:  light = 1.2717, heavy = 0,9629
  //
  //For
  //distance = 1 cm
  //ToF(H) = 0.9152 ns, ToF'(H) = 0 ns, ToF"(H) = 1.04548 ns
  //ToF(L) = 1.232 ns, ToF'(L) = 1.2717 ns, ToF"(L) = 0.76423 ns

  Double_t velocity_light_mean = 1.2327; //cm/ns
  Double_t velocity_light_sigma = 1.2726; //cm/ns

  Double_t velocity_heavy_mean = 0.916; //cm/ns
  Double_t velocity_heavy_sigma = 0.9694; //cm/ns

  Double_t tof_heavy_mean =si_foil_distance/velocity_heavy_mean;
  Double_t tof_heavy_sigma =si_foil_distance/velocity_heavy_sigma;
  Double_t tof_light_mean =si_foil_distance/velocity_light_mean;
  Double_t tof_light_sigma =si_foil_distance/velocity_light_sigma;

  Double_t delta_tof_HL = ( tof_heavy_mean - tof_light_mean);
  //Extract run numbers
  std::size_t pos = input_file_name.find_first_of("_");
  string runNo = input_file_name;
  runNo.erase(0,pos+1);
  pos = runNo.length();
  runNo.erase(4,pos);

  string runNo2 = input_file_name;
  pos = runNo2.find_first_of("_");
  runNo2.erase(0,pos+6);
  pos = runNo2.length();
  runNo2.erase(1,pos);

  UShort_t binx1, binx2; Double_t x1, x2, y1, y2;
  Int_t i1,i2,i3,i4;
  Double_t *par = new Double_t[6];
  Double_t *err_par = new Double_t[6];
  Double_t integral1 =0., err_integral1 =0., integral2 =0., err_integral2 =0., deltax =0., err_deltax =0.;
  Double_t bin_width =0.;
  int w =15;
  //output file
  //structure for optimization
  vector<CFDparameters> data;
  TFile* infile = new TFile(input_file_name.data(),"READ");
  TTree* hTree = (TTree*)infile->Get("header");
  hTree->SetBranchAddress("sed_CFD_frac_start", &sed_cfd_fract_start);
  hTree->SetBranchAddress("sed_CFD_frac_stop", &sed_cfd_fract_stop);
  hTree->SetBranchAddress("sed_CFD_frac_step", &sed_cfd_fract_step);
  hTree->SetBranchAddress("si_CFD_frac_start", &si_cfd_fract_start);
  hTree->SetBranchAddress("si_CFD_frac_stop", &si_cfd_fract_stop);
  hTree->SetBranchAddress("si_CFD_frac_step", &si_cfd_fract_step);
  //delay info
  hTree->SetBranchAddress("sed_CFD_delay_start", &sed_cfd_delay_low);
  hTree->SetBranchAddress("sed_CFD_delay_stop", &sed_cfd_delay_high);
  hTree->SetBranchAddress("sed_CFD_delay_step", &sed_delay_step);
  hTree->SetBranchAddress("si_CFD_delay_start", &si_cfd_delay_low);
  hTree->SetBranchAddress("si_CFD_delay_stop", &si_cfd_delay_high);
  hTree->SetBranchAddress("si_CFD_delay_step", &si_delay_step);
  hTree->GetEntry(0);

  si_n_delay = (si_cfd_delay_high - si_cfd_delay_low)/si_delay_step + 1;
  sed_n_delay = (sed_cfd_delay_high - sed_cfd_delay_low)/sed_delay_step + 1;
  nCFD_si = (Int_t)((si_cfd_fract_stop - si_cfd_fract_start)/si_cfd_fract_step) +1;
  nCFD_sed = (Int_t)((sed_cfd_fract_stop - sed_cfd_fract_start)/sed_cfd_fract_step) + 1;


  std::ofstream ofs ("../../Analysis_plots/optimization_parameters.txt", std::ofstream::out);
  ofs << "sed_delay"<<setw(w)<<"sed_frac"<<setw(w)<<"si_delay"<<setw(w)<<"si_frac"<<setw(w)<<"Intg1"<<setw(w)<<"err(Intg1)"<<setw(w)<<"Intg2"<<setw(w)<<"err(Intg2)"<<setw(w)<<"d(x2-x1)"<<setw(w)<<"err(d(x2-x1))"<<endl;

  TString oFileName;
  runNo = "32";runNo2 = "0";
  oFileName.Form("../../Analysis_plots/optimization_fits_run_%s_r%s-Numexo2_new2.root", runNo.data(), runNo2.data());
  TFile* ofile = new TFile(oFileName,"RECREATE");
  TH2F* hintegrals = new TH2F("hintegrals","peak integrals; heavy; light;", 500,0,20000,500,0,20000);
  TH2F* hsigmas = new TH2F("hsigmas","peak sigma; heavy; light;", 500,0,10,500,0,10);

  TH1F* hdeltatof = new TH1F("hdeltatof","delta tof", 100,0,10);
  TH1F* hdeltatof_cut = new TH1F("hdeltatof_cut","delta tof cut:sigma:integral", 100,0,10);

  TH2F* hintegralsigma_light = new TH2F("hintegralsigmaLight","peak integral vs sigma (light); integral; sigma;", 500,0,20000,500,0,10);
  TH2F* hintegralsigma_heavy = new TH2F("hintegralsigmaHeavy","peak integral vs sigma (heavy); integral; sigma;", 500,0,20000,500,0,10);

  TH2F* hsigmatof_light = new TH2F("hsigmatof_light","sigma vs delta tof (light); sigma; delta tof;", 500,0,10,100,0,10);
  TH2F* hsigmatof_heavy = new TH2F("hsigmatof_heavy","sigma vs delta tof (heavy); sigma; delta tof;", 500,0,10,100,0,10);


  TH2F* hintegraltof_light = new TH2F("hintegraltof_light","integral vs delta tof (light); integral; delta tof;", 500,0,20000,100,0,10);
  TH2F* hintegraltof_heavy = new TH2F("hintegraltof_heavy","integral vs delta tof (heavy); integral; delta tof;", 500,0,20000,100,0,10);

  //si plots
  TH2F*** hsigmaL_sed = new TH2F**[sed_n_delay];
  for(UShort_t d1 = 0; d1 < sed_n_delay;d1++){

    hsigmaL_sed[d1] = new TH2F*[nCFD_sed];
    for(UShort_t n1 = 0; n1 < nCFD_sed;n1++){

      i1 = d1*sed_delay_step + sed_cfd_delay_low;
      i2 = (n1+1)*sampling_period;
      hsigmaL_sed[d1][n1] = new TH2F(Form("sigmaLSed_d%d_f%d",i1,i2),"sigma light; Si CFD frac; Si CFD delay", nCFD_si, si_cfd_fract_start*100,si_cfd_fract_stop*100,si_n_delay, si_cfd_delay_low,si_cfd_delay_high);
    }
  }

  TH2F*** hsigmaH_sed = new TH2F**[sed_n_delay];
  for(UShort_t d1 = 0; d1 < sed_n_delay;d1++){

    hsigmaH_sed[d1] = new TH2F*[nCFD_sed];
    for(UShort_t n1 = 0; n1 < nCFD_sed;n1++){

      i1 = d1*sed_delay_step + sed_cfd_delay_low;
      i2 = (n1+1)*sampling_period;
      hsigmaH_sed[d1][n1] = new TH2F(Form("sigmaHSed_d%d_f%d",i1,i2),"sigma heavy; Si CFD frac; Si CFD delay", nCFD_si, si_cfd_fract_start*100,si_cfd_fract_stop*100,si_n_delay, si_cfd_delay_low,si_cfd_delay_high);
    }
  }


  TH2F*** hdx_sed = new TH2F**[sed_n_delay];
  for(UShort_t d1 = 0; d1 < sed_n_delay;d1++){

    hdx_sed[d1] = new TH2F*[nCFD_sed];
    for(UShort_t n1 = 0; n1 < nCFD_sed;n1++){

      i1 = d1*sed_delay_step + sed_cfd_delay_low;
      i2 = (n1+1)*sampling_period;
      hdx_sed[d1][n1] = new TH2F(Form("dxHLSed_d%d_f%d",i1,i2),"#Delta ToF(heavy - light); Si CFD frac; Si CFD delay", nCFD_si, si_cfd_fract_start*100,si_cfd_fract_stop*100,si_n_delay, si_cfd_delay_low,si_cfd_delay_high);
    }
  }
  //sed plots
  TH2F*** hsigmaL_si = new TH2F**[si_n_delay];
  for(UShort_t d1 = 0; d1 < si_n_delay;d1++){

    hsigmaL_si[d1] = new TH2F*[nCFD_si];
    for(UShort_t n1 = 0; n1 < nCFD_si;n1++){

      i1 = d1*si_delay_step + si_cfd_delay_low;
      i2 = (n1+1)*sampling_period;
      hsigmaL_si[d1][n1] = new TH2F(Form("sigmaLSi_d%d_f%d",i1,i2),"sigma light; Sed CFD frac; Sed CFD delay", nCFD_sed, sed_cfd_fract_start*100,sed_cfd_fract_stop*100,sed_n_delay, sed_cfd_delay_low,sed_cfd_delay_high);
    }
  }

  TH2F*** hsigmaH_si = new TH2F**[si_n_delay];
  for(UShort_t d1 = 0; d1 < si_n_delay;d1++){

    hsigmaH_si[d1] = new TH2F*[nCFD_si];
    for(UShort_t n1 = 0; n1 < nCFD_si;n1++){

      i1 = d1*si_delay_step + si_cfd_delay_low;
      i2 = (n1+1)*sampling_period;
      hsigmaH_si[d1][n1] = new TH2F(Form("sigmaHSi_d%d_f%d",i1,i2),"sigma heavy; Sed CFD frac; Sed CFD delay", nCFD_sed, sed_cfd_fract_start*100,sed_cfd_fract_stop*100,sed_n_delay, sed_cfd_delay_low,sed_cfd_delay_high);
    }
  }


  TH2F*** hdx_si = new TH2F**[si_n_delay];
  for(UShort_t d1 = 0; d1 < si_n_delay;d1++){

    hdx_si[d1] = new TH2F*[nCFD_si];
    for(UShort_t n1 = 0; n1 < nCFD_si;n1++){

      i1 = d1*si_delay_step + si_cfd_delay_low;
      i2 = (n1+1)*sampling_period;
      hdx_si[d1][n1] = new TH2F(Form("dxHLSi_d%d_f%d",i1,i2),"#Delta ToF(heavy - light); Sed CFD frac; Sed CFD delay", nCFD_sed, sed_cfd_fract_start*100,sed_cfd_fract_stop*100,sed_n_delay, sed_cfd_delay_low,sed_cfd_delay_high);
    }
  }




  TCanvas* canvas = new TCanvas;
  TF1* func = new TF1("func","gaus(0)+gaus(3)",0.,200.);
  TF1* func1 = new TF1("func1","gaus",0.,200.);
  func1->SetLineColor(kBlue);
  TF1* func2 = new TF1("func2","gaus",0.,200.);
  func2->SetLineColor(kGreen);
  func->SetParNames("constant1","mean1","sigma1","constant2","mean2","sigma2");
  TH1F* htemp;
  TString histname;
  TSpectrum * spec = new TSpectrum(2);

  Double_t fPositionX[2];
  Double_t fPositionY[2];


  Double_t delta_sigma=0.; double deltax_mt =0.;Double_t integral_ratio =0.;
  Double_t  err_integral_ratio  =0.;
                             
  for(UShort_t d1 = 0; d1 < sed_n_delay  ; d1++){

    i1 = d1*sed_delay_step + sed_cfd_delay_low;

    for(UShort_t n1 = 0; n1 < nCFD_sed     ; n1++){

      i2 = (n1+1)*sampling_period;

      for(UShort_t d2 = 0; d2 < si_n_delay   ; d2++){

        i3 = d2*si_delay_step +si_cfd_delay_low;

        for(UShort_t n2 = 0; n2 < nCFD_si      ; n2++){

          i4 = (n2+1)*sampling_period;

          histname.Form("tof_sed_d%df%d_si_d%df%d",i1,i2,i3,i4);
          htemp = (TH1F*)infile->Get(histname);
          bin_width = htemp->GetBinWidth(1);
          canvas->cd();
          htemp->Draw();
          //Get maxbin

          Int_t nfound = spec->Search(htemp, 1.5, "Markov", 0.1 );
          Float_t *xpeaks = spec->GetPositionX();
          int idx[nfound];
          TMath::Sort(nfound,spec->GetPositionX(),idx,false);
          for (int i = 0; i < nfound; i++) {
            fPositionX[i] = xpeaks[idx[i]];
            Int_t bin = htemp->GetXaxis()->FindBin(fPositionX[i]);
            fPositionY[i] = htemp->GetBinContent(bin);  
          }
          //Set Range for the fit functions
          func->SetRange(fPositionX[0]-10.,fPositionX[1]+10.);
          //set parameters
          func->SetParameters(fPositionY[0],fPositionX[0],.7, fPositionY[1],fPositionX[1],1.);
          //func->SetParLimits(2, 0.1, 2.);
          //func->SetParLimits(5, 0.1, 2.);
            
          /*
            binx1 = htemp->GetMaximumBin();
            x1 = htemp->GetXaxis()->GetBinCenter(binx1);
            y1 = htemp->GetBinContent(binx1);
            //cout<<"binx1 "<<binx1<<"  x1 "<<x1<<" y1 "<<y1<<endl;
            //x2 can be either on the left or right
            //check left first
            x2 = x1 - 3.;
            binx2 = htemp->GetXaxis()->FindBin(x2);
            y2 = htemp->GetBinContent(binx2);
            //cout<<"left binx2 "<<binx2<<"  x2 "<<x2<<" y2 "<<y2<<endl;
            if((y1 -y2) > 100){
            x2 = x1 + 3.;
            binx2 = htemp->GetXaxis()->FindBin(x2);
            y2 = htemp->GetBinContent(binx2);
            //cout<<"right binx2 "<<binx2<<"  x2 "<<x2<<" y2 "<<y2<<endl;
            }
            if(x1 > x2)
            func->SetParameters(y2,x2, 1., y1, x1,1.);
            else
            func->SetParameters(y1,x1, 1., y2, x2,1.);
            htemp->GetXaxis()->SetRangeUser(x1-10.,x1+10.);
            htemp->Fit(func,"MQ","same",x1-10.,x1+10.);
          */

          htemp->GetXaxis()->SetRangeUser(fPositionX[0]-10.,fPositionX[1]+10.);
          htemp->Fit(func,"MEQR","same");
          par = func->GetParameters();
          err_par = func->GetParErrors();
          func1->SetParameters(par);
          func2->SetParameters(&par[3]);
          func1->SetRange(fPositionX[0]-10.,fPositionX[1]+10.);
          func2->SetRange(fPositionX[0]-10.,fPositionX[1]+10.);


	   

          //integral 1 light
          integral1 = par[0]*par[2]*TMath::Sqrt(2.*TMath::Pi()) / bin_width;
          err_integral1 = integral1*TMath::Sqrt((err_par[0]/par[0])*(err_par[0]/par[0]) + (err_par[2]/par[2])*(err_par[2]/par[2]));
          //integral 2 heavy
          integral2 = par[3]*par[5]*TMath::Sqrt(2.*TMath::Pi()) / bin_width;
          err_integral2 = integral1*TMath::Sqrt((err_par[3]/par[3])*(err_par[3]/par[3]) + (err_par[5]/par[5])*(err_par[5]/par[5]));
          //delta x
          deltax = par[4] - par[1];
          err_deltax = TMath::Sqrt(err_par[4]*err_par[4] + err_par[1]*err_par[1]);
          delta_sigma = par[5] - par[2];
          integral_ratio = integral1/integral2;

          err_integral_ratio = integral_ratio*TMath::Sqrt((err_integral1/integral1)*(err_integral1/integral1) +(err_integral2/integral2)*(err_integral2/integral2));
          hdeltatof->Fill(deltax);

          hintegrals->Fill(integral2, integral1);
          hsigmas->Fill(par[5], par[2]);
          hintegralsigma_light->Fill(integral1, par[2]);
          hintegralsigma_heavy->Fill(integral2, par[5]);
          hsigmatof_light->Fill(par[2], deltax);
          hsigmatof_heavy->Fill(par[5], deltax);


          hintegraltof_light->Fill(integral1, par[2]);
          hintegraltof_heavy->Fill(integral2, par[5]);
          htemp->GetListOfFunctions()->Add(func1);
          htemp->GetListOfFunctions()->Add(func2);
  		   

          if(deltax >= 4.17 && deltax <= 4.34){
            if(par[2]>=0.588 && par[2] < 1. && par[5]>= 1.265 && par[5] < 1.6 && par[5] > par[2]) {
              if((integral1 / integral2) <1.1){
                  

                hdeltatof_cut->Fill(deltax);

                // cout<<"intgral1: "<<integral1<<"  err: "<<err_integral1<<"  integral2: "<<integral2<<"  err: "<<err_integral2<<" dx "<<deltax<<"  err: "<<err_deltax<<endl;

                deltax_mt = abs(deltax-delta_tof_HL);//measured -theory
                hsigmaL_sed[d1][n1]->SetBinContent(n2+1,d2+1,par[2]);
                hsigmaH_sed[d1][n1]->SetBinContent(n2+1,d2+1,par[5]);
                hdx_sed[d1][n1]->SetBinContent(n2+1,d2+1,deltax);

                hsigmaL_si[d2][n2]->SetBinContent(n1+1,d1+1,par[2]);
                hsigmaH_si[d2][n2]->SetBinContent(n1+1,d1+1,par[5]);
                hdx_si[d2][n2]->SetBinContent(n1+1,d1+1,deltax);

                CFDparameters a(i1,i2,i3,i4, par[2], err_par[2], par[5],err_par[5], deltax, err_deltax, deltax_mt, integral_ratio,  err_integral_ratio);//light, heavy
                data.push_back(a);
                ofs << i1<<setw(w)<<i2<<setw(w)<<i3<<setw(w)<<i4<<setw(w)<<integral1<<setw(w)<<err_integral1<<setw(w)<<integral2<<setw(w)<<err_integral2<<setw(w)<<deltax<<setw(w)<<err_deltax<<endl;
                ofile->cd();
                htemp->SetName(Form("h_sed_d%df%d_si_d%df%d",i1,i2,i3,i4));
                htemp->Write();
                //canvas->SetName(Form("c_sed_d%df%d_si_d%df%d",i1,i2,i3,i4));
                //canvas->Write();

              }
            }
          }
        }
      }
    }
  }
  delete[] par;
  delete[] err_par;
 
  cout<<"===============light fragments ======================="<<endl;
  //table 1
  //Sort data to find the min value
  std::sort(data.begin(),data.end(), sortFunc_light);//sigma light fisson fragment

  cout<< "sigma light : err : sigma heavy : err : delte Tof : err : integral ratio : err : sed delay : sed frac : si delay : si frac "<<endl;
  for(int i = 0; i < 5; i++){
    printf("i: %d : LF (%f : % f) : HF (%f : %f) : dToF (% f : %f) : ratio (% f : %f) : Sed (%d : %d ): Si (%d : %d)\n",i, data[i].sigma_LF, data[i].err_sigma_LF, data[i].sigma_HF, data[i].err_sigma_HF, data[i].delta_tof, data[i].err_delta_tof, data[i].integral_ratio, data[i].err_integral_ratio, data[i].sed_delay, data[i].sed_frac, data[i].si_delay, data[i].si_frac);
  }
  //Get the electronic file run 37
  Double_t sigma_elec[5] ={0.,0,0,0,0};
  Double_t err_sigma_elec[5] = {0,0,0,0,0};

  Double_t sigma_setup_light[5]= {0,0,0,0,0};
  Double_t err_sigma_setup_light[5] = {0,0,0,0,0};


  Double_t sigma_setup_heavy[5]= {0,0,0,0,0};
  Double_t err_sigma_setup_heavy[5]= {0,0,0,0,0};

  Double_t sigma_intrinsic_light = 0.606;
  Double_t err_sigma_intrinsic_light = sqrt(0.0179*0.0179 + 0.09828*0.09828);

  Double_t sigma_intrinsic_heavy = 1.299;
  Double_t err_sigma_intrinsic_heavy = sqrt(0.03593*0.03359 + 0.1316*0.1316);

  Double_t RMS_geometry_light = 0.255;
  Double_t RMS_geometry_heavy = 0.3456;
  
  TFile* felec= new TFile("../../Analysis_plots/run_0037_r0-Numexo2-SiTM2-SedTM2.root");
  for(int i = 0; i < 5; i++){
    i1 = data[i].sed_delay;
    i2 = data[i].sed_frac;
    i3 = data[i].si_delay;
    i4 = data[i].si_frac;
    cout<< "i1 "<<i1<<"  i2 "<<i2<<"  i3 "<<i3<<"  i4 "<<i4<<endl;
    histname.Form("tof_sed_d%df%d_si_d%df%d",i1,i2,i3,i4);
    htemp = (TH1F*)felec->Get(histname);
    canvas->cd();
    htemp->Draw();
    //Get maxbin
    //set parameters
    fPositionY[0] = htemp->GetMaximum();
    int bin = htemp->GetMaximumBin();
    fPositionX[0] = htemp->GetBinCenter(bin);
    func1->SetParameters(fPositionY[0],fPositionX[0],.3);  

    htemp->GetXaxis()->SetRangeUser(fPositionX[0]-5.,fPositionX[0]+5.);
    func1->SetRange(fPositionX[0]-5.,fPositionX[0]+5.);
    htemp->Fit(func1,"MEQR","same");
    sigma_elec[i] = func1->GetParameter(2);
    err_sigma_elec[i] =  func1->GetParError(2);
    ofile->cd();
    htemp->SetName(Form("h_sed_d%df%d_si_d%df%d_elec_minLight",i1,i2,i3,i4));
    htemp->Write();

    cout<<"electronic  i = "<<i<<" sigma= "<< sigma_elec[i]<<"  err = "<< err_sigma_elec[i]<<endl;
               
  }

  //table 2
  
  for(int i = 0; i < 5; i++){
    //light fission fragment
    sigma_setup_light[i] = TMath::Sqrt(data[i].sigma_LF*data[i].sigma_LF - sigma_elec[i]*sigma_elec[i]);
    //    err_sigma_setup_light[i] = (1./ sigma_setup_light[i])*(data[i].sigma_LF*data[i].err_sigma_LF + sigma_elec[i]*err_sigma_elec[i]);


    err_sigma_setup_light[i] = (1./ sigma_setup_light[i])*TMath::Sqrt((pow(data[i].sigma_LF*data[i].err_sigma_LF,2) +pow(sigma_elec[i]*err_sigma_elec[i],2)));
     
    //heavy fission fragment
    sigma_setup_heavy[i] = TMath::Sqrt(data[i].sigma_HF*data[i].sigma_HF - sigma_elec[i]*sigma_elec[i] );
    //    err_sigma_setup_heavy[i] = (1./ sigma_setup_heavy[i])*(data[i].sigma_HF*data[i].err_sigma_HF + sigma_elec[i]*err_sigma_elec[i]);
    err_sigma_setup_heavy[i] = (1./ sigma_setup_heavy[i])*TMath::Sqrt((pow(data[i].sigma_HF*data[i].err_sigma_HF,2) + pow(sigma_elec[i]*err_sigma_elec[i],2)));
  }


  //print result
  printf("setup Resolution including intrinsic dispersion %s : %s : %s : %s : %s\n","i", "light FF", "error","heavy FF", "error");
  for(int i = 0; i <5; i++){
    printf("i: %d : %f : %f : %f : %f\n", i,  sigma_setup_light[i],  err_sigma_setup_light[i],  sigma_setup_heavy[i], err_sigma_setup_heavy[i] );
  }




  //table 3
  
  for(int i = 0; i < 5; i++){
    //light fission fragment
    sigma_setup_light[i] = TMath::Sqrt( (data[i].sigma_LF*data[i].sigma_LF) - (sigma_elec[i]*sigma_elec[i]) - (sigma_intrinsic_light*sigma_intrinsic_light));
    //    err_sigma_setup_light[i] = (1./ sigma_setup_light[i])*((data[i].sigma_LF*data[i].err_sigma_LF) + (sigma_elec[i]*err_sigma_elec[i])  + (sigma_intrinsic_light* err_sigma_intrinsic_light));


    err_sigma_setup_light[i] = (1./ sigma_setup_light[i])*TMath::Sqrt(pow(data[i].sigma_LF*data[i].err_sigma_LF,2) + pow(sigma_elec[i]*err_sigma_elec[i],2)  + pow(sigma_intrinsic_light* err_sigma_intrinsic_light,2));


    //heavy fission fragment

    sigma_setup_heavy[i] = TMath::Sqrt((data[i].sigma_HF*data[i].sigma_HF) - (sigma_elec[i]*sigma_elec[i]) -( sigma_intrinsic_heavy*sigma_intrinsic_heavy));
    //    err_sigma_setup_heavy[i] = (1./ sigma_setup_heavy[i])*((data[i].sigma_HF*data[i].err_sigma_HF) + (sigma_elec[i]*err_sigma_elec[i])  + (sigma_intrinsic_heavy* err_sigma_intrinsic_heavy));

    
    err_sigma_setup_heavy[i] = (1./ sigma_setup_heavy[i])*TMath::Sqrt(pow(data[i].sigma_HF*data[i].err_sigma_HF,2) + pow(sigma_elec[i]*err_sigma_elec[i],2)  + pow(sigma_intrinsic_heavy* err_sigma_intrinsic_heavy,2));
  }


  //print result
  printf("setup Resolution after electronic and intrinsic %s : %s : %s : %s : %s\n","i", "light FF", "error","heavy FF", "error");
  for(int i = 0; i <5; i++){
    printf("i: %d : %f : %f : %f : %f\n", i,  sigma_setup_light[i],  err_sigma_setup_light[i],  sigma_setup_heavy[i], err_sigma_setup_heavy[i] );
  }

  //table 4 result final

  TMultiGraph *mg1 = new TMultiGraph();
  TMultiGraph *mg1noise = new TMultiGraph();
  TMultiGraph *mg2noise = new TMultiGraph();
  char const *range[2] = {"light","heavy"};
  TH1F* h = new TH1F("h","",2,0.5,8.5);
  h->GetYaxis()->SetTitle("#sigma_{Si-SED} [ns]");
  h->GetXaxis()->SetTitle("Fission Fragments of ^{252}Cf");
  for(int i = 1; i < 3; i++)h->GetXaxis()->SetBinLabel(i, range[i -1]);
  h->SetMaximum(1.5);
  h->SetStats(0);
  mg1->SetNameTitle("MinLFF","minimization with light FFs");
  mg1noise->SetNameTitle("MinLFFElectronic","minimization with light FFs");
  mg2noise->SetNameTitle("MinHFFElectronic","minimization with heavy FFs");
  TGraphErrors** grl = new TGraphErrors*[5];
  int color[5] = {1,2,3,4,6};
  int style[5] = {20,21,29,33,34};

  TGraphErrors** grlnoise = new TGraphErrors*[5];
  TGraphErrors** grhnoise = new TGraphErrors*[5];
  for(int i = 0; i < 5; i++){
    grl[i] = new TGraphErrors();
    i1 = data[i].sed_delay;
    i2 = data[i].sed_frac;
    i3 = data[i].si_delay;
    i4 = data[i].si_frac;
    grl[i]->SetNameTitle(Form("sed_d%d_f%d_si_d%d_f%d", i1,i2,i3,i4), Form("SED(%d, %d) Si(%d, %d)", i1,i2,i3,i4));
    grl[i]->SetMarkerStyle(style[i]);
    grl[i]->SetMarkerColor(color[i]);
    grl[i]->SetLineColor(color[i]);
    grl[i]->SetLineWidth(2);


    grlnoise[i] = new TGraphErrors();
    grlnoise[i]->SetNameTitle(Form("ELFF_sed_d%d_f%d_si_d%d_f%d", i1,i2,i3,i4), Form("SED(%d, %d) Si(%d, %d)", i1,i2,i3,i4));
    grlnoise[i]->SetMarkerStyle(style[i]);
    grlnoise[i]->SetMarkerColor(color[i]);
    grlnoise[i]->SetLineColor(color[i]);
    grlnoise[i]->SetLineWidth(2);


  }
  
  //print result
  printf("Final setup Resolution agetr geomtery %s : %s : %s : %s : %s\n","i", "light FF", "error","heavy FF", "error");
  for(int i = 0; i <5; i++){
    double light = TMath::Sqrt(sigma_setup_light[i]* sigma_setup_light[i] - RMS_geometry_light*RMS_geometry_light);


    double heavy  = TMath::Sqrt( sigma_setup_heavy[i]* sigma_setup_heavy[i] - RMS_geometry_heavy*RMS_geometry_heavy);

    printf("i: %d : %f : %f : %f : %f\n", i, light,  err_sigma_setup_light[i],  heavy, err_sigma_setup_heavy[i] );

    /*
      light = data[i].sigma_LF;
      err_sigma_setup_light[i]= data[i].err_sigma_LF;
      heavy = data[i].sigma_HF;
      err_sigma_setup_heavy[i] =data[i].err_sigma_HF;*/

    grl[i]->SetPoint(0, (i*.5)+1, light);
    grl[i]->SetPointError(0, 0, err_sigma_setup_light[i]);
    grl[i]->SetPoint(1, (i*0.5)+5, heavy);
    grl[i]->SetPointError(1, 0, err_sigma_setup_heavy[i]);
    mg1->Add(grl[i]);


    //electronic
    grlnoise[i]->SetPoint(0, (i*.5)+1, sigma_elec[i]);
    grlnoise[i]->SetPointError(0, 0, err_sigma_elec[i]);
    mg1noise->Add(grlnoise[i]);

    
  }  
  h->Draw("");
  mg1->Draw("CP");
  //sigma heavy fisson fragment
  cout<<"==================heavy fragments ===================="<<endl;
  //table 1

  std::sort(data.begin(),data.end(), sortFunc_heavy);
 
  cout<< "sigma light : err : sigma heavy : err : delte Tof : err : integral ratio : err : sed delay : sed frac : si delay : si frac "<<endl;
  for(int i = 0; i < 5; i++){
    printf("i: %d : LF (%f : % f) : HF (%f : %f) : dToF (% f : %f) : ratio (% f : %f) : Sed (%d : %d ): Si (%d : %d)\n",i, data[i].sigma_LF, data[i].err_sigma_LF, data[i].sigma_HF, data[i].err_sigma_HF, data[i].delta_tof, data[i].err_delta_tof, data[i].integral_ratio, data[i].err_integral_ratio, data[i].sed_delay, data[i].sed_frac, data[i].si_delay, data[i].si_frac);
  }

  for(int i = 0; i < 5; i++){
    i1 = data[i].sed_delay;
    i2 = data[i].sed_frac;
    i3 = data[i].si_delay;
    i4 = data[i].si_frac;
    cout<< "i1 "<<i1<<"  i2 "<<i2<<"  i3 "<<i3<<"  i4 "<<i4<<endl;
    histname.Form("tof_sed_d%df%d_si_d%df%d",i1,i2,i3,i4);
    htemp = (TH1F*)felec->Get(histname);
    canvas->cd();
    htemp->Draw();
    //Get maxbin
    //set parameters
    fPositionY[0] = htemp->GetMaximum();
    int bin = htemp->GetMaximumBin();
    fPositionX[0] = htemp->GetBinCenter(bin);
    func1->SetParameters(fPositionY[0],fPositionX[0],.3);  

    htemp->GetXaxis()->SetRangeUser(fPositionX[0]-5.,fPositionX[0]+5.);
    func1->SetRange(fPositionX[0]-5.,fPositionX[0]+5.);
    htemp->Fit(func1,"MEQR","same");
    sigma_elec[i] = func1->GetParameter(2);
    err_sigma_elec[i] =  func1->GetParError(2);
    ofile->cd();
    htemp->SetName(Form("h_sed_d%df%d_si_d%df%d_elec_minLight",i1,i2,i3,i4));
    htemp->Write();

    cout<<"electronic  i = "<<i<<" sigma= "<< sigma_elec[i]<<"  err = "<< err_sigma_elec[i]<<endl;
               
  }

  //table 2
  
  for(int i = 0; i < 5; i++){
    //light fission fragment
    sigma_setup_light[i] = TMath::Sqrt(data[i].sigma_LF*data[i].sigma_LF - sigma_elec[i]*sigma_elec[i]);
    //err_sigma_setup_light[i] = (1./ sigma_setup_light[i])*(data[i].sigma_LF*data[i].err_sigma_LF + sigma_elec[i]*err_sigma_elec[i]);

    err_sigma_setup_light[i] = (1./ sigma_setup_light[i])*TMath::Sqrt(pow(data[i].sigma_LF*data[i].err_sigma_LF,2) + pow(sigma_elec[i]*err_sigma_elec[i],2));


    //heavy fission fragment
    sigma_setup_heavy[i] = TMath::Sqrt(data[i].sigma_HF*data[i].sigma_HF - sigma_elec[i]*sigma_elec[i] );
    //  err_sigma_setup_heavy[i] = (1./ sigma_setup_heavy[i])*(data[i].sigma_HF*data[i].err_sigma_HF + sigma_elec[i]*err_sigma_elec[i]);

    err_sigma_setup_heavy[i] = (1./ sigma_setup_heavy[i])*TMath::Sqrt((pow(data[i].sigma_HF*data[i].err_sigma_HF,2) + pow(sigma_elec[i]*err_sigma_elec[i],2)));
  }


  //print result
  printf("setup Resolution including intrinsic dispersion %s : %s : %s : %s : %s\n","i", "light FF", "error","heavy FF", "error");
  for(int i = 0; i <5; i++){
    printf("i: %d : %f : %f : %f : %f\n", i,  sigma_setup_light[i],  err_sigma_setup_light[i],  sigma_setup_heavy[i], err_sigma_setup_heavy[i] );
  }

  //table 3
  
  for(int i = 0; i < 5; i++){
    //light fission fragment
    sigma_setup_light[i] = TMath::Sqrt( (data[i].sigma_LF*data[i].sigma_LF) - (sigma_elec[i]*sigma_elec[i]) - (sigma_intrinsic_light*sigma_intrinsic_light));
    //    err_sigma_setup_light[i] = (1./ sigma_setup_light[i])*((data[i].sigma_LF*data[i].err_sigma_LF) + (sigma_elec[i]*err_sigma_elec[i])  + (sigma_intrinsic_light* err_sigma_intrinsic_light));


    err_sigma_setup_light[i] = (1./ sigma_setup_light[i])*TMath::Sqrt(pow(data[i].sigma_LF*data[i].err_sigma_LF,2) + pow(sigma_elec[i]*err_sigma_elec[i],2)  + pow(sigma_intrinsic_light* err_sigma_intrinsic_light,2));

    //heavy fission fragment

    sigma_setup_heavy[i] = TMath::Sqrt((data[i].sigma_HF*data[i].sigma_HF) - (sigma_elec[i]*sigma_elec[i]) -( sigma_intrinsic_heavy*sigma_intrinsic_heavy));
    //    err_sigma_setup_heavy[i] = (1./ sigma_setup_heavy[i])*((data[i].sigma_HF*data[i].err_sigma_HF) + (sigma_elec[i]*err_sigma_elec[i])  + (sigma_intrinsic_heavy* err_sigma_intrinsic_heavy));


    err_sigma_setup_heavy[i] = (1./ sigma_setup_heavy[i])*TMath::Sqrt(pow(data[i].sigma_HF*data[i].err_sigma_HF,2) + pow(sigma_elec[i]*err_sigma_elec[i],2)  + pow(sigma_intrinsic_heavy* err_sigma_intrinsic_heavy,2));
  }


  //print result
  printf("setup Resolution %s : %s : %s : %s : %s\n","i", "light FF", "error","heavy FF", "error");
  for(int i = 0; i <5; i++){
    printf("i: %d : %f : %f : %f : %f\n", i,  sigma_setup_light[i],  err_sigma_setup_light[i],  sigma_setup_heavy[i], err_sigma_setup_heavy[i] );
  }



  //table 4 result final

  TMultiGraph *mg2 = new TMultiGraph();
  mg2->SetNameTitle("MinHFF","minimization with heavy FFs");
  TGraphErrors** grh = new TGraphErrors*[5];

  for(int i = 0; i < 5; i++){
    grh[i] = new TGraphErrors();
    i1 = data[i].sed_delay;
    i2 = data[i].sed_frac;
    i3 = data[i].si_delay;
    i4 = data[i].si_frac;
    grh[i]->SetNameTitle(Form("sed_d%d_f%d_si_d%d_f%d", i1,i2,i3,i4), Form("SED(%d, %d) Si(%d, %d)", i1,i2,i3,i4));
    grh[i]->SetMarkerStyle(style[i]);
    grh[i]->SetMarkerColor(color[i]);
    grh[i]->SetLineColor(color[i]);
    grh[i]->SetLineWidth(2);

    grhnoise[i] = new TGraphErrors();
    grhnoise[i]->SetNameTitle(Form("EHFF_sed_d%d_f%d_si_d%d_f%d", i1,i2,i3,i4), Form("SED(%d, %d) Si(%d, %d)", i1,i2,i3,i4));
    grhnoise[i]->SetMarkerStyle(style[i]);
    grhnoise[i]->SetMarkerColor(color[i]);
    grhnoise[i]->SetLineColor(color[i]);
    grhnoise[i]->SetLineWidth(2);
  }
  
  //print result
  printf("Final setup Resolution %s : %s : %s : %s : %s\n","i", "light FF", "error","heavy FF", "error");
  for(int i = 0; i <5; i++){
    double light = TMath::Sqrt(sigma_setup_light[i]* sigma_setup_light[i] - RMS_geometry_light*RMS_geometry_light);


    double heavy  = TMath::Sqrt( sigma_setup_heavy[i]* sigma_setup_heavy[i] - RMS_geometry_heavy*RMS_geometry_heavy);

    printf("i: %d : %f : %f : %f : %f\n", i, light,  err_sigma_setup_light[i],  heavy, err_sigma_setup_heavy[i] );
    /*light = data[i].sigma_LF;
      err_sigma_setup_light[i]= data[i].err_sigma_LF;
      heavy = data[i].sigma_HF;
      err_sigma_setup_heavy[i]= data[i].err_sigma_HF;*/
    
	
	
    grh[i]->SetPoint(0, (i*.5)+1, light);
    grh[i]->SetPointError(0, 0, err_sigma_setup_light[i]);
    grh[i]->SetPoint(1, (i*.5)+5, heavy);
    grh[i]->SetPointError(1, 0, err_sigma_setup_heavy[i]);
    mg2->Add(grh[i]);

    //electronic
    grhnoise[i]->SetPoint(0, (i*.5)+1, sigma_elec[i]);
    grhnoise[i]->SetPointError(0, 0, err_sigma_elec[i]);
    mg2noise->Add(grhnoise[i]);

  }
  // TCanvas* c22 = new TCanvas;
  h->Draw("");
  mg2->Draw("CP"); 
  
  /* std::sort(data.begin(),data.end(), sortFunc_tof);//sigma delta tofs
     cout<<"================delta tofs ======================"<<endl;
     for(int i = 0; i < 5; i++)
     cout<< "i "<<i<< "smallest delta tofs: "<<data[i].delta_tof_mt<< " sigma light: "<<data[i].sigma_LF<<"  sigma heavy: "<<data[i].sigma_HF<<"  delta tof "<<data[i].delta_tof<<" sed delay: "<<data[i].sed_delay<< "  sed frac: "<<data[i].sed_frac<< "  si delay "<<data[i].si_delay<<"  si frac: "<<data[i].si_frac<<"  integral (light/heavy) "<<data[i].integral_ratio<<endl;
  */


  data.clear();


  //clear memory
  felec->cd();
  felec->Close();
  
  
  ofile->cd();
  h->Write();
  mg1->Write();
  mg2->Write();
  mg1noise->Write();
  mg2noise->Write();
  ofile->Close();
  cout<<"spectra saved.....output file closed"<<endl;
  ofs.close();
  cout<<"output txt file closed"<<endl;

  infile->cd();
  infile->Close();

  cout<<"input file closed"<<endl;

 
  

  cout<<"clearing memory......";
  delete [] grh;
  delete [] grl;
  delete [] grlnoise;
  delete mg1;
  delete mg2;
  delete mg1noise;
  delete infile;
  delete ofile;
  //delete func;
  // delete spec;
  for(UShort_t d1 = 0; d1 < sed_n_delay;d1++){
    delete [] hsigmaL_sed[d1];
    delete [] hsigmaH_sed[d1];
    delete []  hdx_sed[d1]; 
  }
  delete [] hsigmaL_sed;
  delete [] hsigmaH_sed;
  delete []  hdx_sed; 
   
  for(UShort_t d1 = 0; d1 < si_n_delay;d1++){
    delete [] hsigmaL_si[d1];
    delete [] hsigmaH_si[d1];
    delete []  hdx_si[d1]; 
  } 
  delete [] hsigmaL_si;
  delete [] hsigmaH_si;
  delete []  hdx_si; 

  delete canvas; 
  /*
  //get deleted when you close the files
  delete hintegrals;
  delete hsigmas; 

  delete hdeltatof;
  delete hdeltatof_cut;

  delete hintegralsigma_light;
  delete hintegralsigma_heavy;

  delete hsigmatof_light;
  delete hsigmatof_heavy;

  delete func1;
  delete func2;
  delete hintegraltof_light;
  delete hintegraltof_heavy;*/
  
  cout<<"done\n";


}

void CFDOptimization::drawTrace(Int_t evt, std::string input_file_name){
  Int_t CFDmethod =2;
  TGraph* gr_sig = new TGraph();
  gr_sig->SetName("gr_sig");

  hcfd_fracvin = new TGraph();
  hcfd_fracvin->SetNameTitle("fracvin","f*V_{in}");
  hcfd_fracvin->SetLineColor(4);
  hcfd_vin = new TGraph();
  hcfd_vin->SetNameTitle("vin","V_{in}");
  hcfd_vin->SetLineColor(1);
  hcfd_vindelayed = new TGraph();
  hcfd_vindelayed->SetNameTitle("vindelayed","V_{in}(t - t_{d})");
  hcfd_vindelayed->SetLineColor(6);
  hcfd  = new TGraph();
  hcfd->SetNameTitle("vout","V_{out}");
  hcfd->SetLineColor(2);



  TFile* inFile = new TFile(input_file_name.data(), "READ");
  if(!inFile) std::terminate();
  TTree* inTree = (TTree*)inFile->Get("rawDataTree");

  inTree->SetBranchAddress("Time", &ftime);
  inTree->SetBranchAddress("EventNo", &feventNo);
  inTree->SetBranchAddress("Trace", ftrace);
  inTree->SetBranchAddress("Gain", &fgain);
  inTree->SetBranchAddress("BoardID", &fboardID);
  inTree->SetBranchAddress("ChannelID", &fchannelID);

  DataInfo sed_signal;
  sed_signal.initialize_cfd_array(sed_n_delay,nCFD_sed);
  DataInfo si_signal;
  si_signal.initialize_cfd_array(si_n_delay,nCFD_si);


  inTree->GetEntry(evt);
  draw_cfd_components =1;
  if(fchannelID == 0) CFDmethod = Sed_CFDmethod;
  else if(fchannelID == 2) CFDmethod = Si_CFDmethod;
  //  signalInfo signal_info; signal_info.etSignalInfo(ftrace);
  for(UShort_t i = 0;i < TRACE_SIZE;i++)gr_sig->SetPoint(i, (Double_t)i, ftrace[i]);
  /*if(CFDmethod ==1)getCFD_Time1(ftrace, signal_info, si_signal.CFD_Time);
  else if(CFDmethod ==2)getCFD_Time2(ftrace, signal_info, si_signal.CFD_Time);
  else if(CFDmethod ==3)getCFD_Time3(ftrace, signal_info, si_signal.CFD_Time);
*/

  TCanvas * c = new TCanvas;
  c->Divide(2);
  c->cd(1);
  gr_sig->Draw();
  c->cd(2);
  gr->Draw();
  c->Update();

  TCanvas *c2 = new TCanvas;

  hcfd_vin->Draw("");
  hcfd_vindelayed->Draw("same");
  hcfd_fracvin->Draw("same");
  hcfd->Draw("same");

}

void CFDOptimization::test(Int_t evt, std::string input_file_name){

  TH1I* htrace = new TH1I("trace","trace", TRACE_SIZE,0, TRACE_SIZE);
  TH1I* htrigger = new TH1I("trigger","trigger", TRACE_SIZE,0, TRACE_SIZE);

  TH1I* htrace_diff = new TH1I("trace_diff","trace_diff", TRACE_SIZE,0, TRACE_SIZE);

  TH1I* htrace_diff_sum = new TH1I("trace_diff_sum","trace_diff_sum", TRACE_SIZE,0, TRACE_SIZE);

  TFile* inFile = new TFile(input_file_name.data(), "READ");
  if(!inFile) std::terminate();
  TTree* inTree = (TTree*)inFile->Get("rawDataTree");

  inTree->SetBranchAddress("Time", &ftime);
  inTree->SetBranchAddress("EventNo", &feventNo);
  inTree->SetBranchAddress("Trace", ftrace);
  inTree->SetBranchAddress("Gain", &fgain);
  inTree->SetBranchAddress("BoardID", &fboardID);
  inTree->SetBranchAddress("ChannelID", &fchannelID);

  inTree->GetEntry(evt);
  while(fchannelID ==0){evt++;inTree->GetEntry(evt);}
  Int_t diff;
  Double_t baseline = 0., noise =0.;
  //calculate trigger
  Double_t mv_avg =0. ,accu =0.,accu2 =0., mv_sigma =0.;
  for (UShort_t i = 1; i < TRACE_SIZE; i++) {
    accu += (Double_t)ftrace[i];

    mv_avg = accu/(Double_t)i;

    accu2 += ((Double_t)ftrace[i] - mv_avg)*((Double_t)ftrace[i] - mv_avg);
    mv_sigma = TMath::Sqrt(accu2/(Double_t)(i-1));

    if(ftrace[i] > (mv_avg + 5.*mv_sigma)){baseline = mv_avg; noise = mv_sigma;break;}

    htrace->SetBinContent(i+1,ftrace[i]);
    //if(i > 1){
    //diff = ftrace[i] - ftrace[i-1];
    htrace_diff->SetBinContent(i+1,mv_avg);
    //}
    htrace_diff_sum->SetBinContent(i+1,mv_sigma);
  }



  TCanvas *c2 = new TCanvas;
  //c2->Divide(3);
  //c2->cd(1);
  htrace->Draw();
  htrigger->SetLineColor(8);
  htrigger->Draw("same");
  //c2->cd(2);
  htrace_diff->SetLineColor(2);
  htrace_diff->Draw("same");
  //c2->cd(3);
  htrace_diff_sum->SetLineColor(6);
  htrace_diff_sum->Draw("same");
  cout<<"method2 baseline "<<baseline<<" noise "<<noise<<endl;
  //signalInfo sgnal_info; signal_info.GetSignalInfo(ftrace);

  //cout<<"baseline: "<<signal_info[2]<<"  noise: "<<signal_info[3]<<endl;
}

void CFDOptimization::monteCarlo(){
  TCanvas* c = new TCanvas;
  TCanvas* c1 = new TCanvas;
  TCanvas* c2 = new TCanvas;


  TH1F* h_light = new TH1F("h_light","light", 1000,0,20);

  TH1F* h_heavy = new TH1F("h_heavy","heavy", 1000,0,20);

 
  
  double distance_foil_source = 18.;
  double distance_foil_Si = 15.;
  double Si_diameter = 1.;
  double angleL = TMath::ATan(Si_diameter/(distance_foil_source + distance_foil_Si));
  double  angle =0.;
  double height = 0.;
  double distance = 0.;
  double v_light = 1.2351;
  double v_heavy = 0.91455;

  double ToF_light = 0;
  double ToF_heavy = 0;

  TH1F* h_angle = new TH1F("h_angle","angle", 1000,-angleL,angleL);

  // cout<<"angleL : "<<angleL<<"  "<<33.*TMath::Tan(angleL)<<endl;
  for(int i = 0; i < 10000; i++){
    angle = rand->Gaus(0., angleL);
    h_angle->Fill(angle);
    height = distance_foil_source*TMath::Tan(angle);
    distance = distance_foil_Si - height;
    ToF_light = distance / v_light;
    ToF_heavy = distance / v_heavy;
    h_light->Fill(ToF_light);
    h_heavy->Fill(ToF_heavy);

  }

  c->cd();
  h_angle->Draw();
  c1->cd();
  h_light->Draw();
  c2->cd();
  h_heavy->Draw();
}
void CFDOptimization::monteCarlo2(){
  TFile* fout = new TFile("monteCarlo_geometry_gaus.root","recreate");

  /* TCanvas* c = new TCanvas;
     TCanvas* c1 = new TCanvas;
     TCanvas* c2 = new TCanvas;
     TCanvas* c11 = new TCanvas;
     TCanvas* c22 = new TCanvas;
  */
  TH1F* h_light = new TH1F("h_light","ToF light", 1000,0,20);

  TH1F* h_heavy = new TH1F("h_heavy","ToF heavy", 1000,0,20);

  TH1F* h_Elight = new TH1F("h_Elight","E light", 1000,50,100);

  TH1F* h_Eheavy = new TH1F("h_Eheavy","E heavy", 1000,50,100);
 
  TH1F* h_Vlight = new TH1F("h_Vlight","v light", 1000,0,2);

  TH1F* h_Vheavy = new TH1F("h_Vheavy","v heavy", 1000,0,2);
  TH1F * h_distance = new TH1F("distance","distance from foil to Si", 1000,14,16);
  double distance_foil_source = 18.;
  double distance_foil_Si = 15.;
  double Si_diameter = 1.;
  double angleL = TMath::ATan(Si_diameter/(distance_foil_source + distance_foil_Si));
  double  angle =0.;
  double height = 0.;
  double distance = 0.;
  double v_light = 1.2351;
  double v_heavy = 0.91455;
  double E_light_mean = 85.57;
  double E_light_sigma = 0.2259;
  double E_heavy_mean =62.16;
  double E_heavy_sigma = 0.2818;
  double mass_light = 108.39;
  double mass_heavy = 143.61;
  double ToF_light = 0;
  double ToF_heavy = 0;
  double E_light, E_heavy;
  TH1F* h_angle = new TH1F("h_angle","angle", 1000,-angleL,angleL);

  // cout<<"angleL : "<<angleL<<"  "<<33.*TMath::Tan(angleL)<<endl;
  for(int i = 0; i < 10000; i++){
    angle = rand->Gaus(0, angleL);
    h_angle->Fill(angle);
    height = distance_foil_source*TMath::Tan(angle);
    distance = distance_foil_Si - height;
    h_distance->Fill(distance);
    E_light = rand2->Gaus(E_light_mean, E_light_sigma);
    E_heavy = rand3->Gaus(E_heavy_mean, E_heavy_sigma);

    h_Elight->Fill(E_light);
    h_Eheavy->Fill(E_heavy);
   
    v_light = TMath::Sqrt((2.*E_light) / (mass_light*931.5)) *30.;
    v_heavy = TMath::Sqrt((2.*E_heavy) / (mass_heavy*931.5)) *30.;

    h_Vlight->Fill(v_light);
    h_Vheavy->Fill(v_heavy);
    
    ToF_light = distance / v_light;
    ToF_heavy = distance / v_heavy;
    h_light->Fill(ToF_light);
    h_heavy->Fill(ToF_heavy);

  }
  fout->cd();
  fout->Write();
  fout->Close();
  /* c->cd();
     h_angle->Draw();
     c1->cd();
     h_light->Draw();
     c2->cd();
     h_heavy->Draw();
     c11->cd();
     h_Elight->Draw();
     c22->cd();
     h_Eheavy->Draw();

     TCanvas* c3 = new TCanvas;
     h_Vlight->Draw();
     TCanvas* c4 = new TCanvas;
     h_Vheavy->Draw();
     TCanvas* c5 = new TCanvas;
     h_distance->Draw();*/
}
