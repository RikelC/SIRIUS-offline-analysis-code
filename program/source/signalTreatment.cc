#include "signalTreatment.hh"

signalTreatment::signalTreatment()
{
  myFilter = new digitalFilters<UShort_t>();
 rand = new TRandom3;
 
  kStart= 50;
  kStop = 50;
  kStep = 10;
  mStart = 0;
  mStop = 0;
  mStep = 10;
  kTimes = (kStop - kStart)/kStep + 1;
  mTimes = (mStop - mStart)/mStep + 1;
}

signalTreatment::~signalTreatment()
{
  delete myFilter;
    delete rand;
}

template<class T>
void signalTreatment::invert_signal(T* v){
  for(UShort_t i = 0; i < TRACE_SIZE; i++){
    v[i] = 20000 - v[i];
  }
}


void signalTreatment::drawEntry(int entry, std::string input_file_name){
  TH1I* h = new TH1I("hsignal","signal",TRACE_SIZE,0,TRACE_SIZE);
  TH1F* hdk = new TH1F("hdkl","",TRACE_SIZE,0,TRACE_SIZE);
  TH1F* hp = new TH1F("hp","",TRACE_SIZE,0,TRACE_SIZE);
  TH1F* hr = new TH1F("hrec","",TRACE_SIZE,0,TRACE_SIZE);
  TH1F* htrap = new TH1F("htrap","",TRACE_SIZE,0,TRACE_SIZE);
  TFile* inFile = new TFile(input_file_name.data(), "READ");
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
  cout<<"no of entries in the run: "<<inTree->GetEntries()<<endl;
  //for(int i = 0; i < 991;i++)ftrace[i] = ftrace[i+1];
     
  // boardIdx = getboardIndex[fboardID];
  // if(entry < 0)filter->generate_simulted_signal(5,1500,1000.,-500.);
  invert_signal(ftrace);
  
  int m = 50;
  double p1 = 0.;
  for(int n = 0; n < m; n++){
    p1 += ftrace[n];
  }
  double offset1 = p1/m;
        
        
  //Trapezoidal Filtering
  UShort_t  k = 100;
  UShort_t l = 200;
    
  //  Double_t amp =CR_RC4_filter(10, ftrace, 100. );
  // pole_zero_correction(ftrace);
  //offset1 = 0.;
  //
  Double_t amp = myFilter->trapezoidal_filter_algorithm1( k , m , ftrace, offset1);
  //  if and only if  k <= signal max position  k and l can be swapped in algorithms 1 and 2
  //otherwise l > k
  // for algorithms 3 and 4 there is no such constraints
  cout<<"amp "<<amp<<endl;
      
  for(unsigned int n = 1; n < TRACE_SIZE; n++){
    h->SetBinContent(n, ftrace[n]);
    hdk->Fill(n, myFilter->get_array_d()[n]);

    hr->Fill(n, myFilter->get_array_Rectangular()[n]);

    htrap->Fill(n,myFilter->get_array_Trapezoidal()[n]);
  }
  cout<<"signal ampli "<<h->GetMaximum()-offset1<<endl;
  // cout<<" CR amp "<<htrap->GetMaximum()/100.<<endl;
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(4);
  c1->cd(1);
  h->Draw("hist");
  c1->cd(2);
  hdk->Draw("hist");
  c1->cd(3);
  hr->Draw("hist");
  c1->cd(4);
  htrap->Draw("hist");
  c1->Update();
}


void signalTreatment::drawRun(int al, std::string inputFileName){
  std::size_t pos = inputFileName.find_first_of("_");
  string runNo = inputFileName;
  runNo.erase(0,pos+1);
  pos = runNo.length();
  runNo.erase(4,pos);

  string runNo2 = inputFileName;
  pos = runNo2.find_first_of("-");
  runNo2.erase(0,pos-1);
  pos = runNo2.length();
  runNo2.erase(1,pos);
  
  //output file
  TString oFileName;
  oFileName.Form("../results/Spectra_run_%s_r%s.root", runNo.data(), runNo2.data());
    
  TFile* fout= new TFile(oFileName,"RECREATE");
  TH2F*** hsignal = new TH2F**[NBOARDS];
  TH2F*** htraceSum = new TH2F**[NBOARDS];
  TGraph*** gr_baseline = new TGraph**[NBOARDS];
  TGraph*** gr_noise = new TGraph**[NBOARDS];
  // TH1F*** hAmp = new TH1F**[NBOARDS];
  UShort_t board = 0;
  int **counter = new int*[NBOARDS];
  for(int i = 0; i < NBOARDS; i++){
    hsignal[i] = new TH2F*[NCHANNELS];
    htraceSum[i] = new TH2F*[NCHANNELS];
    // hAmp[i] = new TH1F*[NCHANNELS];
    gr_baseline[i] = new TGraph*[NCHANNELS];
    gr_noise[i] = new TGraph*[NCHANNELS];
    UShort_t board = NUMEXO2[i];
    counter[i] = new int[NCHANNELS];
    for(int j = 0; j < NCHANNELS;j++){
      hsignal[i][j] = new TH2F(Form("hsignal_b%d_ch%d",board, j),Form("signal, b%d ch%d",board, j),TRACE_SIZE,0,TRACE_SIZE, 2000,0,20000);
      htraceSum[i][j] = new TH2F(Form("htraceSum_b%d_ch%d",board, j),Form("signal, b%d ch%d",board, j),TRACE_SIZE,0,TRACE_SIZE, 2000,0,20000);
      //    hAmp[i][j] = new TH1F(Form("hAmp_b%d_ch%d",board,j),Form("Amp: b%d ch%d", board, j), 2000,0,20000);
      gr_baseline[i][j] = new TGraph();
      gr_baseline[i][j]->SetNameTitle(Form("gr_baseline_b%d_ch%d",board, j),Form("gr_baseline_b%d_ch%d",board, j));
      gr_noise[i][j] = new TGraph();
      gr_noise[i][j]->SetNameTitle(Form("gr_noise_b%d_ch%d",board, j),Form("gr_noise_b%d_ch%d",board, j));

      counter[i][j] = 0;
    }
  }      
  TH2F* hAmp = new TH2F("hAmpFilter","amp from filter", 2000,0,2000, NSTRIPS,0,NSTRIPS);
  TH2F* hAmpRaw = new TH2F("hAmpRaw","amp Raw", 2000,0,2000, NSTRIPS,0,NSTRIPS);

  TH2F* hbaseline = new TH2F("hbaseline","baseline; ADC ch; strip", 2000,10000,12000, NSTRIPS,0,NSTRIPS);
  TH2F* hnoise = new TH2F("hnoise","noise; ADC ch; strip", 1000,0,100, NSTRIPS,0,NSTRIPS);

  //input file
  TFile* inFile = new TFile(inputFileName.data(), "READ");
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
  // number of entries in the tree
  UShort_t k = 100, m = 10, l = k +m;
  UShort_t max = 0;
  ULong64_t nEntries = inTree->GetEntries();
  cout<<"nEntries "<<nEntries<<endl;
  for(ULong64_t i = 0; i < nEntries;i++){
    inTree->GetEntry(i);
    printProgress (i, inTree );
    boardIdx = getboardIndex(fboardID);
    max = ftrace[400];

    hsignal[boardIdx][fchannelID]->Reset();

    //GetSignalInfo(ftrace);
    counter[boardIdx][fchannelID]++;
    gr_baseline[boardIdx][fchannelID]->SetPoint( counter[boardIdx][fchannelID], counter[boardIdx][fchannelID], info.Baseline);
    gr_noise[boardIdx][fchannelID]->SetPoint( counter[boardIdx][fchannelID], counter[boardIdx][fchannelID], info.Noise);
        
    
    for(UShort_t i = 0; i < 991;i++){
      // ftrace[i] = ftrace[i+1];
      hsignal[boardIdx][fchannelID]->Fill(i+1, ftrace[i+1]);
      htraceSum[boardIdx][fchannelID]->Fill(i+1, ftrace[i+1]);
      /*if(i > 400 && i < 700){
        if(max < ftrace[i])
        max = ftrace[i];
        }*/
    }
   
    //calculate the baseline
    /* double baseline = 0.;
       for(UShort_t n = 0; n < 100; n++){
       baseline += ftrace[n];
       }
       baseline  = baseline/100.;
    */
    // l = k + m;

    //double ampRaw = max - baseline;
    double amplitude = myFilter->trapezoidal_filter_algorithm1(k, m, ftrace, info.Baseline);
    UShort_t sID = (boardIdx *4) + fchannelID;
    hAmp->Fill(amplitude, sID);
    hAmpRaw->Fill(info.Amplitude, sID);
    hbaseline->Fill(info.Baseline, sID);
    hnoise->Fill(info.Noise, sID);

   
    
  }

  cout<<"end of event loop "<<endl;
  inFile->cd();
  inFile->Close();
  cout<<"input file closed"<<endl;


  cout<<"Saving plots......"<<endl;
  fout->cd();
  for(int i = 0; i < NBOARDS; i++){

    for(int j = 0; j < NCHANNELS; j++){
      gr_baseline[i][j]->Write();
      gr_noise[i][j]->Write();
    }
  }
  fout->Write();
    
  fout->Close();
  cout<<"output file closed."<<endl;
    
  cout<<"clear memory.........."<<endl;
     
 
  for(int i = 0; i < NBOARDS; i++){
    delete [] hsignal[i];
    delete [] htraceSum[i];
    delete [] gr_baseline[i];
    delete [] gr_noise[i];
    delete [] counter[i];
    // delete [] hAmp[i];
  }
  delete [] hsignal;
  delete [] htraceSum;
  delete [] gr_baseline;
  delete [] gr_noise;
  delete [] counter;
  //delete [] hAmp;
  delete inFile;
  delete fout;
}

void signalTreatment::perform(){
  info.GetSignalInfo(ftrace);
double amplitude = myFilter->trapezoidal_filter_algorithm1(k, m, ftrace, info.Baseline);

//get CFD time

}
