#include "digitalFilterOptimization.hh"

bool compareFWHM(trap_parameter i,trap_parameter j) { return (i.FWHM < j.FWHM); }


bool compareFWHM2(resolution_strip i,resolution_strip j) { return (i.FWHM < j.FWHM); }


digitalFilterOptimization::digitalFilterOptimization()
{
  filter = new digitalFilters<UShort_t>();
  ktimes = (kStop - kStart)/kStep + 1;
  mtimes = (mStop - mStart)/mStep + 1;
 
  nPeaks = 3;
  fPositionX = new Double_t[nPeaks];
  fPositionY = new Double_t[nPeaks];
  //3 alpha source 5156.6----241Am:  5485.8----244Cm: 5804.8
  Energy[0]= 5156.6;
  Energy[1]= 5804.8;

}

digitalFilterOptimization::~digitalFilterOptimization()
{
  delete filter;
  delete[] fPositionX;
  delete[] fPositionY;
}

void digitalFilterOptimization::setKnM_multiplicity(UShort_t k, UShort_t m){
  ktimes = k;
  mtimes = m;
}

void digitalFilterOptimization::initialize_TObjects()
{
 
  //calibration parameters
  gain = new Double_t***[NBOARDS];
  offset = new Double_t***[NBOARDS];
  for(UShort_t i = 0; i < NBOARDS;i++){
    gain[i] = new Double_t**[NCHANNELS];
    offset[i] = new Double_t**[NCHANNELS];
    for(UShort_t j = 0; j < NCHANNELS;j++){
      gain[i][j] = new Double_t*[ktimes];
      offset[i][j] = new Double_t*[ktimes];
      for(UShort_t k1 =0; k1 < ktimes; k1++){
        gain[i][j][k1] = new Double_t[mtimes];
        offset[i][j][k1] = new Double_t[mtimes];
        for(UShort_t m1 = 0; m1 < mtimes; m1++ ){
          gain[i][j][k1][m1] = 0.;
          offset[i][j][k1][m1] = 0.;
        }
      }
    }
  }
  
  //Fit functions
  fLinear = new TF1("fLinear","pol1",0,10000);
  G0    = new TF1("G0","gaus",0,10000);
  G1    = new TF1("G1","gaus",0,10000);
  G2    = new TF1("G2","gaus",0,10000);
  // graph for calibration
  gr_cal = new TGraph();
  //3 alpha source 5156.6----241Am:  5485.8----244Cm: 5804.8
  const Int_t nbins = 2500;
  Double_t xmin     = 0;
  Double_t xmax     = 2500.;
 
  spectrum = new TSpectrum(nPeaks);
  canvas = new TCanvas;
  canvas2 = new TCanvas;
  TList* list_res = new TList;
  
  hRaw = new TH1F****[NBOARDS];
  hCalib = new TH1F****[NBOARDS];
  hResolution = new TH2F**[NBOARDS];
 
  for(UShort_t i = 0; i < NBOARDS;i++){
    UShort_t board =  NUMEXO2[i];
    hRaw[i] = new TH1F***[NCHANNELS];
    hCalib[i] = new TH1F***[NCHANNELS];
    hResolution[i] = new TH2F*[NCHANNELS];
    for(UShort_t j = 0; j < NCHANNELS;j++){
      hResolution[i][j] = new TH2F(Form("hRes_b%d_ch%d",board, j),Form("hRes_b%d_ch%d ; k; m",board, j), ktimes, kStart, kStop, mtimes, mStart, mStop);
      list_res->Add(hResolution[i][j]);
      hRaw[i][j]= new TH1F**[ktimes];
      hCalib[i][j]= new TH1F**[ktimes];
      for(UShort_t k1 =0; k1 < ktimes; k1++){
        hRaw[i][j][k1]= new TH1F*[mtimes];
        hCalib[i][j][k1]= new TH1F*[mtimes];
        for(UShort_t m1 = 0; m1 < mtimes; m1++ ){
          //convert to  m,k values
          UShort_t k_val = kStart + kStep * k1;
          UShort_t m_val = mStart + mStep * m1;
          // cout<<"i "<<i<<"  j "<<j<<"  k1 "<<k1<<" m1 "<<m1<<endl;
          //  cout<<"k val "<<k_val<<"  m val "<<m_val<<endl;
          hRaw[i][j][k1][m1]= new TH1F(Form("hRaw_b%d_c%d_k%d_m%d", board, j, k_val, m_val),Form("hRaw_boardId%d_channel%d_k%d_m%d", board,j,k_val,m_val),nbins,xmin,xmax);
          hCalib[i][j][k1][m1]= new TH1F(Form("hCalib_b%d_c%d_k%d_m%d", board,j, k_val, m_val),Form("hCalib_boardId%d_channel%d_k%d_m%d", board,j,k_val,m_val),1500,0,6000);
        }
      }
    }
  } 
  //calibrated spectra
  list_calib = new TList;
  hCalibStrip = new TH2F**[ktimes];
  for(UShort_t k1 =0; k1 < ktimes; k1++){
    hCalibStrip[k1]= new TH2F*[mtimes];
    for(UShort_t m1 = 0; m1 < mtimes; m1++ ){
      //convert to  m,k values
      UShort_t k_val = kStart + kStep * k1;
      UShort_t m_val = mStart + mStep * m1;
      hCalibStrip[k1][m1] = new TH2F(Form("hCalib_k%d_m%d", k_val, m_val), Form("hCalib_k%d_m%d", k_val, m_val),1500,0,6000,NSTRIPS,0,NSTRIPS);
      list_calib->Add( hCalibStrip[k1][m1]);
    }
  }

  if(ktimes > 1 && mtimes > 1){
    h_m = new TH2F("h_m", "calib par: m;k;m",ktimes, kStart, kStop, mtimes, mStart, mStop);
    h_c = new TH2F("h_c", "calib par: c;k;m",ktimes, kStart, kStop, mtimes, mStart, mStop);
  }
}

void digitalFilterOptimization::delete_TObjects()
{
  for(UShort_t i = 0; i < NBOARDS;i++){
    for(UShort_t j = 0; j < NCHANNELS;j++){
      delete [] gain[i][j];
      delete [] offset[i][j];
        
    }
  }
  for(UShort_t i = 0; i < NBOARDS;i++){
    delete [] gain[i];
    delete [] offset[i];
  }

  delete [] gain;
  delete [] offset;
      
  delete gr_cal;
  delete spectrum;
  delete canvas;
  delete canvas2;
  delete G0; delete G1; delete G2; delete fLinear;
  list_res->Clear();
  list_calib->Clear();
  for(UShort_t i = 0; i < NBOARDS;i++){
    delete [] hResolution[i];
    for(UShort_t j = 0; j < NCHANNELS;j++){
      for(UShort_t k1 =0; k1 < ktimes; k1++){
        delete [] hRaw[i][j][k1];
        delete [] hCalib[i][j][k1];
      }
    }
  }
  delete [] hResolution;
 
  for(UShort_t i = 0; i < NBOARDS;i++){
    for(UShort_t j = 0; j < NCHANNELS;j++){
      delete [] hRaw[i][j];
      delete [] hCalib[i][j];        
    }
  }
  for(UShort_t i = 0; i < NBOARDS;i++){
    delete [] hRaw[i];
    delete [] hCalib[i];
  }

  delete [] hRaw;
  delete [] hCalib;
  for(UShort_t k1 =0; k1 < ktimes; k1++){
    delete[] hCalibStrip[k1];

  }
  delete[] hCalibStrip;
  delete list_res;
  delete list_calib;
  if(ktimes > 1 && mtimes > 1){
    delete h_m;
    delete h_c;
  }
}
void digitalFilterOptimization::optimize_trapezoidal_filter_parameters(int TrapAlgorithm, vector<string>inFileList)
{
  //Extract run numbers
  std::size_t pos = inFileList[0].find_first_of("_");
  string runNo = inFileList[0];
  runNo.erase(0,pos+1);
  pos = runNo.length();
  runNo.erase(4,pos);

  string runNo2 = inFileList[0];
  pos = runNo2.find_first_of("-");
  runNo2.erase(0,pos-1);
  pos = runNo2.length();
  runNo2.erase(1,pos);

  // open a text file to save the minimization parameters
  TString minParFile;
  minParFile.Form("../results/min_trapezoidal_parameters_run_%s_r%s_algorithm-%d.txt",  runNo.data(), runNo2.data(), TrapAlgorithm);
  ofstream myfile (minParFile.Data());
  if (myfile.is_open())
  {
    myfile << "boardID"<<setw(20)<<"channelID"<<setw(20)<<"k"<<setw(20)<<"m\n";

    //create an output file to save the spectra
    TString oFileName;
    oFileName.Form("../results/Spectra_run_%s_r%s_kStart%dkStop%dkStep%d_mStart%dmStop%dmStep%d-TrapAlg%d.root", runNo.data(), runNo2.data(), kStart,kStop, kStep, mStart,mStop, mStep, TrapAlgorithm);
    
    file_spectra= new TFile(oFileName,"RECREATE");
    //create histograms, functions and lists
    initialize_TObjects();
    
    vector< trap_parameter> myvector[NBOARDS][NCHANNELS];
    //check if trapezoidal amplitude file is already created
    bool file_exist = false;     
    char FileName[200];
    sprintf(FileName,"../results/run_%s_r%s_kStart%dkStop%dkStep%d_mStart%dmStop%dmStep%d-TrapAlg%d.root", runNo.data(), runNo2.data(), kStart,kStop, kStep, mStart,mStop, mStep, TrapAlgorithm);
    vector<string> list(1, FileName);
    file_exist = check_if_files_exist(list)[0];
    if(file_exist)
    {
      file_tree = new TFile(FileName,"READ");
      if(!file_tree->IsOpen()){
        // if file is not open clear memory and exit
        file_spectra->cd();
        file_spectra->Delete("T;*");
        file_spectra->Close();
        remove(oFileName.Data());
        delete_TObjects();         
        delete file_tree;
        delete file_spectra;
        std::cerr << "ERROR: could not open input file\n";
        std::terminate();
      }
      else{
      
        trapTree = (TTree*)file_tree->Get("trapData");
        trapTree->SetBranchAddress("boardId", &fboardID);
        trapTree->SetBranchAddress("boardIdx", &fboardIdx);
        trapTree->SetBranchAddress("channelId", &fchannelID);
        trapTree->SetBranchAddress("trapAmplitude", &trapAmplitude);
        trapTree->SetBranchAddress("k", &k);
        trapTree->SetBranchAddress("m", &m);
        trapTree->SetBranchAddress("kIdx", &kIdx);
        trapTree->SetBranchAddress("mIdx", &mIdx);
        cout<<"Filling raw histograms............."<<endl;
        for(ULong64_t i = 0; i < trapTree->GetEntries(); i++){
          trapTree->GetEntry(i);
          printProgress(i, trapTree);
          hRaw[fboardIdx][fchannelID][kIdx][mIdx]->Fill(trapAmplitude);
        }
      }
    }//file exist condition terminates here
    else{// trap amplitude file does not exist so read the input files
      //create a tree file to save the trapamplitude so that repeated filtering is avoided
      file_tree= new TFile(FileName,"RECREATE");
      trapTree = new TTree("trapData", "trapData");
      trapTree->Branch("boardId", &fboardID,"boardId/s");
      trapTree->Branch("boardIdx", &fboardIdx,"boardIdx/s");
      trapTree->Branch("channelId", &fchannelID,"channelId/s");
      trapTree->Branch("trapAmplitude", &trapAmplitude,"trapAmplitude/D");
      trapTree->Branch("k", &k,"k/s");
      trapTree->Branch("m", &m,"k/s");
      trapTree->Branch("kIdx", &kIdx,"kIdx/s");
      trapTree->Branch("mIdx", &mIdx,"mIdx/s");
      //check if data files exist
      vector<bool> fileStatus =  check_if_files_exist(inFileList);
      for(int f = 0; f <  inFileList.size(); f++){
        if(fileStatus[f]==false){
          file_spectra->cd();
          file_spectra->Delete("T;*");
          file_spectra->Close();
          remove(oFileName.Data());
          delete_TObjects();
          delete file_spectra;
          file_tree->cd();
          file_tree->Delete("T;*");
          file_tree->Close();
          delete trapTree;
          delete file_tree;

          
          fileStatus.clear();
          std::cerr << "ERROR: could not open input file\n";
          std::terminate();
        }
      }
      fileStatus.clear();
      TChain * inTree = new TChain("rawDataTree");
      for(int f =0; f < inFileList.size(); f++){
        cout<<inFileList[f]<<endl;
        inTree->Add(inFileList[f].data());
      }
      inTree->SetBranchAddress("Time", &ftime);
      inTree->SetBranchAddress("EventNo", &feventNo);
      inTree->SetBranchAddress("Trace", ftrace);
      inTree->SetBranchAddress("Gain", &fgain);
      inTree->SetBranchAddress("BoardID", &fboardID);
      inTree->SetBranchAddress("ChannelID", &fchannelID);
      inTree->SetBranchStatus("*",1);
      // number of entries in the tree
      Int_t nEntries = inTree->GetEntries();
      //Start event loops
      // nEntries = 1;
      cout<<"Reading file.............."<<endl;
      cout<<"number of entries in the file: "<<nEntries<<endl;

      for(ULong64_t entry = 0; entry < nEntries; entry++){
        inTree->GetEntry(entry);
        printProgress (entry,inTree);
        // cout<<"e "<<entry<<endl;
        //inTree->Show(entry);
        fboardIdx = getboardIndex(fboardID);
        //cout<<"here"<<endl;
        //calculate the baseline
        double baseline = 0.;
        for(UShort_t n = 0; n < 50; n++){
          baseline += ftrace[n];
        }
        baseline  = baseline/50.;
        //cout<<"baseline "<<baseline<<endl;
        for(k = kStart; k <= kStop; k += kStep){
          for(m= mStart; m <= mStop; m += mStep){
            kIdx = (k - kStart)/kStep;
            mIdx = (m - mStart)/mStep;
            // l = k + m;
            //   if((l + k) > 500) continue;
            //      cout<<"=======> k " <<k << "  l "<<l<<endl;
            if(TrapAlgorithm == 1)
              trapAmplitude = filter->trapezoidal_filter_algorithm1(k, m, ftrace, baseline);
            else if(TrapAlgorithm == 2) trapAmplitude = filter->trapezoidal_filter_algorithm2(k, m, ftrace, baseline);
            else  if(TrapAlgorithm == 3) trapAmplitude = filter->trapezoidal_filter_algorithm3(k, m, ftrace, baseline);
            else  if(TrapAlgorithm == 4) trapAmplitude = filter->trapezoidal_filter_algorithm4(k, m, ftrace, baseline);
            // trapAmplitude += rand->Uniform(-0.5,0.5);
            //           cout<<"k "<<k<<" m "<<m<<" board "<<boardIdx<<" channel "<<fchannelID<<endl;
            hRaw[fboardIdx][fchannelID][kIdx][mIdx]->Fill(trapAmplitude);
            trapTree->Fill();
          }//m loop ends here
        }//k loop ends here
   
      }
      //close input files
      delete inTree;
    }//read raw data files codition ends here
    //----------------------
    // Calibrate
    //-----------------
    for (UShort_t bID = 0; bID < NBOARDS; bID++) {
      for(UShort_t chID = 0; chID < NCHANNELS; chID++) {
        for(UShort_t kId = 0; kId < ktimes; kId++){
          for(UShort_t mId = 0; mId < mtimes; mId++){
            if(hRaw[bID][chID][kId][mId]->GetEntries() < 100) continue;
            hRaw[bID][chID][kId][mId]->GetXaxis()->SetRangeUser(500.,2500.);
            canvas->cd();
            hRaw[bID][chID][kId][mId]->Draw("");
            //Find the peaks
            Int_t nfound = spectrum->Search(hRaw[bID][chID][kId][mId],  2.5, "Markov", 0.2);
            Float_t *xpeaks = spectrum->GetPositionX();
            int idx[nfound];
            TMath::Sort(nfound,spectrum->GetPositionX(),idx,false);
            for (int i = 0; i < nfound; i++) {
              fPositionX[i] = xpeaks[idx[i]];
              Int_t bin = hRaw[bID][chID][kId][mId]->GetXaxis()->FindBin(fPositionX[i]);
              fPositionY[i] = hRaw[bID][chID][kId][mId]->GetBinContent(bin);
              //cout<<"x  "<<fPositionX[i]<<"  y "<<fPositionY[i]<<endl; 
            }

            //fit 1 peak first
            G0->SetRange(fPositionX[0]-10,fPositionX[0]+10);
            G0->SetParameters(fPositionY[0],fPositionX[0],3.);
            hRaw[bID][chID][kId][mId]->Fit(G0,"QLR");
            double sigma = G0->GetParameter(2);
            //Set Range for the fit functions
            G1->SetRange(fPositionX[0]-1*sigma,fPositionX[0]+2*sigma);
            G2->SetRange(fPositionX[2]-1*sigma,fPositionX[2]+2*sigma);
            //set parameters

            //G1->SetParameters(fPositionY[0],fPositionX[0],2.5,1,3, fPositionY[1],fPositionX[1],2.5,1,3);
            //G2->SetParameters(fPositionY[4],fPositionX[4],2.5,1,3, fPositionY[5],fPositionX[5],2.5,1,3);
            G1->SetParameters(fPositionY[0],fPositionX[0],sigma);
            G2->SetParameters(fPositionY[2],fPositionX[2],sigma);
            G1->SetParLimits(2, 1., 20.);
            G2->SetParLimits(2, 1., 20.);
            hRaw[bID][chID][kId][mId]->Fit(G1,"QLR+");
            hRaw[bID][chID][kId][mId]->Fit(G2,"QLR+");

            // Calibration
            fPositionX[0] = G1->GetParameter(1);
            fPositionX[1] = G2->GetParameter(1);
            for (int pt(0); pt<2; pt++) {
              gr_cal->SetPoint(pt,fPositionX[pt], Energy[pt]);
            }
            canvas2->cd();
            gr_cal->Draw();
            gr_cal->Fit(fLinear,"MQ","same");
	   
            //Get Parameters
            gain[bID][chID][kId][mId] =  fLinear->GetParameter(1);
            offset[bID][chID][kId][mId] = fLinear->GetParameter(0);

          }//m loop ends here
        }//k loop ends here
      }//channel loop ends here
    }// board loop ends here

    // fill calibrated spectra
    for(Int_t i = 0; i < trapTree->GetEntries();i++){
      trapTree->GetEntry(i);
      Double_t calibData = gain[fboardIdx][fchannelID][kIdx][mIdx]*trapAmplitude + offset[fboardIdx][fchannelID][kIdx][mIdx];
      hCalib[fboardIdx][fchannelID][kIdx][mIdx]->Fill(calibData);
      // UShort_t sID = (boardIdx *4) + fchannelID;
      // hCalib_tot->Fill(calibData, sID);
    }
 
    //optimize k and m
    for(UShort_t bID = 0; bID < NBOARDS; bID++){
      for(UShort_t chID = 0; chID < NCHANNELS; chID++){
        for(UShort_t kId = 0; kId < ktimes; kId++){
          for(UShort_t mId = 0; mId < mtimes; mId++){
            if( hCalib[bID][chID][kId][mId]->GetEntries() < 100) continue;
            hCalib[bID][chID][kId][mId]->GetXaxis()->SetRangeUser(5700.,6000.);
            canvas->cd();
            hCalib[bID][chID][kId][mId]->Draw();
            double height= hCalib[bID][chID][kId][mId]->GetMaximum();
            G0->SetRange(5800-30.,5800+50.);
            G0->SetParameters(height, 5805,10.);
            hCalib[bID][chID][kId][mId]->Fit(G0,"QLR");
            double sigma = G0->GetParameter(2);
          
            G1->SetRange(5805-1*sigma,5805+3*sigma);
            //int bin =  hCalib[bID][chID][kIdx][mIdx]->GetXaxis()->FindBin(5770);
            //double height1= hCalib[bID][chID][kIdx][mIdx]->GetBinContent(bin);
          
            //G1->SetParameters(height1, 5770, 10,1, 3, height, 5805,10.,1.3);
            G1->SetParameters(height, 5805,sigma);
            hCalib[bID][chID][kId][mId]->Fit(G1,"QLR+");
            k = kId * kStep + kStart;
            m = mId * mStep + mStart;
            double FWHM = G1->GetParameter(2) * 2.35;
            trap_parameter r (k, m, FWHM);
            myvector[bID][chID].push_back(r);
            hResolution[bID][chID]->Fill(k, m,FWHM );
          }
        }   
      }
    }
    //----------------
    // Sort in decrasing order----------------
    //------------
    for(UShort_t bID = 0; bID < NBOARDS; bID++){
      for(UShort_t chID = 0; chID < NCHANNELS; chID++){
        if(myvector[bID][chID].size() != 0){
       
          std::sort (myvector[bID][chID].begin(), myvector[bID][chID].end(),compareFWHM);
          cout<<"bID "<< NUMEXO2[bID]<<"  chId "<<chID<< " FWHM "<<myvector[bID][chID][0].FWHM<<"  k "<<myvector[bID][chID][0].k <<" m "<<myvector[bID][chID][0].m<<endl;
       
        }
      }
    }
    // sort ends here
    cout<<"Saving plots......"<<endl;
    file_tree->cd();
    trapTree->Write();
    file_tree->Close();


    file_spectra->cd();
    list_res->Write();
    list_calib->Write();
    file_spectra->Close();
    cout<<"output file closed."<<endl;
    //clear memory
    cout<<"clearing memory..................";

    for(UShort_t bID = 0; bID < NBOARDS; bID++){
      for(UShort_t chID = 0; chID < NCHANNELS; chID++){
        myvector[bID][chID].clear();
      }
    }
    delete_TObjects();
    // delete inFile;
    delete file_tree;
    delete file_spectra;
    
    cout<<"done\n"<<endl;
    //close 
    myfile.close();
  }//minTrapparameter file condition ends here
    
}

void digitalFilterOptimization::sort_data_with_min_trapParametes(int TrapAlgorithm, std::string  input_file_name){
  setKnM_multiplicity(1,1);
  //Extract run numbers
  std::size_t pos = input_file_name.find_first_of("_");
  string runNo = input_file_name;
  runNo.erase(0,pos+1);
  pos = runNo.length();
  runNo.erase(4,pos);

  string runNo2 = input_file_name;
  pos = runNo2.find_first_of("-");
  runNo2.erase(0,pos-1);
  pos = runNo2.length();
  runNo2.erase(1,pos);
  //get min trapezoidal parameters
  UShort_t **min_k = new UShort_t*[NBOARDS];
  UShort_t **min_m = new UShort_t*[NBOARDS];

  for(UShort_t i = 0; i < NBOARDS; i++){
    min_k[i] = new UShort_t[NCHANNELS];
    min_m[i] = new UShort_t[NCHANNELS];
    for(UShort_t j = 0; j < NCHANNELS; j++){
      min_k[i][j] = 0;
      min_m[i][j] = 0;
    }
  }
  TString minFile;
  minFile.Form("../results/min_trapezoidal_parameters_run_%s_r%s_algorithm-%d.txt",runNo.data(), runNo2.data(), TrapAlgorithm );
  string line; UShort_t boardID, bIdx, channelID, k_val, m_val;
  ifstream myfile (minFile.Data());
  if (myfile.is_open())
  {
    getline (myfile,line);
    while (!myfile.eof() )
    {
      myfile >> boardID >> channelID >> k_val >> m_val;
      bIdx = getboardIndex(boardID);
      min_k[bIdx][channelID] = k_val;
      min_m[bIdx][channelID] = m_val;
    }
    myfile.close();
 
    TString minParFile;
    minParFile.Form("../results/Calibration_Parameters_with_min_trapezoidal_parameters_run_%s_r%s_algorithm-%d.txt",  runNo.data(), runNo2.data(), TrapAlgorithm);
    ofstream myfile (minParFile.Data());
    if (myfile.is_open())
    {
      myfile << "boardID"<<setw(20)<<"channelID"<<setw(20)<<"gain"<<setw(20)<<"offset\n";

      Double_t amplitude;
      //output file
      TString oFileName;
      oFileName.Form("../results/Spectra_run_%s_r%s_minTrapParameters-TrapAlg%d.root", runNo.data(), runNo2.data(), TrapAlgorithm);
    
      TFile* fout= new TFile(oFileName,"RECREATE");
      //Fit functions
      initialize_TObjects();
      TF1 *Cball2    = new TF1("Cball2", fit_function,0,10000,10);
      TH1F* htot = new TH1F("htot","hCalib_tot; E;",1500,0,6000);
 
      //check if the file exists in the directory
      string dir = string("../results");
      vector<string> files = vector<string>();

      DIR *dp;
      struct dirent *dirp;

      if((dp = opendir(dir.c_str())) == NULL)
      {
        cout << "Error opening " << dir << endl;
      }
    
      while ((dirp = readdir(dp)) != NULL)
        files.push_back(string(dirp->d_name));
      closedir(dp);

     
      char inFileName[200];
      sprintf(inFileName,"run_%s_r%s_kStart%dkStop%dkStep%d_mStart%dmStop%dmStep%d-TrapAlg%d.root", runNo.data(), runNo2.data(), kStart,kStop, kStep, mStart,mStop, mStep, TrapAlgorithm);
      //TString inFileName = "../results/run_0038_r0_kStart40kStop250kStep10_mStart0mStop100mStep10-TrapAlg1.root";
      for (unsigned int i = 0;i < files.size();i++){
        if(files[i].compare(inFileName) ==0)
          cout << files[i]<<"  exists in"<<dir.c_str()<< endl;
      }
      files.clear();
        
      char dirpath[100];
  
      // printing current working directory
      printf("Current directory %s\n", getcwd(dirpath, 100));
  
      // using the command
      chdir(dir.c_str());
  
      // printing current working directory
      printf("New current directory %s\n", getcwd(dirpath, 100));
  
      TFile* inFile= new TFile(inFileName,"READ");
      if(!inFile->IsOpen()){
        fout->cd();
        fout->Delete("T;*");
        fout->Close();
        remove(oFileName.Data());
        delete_TObjects();

        for(UShort_t i = 0; i < NBOARDS; i++){
          delete [] min_k[i];
          delete [] min_m[i];
  
        }
        delete [] min_k;
        delete [] min_m;
        delete inFile;
        delete fout;
        std::cerr << "ERROR: could not open input file\n";
        std::terminate();
      }
      TTree* tree = (TTree*)inFile->Get("trapData");
      tree->SetBranchAddress("boardId", &boardID);
      //tree->SetBranchAddress("boardIdx", &bIdx);
      tree->SetBranchAddress("channelId", &channelID);
      tree->SetBranchAddress("trapAmplitude", &amplitude);
      tree->SetBranchAddress("k", &k_val);
      tree->SetBranchAddress("m", &m_val);
      //tree->SetBranchAddress("kIdx", &kIdx);
      // tree->SetBranchAddress("mIdx", &mIdx);
      cout<<"Filling raw histograms............."<<endl;
      for(ULong64_t i = 0; i < tree->GetEntries(); i++){
        tree->GetEntry(i);
        printProgress(i, tree);
        bIdx = getboardIndex(boardID);
        if(k_val == min_k[bIdx][channelID] && m_val == min_m[bIdx][channelID])
          hRaw[bIdx][channelID][0][0]->Fill(amplitude);
      }
    
      //Calibrate here
      for (UShort_t bIdx = 0; bIdx < NBOARDS; bIdx++) {
        for(UShort_t chID = 0; chID < NCHANNELS; chID++) {
   
          //cout<<"here 1"<<endl;
          if( hRaw[bIdx][chID][0][0]->GetEntries() < 100) continue;
          hRaw[bIdx][chID][0][0]->GetXaxis()->SetRangeUser(500.,2500.);
          canvas->cd();
          hRaw[bIdx][chID][0][0]->Draw("");
          //Find the peaks
          Int_t nfound = spectrum->Search(hRaw[bIdx][chID][0][0], 3, "Markov", 0.3); //hist, sigma, option, threshold
          Float_t *xpeaks = spectrum->GetPositionX();
          int idx[nfound];
          TMath::Sort(nfound,spectrum->GetPositionX(),idx,false);
          for (int i = 0; i < nfound; i++) {
            fPositionX[i] = xpeaks[idx[i]];
            Int_t bin = hRaw[bIdx][chID][0][0]->GetXaxis()->FindBin(fPositionX[i]);
            fPositionY[i] = hRaw[bIdx][chID][0][0]->GetBinContent(bin);
            // cout<<"x  "<<fPositionX[i]<<"  y "<<fPositionY[i]<<endl; 
          }
 
          //fit 1 peak first
          G0->SetRange(fPositionX[0]-10,fPositionX[0]+10);
          G0->SetParameters(fPositionY[0],fPositionX[0],3.);
          hRaw[bIdx][chID][0][0]->Fit(G0,"QLR");
          double sigma = G0->GetParameter(2);
          //Set Range for the fit functions
          G1->SetRange(fPositionX[0]-1*sigma,fPositionX[0]+2*sigma);
          G2->SetRange(fPositionX[2]-1*sigma,fPositionX[2]+2*sigma);
          //set parameters

          //G1->SetParameters(fPositionY[0],fPositionX[0],2.5,1,3, fPositionY[1],fPositionX[1],2.5,1,3);
          //G2->SetParameters(fPositionY[4],fPositionX[4],2.5,1,3, fPositionY[5],fPositionX[5],2.5,1,3);
          G1->SetParameters(fPositionY[0],fPositionX[0],sigma);
          G2->SetParameters(fPositionY[2],fPositionX[2],sigma);
          G1->SetParLimits(2, 1., 20.);
          G2->SetParLimits(2, 1., 20.);
          hRaw[bIdx][chID][0][0]->Fit(G1,"QLR+");
          hRaw[bIdx][chID][0][0]->Fit(G2,"QLR+");

          // Calibration
          fPositionX[0] = G1->GetParameter(1);
          fPositionX[1] = G2->GetParameter(1);
          for (int pt(0); pt<2; pt++) {
            gr_cal->SetPoint(pt,fPositionX[pt], Energy[pt]);
          }
          canvas2->cd();
          gr_cal->Draw();
          gr_cal->Fit(fLinear,"MQ","same");
	   
          //Get Parameters
          gain[bIdx][chID][0][0] =  fLinear->GetParameter(1);
          offset[bIdx][chID][0][0] = fLinear->GetParameter(0);
      
          myfile << NUMEXO2[bIdx]<<setw(20)<<chID<<setw(20)<< gain[bIdx][chID][0][0] <<setw(20)<<offset[bIdx][chID][0][0]<<endl;

        }
      }
      //Get Parameters
      cout<<"Calibration started..............."<<endl;

      for(ULong64_t i = 0; i < tree->GetEntries();i++){
        tree->GetEntry(i);
        printProgress(i, tree);
        bIdx = getboardIndex(boardID);
        if(k_val == min_k[bIdx][channelID] && m_val == min_m[bIdx][channelID]){
          Double_t calibData = gain[bIdx][channelID][0][0]*amplitude + offset[bIdx][channelID][0][0];
          hCalib[bIdx][channelID][0][0]->Fill(calibData);
          UShort_t sID = (bIdx *4) + channelID;
          hCalibStrip[0][0]->Fill(calibData, sID);
          htot->Fill(calibData);
        }
      }
      cout<<"Fit 5.8 MeV peak----------------"<<endl;
      vector<resolution_strip>minRes; int bin1, bin2; double height, height1, sigma, temp;
      //optimize k and m
      for(UShort_t bIdx = 0; bIdx < NBOARDS; bIdx++){
        for(UShort_t chID = 0; chID < NCHANNELS; chID++){
          if( hCalib[bIdx][chID][0][0]->GetEntries() < 100) continue;
          hCalib[bIdx][chID][0][0]->GetXaxis()->SetRangeUser(5700.,6000.);
          canvas->cd();
          hCalib[bIdx][chID][0][0]->Draw();
          height= hCalib[bIdx][chID][0][0]->GetMaximum();
          G0->SetRange(5800-30.,5800+50.);
          G0->SetParameters(height, 5805,10.);
          hCalib[bIdx][chID][0][0]->Fit(G0,"QLR");
          sigma = G0->GetParameter(2);
          
          G1->SetRange(5805-1*sigma,5805+3*sigma);
          G1->SetParameters(height, 5805,sigma);
         
          G1->SetParameters(height, 5805,sigma);
          hCalib[bIdx][chID][0][0]->Fit(G1,"QLR");
          //max in the range
          bin1 =  hCalib[bIdx][chID][0][0]->GetXaxis()->FindBin(5750);
          bin2 =  hCalib[bIdx][chID][0][0]->GetXaxis()->FindBin(5780);
          height1 = hCalib[bIdx][chID][0][0]->GetBinContent(bin1);
          for(int i = bin1; i <= bin2;i++){
            temp = hCalib[bIdx][chID][0][0]->GetBinContent(i);
            if(height1 < temp)
              height1 = temp; 
          }

          Cball2->SetRange(5740,5850);
          Cball2->SetParameters(height1, 5762, sigma, 10, 3, height, 5805,sigma, 10, 3);
          hCalib[bIdx][chID][0][0]->Fit(G0,"QLR");
          
          UShort_t board = NUMEXO2[bIdx];
          double FWHM = Cball2->GetParameter(7)*2.35;
          k_val = min_k[bIdx][chID];
          m_val = min_m[bIdx][chID];
          resolution_strip r(board, chID,  k_val, m_val , FWHM);
          minRes.push_back(r);
        }
      }
 
      //get the strip with min Resolution
      std::sort(minRes.begin(), minRes.end(), compareFWHM2);

      cout<<"minimum FWHM: "<<minRes[0].FWHM <<"  board: "<<minRes[0].board<<"  channel: "<<minRes[0].channel<<"  k: "<<minRes[0].k<<"  m: "<<minRes[0].m<<endl;
      minRes.clear();

      ///gettotal res

      bin1 =  htot->GetXaxis()->FindBin(5750);
      bin2 =  htot->GetXaxis()->FindBin(5780);
      height1 = htot->GetBinContent(bin1);
      for(int i = bin1; i <= bin2;i++){
        temp = htot->GetBinContent(i);
        if(height1 < temp)
          height1 = temp; 
      }
      htot->GetXaxis()->SetRangeUser(5700.,6000.);
      canvas->cd();
      htot->Draw();
      height= htot->GetMaximum();
          
      Cball2->SetRange(5740,5850);
      Cball2->SetParameters(height1, 5762, sigma, 10, 3, height, 5805,sigma, 10, 3);
      htot->Fit(Cball2,"QLR");
      sigma = Cball2->GetParameter(7);
      double FWHM = sigma*2.35;
      cout<<"FWHM total: "<<FWHM<<endl;

      //from the function
      Double_t maxf = Cball2->GetMaximum(5805-sigma, 5805+sigma);

      maxf = maxf/2.;

      Double_t x1 = Cball2->GetX(maxf, 5805-2*sigma, 5805);
      Double_t x2 = Cball2->GetX(maxf, 5805, 5805+2*sigma); 

      cout<<"new FWHM "<<x2-x1<<endl;

      cout<<"end of event loop "<<endl;
      inFile->cd();
      inFile->Close();
      cout<<"input file closed"<<endl;

      cout<<"Saving plots......"<<endl;
      fout->cd();
      fout->Write();
      fout->Close();
      cout<<"output file closed."<<endl;
      //clear memory
      cout<<"Clearing memory..................";
      delete_TObjects();

  
      for(UShort_t i = 0; i < NBOARDS; i++){
        delete [] min_k[i];
        delete [] min_m[i];
  
      }
      delete [] min_k;
      delete [] min_m;
      delete inFile;
      delete fout;
      myfile.close();
    }
    cout<<"done\n"<<endl;
    chdir("../build");
  }
  else cout << "Unable to open min trapezoidal parameter file for this run\n"; 
  // printing current working directory
  // printf("New current directory %s\n", getcwd(dirpath, 100));
  

}
