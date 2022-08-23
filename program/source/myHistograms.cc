#include "myHistograms.hh"

myHistograms:: myHistograms()
{
}

myHistograms::~ myHistograms(){}

void myHistograms::rawHistInitialization(TFile* f){
  TString name;
  TString fname;
  int sizehisto = 65536;
  int maxhisto = 65536;
  int sizehisto2D = 992;
  int maxhisto2D = 16384;
  int fReduce_factor = 32;
  //dFile = f->mkdir("raw");
  dFile =gDirectory;// ->cd();// make it current directory
  
  Gain      = new TH1I**[NBOARDS];
  FeedBack  = new TH1I**[NBOARDS];
  Trace     = new TH1I**[NBOARDS];
  Trace_Sum = new TH2I**[NBOARDS];
  hRaw    = new TH1D**[NBOARDS];
  //hTrap = new TH1F**[NBOARDS];

  for (UShort_t iboard = 0;iboard <NBOARDS;iboard++){
    fname.Form("Card_%d",NUMEXO2[iboard]);
    dBoard[iboard] =  f->mkdir(fname.Data());
    dBoard[iboard]->cd();
    Gain[iboard]     = new TH1I*[NCHANNELS];
    FeedBack[iboard] = new TH1I*[NCHANNELS];
    Trace[iboard]    = new TH1I*[NCHANNELS];
    Trace_Sum[iboard]= new TH2I*[NCHANNELS];
    hRaw[iboard]    = new TH1D*[NCHANNELS];
    //hTrap[iboard]   = new TH1F*[NUMEXO_NB_CHANNELS];

    for ( UShort_t channel =0;channel <NCHANNELS;channel++){

      name.Form("Gain_%d_%d",NUMEXO2[iboard],channel);
      Gain[iboard][channel] = new TH1I (name.Data(),name.Data(),sizehisto,0,maxhisto);
      //TH1F *h = (TH1F*)gDirectory->Get(name);
      // h->SetDirectory(dBoard[iboard]);
      name.Form("FeedBack_%d_%d",NUMEXO2[iboard],channel);
      FeedBack[iboard][channel] = new TH1I (name.Data(),name.Data(),sizehisto,0,maxhisto);
      name.Form("Trace_%d_%d",NUMEXO2[iboard],channel);
      Trace[iboard][channel] = new TH1I (name.Data(),name.Data(),sizehisto2D,0,sizehisto2D);
      name.Form("Trace_Sum%d_%d",NUMEXO2[iboard],channel);
      Trace_Sum[iboard][channel] = new TH2I (name.Data(),name.Data(),sizehisto2D,0,sizehisto2D,maxhisto2D/fReduce_factor,0,maxhisto2D);
      name.Form("RawData_%d_%d",NUMEXO2[iboard],channel);
      hRaw[iboard][channel] = new TH1D(name.Data(), name.Data(),10000,0.,10000.);
      /* name.Form("Trap_%d_%d",fBoardList[iboard],channel);
         hTrap[iboard][channel] = new TH1F(name.Data(),name.Data(), SIRIUS_ITEMS_NUMBER, 0, SIRIUS_ITEMS_NUMBER);*/
    }
    dFile->cd();//change current directory to the top directory
  }
}

void myHistograms::deleteRawHistograms(){

  for (UShort_t iboard = 0;iboard <NBOARDS;iboard++){
    delete[] Gain[iboard];
    delete[] FeedBack[iboard];
    delete[] Trace[iboard];
    delete[] Trace_Sum[iboard];
    delete[] hRaw[iboard];
    //    delete[] hTrap[iboard];
  }

  delete[] Gain;
  delete[] FeedBack;
  delete[] Trace;
  delete[] Trace_Sum;
  delete[] hRaw;
  //    delete[] hTrap;
}
