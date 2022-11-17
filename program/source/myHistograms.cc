#include "myHistograms.hh"

myHistograms:: myHistograms()
{
	s1 = myGlobal::getInstance();
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
	dFile =gDirectory;// ->cd();// make it current directory
	dBoard = new TDirectory*[s1->NBOARDS_DSSD];
	Gain      = new TH1I**[s1->NBOARDS_DSSD];
	FeedBack  = new TH1I**[s1->NBOARDS_DSSD];
	Trace     = new TH1I**[s1->NBOARDS_DSSD];
	Trace_Sum = new TH2I**[s1->NBOARDS_DSSD];
	hRaw    = new TH1D**[s1->NBOARDS_DSSD];
	//hTrap = new TH1F**[s1->NBOARDS_DSSD];

	for (UShort_t iboard = 0;iboard <s1->NBOARDS_DSSD;iboard++){
		fname.Form("Card_%d",s1->boardList_DSSD[iboard]);
		dBoard[iboard] =  f->mkdir(fname.Data());
		dBoard[iboard]->cd();
		Gain[iboard]     = new TH1I*[NCHANNELS];
		FeedBack[iboard] = new TH1I*[NCHANNELS];
		Trace[iboard]    = new TH1I*[NCHANNELS];
		Trace_Sum[iboard]= new TH2I*[NCHANNELS];
		hRaw[iboard]    = new TH1D*[NCHANNELS];
		//hTrap[iboard]   = new TH1F*[NUMEXO_NB_CHANNELS];

		for ( UShort_t channel =0;channel <NCHANNELS;channel++){

			name.Form("Gain_%d_%d",s1->boardList_DSSD[iboard],channel);
			Gain[iboard][channel] = new TH1I (name.Data(),name.Data(),sizehisto,0,maxhisto);
			//TH1F *h = (TH1F*)gDirectory->Get(name);
			//h->SetDirectory(dBoard[iboard]);
			name.Form("FeedBack_%d_%d",s1->boardList_DSSD[iboard],channel);
			FeedBack[iboard][channel] = new TH1I (name.Data(),name.Data(),sizehisto,0,maxhisto);
			name.Form("Trace_%d_%d",s1->boardList_DSSD[iboard],channel);
			Trace[iboard][channel] = new TH1I (name.Data(),name.Data(),sizehisto2D,0,sizehisto2D);
			name.Form("Trace_Sum%d_%d",s1->boardList_DSSD[iboard],channel);
			Trace_Sum[iboard][channel] = new TH2I (name.Data(),name.Data(),sizehisto2D,0,sizehisto2D,maxhisto2D/fReduce_factor,0,maxhisto2D);
			name.Form("RawData_%d_%d",s1->boardList_DSSD[iboard],channel);
			hRaw[iboard][channel] = new TH1D(name.Data(), name.Data(),10000,0.,10000.);
			/* name.Form("Trap_%d_%d",fBoardList[iboard],channel);
			   hTrap[iboard][channel] = new TH1F(name.Data(),name.Data(), SIRIUS_ITEMS_NUMBER, 0, SIRIUS_ITEMS_NUMBER);*/
		}
		dFile->cd();//change current directory to the top directory
	}
}

void myHistograms::deleteRawHistograms(){
	for (UShort_t iboard = 0;iboard <s1->NBOARDS_DSSD;iboard++){
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
	delete [] dBoard;
}
