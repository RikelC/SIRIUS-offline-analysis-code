/*!
 * \file myHistograms.cc
 * \author Rikel Chakma
 * \date 18/11/2022
 * \brief A class for creating and resetting histograms.
 * \details The spectra are catagorically saved in specific directories of the output histogram TFile under the name specidied by the user. 
 *
 *
 */

#include "myHistograms.hh"
/*!
 * Constructor
 */
myHistograms:: myHistograms()
{
	s1 = myGlobal::getInstance();
}
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
/*!
 * Destructor
 */
myHistograms::~ myHistograms(){}
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
/*!
 *  Iniatialization of the raw histograms
 */
void myHistograms::rawHistInitialization(TFile* f){
	TString name;
	TString fname;
	TString sfname;
	int sizehisto = 65536;
	int maxhisto = 65536;
	int sizehisto2D = 992;
	int maxhisto2D = 16384;
	int fReduce_factor = 32;
	dFile =gDirectory;// ->cd();// make it current directory

	std::string dssd_boards_char[s1->NBOARDS_DSSD];
	for(int i = 0; i < s1->NBOARDS_DSSD;i++) dssd_boards_char[i] = std::to_string(s1->boardList_DSSD[i]);
	std::string tunnel_boards_char[s1->NBOARDS_TUNNEL];
	for(int i = 0; i < s1->NBOARDS_TUNNEL;i++) tunnel_boards_char[i] = std::to_string(s1->boardList_Tunnel[i]);

	dBoard = new TDirectory*[s1->NBOARDS_DSSD];
	dBoard2 = new TDirectory*[s1->NBOARDS_TUNNEL];

	//Gain      = new TH1I**[s1->NBOARDS_DSSD];
	//FeedBack  = new TH1I**[s1->NBOARDS_DSSD];

	hTrace_Sum     = new TH2F**[s1->NBOARDS_DSSD];
	hRaw           = new TH1F**[s1->NBOARDS_DSSD];
	hCalib         = new TH1F**[s1->NBOARDS_DSSD];
	hTrap          = new TH2F**[s1->NBOARDS_DSSD];
	gr_baseline    = new TGraph**[s1->NBOARDS_DSSD];
	hTrigger      = new TH1F**[s1->NBOARDS_DSSD];
	hBaseline      = new TH1F**[s1->NBOARDS_DSSD];
	hNoise         = new TH1F**[s1->NBOARDS_DSSD];
	hRisetime      = new TH1I**[s1->NBOARDS_DSSD];
	TDirectory *tempDir;
	for (UShort_t iboard = 0;iboard <s1->NBOARDS_DSSD;iboard++){
		fname.Form("Card_%d",s1->boardList_DSSD[iboard]);
		dBoard[iboard] = (TDirectory*) dFile->GetDirectory(fname.Data());
		if(!dBoard[iboard]) dBoard[iboard] = dFile->mkdir(fname.Data());
		dBoard[iboard]->cd();
		//Gain[iboard]     = new TH1I*[NCHANNELS];
		//FeedBack[iboard] = new TH1I*[NCHANNELS];
		hTrace_Sum[iboard]    = new TH2F*[NCHANNELS];
		hRaw[iboard]          = new TH1F*[NCHANNELS];
		hCalib[iboard]        = new TH1F*[NCHANNELS];
		hTrap[iboard]         = new TH2F*[NCHANNELS];
		gr_baseline[iboard]   = new TGraph*[NCHANNELS];
		hTrigger[iboard]     = new TH1F*[NCHANNELS];
		hBaseline[iboard]     = new TH1F*[NCHANNELS];
		hNoise[iboard]        = new TH1F*[NCHANNELS];
		hRisetime[iboard]     = new TH1I*[NCHANNELS];

		for ( UShort_t channel =0;channel <NCHANNELS;channel++){


			/*  name.Form("Gain_%d_%d",boardList[iboard],channel);
			    hGain[iboard][channel] = new TH1I (name.Data(),name.Data(),sizehisto,0,maxhisto);
			    GetSpectra()->AddSpectrum(hGain[iboard][channel],fname);

			    name.Form("FeedBack_%d_%d",boardList[iboard],channel);
			    hFeedBack[iboard][channel] = new TH1I (name.Data(),name.Data(),sizehisto,0,maxhisto);
			    GetSpectra()->AddSpectrum(hFeedBack[iboard][channel],fname);
			    */
			sfname = "Trace";
			tempDir = (TDirectory*) dBoard[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("Trace_Sum%d_%d",s1->boardList_DSSD[iboard],channel);
			hTrace_Sum[iboard][channel] = new TH2F (name.Data(),name.Data(),s1->TRACE_SIZE,0,s1->TRACE_SIZE,1700,-1000,16000);

			sfname = "RawHist";
			tempDir = (TDirectory*) dBoard[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("RawData_%d_%d",s1->boardList_DSSD[iboard],channel);
			hRaw[iboard][channel] = new TH1F(name.Data(), name.Data(),2000,0.,20000.);


			sfname = "CalibHist";
			tempDir = (TDirectory*) dBoard[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("CalibData_%d_%d",s1->boardList_DSSD[iboard],channel);
			hCalib[iboard][channel] = new TH1F(name.Data(), name.Data(),2000,0.,10000.);

			sfname = "Trapezoidal";
			tempDir = (TDirectory*) dBoard[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("Trap_%d_%d",s1->boardList_DSSD[iboard],channel);
			hTrap[iboard][channel] = new TH2F(name.Data(),name.Data(), s1->TRACE_SIZE, 0,s1-> TRACE_SIZE, 4000,-8000,8000);

			sfname = "Trigger";
			tempDir = (TDirectory*) dBoard[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("hTrigger_%d_%d",s1->boardList_DSSD[iboard],channel);
			hTrigger[iboard][channel] = new TH1F(name.Data(),name.Data(), 992,0, 992);

			sfname = "BaselineHist";
			tempDir = (TDirectory*) dBoard[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("hBaseline_%d_%d",s1->boardList_DSSD[iboard],channel);
			hBaseline[iboard][channel] = new TH1F(name.Data(),name.Data(), 4000,0,16000);

			sfname = "BaselineGraph";
			tempDir = (TDirectory*) dBoard[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("Baseline_%d_%d",s1->boardList_DSSD[iboard],channel);
			gr_baseline[iboard][channel] =  new TGraph();
			gr_baseline[iboard][channel]->SetNameTitle(name,name);

			sfname = "Noise";
			tempDir = (TDirectory*) dBoard[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("hNoise%d_%d",s1->boardList_DSSD[iboard],channel);
			hNoise[iboard][channel] = new TH1F(name.Data(),name.Data(), 1000,0,100);

			sfname = "Risetime";
			tempDir = (TDirectory*) dBoard[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("hRisetime%d_%d",s1->boardList_DSSD[iboard],channel);
			hRisetime[iboard][channel] = new TH1I(name.Data(),name.Data(), 1000,0,1000);

		}
		dFile->cd();//change current directory to the top directory
	}
	dFile->cd();
	fname = "DssdSpectra";	
	tempDir = (TDirectory*) dFile->GetDirectory(fname.Data());
	if(!tempDir) tempDir = dFile->mkdir(fname.Data());
	tempDir->cd();

	//- time difference b/w 2 front strips
	h_delT_ff =  new TH1I("h_delT_ff","#Delta T(front1-front2);ns", 1000,-50000,50000);
	//- time difference b/w  front and back strips
	h_delT_fb =  new TH1I("h_delT_fb","#Delta T(front - back);ns",1000,-50000,50000);
	//- time difference b/w back and front strips
	h_delT_bf =  new TH1I("h_delT_bf","#Delta T(back - front);ns",1000,-50000,50000);
	//- time difference b/w 2 back strips
	h_delT_bb =  new TH1I("h_delT_bb","#Delta T(back1 - back2);ns",1000,-50000,50000);
	//-raw energy vs strip number
	h_raw_strip = new TH2F("h_raw_strip",";E (ADC ch);strip number", 2000,0,10000,256,0,256);
	//-calibrated energy vs strip number
	h_calib_strip = new TH2F("h_calib_strip",";E (keV);strip number", 2000,0,10000,256,0,256);
	//-histo for front energy vs back energy
	h_E_frontBack = new TH2F("h_E_frontBack",";frontE;backE",1000,0,10000,1000,0,10000);
	//-hit pattern in the DSSD
	h_DSSD_XY_hit =  new TH2I("h_DSSD_XY_hit","hit pattern ;X;Y",128,0,128,128,0,128);
	//-histo for viewing number of counts in each strip
	h_dssd_count_strip = new TH1I("h_dssd_count_strip",";strip number; counts", s1->NSTRIPS_DSSD,0,s1->NSTRIPS_DSSD);
	//histo for viewing counts per board
	h_dssd_count_board = new TH1I("h_dssd_count_board",";board number; counts", s1->NBOARDS_DSSD,0,s1->NBOARDS_DSSD);
	for(int i = 1; i <= s1->NBOARDS_DSSD; i++)h_dssd_count_board->GetXaxis()->SetBinLabel(i, dssd_boards_char[i-1].c_str());
	//	h_dssd_count_board->GetXaxis()->SetLabelSize(0.1);

	//-event rates in sec in each strip of the DSSD
	//h_dssdStrip_rate = new TH1I("h_dssdStrip_rate",Form(";strip number;counts/%d sec",s1->rate_calcul_time_lapse), s1->NSTRIPS_DSSD,0,s1->NSTRIPS_DSSD);
	//-event rates in sec in each Board of the DSSD
	//h_dssdBoard_rate = new TH1I("h_dssdBoard_rate",Form(";board number;counts/%d sec",s1->rate_calcul_time_lapse), s1->NBOARDS_DSSD,0,s1->NBOARDS_DSSD);
	//for(int i = 1; i <= s1->NBOARDS_DSSD; i++)h_dssdBoard_rate->GetXaxis()->SetBinLabel(i, dssd_boards_char[i-1].c_str());
	//h_dssdBoard_rate->GetXaxis()->SetLabelSize(0.1);	
	//Rates per board graphs
	sfname = "RatePerBoard";
/*	TDirectory * tempDir2 = (TDirectory*) tempDir->GetDirectory(sfname.Data());
	if(!tempDir2) tempDir2 = tempDir->mkdir(sfname.Data());
	tempDir2->cd();
	gr_rate_dssdBoard = new TGraph*[s1->NBOARDS_DSSD];
	for(int iboard = 0;iboard <s1->NBOARDS_DSSD;iboard++){
		name.Form("rate_%d",s1->boardList_DSSD[iboard]);
		gr_rate_dssdBoard[iboard] =  new TGraph();
		gr_rate_dssdBoard[iboard]->SetName(name);
		gr_rate_dssdBoard[iboard]->SetTitle(Form("Rate %d; time(x %d sec);counts/%d sec", s1->boardList_DSSD[iboard], s1->rate_calcul_time_lapse,s1->rate_calcul_time_lapse));
	}	
	tempDir->cd();*/
	dFile->cd();
	//---------tunnel
	h_tunnelRaw = new TH1I**[s1->NBOARDS_TUNNEL];
	h_tunnelCalib = new TH1I**[s1->NBOARDS_TUNNEL];
	for(int iboard = 0;iboard < s1->NBOARDS_TUNNEL;iboard++){
		fname.Form("Card_%d",s1->boardList_Tunnel[iboard]);
		dBoard2[iboard] = (TDirectory*) dFile->GetDirectory(fname.Data());
		if(!dBoard2[iboard]) dBoard2[iboard] = dFile->mkdir(fname.Data());
		dBoard2[iboard]->cd();
		h_tunnelRaw[iboard]    = new TH1I*[NCHANNELS];
		h_tunnelCalib[iboard]    = new TH1I*[NCHANNELS];

		for (int channel =0;channel <NCHANNELS;channel++){
			//Raw data
			sfname = "RawHist";			
			tempDir = (TDirectory*) dBoard2[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard2[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("RawData_%d_%d",s1->boardList_Tunnel[iboard], channel);
			h_tunnelRaw[iboard][channel] = new TH1I (name.Data(),name.Data(),8000,0,20000);
			//Calib data
			sfname = "CalibHist";		
			tempDir = (TDirectory*) dBoard2[iboard]->GetDirectory(sfname);
			if(!tempDir) tempDir = dBoard2[iboard]->mkdir(sfname.Data());
			tempDir->cd();
			name.Form("CalibData_%d_%d",s1->boardList_Tunnel[iboard], channel);
			h_tunnelCalib[iboard][channel] = new TH1I (name.Data(),name.Data(),8000,0,8000);

		}
		dFile->cd();	
	}
	dFile->cd();
	fname = "TunnelSpectra";
	tempDir = (TDirectory*) dFile->GetDirectory(fname.Data());
	if(!tempDir) tempDir = dFile->mkdir(fname.Data());
	tempDir->cd();
	//-raw energy vs strip number
	h_raw_pad = new TH2F("h_raw_pad",";E (ADC ch);pad number", 8000,0,20000,96,0,96);
	//-calibrated energy vs strip number
	h_calib_pad = new TH2F("h_calib_pad",";E (keV);pad number", 2000,0,10000,96,0,96);
	//-hit pattern in the TUNNEL

	h_TUNNEL_XY_hit = new TH2I*[s1->NDETECTOR_TUNNEL];
	for(int i =0; i < s1->NDETECTOR_TUNNEL;i++){
		name.Form("h_tunnel%d_XY_hit",i+1);
		h_TUNNEL_XY_hit[i] =  new TH2I(name,Form("Tunnel %d;X;Y",i+1),8,0,8,8,0,8);
	}
	//-histo for viewing number of counts in each pad
	h_tunnel_count_pad = new TH1I("h_tunnel_count_pad",";pad number; counts", s1->NofMacroPixels,0,s1->NofMacroPixels);
	//histo for viewing counts per board
	h_tunnel_count_board = new TH1I("h_tunnel_count_board",";board number; counts", s1->NBOARDS_TUNNEL,0,s1->NBOARDS_TUNNEL);
	for(int i = 1; i <= s1->NBOARDS_TUNNEL; i++)h_tunnel_count_board->GetXaxis()->SetBinLabel(i, tunnel_boards_char[i-1].c_str());
	//	h_tunnel_count_board->GetXaxis()->SetLabelSize(0.1);


	//-event rates in sec in each strip of the DSSD
	//h_tunnelPad_rate = new TH1I("h_tunnelPad_rate",Form(";pad number;counts/%d sec",s1->rate_calcul_time_lapse), s1->NofMacroPixels,0,s1->NofMacroPixels);

	//-event rates in sec in each Board of the DSSD
	//h_tunnelBoard_rate = new TH1I("h_tunnelBoard_rate",Form(";board number;counts/%d sec",s1->rate_calcul_time_lapse), s1->NBOARDS_TUNNEL,0,s1->NBOARDS_TUNNEL);
	//for(int i = 1; i <= s1->NBOARDS_TUNNEL; i++)h_tunnelBoard_rate->GetXaxis()->SetBinLabel(i, tunnel_boards_char[i-1].c_str());
	//h_tunnelBoard_rate->GetXaxis()->SetLabelSize(0.1);	
	//Rate per board graphs
	sfname = "RatePerBoard";
/*	tempDir2 = (TDirectory*) tempDir->GetDirectory(sfname.Data());
	if(!tempDir2) tempDir2 = tempDir->mkdir(sfname.Data());
	tempDir2->cd();
	gr_rate_tunnelBoard = new TGraph*[s1->NBOARDS_TUNNEL];
	for(int iboard = 0;iboard <s1->NBOARDS_TUNNEL;iboard++){
		name.Form("rate_%d",s1->boardList_Tunnel[iboard]);
		gr_rate_tunnelBoard[iboard] =  new TGraph();
		gr_rate_tunnelBoard[iboard]->SetName(name);
		gr_rate_tunnelBoard[iboard]->SetTitle(Form("Rate %d; time(x %d sec);counts/%d sec", s1->boardList_Tunnel[iboard], s1->rate_calcul_time_lapse,s1->rate_calcul_time_lapse));
	}

*/
	dFile->cd();
}
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
/*!
 * Clearing of memory
 */
void myHistograms::deleteRawHistograms(){
	for (UShort_t i = 0;i <s1->NBOARDS_DSSD;i++){
		//	delete[] Gain[iboard];
		//	delete[] FeedBack[iboard];
		delete [] hTrace_Sum[i];
		delete [] hRaw[i];
		delete [] hCalib[i];
		delete [] hTrap[i];
		delete [] gr_baseline[i];
		delete [] hBaseline[i];
		delete [] hNoise[i];
		delete [] hRisetime[i];
		delete [] hTrigger[i];
	}

	//delete[] Gain;
	//delete[] FeedBack;
	delete [] hTrace_Sum;
	delete [] hRaw;
	delete [] hCalib;
	delete [] hTrap;
	delete [] hBaseline;
	delete [] hNoise;
	delete [] hRisetime;
	delete [] hTrigger;
	delete [] gr_baseline;
	//delete [] gr_rate_dssdBoard;
	//delete [] gr_rate_tunnelBoard;	
	for(int i = 0; i < s1->NBOARDS_TUNNEL; i++){
		delete [] h_tunnelRaw[i];
		delete [] h_tunnelCalib[i];
	}
	delete [] h_tunnelRaw;
	delete [] h_tunnelCalib;
	delete [] h_TUNNEL_XY_hit;
	delete [] dBoard;
	delete [] dBoard2;
}
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
