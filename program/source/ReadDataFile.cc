#include "ReadDataFile.hh"

long get_file_size(const char * filename){
	FILE *fd = fopen(filename, "r");
	if(fd == NULL) {
		cerr << "Error opening file\n";
		return 0;
	}
	fseek(fd, 0, SEEK_END);
	long fileSize = ftell(fd);
	fclose(fd);
	return fileSize;
}

ReadDataFile::ReadDataFile()
{
	s1                             = myGlobal::getInstance();
	fFrame                         = new MFMCommonFrame();
	fCoboframe                     = new MFMCoboFrame();
	fInsideframe                   = new MFMCommonFrame();
	fMergeframe                    = new MFMMergeFrame();
	fSiriusframe                   = new MFMSiriusFrame();
	fGenericframe                  = new MFMReaGenericFrame();
	filter                         = new digitalFilters();
	cfd                            = new digitalCFD();
	dData                          = new dssdData();
	dEvent                         = new dssdEvent();
	tData                          = new tunnelData();
	calib                          = new calibration();
	hist                           = new myHistograms();
	frameCounter                   = new ullint[8];
	tunnel_rate_pad                = new double[s1->NofMacroPixels];
	tunnel_rate_board              = new double[s1->NBOARDS_TUNNEL];
	dssd_rate_strip                = new double[s1->NSTRIPS_DSSD];
	dssd_rate_board                = new double[s1->NBOARDS_DSSD]; 
	rate_counterPoint_dssd         = new ullint[s1->NBOARDS_DSSD];
	rate_counterPoint_tunnel       = new ullint[s1->NBOARDS_TUNNEL];
	dssd_event_counter             = new ullint*[s1->NBOARDS_DSSD];
	dumpsize                       = 0;
	framesize                      = 0;
	maxdump                        = 128;
	channel                        = 0;
	board                          = 0;
	iboard                         = 0;
	value                          = 0;
	trace_size                     = 0;
	rate_counter_dssd              = 0;
	rate_counter_tunnel            = 0;
	dataSet_counter                = 0.;
	rate_calcul_timestamp_lapse    = static_cast<int64_t>(s1->rate_calcul_time_lapse *pow(10,8));
	//--------------
	for(int i = 0; i < 8; i++)frameCounter[i] = 0;
	for(int i = 0; i < s1->NSTRIPS_DSSD; i++){
		dssd_rate_strip[i] =0.;
	}
	for(int i = 0; i < s1->NBOARDS_DSSD; i++){
		rate_counterPoint_dssd[i]=0;	
		dssd_rate_board[i] =0.;
	}
	for(int i = 0; i < s1->NofMacroPixels; i++){
		tunnel_rate_pad[i] =0.;
	}
	for(int i = 0; i < s1->NBOARDS_TUNNEL; i++){
		rate_counterPoint_tunnel[i]=0;	
		tunnel_rate_board[i] =0.;
	}
	for(int i = 0; i < s1->NBOARDS_DSSD; i++){
		dssd_event_counter[i] = new unsigned long long int[NCHANNELS];
	}
	//initialize some preset value
	for(int i = 0; i < s1->NBOARDS_DSSD; i++){
		for(int j = 0; j < NCHANNELS; j++){
			dssd_event_counter[i][j] =0;
		}
	}
	std::string dssd_boards_char[s1->NBOARDS_DSSD];
	for(int i = 0; i < s1->NBOARDS_DSSD;i++) dssd_boards_char[i] = std::to_string(s1->boardList_DSSD[i]);
	std::string tunnel_boards_char[s1->NBOARDS_TUNNEL];
	for(int i = 0; i < s1->NBOARDS_TUNNEL;i++) tunnel_boards_char[i] = std::to_string(s1->boardList_Tunnel[i]);

}

ReadDataFile::~ReadDataFile()
{
	for(int i = 0; i < 8;i++)cout<<"Frame: "<<s1->frameName[i] <<"  no: "<<frameCounter[i]<<endl;
	for(int i = 0; i < s1->NBOARDS_DSSD; i++){
		delete [] dssd_event_counter[i];
	}
	delete [] dssd_event_counter;

	dssdDataVec.clear();
	dssdEventVec.clear();

	delete [] tunnel_rate_pad;
	delete [] tunnel_rate_board;
	delete [] frameCounter;
	delete [] dssd_rate_strip;
	delete [] dssd_rate_board;
	delete [] rate_counterPoint_dssd;
	delete [] rate_counterPoint_tunnel;

	if (fFrame)  delete (fFrame);
	if(fInsideframe) delete fInsideframe;
	if(fMergeframe) delete fMergeframe;
	if(fSiriusframe) delete fSiriusframe;
	if(fGenericframe) delete fGenericframe;
	if(fCoboframe) delete fCoboframe;
	if(filter) delete filter;
	if(cfd) delete cfd;
	if(dEvent) delete dEvent;
	if(calib) delete calib;
	if(tData) delete tData;
	if(dData) delete dData;
	if(hist) delete hist;	

}

void ReadDataFile::initialize_spectra(const char* histFileName){

	oHistFile = new TFile(histFileName, "RECREATE");
	hist->rawHistInitialization(oHistFile);
} 
void ReadDataFile::initialize_ttree(const char* treeFileName, const char* treeName){
	oTreeFile = new TFile(treeFileName, "RECREATE");
	oTree = new TTree(treeName, treeName);
	oTree->Branch("Time", &Time, "Time/l");
	oTree->Branch("EventNo",  &eventnumber, "EventNo/i");
	oTree->Branch("TraceSize",  &trace_size, "trace_size/s");
	oTree->Branch("Trace",  dData->get_trace(), "Trace[trace_size]/s");
	oTree->Branch("Gain",  &gain, "Gain/s");
	oTree->Branch("BoardID",  &board, "BoardID/s");
	oTree->Branch("ChannelID",  &channel, "ChannelID/s");
	oTree->Branch("Energy",  &Energy, "Energy/d");
	oTree->Branch("Baseline",  dData->get_Baseline_address(), "Baseline/d");
	oTree->Branch("Noise",  dData->get_Noise_address(), "Noise/d");	

}
void ReadDataFile::save_spectra(){
	oHistFile->cd();
	oHistFile->Write();
	oHistFile->Close();
	hist->deleteRawHistograms();
	cout<<"Spectra saved..."<<endl;	
}
void ReadDataFile::save_ttree(){
	oTreeFile->cd();
	oTree->Write();
	oTreeFile->Close();
	cout<<"Ttree saved..."<<endl;	
}



void ReadDataFile::initialize(){
	for(int i = 0; i < 8;i++)frameCounter[i] =0;
}
void ReadDataFile::Read(const char * filename, const char * nEvent )
{

	long fileSize = get_file_size(filename);
	long fileSize_kbytes = fileSize/1024.;

	max_nEvents_to_process = atoi(nEvent);
	int fLun = 0; // Logical Unit Number
	fLun = open(filename, (O_RDONLY));
	long long fPtInFile =0;
	long long fNbOfEvent=0;
	char* vector;
	char ** pvector;
	int vectorsize = 8; // minsizeheader		
	vector = (char*) (malloc(vectorsize));
	pvector = &vector;
	int framesize = 0;

	if (fLun <= 0) {
		printf("Error open file : %s", filename);
	}
	cout << endl;
	cout << "-----------------------------------------------------------------------------" << endl;
	cout << "| File : " << filename << endl;
	cout << "-----------------------------------------------------------------------------" << endl;
	cout << endl;

	while (true) {
		framesize = fFrame->ReadInFile(&fLun, pvector, &vectorsize);
		fPtInFile+=framesize;
		cout<<"\rRead: "<<int(double(fPtInFile)/double(fileSize) * 100) << "% of "<<fileSize_kbytes<<" KB" <<std::flush;
		if (framesize <= 0) break;
		if(frameCounter[7] > 1) break;
		ReadUserFrame(fFrame);
		fFrame->SetAttributs();
		if(fNbOfEvent == max_nEvents_to_process) break;
		fNbOfEvent++;
	}

	fLun = close(fLun);
	cout << endl;
	cout << "-------------------End--------with "<<fNbOfEvent<<" frames----------------" << endl;

}

void ReadDataFile::ReadUserFrame(MFMCommonFrame* commonframe) {
	type = commonframe->GetFrameType();
	commonframe->SetAttributs();

	switch (type) {
		case MFM_COBOF_FRAME_TYPE: 
		case MFM_COBO_FRAME_TYPE: 
			frameCounter[0]++;
			ReadCoboFrame(commonframe);
			break;
		case MFM_EBY_EN_FRAME_TYPE:
		case MFM_EBY_TS_FRAME_TYPE:
		case MFM_EBY_EN_TS_FRAME_TYPE: 
			frameCounter[1]++;
			break;
		case MFM_HELLO_FRAME_TYPE: 
			frameCounter[2]++;
			break;
		case MFM_XML_FILE_HEADER_FRAME_TYPE: 
			frameCounter[3]++;
			break;
		case MFM_REA_TRACE_FRAME_TYPE: 
			frameCounter[4]++;
			break;
		case MFM_REA_GENE_FRAME_TYPE: 
			frameCounter[5]++;
			ReadGenericFrame(commonframe);
			break;
		case MFM_SIRIUS_FRAME_TYPE: 
			frameCounter[6]++;
			ReadSiriusFrame(commonframe);
			break;
		case MFM_MERGE_EN_FRAME_TYPE:
		case MFM_MERGE_TS_FRAME_TYPE:
			ReadMergeFrame(commonframe);
			break;
		default: 
			frameCounter[7]++;
			if(frameCounter[7] > 1) break;
			break;
	}// end of switch
	if(oTree)oTree->Fill();
}

void ReadDataFile::ReadDefaultFrame(MFMCommonFrame* commonframe) {
	if (s1->fverbose > 1) {
		commonframe->HeaderDisplay(
				(char*) ((commonframe->WhichFrame()).c_str()));
		if (s1->fverbose > 3) {
			int framesize = commonframe->GetFrameSize();
			int dump = fDumpsize;
			if (framesize < dump)
				dump = framesize;
			if (dump < 4)
				dump = 4;
			commonframe->DumpRaw(dump, 0);
		}
	}
}

//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
//! Treatment of merged frames
/*! The frames within the defined time window are groupped into merged frames. Normally, if everything is ok, the DSSD pixel construction and establishing correlations among the events from different detector are done here. The only catch of enabling the Merger is it consumes a lot of memory.   
*/ 
void ReadDataFile::ReadMergeFrame(MFMCommonFrame* commonframe){
	int i_insframe = 0;
	int nbinsideframe = 0;
	int dumpsize = 16;

	fMergeframe->SetAttributs(commonframe->GetPointHeader());

	nbinsideframe = fMergeframe->GetNbItems();
	//cout<<"------------------- ninside frame = "<<nbinsideframe<<endl;
	framesize= fMergeframe->GetFrameSize();
	if (s1->fverbose >= 3){
		fMergeframe->HeaderDisplay();
	}
	if (s1->fverbose >= 4){
		commonframe->HeaderDisplay();
	}
	if (s1->fverbose >= 5){
		if (framesize < maxdump) dumpsize = framesize;
		else dumpsize = maxdump;
		fMergeframe->DumpRaw(dumpsize, 0);
	}
	fMergeframe->ResetReadInMem();
	//Reset vectors
	dssdDataVec.clear();
	for(i_insframe = 0; i_insframe < nbinsideframe; i_insframe++) {
		//fMergeframe->ReadInFrame(fInsideframe);
		fMergeframe->ReadInFrame(commonframe);
		ReadUserFrame(commonframe);
		//ReadUserFrame(fInsideframe);

	}

	// At this point you can do treatement inter frames


}
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
//! Treatment of Generic frames (for the Tunnel detectors)
/*! The energy and timestamps are obtained from the data. To get the board number, macro pixel number and the detector number, use the methods from the tunnelData class. The macro pixel configurations are specified in the Run.config. The configuration can be changed (for example ABCD to BADC) keeping the macro pixels definitions constant i.e., if pixels 1 and 2 are grouped as A, you should not define it as B. If you do so, remember to update in tunnelData class, otherwise, the correct picture (2D hit pattern) will not be displayed.
 * To calibrate do Calibration::perform(tunnelData* const data).
 * To get the board number use tunnelData::get_tunnelBoardNumber(int* p_board).
 * To get the physical macro pixel number use tunnelData::get_macroPixelPhysicalNumber(int* p_board, int* p_channel). This method also detemines the position of the macro pixel.
 * To get the detector number use  tunnelData::get_tunnelDetectorNumber(int* p_board, int* p_channel).
 */
void ReadDataFile::ReadGenericFrame(MFMCommonFrame* commonframe)
{
	fGenericframe->SetAttributs(commonframe->GetPointHeader());
	framesize=commonframe->GetFrameSize();
	channel =fGenericframe->GetChannelId();
	board =  fGenericframe->GetBoardId();
	iboard =  s1->boardIndex_Tunnel[board];
	value = fGenericframe->GetEnergy();
	trace_size =0;
	timestamp = fGenericframe->GetTimeStamp();
	tData->set_channelID(channel);
	tData->set_boardID(board);
	tData->set_boardIdx(iboard);
	tData->set_timestamp( timestamp);
	tData->set_eventnumber( eventnumber);
	//h_tunnelRaw[iboard][channel]->Fill (value);
	Energy = static_cast<Double_t> (value);
	calibEnergy = calib->perform(tData);
	//h_tunnelCalib[iboard][channel]->Fill (Energy);
	Time = timestamp;
	tunnelBoardNo = tData->get_tunnelBoardNumber(&board); 
	tunnelPadNo = tData->get_macroPixelPhysicalNumber(&board, &channel);//gets the pixels
	tunnelDetectorNo = tData->get_tunnelDetectorNumber(&board, &channel);
	//h_tunnel_count_board->Fill(tunnelBoardNo);
	//h_tunnel_count_pad->Fill(tunnelPadNo);

	//for(unsigned int i = 0; i< tData->get_macroPixel(&board,&channel).pixels.size();i++){
	//h_TUNNEL_XY_hit[tunnelDetectorNo-1]->Fill(tData->get_macroPixel(&board,&channel).pixels[i].get_X(),tData->get_macroPixel(&board,&channel).pixels[i].get_Y());
	//}
	if(s1->fverbose >= 3){
		cout<<"-------Generic Frame-------\n";
		printf (" board = %d , channel =%d , energy=%d , timestamp = %llu\n", board, channel, value, timestamp);
	}
	if (s1->fverbose >= 4){
		commonframe->HeaderDisplay();
	}
	if (s1->fverbose >= 4){
		if (framesize < maxdump) dumpsize = framesize;else dumpsize = maxdump;
		commonframe->DumpRaw(dumpsize, 0);
	}
}
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
//! Treatment of Sirius frames
/*!
 *
 *
 */
void ReadDataFile::ReadSiriusFrame( MFMCommonFrame* commonframe){
	fSiriusframe->SetAttributs(commonframe->GetPointHeader());
	framesize =fSiriusframe->GetFrameSize();
	channel = (int)fSiriusframe->GetChannelId();
	board = (int)fSiriusframe->GetBoardId();
	iboard =  (int)s1->boardIndex_DSSD[board];
	timestamp = fSiriusframe->GetTimeStamp();
	eventnumber = fSiriusframe->GetEventNumber();
	gain = fSiriusframe->GetGain();
	trace_size =s1->TRACE_SIZE;
	dData->set_channelID(channel);
	dData->set_boardID(board);
	dData->set_boardIdx(iboard);
	dData->set_timestamp( timestamp);
	dData->set_eventnumber( eventnumber);
	dData->set_gain(gain);
	stripnumber = dData->get_stripnumber(&board, &channel);
	//h_dssd_count_strip->Fill(stripnumber );
	dssdBoardNo = dData->get_dssdBoardNumber(&board); 
	//h_dssd_count_board->Fill(dssdBoardNo);

	if(s1->fverbose >= 3){
		cout<<"-------Sirius Frame-------\n";
		cout<<"boardId:  "<<board<<" , channel:  "<<channel<<"  strip "<<stripnumber<<" TimeStamp:  "<<timestamp<<" , EventNumber: "<<eventnumber<<endl;
	}
	if (s1->fverbose >= 4){
		commonframe->HeaderDisplay();
	}

	if (s1->fverbose >= 5){
		if (framesize < maxdump) dumpsize = framesize;else dumpsize = maxdump;
		commonframe->DumpRaw(dumpsize, 0);
	}

	// hGain[iboard][channel]->Fill(gain);
	// hFeedBack[iboard][channel]->Fill(fSiriusframe->GetFeedBack(0));
	// NbItems= fSiriusframe->GetNbItems();
	int j =0;	
	for (int i = s1->nStart_trace; i < 992 - s1->nEnd_trace; i++) {
		fSiriusframe->GetParameters(i, &value);
		j = i - s1->nStart_trace;
		//hTrace_Sum[iboard][channel]->Fill(j,value);
		if(i < s1->TRACE_SIZE)
			dData->set_trace_value(j, value);
	}
	//---------------------------
	dData->GetSignalInfo();
	//Fill histograms
	//hBaseline[iboard][channel]->Fill(dData->get_Baseline());
	//hNoise[iboard][channel]->Fill(dData->get_Noise());
	//hRisetime[iboard][channel]->Fill(dData->get_RiseTime());
	//Energy = filter->perform(dData, hTrap[iboard][channel]);
	Energy = filter->perform(dData, NULL);
	//Energy  = dData->get_signalHeight();
	//h_raw_strip->Fill(Energy, stripnumber);
	//hRaw[iboard][channel]->Fill(Energy);
	calibEnergy = calib->perform(dData);//Calibration
	//h_calib_strip->Fill(dData->get_calibrated_energy(),stripnumber);
	//hRaw[iboard][channel]->Fill(dData->get_raw_energy());
	//get CFD time in (ns)
	Time = cfd->perform(dData);

	//Fill baseline graphs
	//gr_baseline[iboard][channel]->SetPoint(dssd_event_counter[iboard][channel],dssd_event_counter[iboard][channel],dData->get_Baseline());
	//	dssd_event_counter[iboard][channel]++;

	//	if(Time < 1E15){
	//----this condition ensures that the time stamp is ok
	//	dPoint.set_time(Time);
	//	dPoint.set_strip(stripnumber);
	//	dPoint.set_energy(dData->get_calibrated_energy());
	//	dssdDataVec.push_back(dPoint);
	//if(s1->data_merged == 0){
	//if(dssdDataVec.size() == s1->buffer_size){
	//dssdEventVec = dEvent->construct(dssdDataVec, h_delT_ff, h_delT_fb, h_delT_bf, h_delT_bb);
	//for (std::vector<dssdPixel>::iterator it = dssdEventVec.begin() ; it != dssdEventVec.end(); ++it){				
	//h_E_frontBack->Fill((*it).get_energyX(), (*it).get_energyY());
	//h_DSSD_XY_hit->Fill((*it).get_X(), (*it).get_Y());
	//}
	//}
	//}
	//}
}
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
//!  Treatment of Cobo frames
/*!
 *
 *
 *
 */
void ReadDataFile::ReadCoboFrame(MFMCommonFrame* commonframe){
	fCoboframe->SetAttributs(commonframe->GetPointHeader());
	framesize=fCoboframe->GetFrameSize();
	eventnumber =fCoboframe->GetEventNumber();
	timestamp = (uint64_t)(fCoboframe->GetTimeStamp());
	trace_size =0;
	if(s1->fverbose >=3){
		cout<<"-------Cobo Frame-------\n";
		cout << " EN = "<<eventnumber <<" TS = " << timestamp <<endl;
	}
	if (s1->fverbose>=4){
		fCoboframe->HeaderDisplay();
	}
	if (s1->fverbose>=5){
		if (framesize < maxdump) dumpsize = framesize;
		else dumpsize = maxdump;
		fCoboframe->DumpRaw(dumpsize, 0);
	}

}
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
