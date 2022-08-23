#include "ReadDataFiles.hh"


ReadDataFiles::ReadDataFiles()
{
  pData = new signalTreatment();
}

ReadDataFiles::~ReadDataFiles()
{
  delete pData;
}
//main
void ReadDataFiles::perform()
{
  TFile* ofile = new TFile("../data/test.root", "RECREATE");
  hist = new myHistograms();
  hist->rawHistInitialization(ofile);
  
  for(int i = 0; i < 8;i++)frameCounter[i] =0;
  fFrame = new MFMCommonFrame();
  fSiriusframe = new MFMSiriusFrame();
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
  cout << "-------------------------------------------------------------" << endl;
  cout << "     test " ;
  cout << "|   READ MFMFrame from file : " << filename << "" << endl;
  cout << "-------------------------------------------------------------" << endl;
  cout << endl;

  while (true) {
    framesize = fFrame->ReadInFile(&fLun, pvector, &vectorsize);
    fPtInFile+=framesize;
    if (framesize <= 0) break;
    if(frameCounter[7] > 1) break;
    ReadUserFrame(fFrame);
    fNbOfEvent++;
    fFrame->SetAttributs();
  }

  fLun = close(fLun);

  for(int i = 0; i < 8;i++)cout<<"Frame: "<<frameName[i] <<"  no: "<<frameCounter[i]<<endl;
  cout << "-------------------End--------with "<<fNbOfEvent<<" frames----------------" << endl;
  if (fFrame)  delete (fFrame);
  if (fSiriusframe)     delete  (fSiriusframe) ;

  ofile->cd();
  ofile->Write();
  ofile->Close();
  delete ofile;
  hist->deleteRawHistograms();
  delete hist;
}

void ReadDataFiles::ReadUserFrame(MFMCommonFrame* commonframe) {
  type = commonframe->GetFrameType();
  commonframe->SetAttributs();

  switch (type) {

    case MFM_COBOF_FRAME_TYPE:
    case MFM_COBO_FRAME_TYPE: {
      // printf(" It is a Cobo frame\n");
      frameCounter[0]++;
      break;
    }
    case MFM_EBY_EN_FRAME_TYPE:
    case MFM_EBY_TS_FRAME_TYPE:
    case MFM_EBY_EN_TS_FRAME_TYPE: {
      // printf(" It is a Eby frame\n");
      frameCounter[1]++;
      break;
    }
    case MFM_HELLO_FRAME_TYPE: {
      // printf(" It is a Hello frame\n");
      frameCounter[2]++;
      break;
    }
    case MFM_XML_FILE_HEADER_FRAME_TYPE: {
      // printf(" It is a Header information frame\n");
      frameCounter[3]++;
      break;
    }

    case MFM_REA_TRACE_FRAME_TYPE: {
      // printf(" It is a Header information frame\n");
      frameCounter[4]++;
      break;
    }
    case MFM_REA_GENE_FRAME_TYPE: {
      //printf(" It is a Header information frame\n");
      frameCounter[5]++;
      //ReadDefaultFrame(commonframe);
      break;
    }
    case MFM_SIRIUS_FRAME_TYPE: {
      frameCounter[6]++;
      ReadSiriusFrame(commonframe);


      break;
    }
	
    default: {
      frameCounter[7]++;
      if(frameCounter[7] > 1) break;
      // printf("\nIt is a other frame\n");    
      break;
    }

  }// end of switch
}

void ReadDataFiles::ReadDefaultFrame(MFMCommonFrame* commonframe) {
  if (fVerbose > 1) {
    commonframe->HeaderDisplay(
        (char*) ((commonframe->WhichFrame()).c_str()));
    if (fVerbose > 3) {
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
//______________________________________________________________________________________________________________________
void ReadDataFiles::ReadSiriusFrame(MFMCommonFrame* commonframe) {
  // executed in case of Sirius Frame
  fSiriusframe->SetAttributs(commonframe->GetPointHeader());
 
  pData->set_channelID((UShort_t) fSiriusframe->GetChannelId());
  pData->set_boardID((UShort_t)fSiriusframe->GetBoardId());
  pData->set_timestamp(fSiriusframe->GetTimeStamp());
  pData->set_eventnumber(fSiriusframe->GetEventNumber());
  pData->set_gain(fSiriusframe->GetGain());
  //  cout<<"boardId:  "<<board<<" , channel:  "<<channel<<" ,    fSiriusframe->GetTimeStamp():  "<<timestamp<<
  //  " ,     fSiriusframe->GetEventNumber(): "<<fSiriusframe->GetEventNumber()<<endl;


  //hist->Trace[iboard][channel]->Reset();// Reset histograms

  if (fVerbose>3){
    fSiriusframe->HeaderDisplay();
  }
  if (fVerbose>5){
    if (fSiriusframe->GetFrameSize()  < maxdump) dumpsize = fSiriusframe->GetFrameSize();
    else dumpsize = maxdump;
    fSiriusframe->DumpRaw(dumpsize, 0);
  }
  // hist->Gain[iboard][channel]->Fill(fSiriusframe->GetGain());
  // hist->FeedBack[iboard][channel]->Fill(fSiriusframe->GetFeedBack(0));
  NbItems= fSiriusframe->GetNbItems();
  for (i = 0; i < NbItems; i++) {
    fSiriusframe->GetParameters(i, &value);
    pData->set_trace_value(i, value);
    //cout << " i = "<< i << " value = " << value << " Value reduced = "<<value/fReduce_factor<<"\n";
    //Trace[iboard][channel]->SetBinContent(i+1,value/fReduce_factor,1);
    //Trace[iboard][channel]->SetBinContent(i+1,value);
    //cout <<"  board = "<<iboard<<"  channel = "<<channel<<endl; 
    // hist->Trace_Sum[iboard][channel]->Fill(i,value);

  }

}
