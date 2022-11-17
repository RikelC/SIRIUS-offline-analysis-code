#ifndef ReadDataFile_h
#define ReadDataFile_h 1
#include "myHistograms.hh"
#include "TFile.h"
#include "TTree.h"
#include "MFMExogamFrame.h"
#include "MFMMergeFrame.h"
#include "MFMCoboFrame.h"
#include "MFMExogamFrame.h"
#include "MFMEbyedatFrame.h"
#include "MFMNumExoFrame.h"
#include "MFMSiriusFrame.h"
#include "MFMAllFrames.h"
#include "MFMReaGenericFrame.h"
#include "Rtypes.h"
#include "TRandom.h"
#include "global.h"
#include "dssdData.h"
#include "tunnelData.h"
#include "digitalFilters.h"
#include "digitalCFD.h"
#include "calibration.h"
#include "dssdDataPoint.h"
#include "dssdEvent.h"

#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>    // std::sort
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <csignal>
#include <stdlib.h>
#include <unistd.h>

#define DSSD_DETECTOR 1
#define TUNNEL_DETECTOR 2

using namespace std;
class ReadDataFile
{
	private:
		MFMCommonFrame            * fFrame ; 
		MFMMergeFrame             * fMergeframe ;/**< MergerFrame type */
		MFMCoboFrame              * fCoboframe ;/**< CoboFrame for the tracker data */
		MFMCommonFrame            * fInsideframe ; /**< MFMCommonFrame type. Used for pointing to frames inside a MFMMergeFrame */
		MFMSiriusFrame            * fSiriusframe ; /**< MFMSiriusFrame type for treating DSSD data.*/
		MFMReaGenericFrame        * fGenericframe;/**< MFMReaGenericFrame type for treating TUnnel data. */

		myHistograms *hist;
		myGlobal                  * s1;/**< pointer to the instance of the Singleton myGlobal class. */
		digitalFilters            * filter;/**< instance for the digital filter algorithms. */
		digitalCFD                * cfd; /**< instance for the digital constant fraction discriminator. */
		calibration               * calib; /**< instance for doing calibration of the detectors. */
		dssdData                  * dData; /*!< dssdData object*/
		dssdEvent                 * dEvent;/*!< dssdEvent object*/
		tunnelData                * tData;/*!< tunnelData object*/
		double     * tunnel_rate_pad; //!
		double     * tunnel_rate_board; //!
		double     * dssd_rate_strip; //!
		double     * dssd_rate_board; //!
		ullint     * rate_counterPoint_dssd;//!
		ullint     * rate_counterPoint_tunnel;//!
		ullint    ** dssd_event_counter;//!
		ullint      * frameCounter; //!

		double       dataSet_counter;
		uint         eventnumber;
		ullint       timestamp;
		ullint       old_timeStamp;
		ullint       rate_counter_dssd;
		ullint       rate_counter_tunnel;
		llint        deltaT_rates;
		ushort       gain;
		ushort       trace_size;
		int          dumpsize;
		int          framesize;
		int          maxdump;
		int          channel;
		int          board;
		int          iboard;
		int          stripnumber;
		int          dssdBoardNo;
		int          tunnelBoardNo;
		int          tunnelPadNo;
		int          tunnelDetectorNo;
		uint16_t     value;
		int64_t      rate_calcul_timestamp_lapse;		
		Double_t     Energy;
		double       calibEnergy;
		ULong64_t    Time;
		int metatype ;
		int blocksize ;
		bool endiness ;
		bool blobness ;
		int fDumpsize =32;
		int datasource;
		int type;
		int rivision;

		int headersize;

		int max_nEvents_to_process;
		TFile* oHistFile;
		TFile* oTreeFile;
		TTree* oTree = NULL;
		dssdDataPoint dPoint;
		std::vector<dssdDataPoint>dssdDataVec;
		std::vector<dssdPixel>dssdEventVec;
		void get_Count_Rates( int &board, int &stripNumber, ullint &time, int type);
		void ReadUserFrame(MFMCommonFrame* commonframe);
		void ReadDefaultFrame(MFMCommonFrame* commonframe) ;
		void ReadSiriusFrame(MFMCommonFrame* commonframe) ;
		void ReadMergeFrame(MFMCommonFrame * commonframe);
		void ReadGenericFrame(MFMCommonFrame * commonframe);
		void ReadCoboFrame(MFMCommonFrame * commonframe);
	public:
		ReadDataFile();
		~ReadDataFile();

		void Read(const char *, const char *);
		void initialize();
		void save_spectra();
		void save_ttree();
		void initialize_spectra(const char*); 
		void initialize_ttree(const char*, const char*);
};

#endif
