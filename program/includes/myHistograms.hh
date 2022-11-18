/*!
 * \file myHistograms.hh
 * \author Rikel Chakma
 * \date 18/11/2022
 * \brief header file of the myHistograms class
 *
 */

#ifndef myHistograms_h
#define myHistograms_h 1

#include "TGraph.h"
#include "global.h"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
class myHistograms
{

	private:
		myGlobal *s1;
		TDirectory *dFile;
		TDirectory **dBoard;   // for dssd
		TDirectory **dBoard2;// for tunnel   

	public:
		myHistograms();
		~ myHistograms();
		//------------DSSD histograms---------------
		//TH1I   *** hGain ; //!
		//TH1I   *** hFeedBack; //!
		TH1F     *** hTrigger;//!
		TH1F     *** hBaseline;//!
		TH1F     *** hNoise;//!
		TH1I     *** hRisetime;//!
		TH2F     *** hTrace_Sum ;  //! -- this silences the error
		TH1F     *** hRaw; //!
		TH1F     *** hCalib; //!
		TH2F     *** hTrap; //!
		TGraph   *** gr_baseline;//!
		//TGraph    ** gr_rate_dssdBoard;//!
		TH1I       * h_delT_ff; //!
		TH1I       * h_delT_fb; //!
		TH1I       * h_delT_bf; //!
		TH1I       * h_delT_bb; //!
		TH2F       * h_E_frontBack; //!
		TH2I       * h_DSSD_XY_hit; //!
		TH1I       * h_dssd_count_strip; //!
		TH1I       * h_dssd_count_board; //!
		TH2F       * h_calib_strip; //!
		TH2F       * h_raw_strip; //!
		//TH1I       * h_dssdStrip_rate; //!
		//TH1I       * h_dssdBoard_rate; //!
		TH2F       * h_dssd_baseline;
		TH2F       * h_dssd_noise;
		//---------TUNNEL Histograms----------
		TH1I     *** h_tunnelRaw;//!
		TH1I     *** h_tunnelCalib;//!
		//TGraph    ** gr_rate_tunnelBoard;//!
		TH2I      ** h_TUNNEL_XY_hit; //!
		TH2F       * h_calib_pad; //!
		TH2F       * h_raw_pad; //!
		//TH1I       * h_tunnelPad_rate; //!
		//TH1I       * h_tunnelBoard_rate; //!
		TH1I       * h_tunnel_count_pad; //!
		TH1I       * h_tunnel_count_board; //!

		void rawHistInitialization(TFile* f);
		void deleteRawHistograms();

};
//---------------ooooooooooooooo---------------ooooooooooooooo---------------ooooooooooooooo---------------
#endif
