// $Id$
//
//    File: JEventProcessor_FCALPpPi.h
// Created: Thu May 28 22:18:51 EDT 2015
// Creator: manlara (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_FCALPpPi_
#define _JEventProcessor_FCALPpPi_

#include <JANA/JEventProcessor.h>

#include <stdint.h>
#include <stdlib.h>
#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include "StringUtilities.h"

#include "TROOT.h"
#include <TString.h>
#include <TH1F.h>
#include <TF1.h>
#include "TMinuit.h"
#include <TVirtualFitter.h>
#include <TTree.h>
#include <TObjString.h>
#include <TObjArray.h>

#include "Pulse.h"

class JEventProcessor_FCALPpPi:public jana::JEventProcessor{
	public:
		JEventProcessor_FCALPpPi();
		~JEventProcessor_FCALPpPi();
		const char* className(void){return "JEventProcessor_FCALPpPi";}
		
		TTree* FCAL_Pulse;
		TTree* FCAL_Pulse2;
    TTree* FCAL_Analysis;
    TTree* FCAL_Analysis2;
    
    TTree* FCAL_Correlation;
    TTree* FCAL_SingleCorrelation;
    
    TTree* FCAL_Waveforms;
    
    std::map<TString, Pulse> pulse_ana;
    std::map<TString, Pulse> pulse_nsp;
    map< pair< pair<int,int>, pair<int,int> >, int > corrmap;
    
    const static int min_nsa = 2;
    const static int max_nsa = 80;
    const static int min_nsp = 0;
    const static int max_nsp = 10;
    
    // key = rocid/slot/channel
    // value = (x,y)
    std::map<TString, pair<int,int> > tranlation_map;
		
		uint32_t rocid;
    uint32_t slot;
    uint32_t channel;
    int x;
    int y;
    uint32_t TC;
    uint32_t nsa;
    uint32_t nsb;
    uint32_t nsp;
    uint32_t start;
    uint32_t end;
    uint32_t threshold;
    uint32_t peak;
    uint32_t integral;
    uint32_t eventnum;
    uint32_t run;
    
    double avg_peak;
    double avg_integral;
    double rms_peak;
    double rms_integral;
    double rms_peak_gaus;
    double rms_integral_gaus;
    
    uint32_t pulse1_rocid   ;
    uint32_t pulse1_slot    ;
    uint32_t pulse1_channel ;
    int pulse1_x            ;
    int pulse1_y            ;
    uint32_t pulse2_rocid   ;
    uint32_t pulse2_slot    ;
    uint32_t pulse2_channel ;
    int pulse2_x;       
    int pulse2_y;       
    double corr;
    double peak_ratio_mean;
    double peak_ratio_rms ;
    double peak_ratio_rms_gaus;
    
    double pulse1_peak      ;
    double pulse1_integral  ;
    double pulse2_peak      ;
    double pulse2_integral  ;
    
    double pulse1_avgPeak;
    double pulse2_avgPeak;

    uint32_t sample;
    uint32_t adc_val;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_FCALPpPi_

