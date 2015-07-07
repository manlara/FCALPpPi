// $Id$
//
//    File: JEventProcessor_FCALPpPi.cc
// Created: Thu May 28 22:18:51 EDT 2015
// Creator: manlara (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_FCALPpPi.h"
using namespace jana;


#include <stdint.h>
// https://halldsvn.jlab.org/repos/branches/sim-recon-commissioning/src/libraries/DAQ/
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseRawData.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulseTime.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250TriggerTime.h>

using namespace std;

// FCAL Hit Objects
// https://halldsvn.jlab.org/repos/branches/sim-recon-commissioning/src/libraries/FCAL/

// Other peoples plugins
// https://halldsvn.jlab.org/repos/branches/sim-recon-commissioning/src/programs/Analysis/plugins/


bool Df250WindowRawData_cmp(const Df250WindowRawData *a,const Df250WindowRawData *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250PulseRawData_cmp(const Df250PulseRawData *a,const Df250PulseRawData *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250PulseIntegral_cmp(const Df250PulseIntegral *a,const Df250PulseIntegral *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250PulseTime_cmp(const Df250PulseTime *a,const Df250PulseTime *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250PulsePedestal_cmp(const Df250PulsePedestal *a,const Df250PulsePedestal *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250TriggerTime_cmp(const Df250TriggerTime *a,const Df250TriggerTime *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	return a->itrigger < b->itrigger;
}

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_FCALPpPi());
}
} // "C"


//------------------
// JEventProcessor_FCALPpPi (Constructor)
//------------------
JEventProcessor_FCALPpPi::JEventProcessor_FCALPpPi()
{

}

//------------------
// ~JEventProcessor_FCALPpPi (Destructor)
//------------------
JEventProcessor_FCALPpPi::~JEventProcessor_FCALPpPi()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FCALPpPi::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//
  
  japp->RootWriteLock();
  
  // Initialize FCAL Channels
  
  nsb = 3;
  threshold = 100 + 160;
  
  // The rocid is between rocid>=11 && rocid<=22
  int max_rocid = 22, min_rocid=11;
  // 3 to 19
  // skip 11 and 12
  int max_slot = 19, min_slot=3;
  // 16 channels starting from 0
  int n_channels = 16;
  for (int roc=min_rocid; roc<=max_rocid; roc++){
    for (int slot=min_slot; slot<=max_slot; slot++){
      for (int channel=0; channel<n_channels; channel++){
        for (int m_nsa=min_nsa; m_nsa<=max_nsa; m_nsa++){
          char name[50];
          sprintf(name,"%02i/%02i/%02i",roc,slot,channel);
          TString key = TString(name)+"_"+StringUtilities::int2TString(m_nsa);
          Pulse::Pulse pulse(key);
          pulse_ana[key] = pulse;
        }
        for (int m_nsp=min_nsp; m_nsp<=max_nsp; m_nsp++){
          char name[50];
          sprintf(name,"%02i/%02i/%02i",roc,slot,channel);
          TString key = TString(name)+"_"+StringUtilities::int2TString(m_nsp);
          Pulse::Pulse pulse(key);
          pulse_nsp[key] = pulse;
        }
      }
    }
  }
  japp->RootUnLock();
  
  // Setup translation table between crate/slot/channel -> x,y
  // loop over DAQCoordinate_Ordered.txt where crate starts from 0
  japp->RootWriteLock();

  string daq_loc;
  int abs_num;
  string det_loc;

  string line;
  std::istringstream lin;
  ifstream inFile("DAQCoordinate_Ordered.txt");
  if(inFile.is_open())
  {
    while( getline (inFile,line) )
    {
      lin.clear();
      lin.str(line);
      std::istringstream iss(line);

      if (iss >> daq_loc >> abs_num >> det_loc){
        vector<TString> parseDAQ = StringUtilities::parseTString(TString(daq_loc),"/");
        vector<TString> parseLOC = StringUtilities::parseTString(TString(det_loc),"/");
        int rocid = parseDAQ[0].Atoi()+min_rocid, slot = parseDAQ[1].Atoi(), channel = parseDAQ[2].Atoi();
        int m_x = parseLOC[0].Atoi(), m_y = parseLOC[1].Atoi();
        char name[50];
        sprintf(name,"%02i/%02i/%02i",rocid,slot,channel);
        tranlation_map[TString(name)] = make_pair(m_x,m_y);
      }
    }
  }
  
  // to avoid double counting of correlation
  // make a map of pulse1,pulse2 and set to 0
  for (int m_x=-30; m_x<=30; m_x++){
    for (int m_y=-30; m_y<=30; m_y++){
      pair<int,int> m_pulse1 = make_pair(m_x,m_y);
      for (int m_x2=-30; m_x2<=30; m_x2++){
        for (int m_y2=-30; m_y2<=30; m_y2++){
          pair<int,int> m_pulse2 = make_pair(m_x2,m_y2);
          corrmap[ make_pair(m_pulse1,m_pulse2) ] = 0;
        }
      }
    }
  }
  
  japp->RootUnLock();
	
	
	FCAL_Pulse = new TTree("FCAL_Pulse","Pulse information");
  FCAL_Pulse->Branch("rocid",    &rocid   , "rocid/i");    // vme crate number
  FCAL_Pulse->Branch("slot",     &slot    , "slot/i");     // board number
  FCAL_Pulse->Branch("channel",  &channel , "channel/i");  // board channel number
  FCAL_Pulse->Branch("x",        &x   ,     "x/I");        // x coordinate where the center of (0,0)
  FCAL_Pulse->Branch("y",        &y   ,     "y/I");        // y coordinate where the center of (0,0)
  FCAL_Pulse->Branch("nsa",    &nsa   , "nsa/i");          // number of samples used in integral after crossing time
  FCAL_Pulse->Branch("nsb",    &nsb   , "nsb/i");          // number of samples used in integral before crossing time
  FCAL_Pulse->Branch("start",    &start   , "start/i");          // start sample
  FCAL_Pulse->Branch("end",    &end   , "end/i");          // end sample
  FCAL_Pulse->Branch("threshold",  &threshold   , "threshold/i");  // threshold
  FCAL_Pulse->Branch("peak", &peak, "peak/i"); // pulse peak in event
  FCAL_Pulse->Branch("integral", &integral, "integral/i");      // sum between samples TC-NSB to TC+NSA
  FCAL_Pulse->Branch("eventnum", &eventnum, "eventnum/i"); // event number
  FCAL_Pulse->Branch("run",      &run,      "run/i");      // run number  
  FCAL_Pulse->Branch("TC",      &TC,      "TC/i");      // time of threshold crossing
  
  FCAL_Pulse2 = new TTree("FCAL_Pulse2","Pulse information 2");
  FCAL_Pulse2->Branch("rocid",    &rocid   , "rocid/i");    // vme crate number
  FCAL_Pulse2->Branch("slot",     &slot    , "slot/i");     // board number
  FCAL_Pulse2->Branch("channel",  &channel , "channel/i");  // board channel number
  FCAL_Pulse2->Branch("x",        &x   ,     "x/I");        // x coordinate where the center of (0,0)
  FCAL_Pulse2->Branch("y",        &y   ,     "y/I");        // y coordinate where the center of (0,0)
  FCAL_Pulse2->Branch("nsp",    &nsp   , "nsp/i");          // number of samples used in integral around the pulse peak
  FCAL_Pulse2->Branch("start",    &start   , "start/i");          // start sample
  FCAL_Pulse2->Branch("end",    &end   , "end/i");          // end sample
  FCAL_Pulse2->Branch("threshold",  &threshold   , "threshold/i");  // threshold
  FCAL_Pulse2->Branch("peak", &peak, "peak/i"); // pulse peak in event
  FCAL_Pulse2->Branch("integral", &integral, "integral/i");      // sum between peakSample +/- NSP
  FCAL_Pulse2->Branch("eventnum", &eventnum, "eventnum/i"); // event number
  FCAL_Pulse2->Branch("run",      &run,      "run/i");      // run number  
  
  // use nsa + nsb for integral
  FCAL_Analysis = new TTree("FCAL_Analysis","Pulse analysis");
  FCAL_Analysis->Branch("rocid",    &rocid   , "rocid/i");    // vme crate number
  FCAL_Analysis->Branch("slot",     &slot    , "slot/i");     // board number
  FCAL_Analysis->Branch("channel",  &channel , "channel/i");  // board channel number
  FCAL_Analysis->Branch("x",        &x   ,     "x/I");        // x coordinate where the center of (0,0)
  FCAL_Analysis->Branch("y",        &y   ,     "y/I");        // y coordinate where the center of (0,0)
  FCAL_Analysis->Branch("nsa",    &nsa   , "nsa/i");          // number of samples used in integral after crossing time
  FCAL_Analysis->Branch("nsb",    &nsb   , "nsb/i");          // number of samples used in integral before crossing time
  FCAL_Analysis->Branch("start",    &start   , "start/i");          // start sample
  FCAL_Analysis->Branch("end",    &end   , "end/i");          // end sample
  FCAL_Analysis->Branch("threshold",  &threshold   , "threshold/i");  // threshold
  FCAL_Analysis->Branch("avg_peak", &avg_peak, "avg_peak/D"); // average of pulse peak over all events
  FCAL_Analysis->Branch("avg_integral",  &avg_integral, "avg_integral/D");      // average of pulse integral over all events
  FCAL_Analysis->Branch("rms_peak", &rms_peak, "rms_peak/D"); // rms of average of pulse peak
  FCAL_Analysis->Branch("rms_integral",&rms_integral, "rms_integral/D");   // rms of average of pulse integral
  FCAL_Analysis->Branch("rms_peak_gaus", &rms_peak_gaus, "rms_peak_gaus/D"); // rms of average of pulse peak
  FCAL_Analysis->Branch("rms_integral_gaus",&rms_integral_gaus, "rms_integral_gaus/D");   // rms of average of pulse integral
  FCAL_Analysis->Branch("eventnum", &eventnum, "eventnum/i"); // event number
  FCAL_Analysis->Branch("run",      &run,      "run/i");      // run number
  FCAL_Analysis->Branch("TC",      &TC,      "TC/i");      // time of threshold crossing
  
  // using the samples around the peak for the integral
  FCAL_Analysis2 = new TTree("FCAL_Analysis2","Pulse analysis 2");
  FCAL_Analysis2->Branch("rocid",    &rocid   , "rocid/i");    // vme crate number
  FCAL_Analysis2->Branch("slot",     &slot    , "slot/i");     // board number
  FCAL_Analysis2->Branch("channel",  &channel , "channel/i");  // board channel number
  FCAL_Analysis2->Branch("x",        &x   ,     "x/I");        // x coordinate where the center of (0,0)
  FCAL_Analysis2->Branch("y",        &y   ,     "y/I");        // y coordinate where the center of (0,0)
  FCAL_Analysis2->Branch("nsp",    &nsp   , "nsp/i");          // number of samples used in integral around the pulse peak
  FCAL_Analysis2->Branch("start",    &start   , "start/i");          // start sample
  FCAL_Analysis2->Branch("end",    &end   , "end/i");          // end sample
  FCAL_Analysis2->Branch("threshold",  &threshold   , "threshold/i");  // threshold
  FCAL_Analysis2->Branch("avg_peak", &avg_peak, "avg_peak/D"); // average of pulse peak over all events
  FCAL_Analysis2->Branch("avg_integral",  &avg_integral, "avg_integral/D");      // average of pulse integral over all events
  FCAL_Analysis2->Branch("rms_peak", &rms_peak, "rms_peak/D"); // rms of average of pulse peak
  FCAL_Analysis2->Branch("rms_integral",&rms_integral, "rms_integral/D");   // rms of average of pulse integral
  FCAL_Analysis2->Branch("rms_peak_gaus", &rms_peak_gaus, "rms_peak_gaus/D"); // rms of average of pulse peak
  FCAL_Analysis2->Branch("rms_integral_gaus",&rms_integral_gaus, "rms_integral_gaus/D");   // rms of average of pulse integral
  FCAL_Analysis2->Branch("eventnum", &eventnum, "eventnum/i"); // event number
  FCAL_Analysis2->Branch("run",      &run,      "run/i");      // run number

  FCAL_Correlation = new TTree("FCAL_Correlation","Pulse correlations");
  FCAL_Correlation->Branch("pulse1_rocid",    &pulse1_rocid   , "pulse1_rocid/i");    // vme crate number
  FCAL_Correlation->Branch("pulse1_slot",     &pulse1_slot    , "pulse1_slot/i");     // board number
  FCAL_Correlation->Branch("pulse1_channel",  &pulse1_channel , "pulse1_channel/i");  // board channel number
  FCAL_Correlation->Branch("pulse1_x",        &pulse1_x       , "pulse1_x/I");        // x coordinate where the center of (0,0)
  FCAL_Correlation->Branch("pulse1_y",        &pulse1_y       , "pulse1_y/I");        // y coordinate where the center of (0,0)
  FCAL_Correlation->Branch("pulse1_avgPeak",  &pulse1_avgPeak , "pulse1_avgPeak/D");        // pulse peak for pulse 1
  FCAL_Correlation->Branch("pulse2_rocid",    &pulse2_rocid   , "pulse2_rocid/i");    // vme crate number
  FCAL_Correlation->Branch("pulse2_slot",     &pulse2_slot    , "pulse2_slot/i");     // board number
  FCAL_Correlation->Branch("pulse2_channel",  &pulse2_channel , "pulse2_channel/i");  // board channel number
  FCAL_Correlation->Branch("pulse2_x",        &pulse2_x       , "pulse2_x/I");        // x coordinate where the center of (0,0)
  FCAL_Correlation->Branch("pulse2_y",        &pulse2_y       , "pulse2_y/I");        // y coordinate where the center of (0,0)
  FCAL_Correlation->Branch("pulse2_avgPeak",  &pulse2_avgPeak , "pulse2_avgPeak/D");        // pulse peak for pulse 1
  FCAL_Correlation->Branch("nsa",    &nsa   , "nsa/i");          // number of samples used in integral after crossing time
  FCAL_Correlation->Branch("nsb",    &nsb   , "nsb/i");          // number of samples used in integral before crossing time
  FCAL_Correlation->Branch("threshold",  &threshold   , "threshold/i");  // threshold
  FCAL_Correlation->Branch("corr", &corr, "corr/D"); // average of pulse peak over all events
  FCAL_Correlation->Branch("peak_ratio_mean", &peak_ratio_mean, "peak_ratio_mean/D"); // mean of pulse peak ratio
  FCAL_Correlation->Branch("peak_ratio_rms", &peak_ratio_rms, "peak_ratio_rms/D"); // rms of pulse peak ratio
  FCAL_Correlation->Branch("peak_ratio_rms_gaus", &peak_ratio_rms_gaus, "peak_ratio_rms_gaus/D"); // rms of pulse peak ratio
  FCAL_Correlation->Branch("eventnum", &eventnum, "eventnum/i"); // event number
  FCAL_Correlation->Branch("run",      &run,      "run/i");      // run number
  FCAL_Correlation->Branch("TC",      &TC,      "TC/i");      // time of threshold crossing
  
  FCAL_SingleCorrelation = new TTree("FCAL_SingleCorrelation","Single Pulse correlation");
  FCAL_SingleCorrelation->Branch("pulse1_rocid",    &pulse1_rocid   , "pulse1_rocid/i");    // vme crate number
  FCAL_SingleCorrelation->Branch("pulse1_slot",     &pulse1_slot    , "pulse1_slot/i");     // board number
  FCAL_SingleCorrelation->Branch("pulse1_channel",  &pulse1_channel , "pulse1_channel/i");  // board channel number
  FCAL_SingleCorrelation->Branch("pulse1_x",        &pulse1_x       , "pulse1_x/I");        // x coordinate where the center of (0,0)
  FCAL_SingleCorrelation->Branch("pulse1_y",        &pulse1_y       , "pulse1_y/I");        // y coordinate where the center of (0,0)
  FCAL_SingleCorrelation->Branch("pulse1_peak",     &pulse1_peak    , "pulse1_peak/D");        // pulse1 peak
  FCAL_SingleCorrelation->Branch("pulse1_integral", &pulse1_integral, "pulse1_integral/D");        // pulse1 integral
  FCAL_SingleCorrelation->Branch("pulse2_rocid",    &pulse2_rocid   , "pulse2_rocid/i");    // vme crate number
  FCAL_SingleCorrelation->Branch("pulse2_slot",     &pulse2_slot    , "pulse2_slot/i");     // board number
  FCAL_SingleCorrelation->Branch("pulse2_channel",  &pulse2_channel , "pulse2_channel/i");  // board channel number
  FCAL_SingleCorrelation->Branch("pulse2_x",        &pulse2_x       , "pulse2_x/I");        // x coordinate where the center of (0,0)
  FCAL_SingleCorrelation->Branch("pulse2_y",        &pulse2_y       , "pulse2_y/I");        // y coordinate where the center of (0,0)
  FCAL_SingleCorrelation->Branch("pulse2_peak",     &pulse2_peak    , "pulse2_peak/D");        // pulse2 peak
  FCAL_SingleCorrelation->Branch("pulse2_integral", &pulse2_integral, "pulse2_integral/D");        // pulse2 integral
  FCAL_SingleCorrelation->Branch("nsa",    &nsa   , "nsa/i");          // number of samples used in integral after crossing time
  FCAL_SingleCorrelation->Branch("nsb",    &nsb   , "nsb/i");          // number of samples used in integral before crossing time
  FCAL_SingleCorrelation->Branch("threshold",  &threshold   , "threshold/i");  // threshold
  FCAL_SingleCorrelation->Branch("corr", &corr, "corr/D"); // average of pulse peak over all events
  FCAL_SingleCorrelation->Branch("eventnum", &eventnum, "eventnum/i"); // event number
  FCAL_SingleCorrelation->Branch("run",      &run,      "run/i");      // run number
  FCAL_SingleCorrelation->Branch("TC",      &TC,      "TC/i");      // time of threshold crossing
  
  FCAL_Waveforms = new TTree("FCAL_Waveforms","Waveforms");
  FCAL_Waveforms->Branch("rocid",    &rocid   , "rocid/i");    // vme crate number
  FCAL_Waveforms->Branch("slot",     &slot    , "slot/i");     // board number
  FCAL_Waveforms->Branch("channel",  &channel , "channel/i");  // board channel number
  FCAL_Waveforms->Branch("x",        &x   ,     "x/I");        // x coordinate where the center of (0,0)
  FCAL_Waveforms->Branch("y",        &y   ,     "y/I");        // y coordinate where the center of (0,0)
  FCAL_Waveforms->Branch("sample",    &sample   , "sample/i");          // sample number
  FCAL_Waveforms->Branch("adc_val",    &adc_val   , "adc_val/i");          // adc value
  FCAL_Waveforms->Branch("eventnum", &eventnum, "eventnum/i"); // event number
  FCAL_Waveforms->Branch("run",      &run,      "run/i");      // run number  

  // Probably will want to add info pertaining to f250 algo performance

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FCALPpPi::brun(JEventLoop *eventLoop, int runnumber)
{
	// This is called whenever the run number changes
	// FreeClear();
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FCALPpPi::evnt(JEventLoop *loop, int eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	
	vector<const Df250WindowRawData*> f250WindowRawData_vec;
  loop->Get(f250WindowRawData_vec);
  sort(f250WindowRawData_vec.begin(), f250WindowRawData_vec.end(), Df250WindowRawData_cmp);
	
	vector<const Df250PulsePedestal*> f250PulsePedestal_vec;
  loop->Get(f250PulsePedestal_vec);
  sort(f250PulsePedestal_vec.begin(), f250PulsePedestal_vec.end(), Df250PulsePedestal_cmp);

  bool get_windowraw_pedestals = f250WindowRawData_vec.size()>0 ? true : false;
  
  if (get_windowraw_pedestals){
    for (int i=0; i<int(f250WindowRawData_vec.size()); i++){
      
      int ROCID   = f250WindowRawData_vec[i]->rocid; 
      int SLOT    = f250WindowRawData_vec[i]->slot; 
      int CHANNEL = f250WindowRawData_vec[i]->channel;
      
      if (ROCID>=11 && ROCID<=23){
        japp->RootWriteLock();
        
        const vector<uint16_t> window = f250WindowRawData_vec[i]->samples;
        int nwindow = window.size();
        
        char name[50];
        sprintf(name,"%02i/%02i/%02i",ROCID,SLOT,CHANNEL);
        
        int m_x = 0;
        int m_y = 0;
        if (tranlation_map.count(name)>0){
          m_x = tranlation_map[TString(name)].first;
          m_y = tranlation_map[TString(name)].second;
        }
        
        // Channel info
        rocid = ROCID;
        slot = SLOT;
        channel = CHANNEL;
        x = m_x;
        y = m_y;
        eventnum = eventnumber;
        run = loop->GetJEvent().GetRunNumber();
        
//********************
// Calculate Pulse Integral using nsa and nsb
//********************
        // find sample that crosses threshold
        int m_TC = 0;
        for (int i=0; i<nwindow; i++){
          //cout << "window: " << window[i] << "\t threshold: " << threshold << endl;
          sample = i;
          adc_val = window[i];
          FCAL_Waveforms->Fill();
          if (window[i]>=threshold && m_TC==0 && window[i]<4096){
            if (i==0) m_TC = 1;
            else m_TC = i;
          }
        }
        if (m_TC==0) {
          japp->RootUnLock();
          continue; // skip channels with no pulse
        }
        for (int m_nsa=min_nsa; m_nsa<nwindow; m_nsa++){
          
          int m_start = (m_TC-nsb)<0 ? 0 : (m_TC-nsb);
          int m_end   = (m_TC+m_nsa)>nwindow ? nwindow+1 : (m_TC+m_nsa);
          if (m_end==nwindow+1) continue;
          
          int m_integral = 0;
          int m_peak = 0;
          for (int i=m_start; i<m_end; i++){
            m_integral += window[i];
            if (window[i]>m_peak) m_peak = window[i];
          }
          // Fill map when varying nsa
          TString key = TString(name)+"_"+StringUtilities::int2TString(m_nsa);
          pulse_ana[key].addPulsePeak(m_peak, eventnumber);
          pulse_ana[key].addPulseIntegral(m_integral);
          pulse_ana[key].setX(m_x);
          pulse_ana[key].setY(m_y);
          pulse_ana[key].setStart(m_start);
          pulse_ana[key].setEnd(m_end);
          pulse_ana[key].setTC(m_TC);
          pulse_ana[key].setEventnum(eventnumber);
          pulse_ana[key].setRun(run);
          
          // Write to TTree
          TC = m_TC;
          nsa = m_nsa;
          peak = m_peak;
          integral = m_integral;
          start = m_start;
          end = m_end;
          FCAL_Pulse->Fill();
        }
        
//********************
// Calculate Pulse Integral using nsp 
//(the number of samples around the peak)
//********************
        int nsp_sample = 0;
        int nsp_peak = 0;
        for (int i=0; i<nwindow; i++){
          if (window[i]>nsp_peak) {
            nsp_peak = window[i];
            nsp_sample = i;
          }
        }
        for (int m_nsp=min_nsp; m_nsp<=max_nsp; m_nsp++){
          int m_integral = 0;
          int m_start = nsp_sample-m_nsp<0 ? 0 : nsp_sample-m_nsp;
          int m_end   = nsp_sample+m_nsp+1>nwindow ? nwindow : nsp_sample+m_nsp+1;
          for (int i=m_start; i<m_end; i++){
            m_integral += window[i];
          }
          // Fill map when varying nsp
          TString key = TString(name)+"_"+StringUtilities::int2TString(m_nsp);
          pulse_nsp[key].addPulsePeak(nsp_peak, eventnumber);
          pulse_nsp[key].addPulseIntegral(m_integral);
          pulse_nsp[key].setX(m_x);
          pulse_nsp[key].setY(m_y);
          pulse_nsp[key].setStart(m_start);
          pulse_nsp[key].setEnd(m_end);
          pulse_nsp[key].setTC(m_TC);
          pulse_nsp[key].setEventnum(eventnumber);
          pulse_nsp[key].setRun(run); 
          
          // Write to TTree
          nsp = m_nsp;
          start = m_start;
          end   = m_end;
          peak = nsp_peak;
          integral = m_integral;
          FCAL_Pulse2->Fill();
        }
        
        japp->RootUnLock();
      }// rocid>=11 && rocid<=23 conditional
    } // End loop over f250WindowRawData_vec
  } // End get_windowraw_pedestals conditional

  // Additional Conditional, f250PulsePedestal
  
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FCALPpPi::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FCALPpPi::fini(void)
{
	// Called before program exit after event processing is finished.
	
	japp->RootWriteLock();
  
//********************
// Using NSA and NSB Get 
// Pulse Integral and Peak rms and mean
//********************
  cout << "Doing NSA + NSB analaysis..." << endl;
  for (std::map<TString, Pulse>::const_iterator iter=pulse_ana.begin(); iter!=pulse_ana.end(); ++iter){
    TString key = iter->first;
    cout << key << endl;
    
    Pulse::Pulse pulse = pulse_ana[key];
    
    if (!pulse.peakExists() || !pulse.integralExists()) continue;
    
    rocid = pulse.getRocid();
    slot = pulse.getSlot();
    channel = pulse.getChannel();
    
    //if ( rocid>11 || slot>5 || channel>10 ) continue;
    
    x = pulse.getX();
    y = pulse.getY();
    TC = pulse.getTC();
    nsa = pulse.getNSA();
    start = pulse.getStart();
    end = pulse.getEnd();
    avg_peak = pulse.getAvgPulsePeak();
    avg_integral = pulse.getAvgPulseIntegral();
    rms_peak = pulse.getRmsPulsePeak();
    rms_integral = pulse.getRmsPulseIntegral();
    rms_peak_gaus = pulse.getRmsPulsePeakFit();
    rms_integral_gaus = pulse.getRmsPulseIntegralFit();
    eventnum = pulse.getEventnum();
    run = pulse.getRun();
    
    vector<double> pp1_pp2_cor;
    if (nsa==15){
      pulse1_rocid = rocid;
      pulse1_slot = slot;
      pulse1_channel = channel;
      pulse1_x = x;
      pulse1_y = y;
      pulse1_avgPeak = pulse.getAvgPulsePeak();
      
      // Get the correlation of this channel with every other channel, for a fixed nsa
      for (std::map<TString, Pulse>::const_iterator iter2=pulse_ana.begin(); iter2!=pulse_ana.end(); ++iter2){
        TString key2 = iter2->first;
        Pulse::Pulse pulse2 = pulse_ana[key2];

        if (!pulse2.peakExists() || !pulse2.integralExists()) continue;
        
        
        uint32_t p2_rocid = pulse2.getRocid();
        uint32_t p2_slot = pulse2.getSlot();
        uint32_t p2_channel = pulse2.getChannel();
        if (rocid==p2_rocid && slot==p2_slot && channel==p2_channel) continue;
        if (nsa!=pulse2.getNSA()) continue;
        
        // avoid double counting correlation
        pair<int,int> coord1 = make_pair(pulse1_x, pulse1_y);
        pair<int,int> coord2 = make_pair(pulse2.getX(), pulse2.getY());
        if (corrmap[ make_pair(coord1, coord2) ]>0 || corrmap[ make_pair(coord2, coord1) ]>0) continue;
        corrmap[ make_pair(coord1, coord2) ]++;
        corrmap[ make_pair(coord2, coord2) ]++;
        
        pulse2_rocid = p2_rocid;
        pulse2_slot = p2_slot;
        pulse2_channel = p2_channel;
        pulse2_x = pulse2.getX();
        pulse2_y = pulse2.getY();
        pulse2_avgPeak = pulse2.getAvgPulsePeak();
        
        map<int,double> p2_map_peak = pulse2.getMapPulsePeak();
        
        vector<double> vec_ratio; 
        if (pulse1_x==3 && pulse1_y==3 && pulse2_x==4 && pulse2_y==3) vec_ratio = pulse.calcPulsePeakRatio( p2_map_peak, true );
        else if (pulse1_x==7 && pulse1_y==13 && pulse2_x==8 && pulse2_y==13) vec_ratio = pulse.calcPulsePeakRatio( p2_map_peak, true );
        else if (pulse1_x==4 && pulse1_y==14 && pulse2_x==5 && pulse2_y==14) vec_ratio = pulse.calcPulsePeakRatio( p2_map_peak, true );
        else vec_ratio = pulse.calcPulsePeakRatio( p2_map_peak );
        peak_ratio_mean = vec_ratio[0];
        peak_ratio_rms  = vec_ratio[1];
        peak_ratio_rms_gaus  = vec_ratio[2];
        
        corr = pulse.calcPulsePeakCorrelation( p2_map_peak );
        if (int(corr)!=101) {
          FCAL_Correlation->Fill();
        }
        
        // These channels fill to the single correlation ttree
        if (nsa==15 && rocid==11 && slot<5 && channel<3){
          // PulsePeak
          vector<double> vec_pulse1PP = pulse.getVecPulsePeak();
          vector<double> vec_pulse2PP = pulse2.getVecPulsePeak();
          // PulseIntegral
          vector<double> vec_pulse1PI = pulse.getVecPulseIntegral();
          vector<double> vec_pulse2PI = pulse2.getVecPulseIntegral();
          // Event numbers
          vector<int> vec_eventnumber = pulse.getVecEventnum();
          int n_eventnumber = vec_eventnumber.size();
        
          int n_pulse1PP = vec_pulse1PP.size();
          int n_pulse2PP = vec_pulse2PP.size();
          int n_pulse1PI = vec_pulse1PI.size();
          int n_pulse2PI = vec_pulse2PI.size();
          if (n_pulse1PP==0 || n_pulse2PP==0 || n_pulse1PI==0 || n_pulse2PI==0) continue;
          int n_pulsePP = n_pulse1PP>n_pulse2PP ? n_pulse2PP : n_pulse1PP;
          int n_pulsePI = n_pulse1PI>n_pulse2PI ? n_pulse2PI : n_pulse1PI;
          int n_pulseP  = n_pulsePP>n_pulsePI ? n_pulsePI : n_pulsePP;
          for (int i=0; i<n_pulseP; i++){
            if (i<n_eventnumber) eventnum = vec_eventnumber[i];
            else eventnum = 0;
            pulse1_peak = vec_pulse1PP[i];
            pulse2_peak = vec_pulse2PP[i];
            pulse1_integral = vec_pulse1PI[i];
            pulse2_integral = vec_pulse2PI[i];
            FCAL_SingleCorrelation->Fill();
          }
        }
        
      }
    }// end nsa == 15 conditional
    
    FCAL_Analysis->Fill();
    // reset
    
  }
	pulse_ana.clear();
	
	// *******
	// Write correlation map to text file for debugging
	fstream corrOut("CorrDebug.txt",fstream::out);
	for (map< pair< pair<int,int>, pair<int,int> >, int >::const_iterator iter=corrmap.begin(); iter!=corrmap.end(); ++iter){
    pair<int,int> m_pulse1 = iter->first.first;
    pair<int,int> m_pulse2 = iter->first.second;
    int exists = iter->second;
    char debug[600];
    sprintf(debug,"Pulse1: (%i,%i), Pulse2: (%i,%i), bool: %i",m_pulse1.first,m_pulse1.second,m_pulse2.first,m_pulse2.second,exists);
    corrOut << debug << endl;
  }
  corrOut.close();
	
//********************
// Using NSP Get 
// Pulse Integral and Peak rms and mean
//********************
  cout << "Doing NSP analaysis..." << endl;
  for (std::map<TString, Pulse>::const_iterator iter=pulse_nsp.begin(); iter!=pulse_nsp.end(); ++iter){
    TString key = iter->first;
    cout << key << endl;
    
    Pulse::Pulse pulse = pulse_nsp[key];
    
    if (!pulse.peakExists() || !pulse.integralExists()) continue;
    
    rocid = pulse.getRocid();
    slot = pulse.getSlot();
    channel = pulse.getChannel();
    x = pulse.getX();
    y = pulse.getY();
    start = pulse.getStart();
    end = pulse.getEnd();
    nsp = pulse.getNSP();
    avg_peak = pulse.getAvgPulsePeak();
    avg_integral = pulse.getAvgPulseIntegral();
    rms_peak = pulse.getRmsPulsePeak();
    rms_integral = pulse.getRmsPulseIntegral();
    rms_peak_gaus = pulse.getRmsPulsePeakFit();
    rms_integral_gaus = pulse.getRmsPulseIntegralFit();
    eventnum = pulse.getEventnum();
    run = pulse.getRun();
    FCAL_Analysis2->Fill();
  }
	
  japp->RootUnLock();
	
	
	return NOERROR;
}
