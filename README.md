# Log into ifarm65

# cd halld_workspace

# Setup env variables

source /w/halld-scifs1a/home/manlara/springgluex_setup.csh

# To compile this plugin do:

scons -u install

# To run this plugin over an evio file do

hd_root -PPLUGINS=FCALPpPi -PTT:SYSTEMS_TO_PARSE=FCAL fcal_n6_132.0 -o test.root

## Plugin information

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

#! Important

Cutting on nEvents and nSamples is important for data analysis. nEvents goes from 1-20 and nSamples goes from 1-20

## Pedestal calculation

The pedestal is calculated by taking the mean of nSamples from the beginning of the channels waveform PLUS nEvents. 
If nEvents = 1, then the pedestal, p_i = Sum_j=1^nSamples q_i,j / nSamples
otherwise, p_k = Sum_i=1^nEvents p_i / nSamples / nEvents

