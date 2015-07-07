additionally to correlation look at 
pulse peak_i / pulse peak_j for every event (would expect about 1.5%) and plot it versus pulse peak_i

Using real beam data look at pulse integral as a function of nsa

Look at pedestal subtracted pulse peak to get an idea of the uncertainty of peak itself

Fold in the uncertainty introduced as NSA gets bigger due to the uncertainty of the pedestal

rms probably should be sqrt(rms_integral**2 + (rms_ped*nsa)**2)
The error as nsa gets bigger is artificially small without this contribution


# Log into ifarm65

# cd halld_workspace

# Setup env variables

source /w/halld-scifs1a/home/manlara/springgluex_setup.csh

# To compile this plugin do:

scons -u install

# To run this plugin over an evio file do

hd_root -PPLUGINS=FCALPedestalPlugin -PTT:SYSTEMS_TO_PARSE=FCAL hd_rawdata_003301_000.evio -o test.root
hd_root -PPLUGINS=FCALPpPi -PTT:SYSTEMS_TO_PARSE=FCAL fcal_n6_132.0 -o test.root

## Plugin information

The FCALPedestalPlugin reads evio and outputs a root file. The root file contains the following variables:
FCAL_Pedestal = new TTree("FCAL_Pedestal","Pedestal information");
FCAL_Pedestal->Branch("rocid",    &rocid   , "rocid/i");    // vme crate number
FCAL_Pedestal->Branch("slot",     &slot    , "slot/i");     // board number
FCAL_Pedestal->Branch("channel",  &channel , "channel/i");  // board channel number
FCAL_Pedestal->Branch("pedestal", &pedestal, "pedestal/D"); // pedestal calculated by averaging nEvents
FCAL_Pedestal->Branch("rms",      &rms,      "rms/D");      // root-mean-square of pedestal average
FCAL_Pedestal->Branch("nSamples", &nSamples, "nSamples/i"); // the number of samples used in pedestal average
FCAL_Pedestal->Branch("nEvents", &nEvents, "nEvents/i");    // the number of events to be averaged in the pedestal calc
FCAL_Pedestal->Branch("eventnum", &eventnum, "eventnum/i"); // event number
FCAL_Pedestal->Branch("run",      &run,      "run/i");      // run number

#! Important

Cutting on nEvents and nSamples is important for data analysis. nEvents goes from 1-20 and nSamples goes from 1-20

## Pedestal calculation

The pedestal is calculated by taking the mean of nSamples from the beginning of the channels waveform PLUS nEvents. 
If nEvents = 1, then the pedestal, p_i = Sum_j=1^nSamples q_i,j / nSamples
otherwise, p_k = Sum_i=1^nEvents p_i / nSamples / nEvents


## Debugging the FCAL Controls

When you turn on or off the FCAL strands some channels may not respond to your command. Power cycling the bases on the bus should be your first attempt. Rebooting the IOC should be your second



FCAL_TestCorrelation->Draw("pulse1_peak:pulse2_peak","pulse1_rocid==11&&pulse1_slot==3&&pulse1_channel==2&&pulse2_rocid==11&&pulse2_slot==3&&pulse2_channel==3")
teh number of entries represents the number of events