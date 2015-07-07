# Log into ifarm65

# cd halld_workspace

# Setup env variables

source /w/halld-scifs1a/home/manlara/springgluex_setup.csh

# To compile this plugin do:

scons -u install

# To run this plugin over an evio file do

hd_root -PPLUGINS=FCALPpPi -PTT:SYSTEMS_TO_PARSE=FCAL fcal_n6_132.0 -o test.root

## Plugin information

### TTrees
##### FCAL_Pulse = new TTree("FCAL_Pulse","Pulse information");
####  Use nsa + nsb for integral
##### FCAL_Analysis = new TTree("FCAL_Analysis","Pulse analysis");

// using the samples around the peak for the integral

FCAL_Analysis2 = new TTree("FCAL_Analysis2","Pulse analysis 2");


FCAL_Correlation = new TTree("FCAL_Correlation","Pulse correlations");


FCAL_SingleCorrelation = new TTree("FCAL_SingleCorrelation","Single Pulse correlation");


FCAL_Waveforms = new TTree("FCAL_Waveforms","Waveforms");


#! Important

Cutting on nEvents and nSamples is important for data analysis. nEvents goes from 1-20 and nSamples goes from 1-20

## Pedestal calculation

The pedestal is calculated by taking the mean of nSamples from the beginning of the channels waveform PLUS nEvents. 
If nEvents = 1, then the pedestal, p_i = Sum_j=1^nSamples q_i,j / nSamples
otherwise, p_k = Sum_i=1^nEvents p_i / nSamples / nEvents

