#include <inttypes.h>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <sstream>

#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TVirtualPad.h"
#include "TFrame.h"
#include "TGraphPainter.h"
#include "TString.h"
#include "TCollection.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TKey.h"
#include "TDatime.h"
#include "TMath.h"
#include "TAxis.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"
#include "TColor.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"


#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>

#include "StringUtilities.h"

const static int min_nsp = 0;
const static int max_nsp = 10;
const static int maxNSA = 27;
const static int minNSA = 2;

uint32_t rocid;
uint32_t slot;
uint32_t channel;
int m_x;
int m_y;
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
double peak_ratio_rms_gaus ;

double pulse1_peak      ;
double pulse1_integral  ;
double pulse2_peak      ;
double pulse2_integral  ;;
double pulse1_avgPeak;

uint32_t sample;
uint32_t adc_val;

const static uint32_t maxADC = 4096;

void SetFCAL_Analysis(TTree* tr){
  tr->SetBranchAddress("rocid",    &rocid  );    // vme crate number
  tr->SetBranchAddress("slot",     &slot   );     // board number
  tr->SetBranchAddress("channel",  &channel);  // board channel number
  tr->SetBranchAddress("x",        &m_x      );        // x coordinate where the center of (0,0)
  tr->SetBranchAddress("y",        &m_y      );        // y coordinate where the center of (0,0)
  tr->SetBranchAddress("nsa",    &nsa      );          // number of samples used in integral after crossing time
  tr->SetBranchAddress("nsb",    &nsb      );          // number of samples used in integral before crossing time
  tr->SetBranchAddress("start",    &start      );          // start sample
  tr->SetBranchAddress("end",    &end      );          // end sample
  tr->SetBranchAddress("threshold",  &threshold   );  // threshold
  tr->SetBranchAddress("avg_peak", &avg_peak); // average of pulse peak over all events
  tr->SetBranchAddress("avg_integral",  &avg_integral);      // average of pulse integral over all events
  tr->SetBranchAddress("rms_peak", &rms_peak); // rms of average of pulse peak
  tr->SetBranchAddress("rms_integral",&rms_integral);   // rms of average of pulse integral
  tr->SetBranchAddress("rms_peak_gaus", &rms_peak_gaus); // rms of average of pulse peak
  tr->SetBranchAddress("rms_integral_gaus",&rms_integral_gaus);   // rms of average of pulse integral
  tr->SetBranchAddress("eventnum", &eventnum); // event number
  tr->SetBranchAddress("run",      &run);      // run number
  tr->Branch("TC",      &TC,      "TC/i");      // time of threshold crossing
}

void SetFCAL_Analysis2(TTree* tr){
  tr->SetBranchAddress("rocid",    &rocid  );    // vme crate number
  tr->SetBranchAddress("slot",     &slot   );     // board number
  tr->SetBranchAddress("channel",  &channel);  // board channel number
  tr->SetBranchAddress("x",        &m_x      );        // x coordinate where the center of (0,0)
  tr->SetBranchAddress("y",        &m_y      );        // y coordinate where the center of (0,0)
  tr->SetBranchAddress("nsp",    &nsp      );          // number of samples used in integral around the pulse peak
  tr->SetBranchAddress("start",    &start      );          // start sample
  tr->SetBranchAddress("end",    &end      );          // end sample
  tr->SetBranchAddress("threshold",  &threshold   );  // threshold
  tr->SetBranchAddress("avg_peak", &avg_peak); // average of pulse peak over all events
  tr->SetBranchAddress("avg_integral",  &avg_integral);      // average of pulse integral over all events
  tr->SetBranchAddress("rms_peak", &rms_peak); // rms of average of pulse peak
  tr->SetBranchAddress("rms_integral",&rms_integral);   // rms of average of pulse integral
  tr->SetBranchAddress("rms_peak_gaus", &rms_peak_gaus); // rms of average of pulse peak
  tr->SetBranchAddress("rms_integral_gaus",&rms_integral_gaus);   // rms of average of pulse integral
  tr->SetBranchAddress("eventnum", &eventnum); // event number
  tr->SetBranchAddress("run",      &run);      // run number
}

void SetFCAL_Correlation(TTree* tr){
  tr->SetBranchAddress("pulse1_rocid",    &pulse1_rocid );    // vme crate number
  tr->SetBranchAddress("pulse1_slot",     &pulse1_slot  );     // board number
  tr->SetBranchAddress("pulse1_channel",  &pulse1_channel);  // board channel number
  tr->SetBranchAddress("pulse1_x",        &pulse1_x      );        // x coordinate where the center of (0,0)
  tr->SetBranchAddress("pulse1_y",        &pulse1_y      );        // y coordinate where the center of (0,0)
  tr->SetBranchAddress("pulse2_rocid",    &pulse2_rocid  );    // vme crate number
  tr->SetBranchAddress("pulse2_slot",     &pulse2_slot   );     // board number
  tr->SetBranchAddress("pulse2_channel",  &pulse2_channel);  // board channel number
  tr->SetBranchAddress("pulse2_x",        &pulse2_x      );        // x coordinate where the center of (0,0)
  tr->SetBranchAddress("pulse2_y",        &pulse2_y      );        // y coordinate where the center of (0,0)
  tr->SetBranchAddress("nsa",    &nsa   );          // number of samples used in integral after crossing time
  tr->SetBranchAddress("nsb",    &nsb   );          // number of samples used in integral before crossing time
  tr->SetBranchAddress("threshold",  &threshold);  // threshold
  tr->SetBranchAddress("corr", &corr); // average of pulse peak over all events
  tr->SetBranchAddress("peak_ratio_mean", &peak_ratio_mean); // mean of pulse peak ratio
  tr->SetBranchAddress("peak_ratio_rms",  &peak_ratio_rms); // rms of pulse peak ratio
  tr->SetBranchAddress("peak_ratio_rms_gaus",  &peak_ratio_rms_gaus); // rms of pulse peak ratio
  tr->SetBranchAddress("eventnum", &eventnum); // event number
  tr->SetBranchAddress("run",      &run);      // run number
  tr->Branch("TC",      &TC,      "TC/i");      // time of threshold crossing
}

void SetFCAL_Waveforms(TTree* tr){
  tr->SetBranchAddress("rocid",    &rocid   );    // vme crate number
  tr->SetBranchAddress("slot",     &slot    );     // board number
  tr->SetBranchAddress("channel",  &channel );  // board channel number
  tr->SetBranchAddress("x",        &m_x   );        // x coordinate where the center of (0,0)
  tr->SetBranchAddress("y",        &m_y   );        // y coordinate where the center of (0,0)
  tr->SetBranchAddress("sample",    &sample);          // sample number
  tr->SetBranchAddress("adc_val",    &adc_val);          // adc value
  tr->SetBranchAddress("eventnum", &eventnum); // event number
  tr->SetBranchAddress("run",      &run);      // run number
}

// Looking at Run 2615
double Weight( double pulsepeak ){
  return TMath::Exp(-2.68073e-02*pulsepeak);
}

void PulseAnalysis(TString filename, bool supportingPlots=false, bool save=false){
  TString dir = "plots/";
  TFile * f = new TFile(filename);
  
  TTree* FCAL_Analysis    = (TTree*)f->Get("FCAL_Analysis");
  SetFCAL_Analysis(FCAL_Analysis);
  
  // nsb is fixed to 5 and the threshold is the same as the Spring run 180 counts above pedestal
  // Histogram: 
  // avg_peak, rms_peak, rms_peak/avg_peak*100
  // rms_integral/avg_integral*100 for each nsa 0-30
  // Graph
  // average of (avg_integral, rms_integral, rms_integral/avg_integral*100) over all channels VS nsa
  
  // Array of avg_peak, rms_peak, error_peak, error_ped_subtracted
  TH1F* h_peak_NoWeight[4];
  TH1F* h_peak_Weight[4];
  
  h_peak_NoWeight[0] = new TH1F("h_nw_avg_peak" ,"Average pulse peak", maxADC/10, 0, maxADC);           h_peak_NoWeight[0]->SetXTitle("Peak");
  h_peak_NoWeight[1] = new TH1F("h_nw_rms_peak" ,"RMS of average pulse peak", maxADC/0.01, 0, maxADC);  h_peak_NoWeight[1]->SetXTitle("RMS");
  h_peak_NoWeight[2] = new TH1F("h_nw_error_peak","Error on pulse peak", 100/0.1, 0, 100);              h_peak_NoWeight[2]->SetXTitle("Error (%)");
  h_peak_NoWeight[3] = new TH1F("h_nw_error_ped_subtracted_peak","Error on pedestal subtracted pulse peak", 100/0.1, 0, 100); h_peak_NoWeight[3]->SetXTitle("Error (%)");
  
  h_peak_Weight[0] = new TH1F("h_w_avg_peak" ,"Average pulse peak (weight)", maxADC/10, 0, maxADC);           h_peak_Weight[0]->SetXTitle("Peak");
  h_peak_Weight[1] = new TH1F("h_w_rms_peak" ,"RMS of average pulse peak (weight)", maxADC/0.01, 0, maxADC);  h_peak_Weight[1]->SetXTitle("RMS");
  h_peak_Weight[2] = new TH1F("h_w_error_peak","Error on pulse peak (weight)", 100/0.1, 0, 100);              h_peak_Weight[2]->SetXTitle("Error (%)");
  h_peak_Weight[3] = new TH1F("h_w_error_ped_subtracted_peak","Error on pedestal subtracted pulse peak (weight)", 100/0.1, 0, 100); h_peak_Weight[3]->SetXTitle("Error (%)");
  
  TH1F* h_integral_NoWeight[4][maxNSA];
  TH1F* h_integral_Weight[4][maxNSA];
  
  for (int i=0; i<maxNSA; i++){
    TString t_nsa = TString::Itoa(i,10);
    h_integral_NoWeight[0][i] = new TH1F("h_nw_avg_integral_"+t_nsa,"Pulse integral mean (nsa="+t_nsa+")", maxADC*100/(10*(i+1)), 0, maxADC*100);
    h_integral_NoWeight[1][i] = new TH1F("h_nw_rms_integral_"+t_nsa,"Pulse integral rms (nsa="+t_nsa+")", maxADC*100/10, 0, maxADC*100);
    h_integral_NoWeight[2][i] = new TH1F("h_nw_error_integral_"+t_nsa,"Pulse integral error (nsa="+t_nsa+")", 100/0.1, 0, 100);    
    h_integral_NoWeight[3][i] = new TH1F("h_nw_error_ped_subtracted_integral_"+t_nsa,"Pedestal subtracted integral error (nsa="+t_nsa+")", 100/0.1, 0, 100);
    
    h_integral_Weight[0][i] = new TH1F("h_w_avg_integral_"+t_nsa,"Pulse integral mean (nsa="+t_nsa+", weight)", maxADC*100/(10*(i+1)), 0, maxADC*100);
    h_integral_Weight[1][i] = new TH1F("h_w_rms_integral_"+t_nsa,"Pulse integral rms (nsa="+t_nsa+", weight)", maxADC*100/10, 0, maxADC*100);
    h_integral_Weight[2][i] = new TH1F("h_w_error_integral_"+t_nsa,"Pulse integral error (nsa="+t_nsa+", weight)", 100/0.1, 0, 100);    
    h_integral_Weight[3][i] = new TH1F("h_w_error_ped_subtracted_integral_"+t_nsa,"Pedestal subtracted integral error (nsa="+t_nsa+", weight)", 100/0.1, 0, 100);
  }
  
  // integral vs nsa
  TGraphErrors* gr_integral_vs_nsa[4];
  for (int i=0; i<4; i++){
    gr_integral_vs_nsa[i] = new TGraphErrors();
    if (i==0){
      gr_integral_vs_nsa[i]->SetName("integral_vs_nsa");
      gr_integral_vs_nsa[i]->SetTitle("Average pulse integral error vs nsa");
    }
    if (i==1){
      gr_integral_vs_nsa[i]->SetName("ped_subtracted_integral_vs_nsa");
      gr_integral_vs_nsa[i]->SetTitle("Pedestal subtracted pulse integral error vs nsa");
    }
    if (i==2){
      gr_integral_vs_nsa[i]->SetName("integral_vs_nsa_weight");
      gr_integral_vs_nsa[i]->SetTitle("Average pulse integral error vs nsa (weighted)");
    }
    if (i==3){
      gr_integral_vs_nsa[i]->SetName("ped_subtracted_integral_vs_nsa_weight");
      gr_integral_vs_nsa[i]->SetTitle("Pedestal subtracted pulse integral error vs nsa (weighted)");
    }
  }
  
  // error vs peak for pulse peak and integral for fixed nsa
  TGraph* gr_error_vs_peak[maxNSA+1];
  for (int i=0; i<maxNSA+1; i++){
    gr_error_vs_peak[i] = new TGraph();
    if (i==0) {
      gr_error_vs_peak[i]->SetName("peakError_vs_peak");
      gr_error_vs_peak[i]->SetTitle("Error on pulse peak vs pulse peak");
    }
    else {
      TString t_nsa = TString::Itoa(i-1,10);
      gr_error_vs_peak[i]->SetName("integralError_vs_peak_"+t_nsa);
      gr_error_vs_peak[i]->SetTitle("Error on pulse integral vs pulse peak (nsa="+t_nsa+")");
    }
  }
  
  // useful for drawing
  double max_peak[4];
  double max_integral[4][maxNSA];
  
  double min_peak[4];
  double min_integral[4][maxNSA];
  
  for (int i=0; i<4; i++){
    max_peak[i] = 0;
    min_peak[i] = maxADC;
    for (int j=0; j<maxNSA; j++) {
      max_integral[i][j] = 0;
      min_integral[i][j] = maxADC*100;
    }
  }
  
  int nEntries = FCAL_Analysis->GetEntries();
  for (int i=0; i<nEntries; i++){
    FCAL_Analysis->GetEntry(i);
    
    // The larger NSA the bigger effect from the pedestal
    // Take the rms of the pedestal to be 0.5 counts
    double ped_subtracted_integral = avg_integral-100.0*double(nsa+nsb);
    
    double error = rms_integral/avg_integral*100.0;
    double ped_subtracted_error    = TMath::Sqrt(TMath::Power(rms_integral/(avg_integral-100.0*double(nsa+nsb)),2)+TMath::Power(0.5/100.0*double(nsa+nsb),2)) * 100.0;
    
    //double error = rms_integral_gaus/avg_integral*100.0;
    //double ped_subtracted_error    = TMath::Sqrt(TMath::Power(rms_integral_gaus/(avg_integral-100.0*double(nsa+nsb)),2)+TMath::Power(0.5/100.0*double(nsa+nsb),2)) * 100.0;
    double weight = Weight(avg_peak);
    
    double peak_error = rms_peak/avg_peak*100.0;
    double peak_pedSubtracted_error = TMath::Sqrt(TMath::Power(rms_peak/(avg_peak-100),2)+TMath::Power(0.5/100.0,2)) * 100.0;
    //double peak_error = rms_peak_gaus/avg_peak*100.0;
    //double peak_pedSubtracted_error = TMath::Sqrt(TMath::Power(rms_peak_gaus/(avg_peak-100),2)+TMath::Power(0.5/100.0,2)) * 100.0;
    
    // Fill Pulse Peak Histograms
    if (nsa==15) {
      // Array of avg_peak, rms_peak, error_peak, error_ped_subtracted
      if (max_peak[0]<avg_peak)             max_peak[0] = avg_peak;
      if (max_peak[1]<rms_peak)             max_peak[1] = rms_peak;
      //if (max_peak[1]<rms_peak_gaus)             max_peak[1] = rms_peak_gaus;
      if (max_peak[2]<error)                max_peak[2] = error;
      if (max_peak[3]<ped_subtracted_error) max_peak[3] = ped_subtracted_error;
      
      if (min_peak[0]>avg_peak)             min_peak[0] = avg_peak;
      if (min_peak[1]>rms_peak)             min_peak[1] = rms_peak;
      //if (min_peak[1]>rms_peak_gaus)             min_peak[1] = rms_peak_gaus;
      if (min_peak[2]>error)                min_peak[2] = error;
      if (min_peak[3]>ped_subtracted_error) min_peak[3] = ped_subtracted_error;
      
      h_peak_NoWeight[0]->Fill(avg_peak);
      h_peak_NoWeight[1]->Fill(rms_peak);
      //h_peak_NoWeight[1]->Fill(rms_peak_gaus);
      h_peak_NoWeight[2]->Fill(error);
      h_peak_NoWeight[3]->Fill(ped_subtracted_error);
      
      h_peak_Weight[0]->Fill(avg_peak, weight);
      h_peak_Weight[1]->Fill(rms_peak, weight);
      //h_peak_Weight[1]->Fill(rms_peak_gaus, weight);
      h_peak_Weight[2]->Fill(peak_error, weight);
      h_peak_Weight[3]->Fill(peak_pedSubtracted_error, weight);
      
      gr_error_vs_peak[0]->SetPoint(gr_error_vs_peak[0]->GetN(), avg_peak, peak_error);
    }
    
    // Fill Pulse Integral Histograms
    if (nsa<maxNSA) {
      if (max_integral[0][nsa]<avg_integral)          max_integral[0][nsa] = avg_integral;
      if (max_integral[1][nsa]<rms_integral)          max_integral[1][nsa] = rms_integral;
      //if (max_integral[1][nsa]<rms_integral_gaus)          max_integral[1][nsa] = rms_integral_gaus;
      if (max_integral[2][nsa]<error)                 max_integral[2][nsa] = error;
      if (max_integral[3][nsa]<ped_subtracted_error)  max_integral[3][nsa] = ped_subtracted_error;
      
      if (min_integral[0][nsa]>avg_integral)          min_integral[0][nsa] = avg_integral;
      if (min_integral[1][nsa]>rms_integral)          min_integral[1][nsa] = rms_integral;
      //if (min_integral[1][nsa]>rms_integral_gaus)          min_integral[1][nsa] = rms_integral_gaus;
      if (min_integral[2][nsa]>error)                 min_integral[2][nsa] = error;
      if (min_integral[3][nsa]>ped_subtracted_error)  min_integral[3][nsa] = ped_subtracted_error;
      
      h_integral_NoWeight[0][nsa]->Fill(avg_integral);
      h_integral_NoWeight[1][nsa]->Fill(rms_integral);
      //h_integral_NoWeight[1][nsa]->Fill(rms_integral_gaus);
      h_integral_NoWeight[2][nsa]->Fill(error);
      h_integral_NoWeight[3][nsa]->Fill(ped_subtracted_error);
      
      h_integral_Weight[0][nsa]->Fill(avg_integral, weight);
      h_integral_Weight[1][nsa]->Fill(rms_integral, weight);
      //h_integral_Weight[1][nsa]->Fill(rms_integral_gaus, weight);
      h_integral_Weight[2][nsa]->Fill(error, weight);
      h_integral_Weight[3][nsa]->Fill(ped_subtracted_error, weight);
      
      if (nsa>0) gr_error_vs_peak[nsa]->SetPoint(gr_error_vs_peak[nsa]->GetN(), avg_peak, error);
    }
  }
  
  for (int i=minNSA; i<maxNSA; i++){
    TString t_nsa = TString::Itoa(i,10);
    
    gr_integral_vs_nsa[0]->SetPoint(i-minNSA,i,h_integral_NoWeight[2][i]->GetMean());
    gr_integral_vs_nsa[0]->SetPointError(i-minNSA,0,h_integral_NoWeight[2][i]->GetRMS());
    
    gr_integral_vs_nsa[1]->SetPoint(i-minNSA,i,h_integral_NoWeight[3][i]->GetMean());
    gr_integral_vs_nsa[1]->SetPointError(i-minNSA,0,h_integral_NoWeight[3][i]->GetRMS());
    
    gr_integral_vs_nsa[2]->SetPoint(i-minNSA,i,h_integral_Weight[2][i]->GetMean());
    gr_integral_vs_nsa[2]->SetPointError(i-minNSA,0,h_integral_Weight[2][i]->GetRMS());
    
    gr_integral_vs_nsa[3]->SetPoint(i-minNSA,i,h_integral_Weight[3][i]->GetMean());
    gr_integral_vs_nsa[3]->SetPointError(i-minNSA,0,h_integral_Weight[3][i]->GetRMS());
  }
  
  if (supportingPlots){
    // Pulse Peak and Integral Canvas'
    TCanvas* c_peak_NoWeight[4];
    TCanvas* c_peak_Weight[4];
    TCanvas* c_integral_NoWeight[4][maxNSA];
    TCanvas* c_integral_Weight[4][maxNSA];
    for (int i=0; i<4; i++){
      TString i_index = TString::Itoa(i,10);
      max_peak[i] = max_peak[i]*1.05;
      min_peak[i] = min_peak[i]*0.9<1 ? 0 : min_peak[i]*0.9;
      c_peak_NoWeight[i] = new TCanvas("c_nw_"+i_index,"c_nw_"+i_index,900,900);
      c_peak_Weight[i]   = new TCanvas("c_w_"+i_index, "c_w_"+i_index,900,900);
      for (int j=minNSA; j<maxNSA; j++){
        TString j_index = TString::Itoa(j,10);
        max_integral[i][j] = max_integral[i][j]*1.05;
        min_integral[i][j] = min_integral[i][j]*0.9<1 ? 0 : min_integral[i][j]*0.9;
        c_integral_NoWeight[i][j] = new TCanvas("c_nw_"+i_index+"-"+j_index,"c_nw_"+i_index+"-"+j_index,900,900);
        c_integral_Weight[i][j]   = new TCanvas("c_w_"+i_index+"-"+j_index, "c_w_"+i_index+"-"+j_index,900,900);
      }
    }
    // Draw peak and integral histograms
    for (int i=0; i<4; i++){
      h_peak_NoWeight[i]->GetXaxis()->SetRangeUser(min_peak[i],max_peak[i]);
      h_peak_Weight[i]->GetXaxis()->SetRangeUser(min_peak[i],max_peak[i]);
    
      c_peak_NoWeight[i]->cd(); h_peak_NoWeight[i]->Draw();
      c_peak_Weight[i]->cd();   h_peak_Weight[i]->Draw();
      for (int j=minNSA; j<maxNSA; j++){
        h_integral_NoWeight[i][j]->GetXaxis()->SetRangeUser(min_integral[i][j],max_integral[i][j]);
        h_integral_Weight[i][j]->GetXaxis()->SetRangeUser(min_integral[i][j],max_integral[i][j]);
      
        c_integral_NoWeight[i][j]->cd(); h_integral_NoWeight[i][j]->Draw();
        c_integral_Weight[i][j]->cd();   h_integral_Weight[i][j]->Draw();
      }
    }
    TCanvas* c_error_vs_peak[maxNSA+1];
    // Draw error vs peak for pulse peak and integral
    for (int i=0; i<maxNSA+1; i++){
      c_error_vs_peak[i] = new TCanvas("c_error_vs_peak_"+TString::Itoa(i,10), "c_error_vs_peak_"+TString::Itoa(i,10), 900, 900);
      c_error_vs_peak[i]->cd();
      gr_error_vs_peak[i]->Draw("A");
      gr_error_vs_peak[i]->GetXaxis()->SetTitle("peak");
      gr_error_vs_peak[i]->GetYaxis()->SetTitle("%");
      gr_error_vs_peak[i]->Draw("AP");
    }
    if (save){
      for (int i=0; i<4; i++){
        if (i==0){
          c_peak_NoWeight[i]->SaveAs(dir+"NoWeights/Avg_Peak.pdf");
          c_peak_Weight[i]->SaveAs(dir+"Weights/Avg_Peak.pdf");
        }
        if (i==1){
          c_peak_NoWeight[i]->SaveAs(dir+"NoWeights/RMS_Peak.pdf");
          c_peak_Weight[i]->SaveAs(dir+"Weights/RMS_Peak.pdf");
        }
        if (i==2){
          c_peak_NoWeight[i]->SaveAs(dir+"NoWeights/Error_Peak.pdf");
          c_peak_Weight[i]->SaveAs(dir+"Weights/Error_Peak.pdf");
        }
        if (i==3){
          c_peak_NoWeight[i]->SaveAs(dir+"NoWeights/PedSubtracted_Error_Peak.pdf");
          c_peak_Weight[i]->SaveAs(dir+"Weights/PedSubtracted_Error_Peak.pdf");
        }
        for (int j=minNSA; j<maxNSA; j++){
          TString t_nsa = TString::Itoa(j,10);
          if (i==0){
            c_integral_NoWeight[i][j]->SaveAs(dir+"NoWeights/Avg_Integral_"+t_nsa+".pdf");
            c_integral_Weight[i][j]->SaveAs(dir+"Weights/Avg_Integral_"+t_nsa+".pdf");
          }
          if (i==1){
            c_integral_NoWeight[i][j]->SaveAs(dir+"NoWeights/Rms_Integral_"+t_nsa+".pdf");
            c_integral_Weight[i][j]->SaveAs(dir+"Weights/Rms_Integral_"+t_nsa+".pdf");
          }
          if (i==2){
            c_integral_NoWeight[i][j]->SaveAs(dir+"NoWeights/Error_Integral_"+t_nsa+".pdf");
            c_integral_Weight[i][j]->SaveAs(dir+"Weights/Error_Integral_"+t_nsa+".pdf");
          }
          if (i==3){
            c_integral_NoWeight[i][j]->SaveAs(dir+"NoWeights/PedSubtracted_Error_Integral_"+t_nsa+".pdf");
            c_integral_Weight[i][j]->SaveAs(dir+"Weights/PedSubtracted_Error_Integral_"+t_nsa+".pdf");
          }
        }
      }
      for (int i=0; i<maxNSA+1; i++){
        TString t_nsa = TString::Itoa(i,10);
        if (i==0) c_error_vs_peak[i]->SaveAs(dir+"Errors_VS_Peak/PulsePeak.pdf");
        else c_error_vs_peak[i]->SaveAs(dir+"Errors_VS_Peak/PulseIntegral_"+t_nsa+".pdf");
      }
    }
  }
  
  TCanvas* c_integral_vs_nsa[4];
  for (int i=0; i<4; i++){
    TString i_index = TString::Itoa(i,10);
    c_integral_vs_nsa[i] = new TCanvas("c_integral_vs_nsa_"+i_index,"c_integral_vs_nsa_"+i_index,900,900);
    c_integral_vs_nsa[i]->cd();
    gr_integral_vs_nsa[i]->Draw("A");
    gr_integral_vs_nsa[i]->GetXaxis()->SetTitle("nsa");
    gr_integral_vs_nsa[i]->GetYaxis()->SetTitle("error (%)");
    gr_integral_vs_nsa[i]->SetMarkerSize(1.4);
    gr_integral_vs_nsa[i]->Draw("AP");
  }
  
  if (save){
    for (int i=0; i<4; i++){
      //if (i==0) c_integral_vs_nsa[i]->SaveAs(dir+"Integral_VS_NSA/IntegralVsNSA_gaus.pdf");
      //if (i==1) c_integral_vs_nsa[i]->SaveAs(dir+"Integral_VS_NSA/PedSubtracted_IntegralVsNSA_gaus.pdf");
      //if (i==2) c_integral_vs_nsa[i]->SaveAs(dir+"Integral_VS_NSA/IntegralVsNSA_weighted_gaus.pdf");
      //if (i==3) c_integral_vs_nsa[i]->SaveAs(dir+"Integral_VS_NSA/PedSubtracted_IntegralVsNSA_weighted_gaus.pdf");
      if (i==0) c_integral_vs_nsa[i]->SaveAs(dir+"Integral_VS_NSA/IntegralVsNSA.pdf");
      if (i==1) c_integral_vs_nsa[i]->SaveAs(dir+"Integral_VS_NSA/PedSubtracted_IntegralVsNSA.pdf");
      if (i==2) c_integral_vs_nsa[i]->SaveAs(dir+"Integral_VS_NSA/IntegralVsNSA_weighted.pdf");
      if (i==3) c_integral_vs_nsa[i]->SaveAs(dir+"Integral_VS_NSA/PedSubtracted_IntegralVsNSA_weighted.pdf");
    }
  }
}


void PulseAnalysis_Waveform(TString filename){
  
  TFile * f = new TFile(filename);
  
  TTree* FCAL_Waveforms    = (TTree*)f->Get("FCAL_Waveforms");
  SetFCAL_Waveforms(FCAL_Waveforms);
  
  
  TGraph* gr_all = new TGraph();
  gr_all->SetName("gr_all");
  gr_all->SetTitle("Waveform for every channel and every event");
  
  int nEntries = FCAL_Waveforms->GetEntries();
  for (int i=0; i<nEntries; i++){
    FCAL_Waveforms->GetEntry(i);
    gr_all->SetPoint(i,sample+1,adc_val);
  }
  TCanvas* c_all = new TCanvas("c_all","c_all",900,900);
  c_all->cd();
  gr_all->Draw("A");
  gr_all->GetXaxis()->SetTitle("sample");
  gr_all->Draw("AP");
}

void PulseAnalysis_PeakRatio(TString filename, bool save=false){
  
  TFile * f = new TFile(filename);
  
  TTree* FCAL_Correlation    = (TTree*)f->Get("FCAL_Correlation");
  SetFCAL_Correlation(FCAL_Correlation);
  
  
  double xmin=-29.5; double xmax=29.5; double ymin=-29.5; double ymax=29.5; double zmin=0; double zmax=-1;
  TH2F* boxBlock4_mean = new TH2F("boxBlock4_mean","(mean) pulse peak ratio, P_{i}/P_{j}, where 5x5 block centered at i of pulse peak ratio mean",xmax-xmin,xmin,xmax,ymax-ymin,ymin,ymax);
  if (zmax > zmin)  boxBlock4_mean->SetMaximum(zmax);
  
  TGraph* gr_peakRatioRms = new TGraph();
  gr_peakRatioRms->SetName("gr_peakRatioRms");
  gr_peakRatioRms->SetTitle("Pulse peak ratio rms vs pulse peak");
  TGraph* gr_peakRatio = new TGraph();
  gr_peakRatio->SetName("gr_peakRatio");
  gr_peakRatio->SetTitle("Pulse peak ratio vs pulse peak");
  
  TH1F* h_mean_ratio  = new TH1F("h_mean_ratio" ,"Average pulse peak ratio", maxADC/0.10, 0, maxADC);
  TH1F* h_rms_ratio  = new TH1F("h_rms_ratio" ,"RMS of pulse peak ratio", maxADC/0.0010, 0, maxADC);
  TH1F* h_error_ratio  = new TH1F("h_error_ratio" ,"Error of pulse peak ratio", 1000/0.10, 0, 1000);
  double max_mean_ratio = 0;
  double max_rms_ratio = 0;
  double max_error_ratio = 0;
  
  map< pair<int,int>, int > theChosen;
  map< pair<int,int>, double> peak_ratio_map;
  map< pair<int,int>, double> peak_counter_map;
  map< pair<int,int>, double>::const_iterator iter;
  // initialize val to 0 for all blocks in 60x60 blocks or less
  for(int index_x=-30; index_x<30; index_x++){
    for(int index_y=-30; index_y<30; index_y++){
      if (abs(index_x)%3==0 && abs(index_y)%3==0) {
        theChosen[ make_pair(index_x,index_y) ] = 1;
      } else{
        theChosen[ make_pair(index_x,index_y) ] = 1;
      }
      peak_counter_map[ make_pair(index_x,index_y) ] = 0;
      peak_ratio_map[ make_pair(index_x,index_y) ] = 0;
    }
  }
  
  
  int nEntries = FCAL_Correlation->GetEntries();
  
  for (int i=0; i<nEntries; i++){
    FCAL_Correlation->GetEntry(i);
    
    if (abs(pulse1_x-pulse2_x)<=4&&abs(pulse1_y-pulse2_y)<=4){
      if (theChosen[ make_pair(pulse1_x,pulse1_y) ]==1){
        peak_counter_map[ make_pair(pulse1_x,pulse1_y) ]++;
        peak_ratio_map[ make_pair(pulse1_x,pulse1_y) ] += peak_ratio_mean;
      }
    }
    
    //if (abs(pulse1_x-pulse2_x)<=3&&abs(pulse1_y-pulse2_y)<=3){
    //cout << peak_ratio_mean << endl;
    //if (peak_ratio_mean<4096 && peak_ratio_mean>0&&abs(pulse1_x-pulse2_x)<=3&&abs(pulse1_y-pulse2_y)<=3){
    //if (pulse1_x==(pulse2_x-1)&&pulse1_y==pulse2_y){
    if (abs(pulse1_x-pulse2_x)<=2&&abs(pulse1_y-pulse2_y)<=2){
      //double error = peak_ratio_rms_gaus/peak_ratio_mean*100.0;
      double error = peak_ratio_rms/peak_ratio_mean*100.0;
      if (error>100) {
        cout << "Pulse1: " << "(" << pulse1_x << "," << pulse1_y << ")" << endl;
        cout << "Pulse2: " << "(" << pulse2_x << "," << pulse2_y << ")" << endl;
      }
      //gr_peakRatioRms->SetPoint(gr_peakRatioRms->GetN(),pulse1_avgPeak,peak_ratio_rms_gaus);
      gr_peakRatioRms->SetPoint(gr_peakRatioRms->GetN(),pulse1_avgPeak,peak_ratio_rms);
      gr_peakRatio   ->SetPoint(gr_peakRatio->GetN(),pulse1_avgPeak,peak_ratio_mean);
      h_mean_ratio   ->Fill(peak_ratio_mean);
      h_rms_ratio    ->Fill(peak_ratio_rms);
      h_error_ratio  ->Fill(error);
      if (max_rms_ratio<peak_ratio_rms) max_rms_ratio = peak_ratio_rms;
      if (max_mean_ratio<peak_ratio_mean) max_mean_ratio = peak_ratio_mean;
      if (max_error_ratio<error) max_error_ratio = error;
    }
    
  }
  
  TCanvas* c_2D_4x4_mean = new TCanvas("c_2D_4x4_mean","c_2D_4x4_mean",900,900);
  c_2D_4x4_mean->cd();
  boxBlock4_mean->SetXTitle("x");
  boxBlock4_mean->SetYTitle("y");
  boxBlock4_mean->SetStats(false);
  for (iter=peak_ratio_map.begin(); iter!=peak_ratio_map.end(); ++iter){
    if (iter->second>0){
      boxBlock4_mean->Fill(iter->first.first,iter->first.second,iter->second/peak_counter_map[iter->first]);
    }
  }
  boxBlock4_mean->Draw("colz");
  
  TCanvas* c_peakRatioRms = new TCanvas("c_peakRatioRms","c_peakRatioRms",900,900);
  c_peakRatioRms->cd();
  gr_peakRatioRms->SetMarkerStyle(3);
  gr_peakRatioRms->SetMarkerSize(1);
  gr_peakRatioRms->Draw("A");
  gr_peakRatioRms->GetXaxis()->SetTitle("pulse peak");
  gr_peakRatioRms->Draw("AP");
  cout << "N = " << gr_peakRatioRms->GetN() << endl;
  
  TCanvas* c_peakRatio = new TCanvas("c_peakRatio","c_peakRatio",900,900);
  c_peakRatio->cd();
  gr_peakRatio->SetMarkerStyle(3);
  gr_peakRatio->SetMarkerSize(1);
  gr_peakRatio->Draw("A");
  gr_peakRatio->GetXaxis()->SetTitle("pulse peak");
  gr_peakRatio->Draw("AP");
  
  TCanvas* c_h_peakRatio = new TCanvas("c_h_peakRatio","c_h_peakRatio",900,900);
  c_h_peakRatio->cd();
  h_mean_ratio->GetXaxis()->SetRangeUser(0,max_mean_ratio);
  h_mean_ratio->Draw();
  
  TCanvas* c_h_peakRatioRms = new TCanvas("c_h_peakRatioRms","c_h_peakRatioRms",900,900);
  c_h_peakRatioRms->cd();
  h_rms_ratio->GetXaxis()->SetRangeUser(0,max_rms_ratio);
  h_rms_ratio->Draw();
  
  TCanvas* c_h_peakRatioError = new TCanvas("c_h_peakRatioError","c_h_peakRatioError",900,900);
  c_h_peakRatioError->cd();
  h_error_ratio->GetXaxis()->SetRangeUser(0,max_error_ratio);
  h_error_ratio->Draw();
  
  // save plots
  if (save){
    TString dir="plots/PulsePeakRatioAnalysis/";
    c_2D_4x4_mean->SaveAs(dir+"hist2d_4x4PulsePeakRatioMean.jpg");
    c_peakRatioRms->SaveAs(dir+"gr_PeakRatioRms.jpg");
    c_peakRatio->SaveAs(dir+"gr_PeakRatioMean.jpg");
    c_h_peakRatio->SaveAs(dir+"hist_PeakRatioMean.jpg");
    c_h_peakRatioRms->SaveAs(dir+"hist_PeakRatioRms.jpg");
    c_h_peakRatioError->SaveAs(dir+"hist_PeakRatioError.jpg");
  }
}

void PulseAnalysis_Correlation(TString filename, bool save=false){
  
  TFile * f = new TFile(filename);
  
  TTree* FCAL_Correlation    = (TTree*)f->Get("FCAL_Correlation");
  SetFCAL_Correlation(FCAL_Correlation);
  
  double xmin=-29.5; double xmax=29.5; double ymin=-29.5; double ymax=29.5; double zmin=0; double zmax=-1;
  TH2F* boxAvgCorr5x5 = new TH2F("boxAvgCorr5x5","Avg Pulse peak correlation (5x5)",xmax-xmin,xmin,xmax,ymax-ymin,ymin,ymax);
  TH2F* boxAvgCorr3x3 = new TH2F("boxAvgCorr3x3","Avg Pulse peak correlation (3x3)",xmax-xmin,xmin,xmax,ymax-ymin,ymin,ymax);
  TH2F* boxAvgCorr = new TH2F("boxAvgCorr","Avg Pulse peak correlation of all channels against each other",xmax-xmin,xmin,xmax,ymax-ymin,ymin,ymax);
  TH2F* boxMaxCorr5x5 = new TH2F("boxMaxCorr5x5","Max Pulse peak correlation (5x5)",xmax-xmin,xmin,xmax,ymax-ymin,ymin,ymax);
  TH2F* boxMaxCorr3x3 = new TH2F("boxMaxCorr3x3","Max Pulse peak correlation (3x3)",xmax-xmin,xmin,xmax,ymax-ymin,ymin,ymax);
  TH2F* boxMaxCorr = new TH2F("boxMaxCorr","Max Pulse peak correlation of all channels against each other",xmax-xmin,xmin,xmax,ymax-ymin,ymin,ymax);
  boxAvgCorr5x5->SetXTitle("x"); boxAvgCorr5x5->SetYTitle("y"); boxAvgCorr5x5->SetStats(false);
  boxAvgCorr3x3->SetXTitle("x"); boxAvgCorr3x3->SetYTitle("y"); boxAvgCorr3x3->SetStats(false);
  boxAvgCorr->SetXTitle("x");    boxAvgCorr->SetYTitle("y");    boxAvgCorr->SetStats(false);
  boxMaxCorr5x5->SetXTitle("x"); boxMaxCorr5x5->SetYTitle("y"); boxMaxCorr5x5->SetStats(false);
  boxMaxCorr3x3->SetXTitle("x"); boxMaxCorr3x3->SetYTitle("y"); boxMaxCorr3x3->SetStats(false);
  boxMaxCorr->SetXTitle("x");    boxMaxCorr->SetYTitle("y");    boxMaxCorr->SetStats(false);
  //if (zmax > zmin)  boxBlock4_mean->SetMaximum(zmax);
  
  TH1F* h_corr_5x5  = new TH1F("h_corr_5x5" ,"Pulse peak correlation (5x5)", 2.0/0.02, -1, 1);
  TH1F* h_corr_3x3  = new TH1F("h_corr_3x3" ,"Pulse peak correlation (3x3)", 2.0/0.02, -1, 1);
  TH1F* h_corr  = new TH1F("h_corr" ,"Pulse peak correlation of all channels against each other", 2.0/0.02, -1, 1);

  double max_corr_5x5 = -1;
  double min_corr_5x5 = 1;
  double max_corr_3x3 = -1;
  double min_corr_3x3 = 1;
  double max_corr     = -1;
  double min_corr     = 1;
  
  int nEntries = FCAL_Correlation->GetEntries();
  
  //coordinate, sum,entries
  map< pair<int,int>, pair<double,double> >::const_iterator iter_avg;
  map< pair<int,int>, pair<double,double> > corr_avg_map_5x5;
  map< pair<int,int>, pair<double,double> > corr_avg_map_3x3;
  map< pair<int,int>, pair<double,double> > corr_avg_map;
  map< pair<int,int>, double >::const_iterator iter_max;
  map< pair<int,int>, double > corr_max_map_5x5;
  map< pair<int,int>, double > corr_max_map_3x3;
  map< pair<int,int>, double > corr_max_map;
  // initialize val to 0 for all blocks in 60x60 blocks or less
  for(int index_x=-30; index_x<30; index_x++){
    for(int index_y=-30; index_y<30; index_y++){
      corr_max_map_5x5[ make_pair(index_x,index_y) ] = -1.1;
      corr_max_map_3x3[ make_pair(index_x,index_y) ] = -1.1;
      corr_max_map[ make_pair(index_x,index_y) ] = -1.1;
      corr_avg_map_5x5[ make_pair(index_x,index_y) ] = make_pair(0,0);
      corr_avg_map_3x3[ make_pair(index_x,index_y) ] = make_pair(0,0);
      corr_avg_map[ make_pair(index_x,index_y) ] = make_pair(0,0);
    }
  }
  
  for (int i=0; i<nEntries; i++){
    FCAL_Correlation->GetEntry(i);
    
    if (abs(pulse1_x-pulse2_x)<=2&&abs(pulse1_y-pulse2_y)<=2){
      if (max_corr_5x5<corr) max_corr_5x5 = corr;
      if (min_corr_5x5>corr) min_corr_5x5 = corr;
      h_corr_5x5->Fill(corr);
      // get max
      if ( corr_max_map_5x5[make_pair(pulse1_x,pulse1_y)]<corr ) corr_max_map_5x5[make_pair(pulse1_x,pulse1_y)] = corr;
      if ( corr_max_map_5x5[make_pair(pulse2_x,pulse2_y)]<corr ) corr_max_map_5x5[make_pair(pulse2_x,pulse2_y)] = corr;
      // get avg
      pair<double,double> avg_pair = corr_avg_map_5x5[make_pair(pulse1_x,pulse1_y)];
      corr_avg_map_5x5[make_pair(pulse1_x,pulse1_y)] = make_pair(avg_pair.first+corr,avg_pair.second+1);
      pair<double,double> avg_pair2 = corr_avg_map_5x5[make_pair(pulse2_x,pulse2_y)];
      corr_avg_map_5x5[make_pair(pulse2_x,pulse2_y)] = make_pair(avg_pair2.first+corr,avg_pair2.second+1);
    }
    if (abs(pulse1_x-pulse2_x)<=1&&abs(pulse1_y-pulse2_y)<=1){
      if (max_corr_3x3<corr) max_corr_3x3 = corr;
      if (min_corr_3x3>corr) min_corr_3x3 = corr;
      h_corr_3x3->Fill(corr);
      // get max
      if ( corr_max_map_3x3[make_pair(pulse1_x,pulse1_y)]<corr ) corr_max_map_3x3[make_pair(pulse1_x,pulse1_y)] = corr;
      if ( corr_max_map_3x3[make_pair(pulse2_x,pulse2_y)]<corr ) corr_max_map_3x3[make_pair(pulse2_x,pulse2_y)] = corr;
      // get avg
      pair<double,double> avg_pair = corr_avg_map_3x3[make_pair(pulse1_x,pulse1_y)];
      corr_avg_map_3x3[make_pair(pulse1_x,pulse1_y)] = make_pair(avg_pair.first+corr,avg_pair.second+1);
      pair<double,double> avg_pair2 = corr_avg_map_3x3[make_pair(pulse2_x,pulse2_y)];
      corr_avg_map_3x3[make_pair(pulse2_x,pulse2_y)] = make_pair(avg_pair2.first+corr,avg_pair2.second+1);
    }
    if (max_corr<corr) max_corr = corr;
    if (min_corr>corr) min_corr = corr;
    h_corr->Fill(corr);
    // get max
    if ( corr_max_map[make_pair(pulse1_x,pulse1_y)]<corr ) corr_max_map[make_pair(pulse1_x,pulse1_y)] = corr;
    if ( corr_max_map[make_pair(pulse2_x,pulse2_y)]<corr ) corr_max_map[make_pair(pulse2_x,pulse2_y)] = corr;
    // get avg
    pair<double,double> avg_pair = corr_avg_map[make_pair(pulse1_x,pulse1_y)];
    corr_avg_map[make_pair(pulse1_x,pulse1_y)] = make_pair(avg_pair.first+corr,avg_pair.second+1);
    pair<double,double> avg_pair2 = corr_avg_map[make_pair(pulse2_x,pulse2_y)];
    corr_avg_map[make_pair(pulse2_x,pulse2_y)] = make_pair(avg_pair2.first+corr,avg_pair2.second+1);
  }
  
  for (iter_max=corr_max_map_5x5.begin(); iter_max!=corr_max_map_5x5.end(); ++iter_max){
    if (iter_max->second>-1.1){
      boxMaxCorr5x5->Fill(iter_max->first.first,iter_max->first.second,iter_max->second);
    }
  }
  for (iter_max=corr_max_map_3x3.begin(); iter_max!=corr_max_map_3x3.end(); ++iter_max){
    if (iter_max->second>-1.1){
      boxMaxCorr3x3->Fill(iter_max->first.first,iter_max->first.second,iter_max->second);
    }
  }
  for (iter_max=corr_max_map.begin(); iter_max!=corr_max_map.end(); ++iter_max){
    if (iter_max->second>-1.1){
      boxMaxCorr->Fill(iter_max->first.first,iter_max->first.second,iter_max->second);
    }
  }
  
  for (iter_avg=corr_avg_map_5x5.begin(); iter_avg!=corr_avg_map_5x5.end(); ++iter_avg){
    boxAvgCorr5x5->Fill(iter_avg->first.first,iter_avg->first.second,iter_avg->second.first/iter_avg->second.second);
  }
  for (iter_avg=corr_avg_map_3x3.begin(); iter_avg!=corr_avg_map_3x3.end(); ++iter_avg){
    boxAvgCorr3x3->Fill(iter_avg->first.first,iter_avg->first.second,iter_avg->second.first/iter_avg->second.second);
  }
  for (iter_avg=corr_avg_map.begin(); iter_avg!=corr_avg_map.end(); ++iter_avg){
    boxAvgCorr->Fill(iter_avg->first.first,iter_avg->first.second,iter_avg->second.first/iter_avg->second.second);
  }
  
  TCanvas* c_2d_max_corr[3];
  c_2d_max_corr[0] = new TCanvas("c_2D_max_corr_5x5","c_2D_max_corr_5x5",1200,900);
  c_2d_max_corr[1] = new TCanvas("c_2D_max_corr_3x3","c_2D_max_corr_3x3",1200,900);
  c_2d_max_corr[2] = new TCanvas("c_2D_max_corr","c_2D_max_corr",1200,900);
  TCanvas* c_2d_avg_corr[3];
  c_2d_avg_corr[0] = new TCanvas("c_2D_avg_corr_5x5","c_2D_avg_corr_5x5",1200,900);
  c_2d_avg_corr[1] = new TCanvas("c_2D_avg_corr_3x3","c_2D_avg_corr_3x3",1200,900);
  c_2d_avg_corr[2] = new TCanvas("c_2D_avg_corr","c_2D_avg_corr",1200,900);
  
  for (int i=0; i<3; i++){
    c_2d_max_corr[i]->cd();
    if (i==0) boxMaxCorr5x5->Draw("colz");
    if (i==1) boxMaxCorr3x3->Draw("colz");
    if (i==2) boxMaxCorr->Draw("colz");
  }
  for (int i=0; i<3; i++){
    c_2d_avg_corr[i]->cd();
    if (i==0) boxAvgCorr5x5->Draw("colz");
    if (i==1) boxAvgCorr3x3->Draw("colz");
    if (i==2) boxAvgCorr->Draw("colz");
  }
  
  TCanvas* c_h_corr[3];
  c_h_corr[0] = new TCanvas("c_h_corr_5x5","c_h_corr_5x5",900,900);
  c_h_corr[1] = new TCanvas("c_h_corr_3x3","c_h_corr_3x3",900,900);
  c_h_corr[2] = new TCanvas("c_h_corr","c_h_corr",900,900);
  
  for (int i=0; i<3; i++){
    c_h_corr[i]->cd();
    if (i==0) {h_corr_5x5->GetXaxis()->SetRangeUser(min_corr_5x5,max_corr_5x5); h_corr_5x5->Draw();}
    if (i==1) {h_corr_3x3->GetXaxis()->SetRangeUser(min_corr_3x3,max_corr_3x3); h_corr_3x3->Draw();}
    if (i==2) {h_corr->GetXaxis()->SetRangeUser(min_corr,max_corr); h_corr->Draw();}
  }
  
  // save plots
  if (save){
    TString dir="plots/PulsePeakCorrelationAnalysis/";
    cout << "plot" << endl;
    for (int i=0; i<3; i++){
      if (i==0) {
        c_2d_max_corr[i]->SaveAs(dir+"hist2d_max_corr_5x5.jpg");
        c_2d_avg_corr[i]->SaveAs(dir+"hist2d_avg_corr_5x5.jpg");
        c_h_corr[i]->SaveAs(dir+"hist_corr_5x5.jpg");
      }
      if (i==1) {
        c_2d_max_corr[i]->SaveAs(dir+"hist2d_max_corr_3x3.jpg");
        c_2d_avg_corr[i]->SaveAs(dir+"hist2d_avg_corr_3x3.jpg");
        c_h_corr[i]->SaveAs(dir+"hist_corr_3x3.jpg");
      }
      if (i==2) {
        c_2d_max_corr[i]->SaveAs(dir+"hist2d_max_corr.jpg");
        c_2d_avg_corr[i]->SaveAs(dir+"hist2d_avg_corr.jpg");
        c_h_corr[i]->SaveAs(dir+"hist_corr.jpg");
      }
    }
  }
}

// ****************************
// z + dz = a + da / b + db => dz/z = sqrt(da/a^2 + db/b^2) if and only if correlation is 0
// plot dz/z - sqrt(da/a^2 + db/b^2)
// ****************************
void PulseAnalysis_ErrorAnalysis(TString filename, bool save=false){
  TFile * f = new TFile(filename);
  
  TTree* FCAL_Correlation    = (TTree*)f->Get("FCAL_Correlation");
  SetFCAL_Correlation(FCAL_Correlation);
  TTree* FCAL_Analysis    = (TTree*)f->Get("FCAL_Analysis");
  SetFCAL_Analysis(FCAL_Analysis);
  
  int nCorr = FCAL_Correlation->GetEntries();
  int nAna  = FCAL_Analysis->GetEntries();
  
  // fractional error for ratio and peaks
  map< pair< pair<int,int>, pair<int,int> >, double> ratio_map;
  map< pair< pair<int,int>, pair<int,int> >, double> ratio_map_3x3;
  map< pair< pair<int,int>, pair<int,int> >, double> ratio_map_5x5;
  //map< pair< pair<int,int>, pair<int,int> >, int> ratio_count;
  map< pair<int,int>, double> peak_map;
  
  for (int i=0; i<nCorr; i++){
    FCAL_Correlation->GetEntry(i);
    //corr_map[make_pair(pulse1_x,pulse1_y)] = peak_ratio_rms_gaus/peak_ratio_mean;
    //corr_map[make_pair(pulse2_x,pulse2_y)] = peak_ratio_rms_gaus/peak_ratio_mean;
    pair<int,int> coord1 = make_pair(pulse1_x,pulse1_y);
    pair<int,int> coord2 = make_pair(pulse2_x,pulse2_y);
    //ratio_map[make_pair(coord1,coord2)] = peak_ratio_rms_gaus/peak_ratio_mean;
    ratio_map[make_pair(coord1,coord2)] = peak_ratio_rms/peak_ratio_mean;
    if (abs(pulse1_x-pulse2_x)<=1&&abs(pulse1_y-pulse2_y)<=1){
      //ratio_map_3x3[make_pair(coord1,coord2)] = peak_ratio_rms_gaus/peak_ratio_mean;
      ratio_map_3x3[make_pair(coord1,coord2)] = peak_ratio_rms/peak_ratio_mean;
    }
    if (abs(pulse1_x-pulse2_x)<=2&&abs(pulse1_y-pulse2_y)<=2){
      //ratio_map_5x5[make_pair(coord1,coord2)] = peak_ratio_rms_gaus/peak_ratio_mean;
      ratio_map_5x5[make_pair(coord1,coord2)] = peak_ratio_rms/peak_ratio_mean;
    }
    
    //ratio_count[make_pair(coord1,coord2)]++;
  }
  for (int i=0; i<nAna; i++){
    FCAL_Analysis->GetEntry(i);
    //peak_map[make_pair(m_x,m_y)] = rms_peak_gaus/avg_peak;
    peak_map[make_pair(m_x,m_y)] = rms_peak/avg_peak;
  }

  TH1F* h_fracErrorDiff = new TH1F("h_fracErrorDiff","peakRatio_fracError - sqrt(peak1_fracError**2+peak2_fracError**2))",400,-1,3);
  TH1F* h_fracErrorDiff_3x3 = new TH1F("h_fracErrorDiff_3x3","peakRatio_fracError - sqrt(peak1_fracError**2+peak2_fracError**2)) (3x3)",400,-1,3);
  TH1F* h_fracErrorDiff_5x5 = new TH1F("h_fracErrorDiff_5x5","peakRatio_fracError - sqrt(peak1_fracError**2+peak2_fracError**2)) (5x5)",400,-1,3);
  for (map< pair< pair<int,int>, pair<int,int> >, double>::const_iterator iter=ratio_map.begin(); iter!=ratio_map.end(); ++iter){
    pair<int,int> coord1 = iter->first.first;
    pair<int,int> coord2 = iter->first.second;
    double ratioFracError = iter->second;
    double peak1FracError = peak_map[coord1];
    double peak2FracError = peak_map[coord2];
    h_fracErrorDiff->Fill(ratioFracError- sqrt(peak1FracError*peak1FracError+peak2FracError*peak2FracError));
    //if (iter->second>1) printf("Pulse1: (%i,%i), Pulse2: (%i,%i) = %i\n",coord1.first,coord1.second,coord2.first,coord2.second,iter->second);
  }
  for (map< pair< pair<int,int>, pair<int,int> >, double>::const_iterator iter=ratio_map_3x3.begin(); iter!=ratio_map_3x3.end(); ++iter){
    pair<int,int> coord1 = iter->first.first;
    pair<int,int> coord2 = iter->first.second;
    double ratioFracError = iter->second;
    double peak1FracError = peak_map[coord1];
    double peak2FracError = peak_map[coord2];
    h_fracErrorDiff_3x3->Fill(ratioFracError- sqrt(peak1FracError*peak1FracError+peak2FracError*peak2FracError));
  }
  for (map< pair< pair<int,int>, pair<int,int> >, double>::const_iterator iter=ratio_map_5x5.begin(); iter!=ratio_map_5x5.end(); ++iter){
    pair<int,int> coord1 = iter->first.first;
    pair<int,int> coord2 = iter->first.second;
    double ratioFracError = iter->second;
    double peak1FracError = peak_map[coord1];
    double peak2FracError = peak_map[coord2];
    h_fracErrorDiff_5x5->Fill(ratioFracError- sqrt(peak1FracError*peak1FracError+peak2FracError*peak2FracError));
  }
  
  TCanvas* c_fracErrorDiff = new TCanvas("c_fracErrorDiff","c_fracErrorDiff",900,900);
  c_fracErrorDiff->cd();
  h_fracErrorDiff->Draw();
  
  TCanvas* c_fracErrorDiff_3x3 = new TCanvas("c_fracErrorDiff_3x3","c_fracErrorDiff_3x3",900,900);
  c_fracErrorDiff_3x3->cd();
  h_fracErrorDiff_3x3->Draw();
  
  TCanvas* c_fracErrorDiff_5x5 = new TCanvas("c_fracErrorDiff_5x5","c_fracErrorDiff_5x5",900,900);
  c_fracErrorDiff_5x5->cd();
  h_fracErrorDiff_5x5->Draw();
  
  if (save){
    TString dir = "plots/PulsePeakRatioAnalysis/";
    c_fracErrorDiff->SaveAs(dir+"FracErrorDiff.pdf");
    c_fracErrorDiff_3x3->SaveAs(dir+"FracErrorDiff_3x3.pdf");
    c_fracErrorDiff_5x5->SaveAs(dir+"FracErrorDiff_5x5.pdf");
  }
  
}

// Calculate pulse integral from pulse peak
void PulseAnalysis_NSP(TString filename, bool save=false, bool supporting=false){

  TFile * f = new TFile(filename);
  
  TTree* FCAL_Analysis2    = (TTree*)f->Get("FCAL_Analysis2");
  SetFCAL_Analysis2(FCAL_Analysis2);
  
  int n_nsp = max_nsp-min_nsp+1;
  //TH1F* h_nsp_avgIntegral[n_nsp];
  //TH1F* h_nsp_rmsIntegral[n_nsp];
  TH1F* h_nsp_errorIntegral_W[n_nsp];
  TH1F* h_nsp_errorIntegral_NW[n_nsp];
  TH1F* h_nsp_PedSub_errorIntegral_W[n_nsp];
  TH1F* h_nsp_PedSub_errorIntegral_NW[n_nsp];
  
  double max_nsp_errorIntegral[n_nsp];
  double min_nsp_errorIntegral[n_nsp];
  double max_nsp_PedSub_errorIntegral[n_nsp];
  double min_nsp_PedSub_errorIntegral[n_nsp];
  
  for (int i=0; i<n_nsp; i++){
    max_nsp_errorIntegral[i] = 0;
    min_nsp_errorIntegral[i] = 100;
    max_nsp_PedSub_errorIntegral[i] = 0;
    min_nsp_PedSub_errorIntegral[i] = 100;
  }
  
  for (int i=0; i<n_nsp; i++){
    int _nsp = i+min_nsp;
    TString t_nsp = TString::Itoa(_nsp,10);
    //h_nsp_avgIntegral[i] = new TH1F("h_nsp_avgIntegral_"+t_nsp,"Avg Pulse Integral (nsp="+t_nsp+")",4096*100,0,4096*100);
    //h_nsp_rmsIntegral[i] = new TH1F("h_nsp_rmsIntegral_"+t_nsp,"Rms Pulse Integral (nsp="+t_nsp+")",4096*100,0,4096*100);
    h_nsp_errorIntegral_W[i] = new TH1F("h_nsp_errorIntegral_w_"+t_nsp,"Error Pulse Integral (nsp="+t_nsp+", weight)",1000,0,100); 
    h_nsp_errorIntegral_W[i]->SetXTitle("(%)");
    h_nsp_errorIntegral_NW[i] = new TH1F("h_nsp_errorIntegral_nw_"+t_nsp,"Error Pulse Integral (nsp="+t_nsp+", no weight)",1000,0,100);
    h_nsp_errorIntegral_NW[i]->SetXTitle("(%)");
    h_nsp_PedSub_errorIntegral_W[i] =  new TH1F("h_nsp_pedsub_errorIntegral_w_"+t_nsp,"Error Pedestal Subtracted Pulse Integral (nsp="+t_nsp+", weight)",1000,0,100);
    h_nsp_PedSub_errorIntegral_W[i]->SetXTitle("(%)");
    h_nsp_PedSub_errorIntegral_NW[i] = new TH1F("h_nsp_pedsub_errorIntegral_nw_"+t_nsp,"Error Pedestal Subtracted Pulse Integral (nsp="+t_nsp+", no weight)",1000,0,100);
    h_nsp_PedSub_errorIntegral_NW[i]->SetXTitle("(%)");
  }

  int nAna  = FCAL_Analysis2->GetEntries();  
  for (int i=0; i<nAna; i++){
    FCAL_Analysis2->GetEntry(i);
    
    int _nsp = nsp-min_nsp;
    
    double weight = Weight( avg_peak );
    //double errorIntegral = rms_integral_gaus/avg_integral*100;
    double errorIntegral = rms_integral/avg_integral*100;
    //double PedSub_errorIntegral = TMath::Sqrt(TMath::Power(rms_integral_gaus/(avg_integral-100.0*(2.0*double(nsp)+1.0)),2)+ TMath::Power(0.5*2*double(nsp)/100.0,2))*100;
    double PedSub_errorIntegral = TMath::Sqrt(TMath::Power(rms_integral/(avg_integral-100.0*(2.0*double(nsp)+1.0)),2)+ TMath::Power(0.5*2*double(nsp)/100.0,2))*100;
    h_nsp_errorIntegral_NW[_nsp]->Fill(errorIntegral);
    h_nsp_errorIntegral_W[_nsp]->Fill(errorIntegral,weight);
    h_nsp_PedSub_errorIntegral_NW[_nsp]->Fill(PedSub_errorIntegral);
    h_nsp_PedSub_errorIntegral_W[_nsp] ->Fill(PedSub_errorIntegral,weight);
    
    if (max_nsp_errorIntegral[_nsp]<errorIntegral) max_nsp_errorIntegral[_nsp] = errorIntegral;
    if (min_nsp_errorIntegral[_nsp]>errorIntegral) min_nsp_errorIntegral[_nsp] = errorIntegral;
    if (max_nsp_PedSub_errorIntegral[_nsp]<PedSub_errorIntegral) max_nsp_PedSub_errorIntegral[_nsp] = PedSub_errorIntegral;
    if (min_nsp_PedSub_errorIntegral[_nsp]>PedSub_errorIntegral) min_nsp_PedSub_errorIntegral[_nsp] = PedSub_errorIntegral;
  }
  TGraphErrors* gr_errorIntegral_NW = new TGraphErrors();
  gr_errorIntegral_NW->SetName("nw"); gr_errorIntegral_NW->SetTitle("Avg pulse integral error vs nsp");
  TGraphErrors* gr_errorIntegral_W  = new TGraphErrors();
  gr_errorIntegral_W->SetName("w"); gr_errorIntegral_W->SetTitle("Avg pulse integral error vs nsp (weighted)");
  
  TGraphErrors* gr_pedsub_errorIntegral_NW = new TGraphErrors();
  gr_pedsub_errorIntegral_NW->SetName("pednw"); gr_pedsub_errorIntegral_NW->SetTitle("Avg pedestal subtracted pulse integral error vs nsp");
  TGraphErrors* gr_pedsub_errorIntegral_W  = new TGraphErrors();
  gr_pedsub_errorIntegral_W->SetName("pedw"); gr_pedsub_errorIntegral_W->SetTitle("Avg pedestal subtracted pulse integral error vs nsp (weighted)");
  for (int i=0; i<n_nsp; i++){
    gr_errorIntegral_NW->SetPoint(i,i+min_nsp,h_nsp_errorIntegral_NW[i]->GetMean());
    gr_errorIntegral_NW->SetPointError(i,0,h_nsp_errorIntegral_NW[i]->GetRMS());
    gr_errorIntegral_W->SetPoint(i,i+min_nsp,h_nsp_errorIntegral_W[i]->GetMean());
    gr_errorIntegral_W->SetPointError(i,0,h_nsp_errorIntegral_W[i]->GetRMS());
    
    gr_pedsub_errorIntegral_NW->SetPoint(i,i+min_nsp,h_nsp_PedSub_errorIntegral_NW[i]->GetMean());
    gr_pedsub_errorIntegral_NW->SetPointError(i,0,h_nsp_PedSub_errorIntegral_NW[i]->GetRMS());
    gr_pedsub_errorIntegral_W->SetPoint(i,i+min_nsp,h_nsp_PedSub_errorIntegral_W[i]->GetMean());
    gr_pedsub_errorIntegral_W->SetPointError(i,0,h_nsp_PedSub_errorIntegral_W[i]->GetRMS());
  }
  TCanvas* c_errorIntegral_NW = new TCanvas("c_errorIntegral_NW","c_errorIntegral_NW",900,900);
  c_errorIntegral_NW->cd();
  gr_errorIntegral_NW->GetYaxis()->SetTitle("%");
  gr_errorIntegral_NW->GetXaxis()->SetTitle("nsp");
  gr_errorIntegral_NW->GetXaxis()->SetRangeUser(-1,max_nsp+1);
  gr_errorIntegral_NW->Draw("AP");
  TCanvas* c_errorIntegral_W = new TCanvas("c_errorIntegral_W","c_errorIntegral_W",900,900);
  c_errorIntegral_W->cd();
  gr_errorIntegral_W->GetYaxis()->SetTitle("%");
  gr_errorIntegral_W->GetXaxis()->SetTitle("nsp");
  gr_errorIntegral_W->GetXaxis()->SetRangeUser(-1,max_nsp+1);
  gr_errorIntegral_W->Draw("AP");
  
  TCanvas* c_pedsub_errorIntegral_NW = new TCanvas("c_pedsub_errorIntegral_NW","c_pedsub_errorIntegral_NW",900,900);
  c_pedsub_errorIntegral_NW->cd();
  gr_pedsub_errorIntegral_NW->GetYaxis()->SetTitle("%");
  gr_pedsub_errorIntegral_NW->GetXaxis()->SetTitle("nsp");
  gr_pedsub_errorIntegral_NW->GetXaxis()->SetRangeUser(-1,max_nsp+1);
  gr_pedsub_errorIntegral_NW->Draw("AP");
  TCanvas* c_pedsub_errorIntegral_W = new TCanvas("c_pedsub_errorIntegral_W","c_pedsub_errorIntegral_W",900,900);
  c_pedsub_errorIntegral_W->cd();
  gr_pedsub_errorIntegral_W->GetYaxis()->SetTitle("%");
  gr_pedsub_errorIntegral_W->GetXaxis()->SetTitle("nsp");
  gr_pedsub_errorIntegral_W->GetXaxis()->SetRangeUser(-1,max_nsp+1);
  gr_pedsub_errorIntegral_W->Draw("AP");
  
  if (supporting){
    TCanvas* c_nsp_errorIntegral_W[n_nsp];
    TCanvas* c_nsp_errorIntegral_NW[n_nsp];
    TCanvas* c_nsp_PedSub_errorIntegral_W[n_nsp];
    TCanvas* c_nsp_PedSub_errorIntegral_NW[n_nsp];
    for (int i=0; i<n_nsp; i++){
      TString t_nsp = TString::Itoa(i+min_nsp,10);
      min_nsp_errorIntegral[i] = min_nsp_errorIntegral[i]*0.9;
      max_nsp_errorIntegral[i] = max_nsp_errorIntegral[i]*1.05;
      min_nsp_PedSub_errorIntegral[i] = min_nsp_PedSub_errorIntegral[i]*0.9;
      max_nsp_PedSub_errorIntegral[i] = max_nsp_PedSub_errorIntegral[i]*1.05;
      
      c_nsp_errorIntegral_W[i] = new TCanvas("c_nsp_errorIntegral_W_"+t_nsp,"c_nsp_errorIntegral_W_"+t_nsp,900,900);
      c_nsp_errorIntegral_W[i]->cd();
      h_nsp_errorIntegral_W[i]->GetXaxis()->SetRangeUser(min_nsp_errorIntegral[i],max_nsp_errorIntegral[i]);
      h_nsp_errorIntegral_W[i]->Draw();
      
      c_nsp_errorIntegral_NW[i] = new TCanvas("c_nsp_errorIntegral_NW_"+t_nsp,"c_nsp_errorIntegral_NW_"+t_nsp,900,900);
      c_nsp_errorIntegral_NW[i]->cd();
      h_nsp_errorIntegral_NW[i]->GetXaxis()->SetRangeUser(min_nsp_errorIntegral[i],max_nsp_errorIntegral[i]);
      h_nsp_errorIntegral_NW[i]->Draw();
      
      c_nsp_PedSub_errorIntegral_W[i] = new TCanvas("c_nsp_PedSub_errorIntegral_W_"+t_nsp,"c_nsp_PedSub_errorIntegral_W_"+t_nsp,900,900);
      c_nsp_PedSub_errorIntegral_W[i]->cd();
      h_nsp_PedSub_errorIntegral_W[i]->GetXaxis()->SetRangeUser(min_nsp_PedSub_errorIntegral[i],max_nsp_PedSub_errorIntegral[i]);
      h_nsp_PedSub_errorIntegral_W[i]->Draw();
      
      c_nsp_PedSub_errorIntegral_NW[i] = new TCanvas("c_nsp_PedSub_errorIntegral_NW_"+t_nsp,"c_nsp_PedSub_errorIntegral_NW_"+t_nsp,900,900);
      c_nsp_PedSub_errorIntegral_NW[i]->cd();
      h_nsp_PedSub_errorIntegral_NW[i]->GetXaxis()->SetRangeUser(min_nsp_PedSub_errorIntegral[i],max_nsp_PedSub_errorIntegral[i]);
      h_nsp_PedSub_errorIntegral_NW[i]->Draw();
      
      if (save){
        TString dir = "plots/NSP/";
        c_nsp_errorIntegral_W[i]->SaveAs(dir+"Weight/IntegralError_NSP-"+t_nsp+".jpg");
        c_nsp_errorIntegral_NW[i]->SaveAs(dir+"NoWeight/IntegralError_NSP-"+t_nsp+".jpg");
        c_nsp_PedSub_errorIntegral_W[i]->SaveAs(dir+"Weight/PedSubIntegralError_NSP-"+t_nsp+".jpg");
        c_nsp_PedSub_errorIntegral_NW[i]->SaveAs(dir+"NoWeight/PedSubIntegralError_NSP-"+t_nsp+".jpg");
      }
    }
  }
  
  if (save){
    TString dir = "plots/Integral_VS_NSP/";
    c_errorIntegral_NW->SaveAs(dir+"IntegralVsNSP.pdf");
    c_errorIntegral_W ->SaveAs(dir+"IntegralVsNSP_weighted.pdf");
    c_pedsub_errorIntegral_NW->SaveAs(dir+"PedSubtracted_IntegralVsNSP.pdf");
    c_pedsub_errorIntegral_W ->SaveAs(dir+"PedSubtracted_IntegralVsNSP_weighted.pdf");
  }
  
}


// Pulse integral vs pulse peak as a function of nsa and nsp
void PulseAnalysis_PulseIntegralVsPulsePeak(TString filename, bool save=false){
  TFile * f = new TFile(filename);
  
  TTree* FCAL_Analysis    = (TTree*)f->Get("FCAL_Analysis");
  SetFCAL_Analysis(FCAL_Analysis);
  TTree* FCAL_Analysis2    = (TTree*)f->Get("FCAL_Analysis2");
  SetFCAL_Analysis2(FCAL_Analysis2);
  
  int nAna  = FCAL_Analysis->GetEntries();
  int nAna2 = FCAL_Analysis2->GetEntries();
  
  int n_nsa = maxNSA-minNSA+1;
  int n_nsp = max_nsp-min_nsp+1;
  TGraph* gr_IntegralVsNSA[n_nsa];
  TGraph* gr_IntegralVsNSP[n_nsp];
  TGraph* gr_IntegralErrorVsNSA[n_nsa];
  TGraph* gr_IntegralErrorVsNSP[n_nsp];
  TGraph* gr_PedSubIntegralErrorVsNSA[n_nsa];
  TGraph* gr_PedSubIntegralErrorVsNSP[n_nsp];
  
  for (int i=0; i<n_nsa; i++){
    TString t_nsa = TString::Itoa(minNSA+i,10);
    gr_IntegralVsNSA[i] = new TGraph();
    gr_IntegralVsNSA[i]->SetName("int_vs_peak_nsa-"+t_nsa);
    gr_IntegralVsNSA[i]->SetTitle("Integral vs Peak (nsa="+t_nsa+")");
    gr_IntegralErrorVsNSA[i] = new TGraph();
    gr_IntegralErrorVsNSA[i]->SetName("interror_vs_peak_nsa-"+t_nsa);
    gr_IntegralErrorVsNSA[i]->SetTitle("Integral rms/avg*100 vs Peak (nsa="+t_nsa+")");
    gr_PedSubIntegralErrorVsNSA[i] = new TGraph();
    gr_PedSubIntegralErrorVsNSA[i]->SetName("interrorPedSub_vs_peak_nsa-"+t_nsa);
    gr_PedSubIntegralErrorVsNSA[i]->SetTitle("Pedestal Subtracted Integral Rel. Error vs Peak (nsa="+t_nsa+")");
  }
  
  for (int i=0; i<n_nsp; i++){
    TString t_nsp = TString::Itoa(min_nsp+i,10);
    gr_IntegralVsNSP[i] = new TGraph();
    gr_IntegralVsNSP[i]->SetName("int_vs_peak_nsp-"+t_nsp);
    gr_IntegralVsNSP[i]->SetTitle("Integral vs Peak (nsp="+t_nsp+")");
    gr_IntegralErrorVsNSP[i] = new TGraph();
    gr_IntegralErrorVsNSP[i]->SetName("interror_vs_peak_nsp-"+t_nsp);
    gr_IntegralErrorVsNSP[i]->SetTitle("Integral rms/avg*100 vs Peak (nsp="+t_nsp+")");
    gr_PedSubIntegralErrorVsNSP[i] = new TGraph();
    gr_PedSubIntegralErrorVsNSP[i]->SetName("interrorPedSub_vs_peak_nsp-"+t_nsp);
    gr_PedSubIntegralErrorVsNSP[i]->SetTitle("Pedestal Subtracted Integral Rel. Error vs Peak (nsp="+t_nsp+")");
  }
  
  for (int i=0; i<nAna; i++){
    FCAL_Analysis->GetEntry(i);
    int index_nsa = nsa-minNSA;
    if (index_nsa>=n_nsa || index_nsa<0) continue;
    gr_IntegralVsNSA[index_nsa]->SetPoint(gr_IntegralVsNSA[index_nsa]->GetN(), avg_peak, avg_integral);
    //gr_IntegralErrorVsNSA[index_nsa]->SetPoint(gr_IntegralErrorVsNSA[index_nsa]->GetN(), avg_peak, rms_integral_gaus/avg_integral*100);
    gr_IntegralErrorVsNSA[index_nsa]->SetPoint(gr_IntegralErrorVsNSA[index_nsa]->GetN(), avg_peak, rms_integral/avg_integral*100);
    //gr_PedSubIntegralErrorVsNSA[index_nsa]->SetPoint(gr_IntegralErrorVsNSA[index_nsa]->GetN(), avg_peak, TMath::Sqrt(TMath::Power(rms_integral_gaus/(avg_integral-100.0*double(nsa+nsb)),2)+ TMath::Power(0.5*double(nsa+nsb)/100.0,2))*100);
    gr_PedSubIntegralErrorVsNSA[index_nsa]->SetPoint(gr_IntegralErrorVsNSA[index_nsa]->GetN(), avg_peak, TMath::Sqrt(TMath::Power(rms_integral/(avg_integral-100.0*double(nsa+nsb)),2)+ TMath::Power(0.5*double(nsa+nsb)/100.0,2))*100);
  }
  
  for (int i=0; i<nAna2; i++){
    FCAL_Analysis2->GetEntry(i);
    int index_nsp = nsp-min_nsp;
    if (index_nsp>=n_nsp || index_nsp<0) continue;
    gr_IntegralVsNSP[index_nsp]->SetPoint(gr_IntegralVsNSP[index_nsp]->GetN(), avg_peak, avg_integral);
    //gr_IntegralErrorVsNSP[index_nsp]->SetPoint(gr_IntegralErrorVsNSP[index_nsp]->GetN(), avg_peak, rms_integral_gaus/avg_integral*100);
    //gr_PedSubIntegralErrorVsNSP[index_nsp]->SetPoint(gr_IntegralErrorVsNSP[index_nsp]->GetN(), avg_peak, TMath::Sqrt(TMath::Power(rms_integral_gaus/(avg_integral-100.0*(2.0*double(nsp)+1.0)),2)+ TMath::Power(0.5*2*double(nsp)/100.0,2))*100);
    gr_IntegralErrorVsNSP[index_nsp]->SetPoint(gr_IntegralErrorVsNSP[index_nsp]->GetN(), avg_peak, rms_integral/avg_integral*100);
    gr_PedSubIntegralErrorVsNSP[index_nsp]->SetPoint(gr_IntegralErrorVsNSP[index_nsp]->GetN(), avg_peak, TMath::Sqrt(TMath::Power(rms_integral/(avg_integral-100.0*(2.0*double(nsp)+1.0)),2)+ TMath::Power(0.5*2*double(nsp)/100.0,2))*100);
  }

  TCanvas* c_IntegralVsNSA[n_nsa];
  TCanvas* c_IntegralVsNSP[n_nsp];
  TCanvas* c_IntegralErrorVsNSA[n_nsa];
  TCanvas* c_IntegralErrorVsNSP[n_nsp];
  TCanvas* c_PedSubIntegralErrorVsNSA[n_nsa];
  TCanvas* c_PedSubIntegralErrorVsNSP[n_nsp];
  for (int i=0; i<n_nsa; i++){
    c_IntegralVsNSA[i] = new TCanvas("c_IntegralVsNSA_"+TString::Itoa(i+minNSA,10),"c_IntegralVsNSA_"+TString::Itoa(i+minNSA,10),900,900);
    c_IntegralVsNSA[i]->cd();
    gr_IntegralVsNSA[i]->Draw("A");
    gr_IntegralVsNSA[i]->GetXaxis()->SetTitle("pulse peak");
    gr_IntegralVsNSA[i]->GetYaxis()->SetTitle("pulse integral");
    gr_IntegralVsNSA[i]->Draw("AP");
    c_IntegralErrorVsNSA[i] = new TCanvas("c_IntegralErrorVsNSA_"+TString::Itoa(i+minNSA,10),"c_IntegralErrorVsNSA_"+TString::Itoa(i+minNSA,10),900,900);
    c_IntegralErrorVsNSA[i]->cd();
    gr_IntegralErrorVsNSA[i]->Draw("A");
    gr_IntegralErrorVsNSA[i]->GetXaxis()->SetTitle("pulse peak");
    gr_IntegralErrorVsNSA[i]->GetYaxis()->SetTitle("relative error (%)");
    gr_IntegralErrorVsNSA[i]->Draw("AP");
    c_PedSubIntegralErrorVsNSA[i] = new TCanvas("c_PedSubIntegralErrorVsNSA_"+TString::Itoa(i+minNSA,10),"c_PedSubIntegralErrorVsNSA_"+TString::Itoa(i+minNSA,10),900,900);
    c_PedSubIntegralErrorVsNSA[i]->cd();
    gr_PedSubIntegralErrorVsNSA[i]->Draw("A");
    gr_PedSubIntegralErrorVsNSA[i]->GetXaxis()->SetTitle("pulse peak");
    gr_PedSubIntegralErrorVsNSA[i]->GetYaxis()->SetTitle("relative error (%)");
    gr_PedSubIntegralErrorVsNSA[i]->Draw("AP");
  }
  
  for (int i=0; i<n_nsp; i++){
    c_IntegralVsNSP[i] = new TCanvas("c_IntegralVsNSP_"+TString::Itoa(i+min_nsp,10),"c_IntegralVsNSP_"+TString::Itoa(i+min_nsp,10),900,900);
    c_IntegralVsNSP[i]->cd();
    gr_IntegralVsNSP[i]->Draw("A");
    gr_IntegralVsNSP[i]->GetXaxis()->SetTitle("pulse peak");
    gr_IntegralVsNSP[i]->GetYaxis()->SetTitle("pulse integral");
    gr_IntegralVsNSP[i]->Draw("AP");
    c_IntegralErrorVsNSP[i] = new TCanvas("c_IntegralErrorVsNSP_"+TString::Itoa(i+min_nsp,10),"c_IntegralErrorVsNSP_"+TString::Itoa(i+min_nsp,10),900,900);
    c_IntegralErrorVsNSP[i]->cd();
    gr_IntegralErrorVsNSP[i]->Draw("A");
    gr_IntegralErrorVsNSP[i]->GetXaxis()->SetTitle("pulse peak");
    gr_IntegralErrorVsNSP[i]->GetYaxis()->SetTitle("relative error (%)");
    gr_IntegralErrorVsNSP[i]->Draw("AP");
    c_PedSubIntegralErrorVsNSP[i] = new TCanvas("c_PedSubIntegralErrorVsNSP_"+TString::Itoa(i+min_nsp,10),"c_PedSubIntegralErrorVsNSP_"+TString::Itoa(i+min_nsp,10),900,900);
    c_PedSubIntegralErrorVsNSP[i]->cd();
    gr_PedSubIntegralErrorVsNSP[i]->Draw("A");
    gr_PedSubIntegralErrorVsNSP[i]->GetXaxis()->SetTitle("pulse peak");
    gr_PedSubIntegralErrorVsNSP[i]->GetYaxis()->SetTitle("relative error (%)");
    gr_PedSubIntegralErrorVsNSP[i]->Draw("AP");
  }
  
  if (save){
    TString dir = "plots/Integral_VS_Peak/";
    for (int i=0; i<n_nsa; i++){
      TString t_nsa = TString::Itoa(minNSA+i,10);
      c_IntegralVsNSA[i]->SaveAs(dir+"IntegralVsPeak_NSA-"+t_nsa+".jpg");
      c_IntegralErrorVsNSA[i]->SaveAs(dir+"IntegralErrorVsPeak_NSA-"+t_nsa+".jpg");
      c_PedSubIntegralErrorVsNSA[i]->SaveAs(dir+"PedSubIntegralErrorVsPeak_NSA-"+t_nsa+".jpg");
    }
    for (int i=0; i<n_nsp; i++){
      TString t_nsp = TString::Itoa(min_nsp+i,10);
      c_IntegralVsNSP[i]->SaveAs(dir+"IntegralVsPeak_NSP-"+t_nsp+".jpg");
      c_IntegralErrorVsNSP[i]->SaveAs(dir+"IntegralErrorVsPeak_NSP-"+t_nsp+".jpg");
      c_PedSubIntegralErrorVsNSP[i]->SaveAs(dir+"PedSubIntegralErrorVsPeak_NSP-"+t_nsp+".jpg");
    }
  }
  
}

