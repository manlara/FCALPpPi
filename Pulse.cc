#include "Pulse.h"

Pulse::~Pulse(){
  
}
Pulse::Pulse(){
  
}
Pulse::Pulse( TString name )
{
  m_existsP = false;
  m_existsI = false;
  //cout << "Initializing Constructor, by name" << endl;
  m_name = name;
  // crate/slot/channel_nsa
  vector<TString> parse_key = StringUtilities::parseTString(m_name,"_");
  vector<TString> parse_name = StringUtilities::parseTString(parse_key[0],"/");
  for (int i=0; i<int(parse_name.size()); i++){
    if (i==0) m_rocid   = StringUtilities::TString2int(parse_name[i]);
    if (i==1) m_slot    = StringUtilities::TString2int(parse_name[i]);
    if (i==2) m_channel = StringUtilities::TString2int(parse_name[i]);
  }
  if (int(parse_key.size())==2){
    m_nsa = StringUtilities::TString2int(parse_key[1]);
  }
}
Pulse::Pulse( int rocid, int slot, int channel, int nsa )
{
  m_existsP = false;
  m_existsI = false;
  //cout << "Initializing Constructor, by roc/slot/chan" << endl;
  char name[50];
  sprintf(name, "%02i/%02i/%02i", m_rocid, m_slot, m_channel);
  m_name = name;
  m_rocid = rocid;
  m_slot = slot;
  m_channel = channel;
  m_nsa = nsa;
}

TString Pulse::getName( void ){ return m_name;    }

int Pulse::getRocid( void )   { return m_rocid;   }
int Pulse::getSlot( void )    { return m_slot;    }
int Pulse::getChannel( void ) { return m_channel; }
int Pulse::getNSA( void )     { return m_nsa; }
int Pulse::getNSP( void )     { return m_nsa; }

int Pulse::getStart( void )       { return m_start; }
int Pulse::getEnd( void )         { return m_end; }
void Pulse::setStart( int start ) { m_start = start; }
void Pulse::setEnd( int end )     { m_end = end; }

void Pulse::addPulsePeak( double peak, int eventnum ) { 
  m_existsP = true;
  m_vec_peak.push_back(peak);
  m_map_peak[eventnum] = peak;
}
double Pulse::getAvgPulsePeak( void ){ 
  int n = int(m_vec_peak.size());
  if (n==0) return 0.0;
  return TMath::Mean(n,&m_vec_peak[0]);
}
double Pulse::getRmsPulsePeak( void ){ 
  int n = int(m_vec_peak.size());
  if (n==0) return 0.0;
  return TMath::RMS(n,&m_vec_peak[0]);
}
double Pulse::getRmsPulsePeakFit( void ){
  // convert vector into histogram
  int n = int(m_vec_peak.size());
  if (n==0) return 0.0;
  TH1::AddDirectory(kFALSE);
  TH1F* hist = new TH1F("hist","hist",4096, 0, 4096);
  for (int i=0; i<n; i++){
    hist->Fill(m_vec_peak[i]);
  }
  hist->Fit("gaus","0Q");
  TF1* fit = hist->GetFunction("gaus"); // par 1=mean, par 2=sigma
  double answer = fit->GetParameter(2);
  delete hist;
  return answer;
}

void Pulse::addPulseIntegral( double integral ) { 
  m_existsI = true;
  m_vec_integral.push_back(integral);
}
double Pulse::getAvgPulseIntegral( void ){ 
  int n = int(m_vec_integral.size());
  if (n==0) return 0.0;
  return TMath::Mean(n,&m_vec_integral[0]);
}
double Pulse::getRmsPulseIntegral( void ){ 
  int n = int(m_vec_integral.size());
  if (n==0) return 0.0;
  return TMath::RMS(n,&m_vec_integral[0]);
}
double Pulse::getRmsPulseIntegralFit( void ){ 
  // convert vector into histogram
  int n = int(m_vec_integral.size());
  if (n==0) return 0.0;
  TH1::AddDirectory(kFALSE);
  TH1F* hist = new TH1F("hist","hist",4096*100, 0, 4096*100);
  for (int i=0; i<n; i++){
    hist->Fill(m_vec_integral[i]);
  }
  hist->Fit("gaus","0Q");
  TF1* fit = hist->GetFunction("gaus"); // par 1=mean, par 2=sigma
  double answer = fit->GetParameter(2);
  delete hist;
  return answer;
}

vector<double> Pulse::getVecPulsePeak( void ){ return m_vec_peak; }
void Pulse::clearVecPulsePeak( void ){ m_vec_peak.clear(); }

map<int,double> Pulse::getMapPulsePeak( void ){ return m_map_peak; }
void Pulse::clearMapPulsePeak( void ){ m_map_peak.clear(); }

vector<double> Pulse::getVecPulseIntegral( void ){ return m_vec_integral; }
void Pulse::clearVecPulseIntegral( void ){ m_vec_integral.clear(); }

vector<int> Pulse::getVecEventnum( void ){ return m_vec_eventnum; }
void Pulse::clearVecEventnum( void ){ m_vec_eventnum.clear(); }

void Pulse::setEventnum( int eventnum ){ m_eventnum = eventnum; m_vec_eventnum.push_back(eventnum);}
int Pulse::getEventnum( void ){return m_eventnum;}

void Pulse::setRun( int run ){m_run = run;}
int Pulse::getRun( void ){return m_run;}

void Pulse::setX( int x ){m_x = x;}
int Pulse::getX( void ){return m_x;}

void Pulse::setY( int y ){m_y = y;}
int Pulse::getY( void ){return m_y;}

void Pulse::setTC( int TC ){m_TC = TC;}
int Pulse::getTC( void ){return m_TC;}

bool Pulse::peakExists( void ){ return m_existsP;}
bool Pulse::integralExists( void ){ return m_existsI;}

// This calculated by looking at the orrelation between this channel
// and the input channel
// Will probably need to make sure you don't calculate the correlation between or on an empty channel
// **************************
// Correlation:
// cov(P1,P2) / (std(P1) std(P2)) = sum_Y sum_X (Xi-X)/Nx (Yi-Y)/Ny
double Pulse::calcPulsePeakCorrelation( map<int,double> pulse2_map ){
  
  vector<double> pulse1;
  vector<double> pulse2;
  
  for (map<int,double>::const_iterator iter=m_map_peak.begin(); iter!=m_map_peak.end(); ++iter){
    int event = iter->first;
    if (pulse2_map.count(event)>0){
      pulse1.push_back(iter->second);
      pulse2.push_back(pulse2_map[event]);
    }
    
  }
  
  int n_pulse1 = pulse1.size();
  int n_pulse2 = pulse2.size();
  if (n_pulse1==0 || n_pulse2==0) {
    cout << "No elements in either pulse1 or pulse2... \nPulse1 = " << n_pulse1 << "\tPulse2 = " << n_pulse2 << endl;
    return 101;
  }
  if (n_pulse1!=n_pulse2){
    cout << "Pulses not the same size:\n\tPulse1: " << n_pulse1 << "\tPulse2: " << n_pulse2 << endl;
    return 101;
  }
  // pedestal subtract
  for (int i=0; i<n_pulse1; i++){
    pulse1[i] = pulse1[i]-100.0;
    pulse2[i] = pulse2[i]-100.0;
  }
  
  double avg_pulse1 = TMath::Mean(n_pulse1, &pulse1[0]);
  double avg_pulse2 = TMath::Mean(n_pulse2, &pulse2[0]);
  double cov = 0;
  double std_pulse1 = 0;
  double std_pulse2 = 0;
  for (int i=0; i<n_pulse1; i++){
    std_pulse1 += (pulse1[i]-avg_pulse1)*(pulse1[i]-avg_pulse1);
    std_pulse2 += (pulse2[i]-avg_pulse2)*(pulse2[i]-avg_pulse2);
    cov += (pulse1[i]-avg_pulse1) * (pulse2[i]-avg_pulse2);
  }
  cov *= 1.0 / double(n_pulse1);
  std_pulse1 = TMath::Sqrt(std_pulse1/double(n_pulse1));
  std_pulse2 = TMath::Sqrt(std_pulse2/double(n_pulse2));
  return cov / (std_pulse1 * std_pulse2);
}

vector<double> Pulse::calcPulsePeakRatio( map<int,double> pulse2_map, bool print ){
  
  vector<double> vec_ratio;
  
  for (map<int,double>::const_iterator iter=m_map_peak.begin(); iter!=m_map_peak.end(); ++iter){
    int event = iter->first;
    
    if (pulse2_map.count(event)>0){
      vec_ratio.push_back( (m_map_peak[event]-100.0)/(pulse2_map[event]-100.0) ); 
    }
    
  }
  vector<double> vec_mean_rms;
  vec_mean_rms.push_back(TMath::Mean(vec_ratio.size(), &vec_ratio[0]));
  vec_mean_rms.push_back(TMath::RMS(vec_ratio.size(),  &vec_ratio[0]));
  
  // convert vector into histogram
  TH1::AddDirectory(kFALSE);
  TH1F* hist = new TH1F("hist","hist",4096, 0, 4096);
  for (int i=0; i<int(vec_ratio.size()); i++){
    hist->Fill(vec_ratio[i]);
  }
  hist->Fit("gaus","0Q");
  TF1* fit = hist->GetFunction("gaus"); // par 1=mean, par 2=sigma
  double answer = fit->GetParameter(2);
  delete hist;
  vec_mean_rms.push_back(answer);
  
  return vec_mean_rms;
}

bool operator== (Pulse &p1, Pulse &p2){
  return (p1.m_rocid == p2.m_rocid && p1.m_slot == p2.m_slot && p1.m_channel == p2.m_channel && p1.m_eventnum == p2.m_eventnum && p1.m_run == p2.m_run && (p1.m_nsa == p2.m_nsa || p1.m_nsp == p2.m_nsp));
}
bool operator!= (Pulse &p1, Pulse &p2){
  return !(p1 == p2);
}
