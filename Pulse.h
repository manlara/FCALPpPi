#if !defined(PULSE_H)
#define PULSE_H

#include <stdint.h>
#include <stdlib.h>
#include <map>
#include <iterator>
#include <vector>

#include "StringUtilities.h"

#include "TROOT.h"
#include <TString.h>
#include "TMath.h"
#include <TH1F.h>
#include <TF1.h>
#include "TMinuit.h"
#include <TVirtualFitter.h>
#include <TTree.h>
#include <TObjString.h>
#include <TObjArray.h>


class Pulse{
  
  public:
    // Construct class roc/slot/channel
    Pulse();
    Pulse(TString name);
    Pulse(int rocid, int slot, int channel, int nsa );
    ~Pulse();
    
    bool peakExists( void );
    bool integralExists( void );
    
    TString getName( void );
    int getRocid( void );
    int getSlot( void );
    int getChannel( void );
    int getNSA( void );
    int getNSP( void );
    
    int getEventnum( void );
    void setEventnum( int eventnum );
    
    int getRun( void );
    void setRun( int run );
    
    int getX( void );
    void setX( int x );
    
    int getY( void );
    void setY( int y );
    
    int getTC( void );
    void setTC( int TC );
    
    void setStart( int start );
    void setEnd( int end );
    int getStart( void );
    int getEnd( void );
    
    void addPulsePeak( double peak, int eventnum=0 );
    void addPulseIntegral( double integral );
    double getAvgPulsePeak( void );
    double getRmsPulsePeak( void );
    double getRmsPulsePeakFit( void );
    double getAvgPulseIntegral( void );
    double getRmsPulseIntegral( void );
    double getRmsPulseIntegralFit( void );
    
    vector<double> getVecPulsePeak( void );
    void clearVecPulsePeak( void );
    
    map<int,double> getMapPulsePeak( void );
    void clearMapPulsePeak( void );
    
    vector<double> getVecPulseIntegral( void );
    void clearVecPulseIntegral( void );
    
    vector<int> getVecEventnum( void );
    void clearVecEventnum( void );
    
    double calcPulsePeakCorrelation( map<int,double> pulse2_map );
    
    vector<double> calcPulsePeakRatio( map<int,double> pulse2_map, bool print=false );
    
    // overload == and !=
    friend bool operator== (Pulse &p1, Pulse &p2);
    friend bool operator!= (Pulse &p1, Pulse &p2);
    
  private:
    TString m_name;
    int m_rocid, m_slot, m_channel, m_nsa, m_nsp, m_eventnum, m_run, m_start, m_end;
    bool m_existsI, m_existsP;
    int m_x, m_y, m_TC;
    double m_peak, m_integral;
    //TH1F* m_hist_ped;
    vector<double> m_vec_peak;
    map<int,double> m_map_peak; // key=eventnum
    vector<double> m_vec_integral;
    vector<int> m_vec_eventnum;
  
};

#endif



