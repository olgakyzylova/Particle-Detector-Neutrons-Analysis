/******************************************************
 * study the neutron number change with the reactor 
 * power variation.
 * Include spatial distribution
 * Old code; can be improved
 ******************************************************/
#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TObject.h>

int getGlobalChannel(int run, int boardID, int adccha) {
  int globalChannel = -1; 
  if (boardID == 1) globalChannel = 16 + adccha;
  else {
    if (adccha<=15) globalChannel = adccha; 
    else globalChannel = 48 + adccha;
  }  
  return globalChannel;
}

// Setup a created tree
struct structNeutronEvent {
  Double_t cubeX;
  Double_t cubeY;
  Double_t cubeZ;
  Double_t time;
  Double_t chi2n;
  Double_t chi2gamma;
};

void setupNeutronTreeFill(TTree &TNeu, structNeutronEvent& iNeu){
  TNeu.Branch("cubeX", &iNeu.cubeX, "cubeX/D");
  TNeu.Branch("cubeY", &iNeu.cubeY, "cubeY/D");
  TNeu.Branch("cubeZ", &iNeu.cubeZ, "cubeZ/D");
  TNeu.Branch("time", &iNeu.time, "time/D");
  TNeu.Branch("chi2n", &iNeu.chi2n, "chi2n/D");
  TNeu.Branch("chi2gamma", &iNeu.chi2gamma, "chi2gamma/D");
}


//chi-square function - a void function that does not return anything but changes the values of chi[2] from 0 to calculated ones
void get_chi(short int *adcval, double baseline, double chi[2]){
  double n_tem[8] = {0.38, 0.20, 0.12, 0.08, 0.07, 0.06, 0.05, 0.04}; //template
  double err_n_tem[8] = {0.03, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01}; //template error
  double g_tem[8] = {1};
  int start_pos = 0;
  for (int j=0; j<129; j++){
    if(adcval[j] > 200 + 12){ //hardcode 200 as baseline
      start_pos = j;
      break;
    }
  }

  //calculate baseline if it is zero
  if(baseline == 0){
    for (int j=0; j<30; j++){
      baseline += adcval[j];
    }
    baseline = baseline/30;
    for (int j=30; j<120; j++){
      if (adcval[j] > baseline + 12){
        start_pos = j;
        break;
      }
    }
  }

  double n_chi = 0;
  double g_chi = 0;
  double data[12] = {0};

  int k = 0;
  double sum_data = 0;
  for(int j=0; j+start_pos<129; j++){
    data[k] += adcval[j+start_pos] - baseline;
    sum_data += adcval[j+start_pos] - baseline;
    if(j>1 && j%11 == 0) k+= 1;
  }

  double norm_data[8] = {0};
  double error[8] = {0};
  int ndf = 0;
  for (int j=0; j<8; j++){
    if (data[j] == 0) continue;
    ++ndf;
    norm_data[j] = data[j]/sum_data;
    error[j] = sqrt(6*abs(data[j]))/sum_data;
    n_chi += (norm_data[j] - n_tem[j])*(norm_data[j] - n_tem[j])/(pow((error[j]+err_n_tem[j])/2,2));
    g_chi += (norm_data[j] - g_tem[j])*(norm_data[j] - g_tem[j])/(pow((error[j]+err_n_tem[j])/2,2));
  }
  if (ndf > 2){
    g_chi = g_chi/ndf;
    n_chi = n_chi/ndf;
  }
  else {
    g_chi = 0;
    n_chi = 0;
  }
  chi[0] = n_chi;
  chi[1] = g_chi;
  if (start_pos < 30 || start_pos > 60){
    chi[0] = 0;
    chi[1] = 0;
  }
}

int main(int argc, char *argv[]) {

  int run {0};
  int last_run {0};
  if (argc<3) {
    std::cout << "Usage: " << argv[0] << " [run]" << "[last run]" << std::endl;
    return -1;
  }
  sscanf(argv[1], "%i", &run);
  sscanf(argv[2], "%i", &last_run);
  std::cout<<"run: "<< run << ", last_run: "<< last_run << std::endl;

//  TH1D *h_wf = new TH1D("h_wf", "Neutron Waveform", 129, 0, 129);
//  TH1D *h_adc1 = new TH1D("h_adc1", "ADC1", 25, 200, 225);
//  TH1D *h_adc2 = new TH1D("h_adc2", "ADC2", 25, 200, 225);

  int ADCSIZE = 129;
  int n_threshold = 14;
  int gl_threshold = 30;
  int single_g_thresh = 15;
  int single_n_thresh = 20;

  std::vector <int> x_run, y_neutron;
  std::vector<std::pair<double,double>> vetoTime011, vetoTime021, vetoTime031, vetoTime041, vetoTime051, vetoTime061, vetoTime071, vetoTime081;
  std::vector<std::pair<double,double>> vetoTime101, vetoTime201, vetoTime301, vetoTime401, vetoTime501, vetoTime601, vetoTime701, vetoTime801;
  std::vector<std::pair<double,double>> vetoTime012, vetoTime022, vetoTime032, vetoTime042, vetoTime052, vetoTime062, vetoTime072, vetoTime082;
  std::vector<std::pair<double,double>> vetoTime102, vetoTime202, vetoTime302, vetoTime402, vetoTime502, vetoTime602, vetoTime702, vetoTime802;
  std::vector<std::pair<double,double>> vetoTime013, vetoTime023, vetoTime033, vetoTime043, vetoTime053, vetoTime063, vetoTime073, vetoTime083;
  std::vector<std::pair<double,double>> vetoTime103, vetoTime203, vetoTime303, vetoTime403, vetoTime503, vetoTime603, vetoTime703, vetoTime803;
  std::vector<std::pair<double,double>> vetoTime014, vetoTime024, vetoTime034, vetoTime044, vetoTime054, vetoTime064, vetoTime074, vetoTime084;
  std::vector<std::pair<double,double>> vetoTime104, vetoTime204, vetoTime304, vetoTime404, vetoTime504, vetoTime604, vetoTime704, vetoTime804;
  std::vector<std::pair<double,double>> vetoTime015, vetoTime025, vetoTime035, vetoTime045, vetoTime055, vetoTime065, vetoTime075, vetoTime085;
  std::vector<std::pair<double,double>> vetoTime105, vetoTime205, vetoTime305, vetoTime405, vetoTime505, vetoTime605, vetoTime705, vetoTime805;

//  std::ofstream myfile;
//  myfile.open("Vetoed_waveforms.csv", std::ios_base::app);

  for (int iRun = run; iRun <= last_run ; iRun++){
    std::cout <<"run: "<< iRun << std::endl;
    int n_neutrons = 0;

    TFile *f;
    if(gSystem->AccessPathName(Form("/e/h.1/Chandler_BL_corr_data_v3/run%i/run%i_bl_corr_merged_evt_v3.root", iRun, iRun))) continue;
    else f = new TFile(Form("/e/h.1/Chandler_BL_corr_data_v3/run%i/run%i_bl_corr_merged_evt_v3.root", iRun, iRun));

    int globalChannelToIndex[80]; 

    // channel map
    int channel2cubeX[80];
    int channel2cubeY[80];
    int channel2cubeZ[80];

    for (int i=0; i<80; i++) {
      channel2cubeX[i] = -1;
      channel2cubeY[i] = -1;
      channel2cubeZ[i] = -1;
    }

    TTree *t_ch_map = (TTree *) f->Get("ch_map");

    {
      // Declaration of leaf types
      UChar_t         nChannels;
      UChar_t         boardID[80];   //[nChannels]
      UChar_t         adccha[80];   //[nChannels]
      Char_t          _cubeX[80];   //[nChannels]
      Char_t          _cubeY[80];   //[nChannels]
      Char_t          _cubeZ[80];   //[nChannels]

      // List of branches
      TBranch        *b_nChannels;   //!
      TBranch        *b_boardID;   //!
      TBranch        *b_adccha;   //!
      TBranch        *b_cubeX;   //!
      TBranch        *b_cubeY;   //!
      TBranch        *b_cubeZ;   //!

      t_ch_map->SetBranchAddress("nChannels", &nChannels, &b_nChannels);
      t_ch_map->SetBranchAddress("boardID", boardID, &b_boardID);
      t_ch_map->SetBranchAddress("adccha", adccha, &b_adccha);
      t_ch_map->SetBranchAddress("cubeX", _cubeX, &b_cubeX);
      t_ch_map->SetBranchAddress("cubeY", _cubeY, &b_cubeY);
      t_ch_map->SetBranchAddress("cubeZ", _cubeZ, &b_cubeZ);

      b_nChannels->GetEntry(0);
      b_boardID->GetEntry(0);
      b_adccha->GetEntry(0);
      b_cubeX->GetEntry(0);
      b_cubeY->GetEntry(0);
      b_cubeZ->GetEntry(0);

      for (int i=0; i<(int) nChannels; i++) {
        int globalChannel = getGlobalChannel(iRun, (int) boardID[i], (int) adccha[i]);
        channel2cubeX[globalChannel] = (int) _cubeX[i];
        channel2cubeY[globalChannel] = (int) _cubeY[i];
        channel2cubeZ[globalChannel] = (int) _cubeZ[i];
      }
    }

    TTree *t = (TTree *) f->Get("wf");

    // Declaration of leaf types
    Int_t           eventID;
    UChar_t         nWaveforms;
    UChar_t         boardID[80];   //[nWaveforms]
    UChar_t         adccha[80];   //[nWaveforms]
    Short_t         adcval[80][129];   //[nWaveforms]
    Float_t         baseline[80];   //[nWaveforms]
    Float_t         baseline_RMS[80];   //[nWaveforms]
    Short_t         trigBin[80];   //[nWaveforms]
    Float_t         pulseH[80];   //[nWaveforms]
    Float_t         area[80];   //[nWaveforms]
    Float_t         PSD[80];   //[nWaveforms]
    UChar_t         nBoards;
    UChar_t         board[2];   //[nBoards]
    ULong64_t       adcTime[2];   //[nBoards]
    ULong64_t       globalTime[2];   //[nBoards]
    Int_t           NofFlip[2];   //[nBoards]
    Int_t           board_time_offset[2];   //[nBoards]
    Short_t         trig_pattern[2];   //[nBoards]
    Short_t         adcSize[2];   //[nBoards]
    Int_t           old_eventID[2];   //[nBoards]
    Bool_t          badBaseline[80];//[nWaveforms]
    Bool_t          badCorrBaseline[80];//[nWaveforms]
    Bool_t          anySaturation;

    // List of branches
    TBranch        *b_eventID;
    TBranch        *b_nWaveforms;
    TBranch        *b_boardID;
    TBranch        *b_adccha;
    TBranch        *b_adcval;
    TBranch        *b_baseline;
    TBranch        *b_baseline_RMS;
    TBranch        *b_trigBin;
    TBranch        *b_pulseH;
    TBranch        *b_area;
    TBranch        *b_PSD;
    TBranch        *b_nBoards;
    TBranch        *b_board;
    TBranch        *b_adcTime;
    TBranch        *b_globalTime;
    TBranch        *b_NofFlip;
    TBranch        *b_board_time_offset;
    TBranch        *b_trig_pattern;
    TBranch        *b_adcSize;
    TBranch        *b_old_eventID;
    TBranch        *b_badBaseline;
    TBranch        *b_badCorrBaseline;
    TBranch        *b_anySaturation;

    t->SetBranchAddress("eventID", &eventID, &b_eventID);
    t->SetBranchAddress("nWaveforms", &nWaveforms, &b_nWaveforms);
    t->SetBranchAddress("boardID", boardID, &b_boardID);
    t->SetBranchAddress("adccha", adccha, &b_adccha);
    t->SetBranchAddress("adcval", adcval, &b_adcval);
    t->SetBranchAddress("baseline", baseline, &b_baseline);
    t->SetBranchAddress("baseline_RMS", baseline_RMS, &b_baseline_RMS);
    t->SetBranchAddress("trigBin", trigBin, &b_trigBin);
    t->SetBranchAddress("pulseH", pulseH, &b_pulseH);
    t->SetBranchAddress("area", area, &b_area);
    t->SetBranchAddress("PSD", PSD, &b_PSD);
    t->SetBranchAddress("nBoards", &nBoards, &b_nBoards);
    t->SetBranchAddress("board", board, &b_board);
    t->SetBranchAddress("adcTime", adcTime, &b_adcTime);
    t->SetBranchAddress("globalTime", globalTime, &b_globalTime);
    t->SetBranchAddress("NofFlip", NofFlip, &b_NofFlip);
    t->SetBranchAddress("board_time_offset", board_time_offset, &b_board_time_offset);
    t->SetBranchAddress("trig_pattern", trig_pattern, &b_trig_pattern);
    t->SetBranchAddress("adcSize", adcSize, &b_adcSize);
    t->SetBranchAddress("old_eventID", old_eventID, &b_old_eventID);
    t->SetBranchAddress("badBaseline", &badBaseline, &b_badBaseline);
    t->SetBranchAddress("badCorrBaseline", &badCorrBaseline, &b_badCorrBaseline);
    t->SetBranchAddress("anySaturation", &anySaturation, &b_anySaturation);
    
    std::vector<double> vEventID, vCubeX, vCubeY, vCubeZ, vTime, vChi2n, vChi2gamma;

// ++++++++++++++++++++++++++++ Looping over entries 1. Finding the waveforms to veto and creating the database ++++++++++++++++++++++++++++++++
    for (int iEnt=0; iEnt<t->GetEntries(); iEnt++) {
//      if(iEnt > 700000) break;
      Long64_t tentry = t->LoadTree(iEnt);
      b_nWaveforms->GetEntry(tentry);
      b_boardID->GetEntry(tentry);
      b_adccha->GetEntry(tentry);
      b_adcval->GetEntry(tentry);
      b_pulseH->GetEntry(tentry);
      b_nBoards->GetEntry(tentry);
      b_globalTime->GetEntry(tentry);
      b_trig_pattern->GetEntry(tentry);
      b_adcSize->GetEntry(tentry);
      double time = 0;
      // Is it a shorter adcSize
      bool one_of_boards_shorter_adcSize = false;
      for (int i=0; i<(int) nBoards; i++) {
        if (adcSize[i]!=ADCSIZE) one_of_boards_shorter_adcSize = true;
      }
      // Is it an external trigger
      bool one_of_boards_external_trigger = false;
      for (int i=0; i<(int) nBoards; i++) {
        if ((trig_pattern[i] & 0x200) == 0x200) one_of_boards_external_trigger = true;
      }
      if (one_of_boards_external_trigger) continue;

      // +++++++++++++ Looping over waveforms I. Finding the waveforms to veto and creating the database +++++++++++++++++++++++++++++++++++++++
      for(int i=0; i<(int) nWaveforms ; i++){  // loop over waveforms
        int globalChannel = getGlobalChannel(run, (int) boardID[i], (int) adccha[i]);
        int cubeX = channel2cubeX[globalChannel];
        int cubeY = channel2cubeY[globalChannel];
        int cubeZ = channel2cubeZ[globalChannel];
        time = globalTime[0]*8/1000; // time in microseconds from the beginning of the run

        // Vetoing the time around the channel with pulse height > 500 ADC
        if (pulseH[i] > 500) {
          if (cubeX==-1) {
            if (cubeY==1) {
              if (cubeZ==1) vetoTime011.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime012.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime013.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime014.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime015.push_back(std::make_pair(time, time+40));
            }
            if (cubeY==2) {
              if (cubeZ==1) vetoTime021.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime022.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime023.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime024.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime025.push_back(std::make_pair(time, time+40));
            }
            if (cubeY==3) {
              if (cubeZ==1) vetoTime031.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime032.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime033.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime034.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime035.push_back(std::make_pair(time, time+40));
            }
            if (cubeY==4) {
              if (cubeZ==1) vetoTime041.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime042.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime043.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime044.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime045.push_back(std::make_pair(time, time+40));
            }
            if (cubeY==5) {
              if (cubeZ==1) vetoTime051.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime052.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime053.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime054.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime055.push_back(std::make_pair(time, time+40));
            }
            if (cubeY==6) {
              if (cubeZ==1) vetoTime061.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime062.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime063.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime064.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime065.push_back(std::make_pair(time, time+40));
            }
            if (cubeY==7) {
              if (cubeZ==1) vetoTime071.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime072.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime073.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime074.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime075.push_back(std::make_pair(time, time+40));
            }
            if (cubeY==8) {
              if (cubeZ==1) vetoTime081.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime082.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime083.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime084.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime085.push_back(std::make_pair(time, time+40));
            }
          }
          if (cubeY==-1) {
            if (cubeX==1) {
              if (cubeZ==1) vetoTime101.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime102.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime103.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime104.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime105.push_back(std::make_pair(time, time+40));
            }
            if (cubeX==2) {
              if (cubeZ==1) vetoTime201.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime202.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime203.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime204.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime205.push_back(std::make_pair(time, time+40));
            }
            if (cubeX==3) {
              if (cubeZ==1) vetoTime301.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime302.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime303.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime304.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime305.push_back(std::make_pair(time, time+40));
            }
            if (cubeX==4) {
              if (cubeZ==1) vetoTime401.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime402.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime403.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime404.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime405.push_back(std::make_pair(time, time+40));
            }
            if (cubeX==5) {
              if (cubeZ==1) vetoTime501.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime502.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime503.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime504.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime505.push_back(std::make_pair(time, time+40));
            }
            if (cubeX==6) {
              if (cubeZ==1) vetoTime601.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime602.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime603.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime604.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime605.push_back(std::make_pair(time, time+40));
            }
            if (cubeX==7) {
              if (cubeZ==1) vetoTime701.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime702.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime703.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime704.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime705.push_back(std::make_pair(time, time+40));
            }
            if (cubeX==8) {
              if (cubeZ==1) vetoTime801.push_back(std::make_pair(time, time+40));
              if (cubeZ==2) vetoTime802.push_back(std::make_pair(time, time+40));
              if (cubeZ==3) vetoTime803.push_back(std::make_pair(time, time+40));
              if (cubeZ==4) vetoTime804.push_back(std::make_pair(time, time+40));
              if (cubeZ==5) vetoTime805.push_back(std::make_pair(time, time+40));
            }
          }
        }
      }
    }

    std::cout << "We finished here	"<< std::endl;

    // ++++++++++++++++++++++++++++ Looping over entries for neutron selection ++++++++++++++++++++++++++++++++
    for (int iEnt=0; iEnt<t->GetEntries(); iEnt++) {
//      if(iEnt > 700000) break;

      Long64_t tentry = t->LoadTree(iEnt);

      b_eventID->GetEntry(tentry);
      b_nWaveforms->GetEntry(tentry);
      b_boardID->GetEntry(tentry);
      b_adccha->GetEntry(tentry);
      b_adcval->GetEntry(tentry);
      b_baseline->GetEntry(tentry);
      b_pulseH->GetEntry(tentry);
      b_area->GetEntry(tentry);
      b_PSD->GetEntry(tentry);
      b_nBoards->GetEntry(tentry);
      b_board->GetEntry(tentry);
      b_adcTime->GetEntry(tentry);
      b_globalTime->GetEntry(tentry);
      b_trig_pattern->GetEntry(tentry);
      b_adcSize->GetEntry(tentry);
      b_old_eventID->GetEntry(tentry);
      b_badBaseline->GetEntry(tentry);
      b_badCorrBaseline->GetEntry(tentry);
      b_anySaturation->GetEntry(tentry);

      // Is it a shorter adcSize
      bool one_of_boards_shorter_adcSize = false;
      for (int i=0; i<(int) nBoards; i++) {
        if (adcSize[i]!=ADCSIZE) one_of_boards_shorter_adcSize = true;
      }

      // Is it an external trigger
      bool one_of_boards_external_trigger = false;
      for (int i=0; i<(int) nBoards; i++) {
        if ((trig_pattern[i] & 0x200) == 0x200) one_of_boards_external_trigger = true;
      }

      if (one_of_boards_external_trigger) continue;

      std::vector<std::array<int, 129>> neutrons;
      bool doesXexist = false;
      bool doesYexist = false;
      double total_nchi = 0;
      double total_gchi = 0;
      double highest_PH = 0;
      double highest_PHnchi = -1;
      double highest_PHgchi = -1;

      double max_cubeX = 0;
      double max_cubeY = 0;
      double cubeX_PH = 0;
      double cubeY_PH = 0;
      double time = 0;

      //get neutron like waveforms and save in neutrons vector
      gStyle->SetOptStat(0);

      bool cubeZone = true;
      bool cubeZtwo = true;
      bool cubeZthree = true;
      bool cubeZfour = true;
      bool cubeZfive = true;

      bool sheetZonehalf = true;
      bool sheetZtwohalf = true;
      bool sheetZthreehalf = true;
      bool sheetZfourhalf = true;

      bool studyFlag = false;

      short int add_neutrons_1[129] = {0};
      short int add_neutrons_2[129] = {0};
      short int add_neutrons_3[129] = {0};
      short int add_neutrons_4[129] = {0};
      short int add_neutrons_5[129] = {0};

      // Implementing trigger threshold for boards 1 and 2
      bool trigger_threshold_1 = false;
      bool trigger_threshold_2 = false;

      for(int i=0; i<(int) nWaveforms ; i++){  // loop over waveforms
        double adc_threshold = pulseH[i] + baseline[i]; // maximum adc value in a waveform
        if ((int) boardID[i] == 1) {
          if (adc_threshold > 214) trigger_threshold_1 = true; // trigger for board1
        }
        if ((int) boardID[i] == 2) {
          if (adc_threshold > 214) trigger_threshold_2 = true; // trigger for board2
        }
      }

      if (!(trigger_threshold_1) && !(trigger_threshold_2)) continue; // if neither of the boards is triggered

// ++++++++++++++++++++++++++++ Looping over waveforms 2. Selecting neutrons ++++++++++++++++++++++++++++++++

      for(int i=0; i<(int) nWaveforms ; i++){  // loop over waveforms

        if ((int) boardID[i] == 1 && !(trigger_threshold_1)) continue; // if this wf is from board1, and board1 was not triggered, skip
        if ((int) boardID[i] == 2 && !(trigger_threshold_2)) continue; // if this wf is from board2, and board2 was not triggered, skip

        // Continuing only with the waveforms that were on triggered board
        double adc_threshold = pulseH[i] + baseline[i];
        if (adc_threshold < 212) continue; // Zero suppression threshold

        int globalChannel = getGlobalChannel(run, (int) boardID[i], (int) adccha[i]);        
//        if ((int) boardID[i] == 1) h_adc1->Fill(pulseH[i]+baseline[i]); 
//        if ((int) boardID[i] == 2) h_adc2->Fill(pulseH[i]+baseline[i]);        
        
        int cubeX = channel2cubeX[globalChannel];
        int cubeY = channel2cubeY[globalChannel];
        int cubeZ = channel2cubeZ[globalChannel];
        time = globalTime[0]*8/1000; // time in microseconds from the beginning of the run
        
        std::array<int, 129> n_array = {0};
        
        // Conditions to assign correct Z for sheet and cube neutrons        
        if (cubeZ != 1) cubeZone = false;
        if (cubeZ != 2) cubeZtwo = false;
        if (cubeZ != 3) cubeZthree = false;
        if (cubeZ != 4) cubeZfour = false;
        if (cubeZ != 5) cubeZfive = false;
        
        if (cubeZ != 1 && cubeZ != 2) sheetZonehalf = false;
        if (cubeZ != 2 && cubeZ != 3) sheetZtwohalf = false;
        if (cubeZ != 3 && cubeZ != 4) sheetZthreehalf = false;
        if (cubeZ != 4 && cubeZ != 5) sheetZfourhalf = false;
        
        double chi[2] = {0};  // Set up the chi_n and chi_g values to 0
        get_chi(adcval[i], baseline[i], chi); // Calculate chis for this waveform and fill the chi variable
        double n_chi = chi[0]; // Now this are double single values
        double g_chi = chi[1];
        
        if (highest_PH < pulseH[i]){ // Finding the highest pulse among looped waveforms
          highest_PH = pulseH[i];
          highest_PHnchi = n_chi; // Assigning the chi2s of the highest waveform in the cluster of the event
          highest_PHgchi = g_chi;
        }

        // Locating X and Y of the event. Based on the maximum pulse of neutron-passing waveforms only
        if (!(g_chi > single_g_thresh && n_chi < single_n_thresh )) continue;

        // Checking if this waveform is in muted (vetoed) PMT channel
        bool tag = false;
        if (cubeX==-1) {
          if (cubeY==1) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime011.size(); m++) {
                if (vetoTime011.at(m).first<time && time<vetoTime011.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime012.size(); m++) {
                if (vetoTime012.at(m).first<time && time<vetoTime012.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime013.size(); m++) {
                if (vetoTime013.at(m).first<time && time<vetoTime013.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime014.size(); m++) {
                if (vetoTime014.at(m).first<time && time<vetoTime014.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime015.size(); m++) {
                if (vetoTime015.at(m).first<time && time<vetoTime015.at(m).second) tag = true;
              }
            }
          }
          if (cubeY==2) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime021.size(); m++) {
                if (vetoTime021.at(m).first<time && time<vetoTime021.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime022.size(); m++) {
                if (vetoTime022.at(m).first<time && time<vetoTime022.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime023.size(); m++) {
                if (vetoTime023.at(m).first<time && time<vetoTime023.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime024.size(); m++) {
                if (vetoTime024.at(m).first<time && time<vetoTime024.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime025.size(); m++) {
                if (vetoTime025.at(m).first<time && time<vetoTime025.at(m).second) tag = true;
              }
            }
          }
          if (cubeY==3) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime031.size(); m++) {
                if (vetoTime031.at(m).first<time && time<vetoTime031.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime032.size(); m++) {
                if (vetoTime032.at(m).first<time && time<vetoTime032.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime033.size(); m++) {
                if (vetoTime033.at(m).first<time && time<vetoTime033.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime034.size(); m++) {
                if (vetoTime034.at(m).first<time && time<vetoTime034.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime035.size(); m++) {
                if (vetoTime035.at(m).first<time && time<vetoTime035.at(m).second) tag = true;
              }
            }
          }
          if (cubeY==4) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime041.size(); m++) {
                if (vetoTime041.at(m).first<time && time<vetoTime041.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime042.size(); m++) {
                if (vetoTime042.at(m).first<time && time<vetoTime042.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime043.size(); m++) {
                if (vetoTime043.at(m).first<time && time<vetoTime043.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime044.size(); m++) {
                if (vetoTime044.at(m).first<time && time<vetoTime044.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime045.size(); m++) {
                if (vetoTime045.at(m).first<time && time<vetoTime045.at(m).second) tag = true;
              }
            }
          }
          if (cubeY==5) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime051.size(); m++) {
                if (vetoTime051.at(m).first<time && time<vetoTime051.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime052.size(); m++) {
                if (vetoTime052.at(m).first<time && time<vetoTime052.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime053.size(); m++) {
                if (vetoTime053.at(m).first<time && time<vetoTime053.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime054.size(); m++) {
                if (vetoTime054.at(m).first<time && time<vetoTime054.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime055.size(); m++) {
                if (vetoTime055.at(m).first<time && time<vetoTime055.at(m).second) tag = true;
              }
            }
          }
          if (cubeY==6) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime061.size(); m++) {
                if (vetoTime061.at(m).first<time && time<vetoTime061.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime062.size(); m++) {
                if (vetoTime062.at(m).first<time && time<vetoTime062.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime063.size(); m++) {
                if (vetoTime063.at(m).first<time && time<vetoTime063.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime064.size(); m++) {
                if (vetoTime064.at(m).first<time && time<vetoTime064.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime065.size(); m++) {
                if (vetoTime065.at(m).first<time && time<vetoTime065.at(m).second) tag = true;
              }
            }
          }
          if (cubeY==7) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime071.size(); m++) {
                if (vetoTime071.at(m).first<time && time<vetoTime071.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime072.size(); m++) {
                if (vetoTime072.at(m).first<time && time<vetoTime072.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime073.size(); m++) {
                if (vetoTime073.at(m).first<time && time<vetoTime073.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime074.size(); m++) {
                if (vetoTime074.at(m).first<time && time<vetoTime074.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime075.size(); m++) {
                if (vetoTime075.at(m).first<time && time<vetoTime075.at(m).second) tag = true;
              }
            }
          }
          if (cubeY==8) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime081.size(); m++) {
                if (vetoTime081.at(m).first<time && time<vetoTime081.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime082.size(); m++) {
                if (vetoTime082.at(m).first<time && time<vetoTime082.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime083.size(); m++) {
                if (vetoTime083.at(m).first<time && time<vetoTime083.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime084.size(); m++) {
                if (vetoTime084.at(m).first<time && time<vetoTime084.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime085.size(); m++) {
                if (vetoTime085.at(m).first<time && time<vetoTime085.at(m).second) tag = true;
              }
            }
          }
        }
        if (cubeY==-1) {
          if (cubeX==1) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime101.size(); m++) {
                if (vetoTime101.at(m).first<time && time<vetoTime101.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime102.size(); m++) {
                if (vetoTime102.at(m).first<time && time<vetoTime102.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime103.size(); m++) {
                if (vetoTime103.at(m).first<time && time<vetoTime103.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime104.size(); m++) {
                if (vetoTime104.at(m).first<time && time<vetoTime104.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime105.size(); m++) {
                if (vetoTime105.at(m).first<time && time<vetoTime105.at(m).second) tag = true;
              }
            }
          }
          if (cubeX==2) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime201.size(); m++) {
                if (vetoTime201.at(m).first<time && time<vetoTime201.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime202.size(); m++) {
                if (vetoTime202.at(m).first<time && time<vetoTime202.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime203.size(); m++) {
                if (vetoTime203.at(m).first<time && time<vetoTime203.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime204.size(); m++) {
                if (vetoTime204.at(m).first<time && time<vetoTime204.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime205.size(); m++) {
                if (vetoTime205.at(m).first<time && time<vetoTime205.at(m).second) tag = true;
              }
            }
          }
          if (cubeX==3) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime301.size(); m++) {
                if (vetoTime301.at(m).first<time && time<vetoTime301.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime302.size(); m++) {
                if (vetoTime302.at(m).first<time && time<vetoTime302.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime303.size(); m++) {
                if (vetoTime303.at(m).first<time && time<vetoTime303.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime304.size(); m++) {
                if (vetoTime304.at(m).first<time && time<vetoTime304.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime305.size(); m++) {
                if (vetoTime305.at(m).first<time && time<vetoTime305.at(m).second) tag = true;
              }
            }
          }
          if (cubeX==4) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime401.size(); m++) {
                if (vetoTime401.at(m).first<time && time<vetoTime401.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime402.size(); m++) {
                if (vetoTime402.at(m).first<time && time<vetoTime402.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime403.size(); m++) {
                if (vetoTime403.at(m).first<time && time<vetoTime403.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime404.size(); m++) {
                if (vetoTime404.at(m).first<time && time<vetoTime404.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime405.size(); m++) {
                if (vetoTime405.at(m).first<time && time<vetoTime405.at(m).second) tag = true;
              }
            }
          }
          if (cubeX==5) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime501.size(); m++) {
                if (vetoTime501.at(m).first<time && time<vetoTime501.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime502.size(); m++) {
                if (vetoTime502.at(m).first<time && time<vetoTime502.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime503.size(); m++) {
                if (vetoTime503.at(m).first<time && time<vetoTime503.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime504.size(); m++) {
                if (vetoTime504.at(m).first<time && time<vetoTime504.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime505.size(); m++) {
                if (vetoTime505.at(m).first<time && time<vetoTime505.at(m).second) tag = true;
              }
            }
          }
          if (cubeX==6) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime601.size(); m++) {
                if (vetoTime601.at(m).first<time && time<vetoTime601.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime602.size(); m++) {
                if (vetoTime602.at(m).first<time && time<vetoTime602.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime603.size(); m++) {
                if (vetoTime603.at(m).first<time && time<vetoTime603.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime604.size(); m++) {
                if (vetoTime604.at(m).first<time && time<vetoTime604.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime605.size(); m++) {
                if (vetoTime605.at(m).first<time && time<vetoTime605.at(m).second) tag = true;
              }
            }
          }
          if (cubeX==7) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime701.size(); m++) {
                if (vetoTime701.at(m).first<time && time<vetoTime701.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime702.size(); m++) {
                if (vetoTime702.at(m).first<time && time<vetoTime702.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime703.size(); m++) {
                if (vetoTime703.at(m).first<time && time<vetoTime703.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime704.size(); m++) {
                if (vetoTime704.at(m).first<time && time<vetoTime704.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime705.size(); m++) {
                if (vetoTime705.at(m).first<time && time<vetoTime705.at(m).second) tag = true;
              }
            }
          }
          if (cubeX==8) {
            if (cubeZ==1) {
              for (int m=0; m<vetoTime801.size(); m++) {
                if (vetoTime801.at(m).first<time && time<vetoTime801.at(m).second) tag = true;
              }
            }
            if (cubeZ==2) {
              for (int m=0; m<vetoTime802.size(); m++) {
                if (vetoTime802.at(m).first<time && time<vetoTime802.at(m).second) tag = true;
              }
            }
            if (cubeZ==3) {
              for (int m=0; m<vetoTime803.size(); m++) {
                if (vetoTime803.at(m).first<time && time<vetoTime803.at(m).second) tag = true;
              }
            }
            if (cubeZ==4) {
              for (int m=0; m<vetoTime804.size(); m++) {
                if (vetoTime804.at(m).first<time && time<vetoTime804.at(m).second) tag = true;
              }
            }
            if (cubeZ==5) {
              for (int m=0; m<vetoTime805.size(); m++) {
                if (vetoTime805.at(m).first<time && time<vetoTime805.at(m).second) tag = true;
              }
            }
          }
        }  
        if (tag) continue;

        if (cubeX == -1){
          doesYexist = true;
          if (cubeY_PH < pulseH[i]){
            cubeY_PH = pulseH[i];
            max_cubeY = cubeY;
          }
        }
        else if(cubeY == -1){
          doesXexist = true;
          if(cubeX_PH < pulseH[i]){
            cubeX_PH = pulseH[i];
            max_cubeX = cubeX;
          }
        }

        for (int j = 0; j<129; j++){
          n_array[j] = adcval[i][j]; // For each waveform we fill the vector with 129 adc values
        }
        neutrons.push_back(n_array); // For each waveform we save adc values in the vector of arrays
      } // waveform loop end

      // Reject event if highest PH is not neutron wf 
      if(highest_PHnchi > single_n_thresh || highest_PHgchi < single_g_thresh) continue;
            
      //add neutron like waveforms to get a single waveform
      short int add_neutrons[129] = {0};
      if(!(neutrons.size() > 1 && (doesXexist && doesYexist))) continue;

      for(int i=0; i< neutrons.size(); i++){
        for(int j=0; j< 129; j++) add_neutrons[j] += neutrons[i][j]; // Adding adc values of all waveforms
      }

      //calculate chi-square of added neutrons
      double total_chi[2] = {0};
      get_chi(add_neutrons, 0, total_chi);

      // Selection based on comparison of chi2 with templates
      if (!(total_chi[1] > gl_threshold && total_chi[0] > 0 && total_chi[0] < n_threshold)) continue;

      ++n_neutrons;
//      std::cout << n_neutrons << std::endl;
      vEventID.push_back(iEnt);
      vCubeX.push_back(max_cubeX);
      vCubeY.push_back(max_cubeY);
      vTime.push_back(time);
      vChi2n.push_back(total_chi[0]);
      vChi2gamma.push_back(total_chi[1]);

//-------------------Choosing Zs for the event - the most straightforward cases--------------------------------------------------
      if (cubeZone) vCubeZ.push_back(1);
      else if (sheetZonehalf && !(cubeZone) && !(cubeZtwo)) vCubeZ.push_back(1.5);
      else if (cubeZtwo) vCubeZ.push_back(2);
      else if (sheetZtwohalf && !(cubeZtwo) && !(cubeZthree)) vCubeZ.push_back(2.5);        
      else if (cubeZthree) vCubeZ.push_back(3);
      else if (sheetZthreehalf && !(cubeZthree) && !(cubeZfour)) vCubeZ.push_back(3.5);
      else if (cubeZfour) vCubeZ.push_back(4);
      else if (sheetZfourhalf && !(cubeZfour) && !(cubeZfive)) vCubeZ.push_back(4.5);
      else if (cubeZfive) vCubeZ.push_back(5);
//-----------------------------------More complicated cases----------------------------------------------------------------------
      else {
        for(int i=0; i<(int) nWaveforms; i++){  // loop over waveforms again but for selected neutrons
          int globalChannel = getGlobalChannel(run, (int) boardID[i], (int) adccha[i]);
          int cubeZ = channel2cubeZ[globalChannel];

          //add neutron like waveforms to get a single waveform for each cubeZ layer
          if (cubeZ == 1) { for(int j=0; j<129; j++) add_neutrons_1[j] += adcval[i][j]; }
          if (cubeZ == 2) { for(int j=0; j<129; j++) add_neutrons_2[j] += adcval[i][j]; }
          if (cubeZ == 3) { for(int j=0; j<129; j++) add_neutrons_3[j] += adcval[i][j]; }
          if (cubeZ == 4) { for(int j=0; j<129; j++) add_neutrons_4[j] += adcval[i][j]; }
          if (cubeZ == 5) { for(int j=0; j<129; j++) add_neutrons_5[j] += adcval[i][j]; }          
        } // waveform looping

        double chi1[2] = {0};
        double chi2[2] = {0};
        double chi3[2] = {0};
        double chi4[2] = {0};
        double chi5[2] = {0};

        bool flag1 = false;
        bool flag2 = false;
        bool flag3 = false;
        bool flag4 = false;
        bool flag5 = false;

        get_chi(add_neutrons_1, 0, chi1);
        get_chi(add_neutrons_2, 0, chi2);
        get_chi(add_neutrons_3, 0, chi3);
        get_chi(add_neutrons_4, 0, chi4);
        get_chi(add_neutrons_5, 0, chi5);

        if (chi1[1] > gl_threshold && chi1[0] > 0 && chi1[0] < n_threshold) flag1 = true;
        if (chi2[1] > gl_threshold && chi2[0] > 0 && chi2[0] < n_threshold) flag2 = true;
        if (chi3[1] > gl_threshold && chi3[0] > 0 && chi3[0] < n_threshold) flag3 = true;
        if (chi4[1] > gl_threshold && chi4[0] > 0 && chi4[0] < n_threshold) flag4 = true;
        if (chi5[1] > gl_threshold && chi5[0] > 0 && chi5[0] < n_threshold) flag5 = true;

        if (flag1) {
          if (flag2 || chi2[0] > 0) vCubeZ.push_back(1.5);
          else if (chi2[0] == 0) vCubeZ.push_back(1);
        }
        else if (flag2 && !(flag1)) {
          if (flag3) vCubeZ.push_back(2.5);
          else if (chi1[0] > 0 && chi3[0] == 0) vCubeZ.push_back(1.5);
          else if (chi1[0] == 0 && chi3[0] > 0) vCubeZ.push_back(2.5);
          else if (chi1[0] == 0 && chi3[0] == 0) vCubeZ.push_back(2);
          else if (chi1[0] > 0 && chi3[0] > 0) {
            if (chi1[0] > chi3[0]) vCubeZ.push_back(2.5);
            else if (chi1[0] < chi3[0]) vCubeZ.push_back(1.5);
          }
        }
        else if (flag3 && !(flag2)) {
          if (flag4) vCubeZ.push_back(3.5);
          else if (chi2[0] > 0 && chi4[0] == 0) vCubeZ.push_back(2.5);
          else if (chi2[0] == 0 && chi4[0] > 0) vCubeZ.push_back(3.5);
          else if (chi2[0] == 0 && chi4[0] == 0) vCubeZ.push_back(3);
          else if (chi2[0] > 0 && chi4[0] > 0) {
            if (chi2[0] > chi4[0]) vCubeZ.push_back(3.5);
            else if (chi2[0] < chi4[0]) vCubeZ.push_back(2.5);
          }
        }
        else if (flag4 && !(flag3)) {
          if (flag5) vCubeZ.push_back(4.5);
          else if (chi3[0] > 0 && chi5[0] == 0) vCubeZ.push_back(3.5);
          else if (chi3[0] == 0 && chi5[0] > 0) vCubeZ.push_back(4.5);
          else if (chi3[0] == 0 && chi5[0] == 0) vCubeZ.push_back(4);
          else if (chi3[0] > 0 && chi5[0] > 0) {
            if (chi3[0] > chi5[0]) vCubeZ.push_back(4.5);
            else if (chi3[0] < chi5[0]) vCubeZ.push_back(3.5);
          }
        }
        else if (flag5 && !(flag4)) {
          if (chi4[0] > 0) vCubeZ.push_back(4.5);
          else if (chi4[0] == 0) vCubeZ.push_back(5);
        }
//------------------------------------Most Complicated cases---------------------------------------------------------------------
        else if (!(flag1) && !(flag2) && !(flag3) && !(flag4) && !(flag5)) {          
          std::vector<double> vChooseCubeZ;          
          for (int i=0; i<(int) nWaveforms ; i++){  // loop over waveforms            
            int globalChannel = getGlobalChannel(run, (int) boardID[i], (int) adccha[i]);
            int cubeZ = channel2cubeZ[globalChannel];
            double chi[2] = {0};
            get_chi(adcval[i], baseline[i], chi); // Calculate chis for this waveform and fill the chi variable            

            if (!(chi[1] > single_g_thresh && 0 < chi[0] && chi[0] < single_n_thresh)) continue;            
            vChooseCubeZ.push_back(cubeZ);
          } // waveform loop end

          double chooseCubeZ = 0;
          for (int i=0; i<vChooseCubeZ.size(); i++) {
            chooseCubeZ = vChooseCubeZ.at(i);
            if (i>0) { if (vChooseCubeZ.at(i) != vChooseCubeZ.at(i-1)) chooseCubeZ = 0; }
          }
          vCubeZ.push_back(chooseCubeZ);
        }
//-------------------------------------------------------------------------------------------------------------------------------
        else vCubeZ.push_back(0); // All the other cases if left and not covered
      }
    } // event loop
    
    x_run.push_back(iRun);
    y_neutron.push_back(n_neutrons);
    
    f->Close();
//    myfile.close();
    
    TFile* OutputFile = new TFile(Form("/home/NeutronTree_%i.root", iRun), "RECREATE");
    
    structNeutronEvent iNeu;
    TTree *TNeu = new TTree("TNeu","Neutron Events");
    setupNeutronTreeFill(*TNeu, iNeu);
    
    for (int i=0; i<vTime.size(); i++) {
      iNeu.cubeX = vCubeX.at(i);
      iNeu.cubeY = vCubeY.at(i);
      iNeu.cubeZ = vCubeZ.at(i);
      iNeu.time = vTime.at(i);
      iNeu.chi2n = vChi2n.at(i);
      iNeu.chi2gamma = vChi2gamma.at(i);
      TNeu->Fill();
    }
    
    TNeu->Write();
    OutputFile->Write();
    OutputFile->Close();
  }
/*
  h_adc1->Draw();
  std::cout << "Draw" << std::endl;
  gPad->Print(Form("draw%i.png", 1));

  h_adc2->Draw();
  gPad->Print(Form("draw%i.png", 2));
*/
  for(int i=0; i<x_run.size(); i++){
      std::cout << x_run[i] << ","<< y_neutron[i] << std::endl;
  }
  return 0;
}
