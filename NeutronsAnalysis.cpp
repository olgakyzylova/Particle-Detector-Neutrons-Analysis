/******************************************************
 * study the neutron number change with the reactor 
 * power variation.
 * Include spatial distribution
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


using std::cout; 
using std::endl;
int getGlobalChannel(int run, int boardID, int adccha) {
    int globalChannel = -1; 
    if (run < 10000) {
        globalChannel = 48*boardID + adccha;
    } else {
        if (boardID == 1) {
            globalChannel = 16 + adccha;
        } else {
            if (adccha<=15) {
                globalChannel = adccha;
            } else {
                globalChannel = 48 + adccha;
            }
        }
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


//chi-square function
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
            if(adcval[j] > baseline + 12){
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
        if(j>1 && j%11 == 0){
            k+= 1;
        }
    }

    double norm_data[8] = {0};
    double error[8] = {0};
    int ndf = 0;
    for (int j=0; j<8; j++){
        if(data[j] == 0) continue;
        ++ndf;
        norm_data[j] = data[j]/sum_data;
        error[j] = sqrt(6*abs(data[j]))/sum_data;
        n_chi += (norm_data[j] - n_tem[j])*(norm_data[j] - n_tem[j])/(pow((error[j]+err_n_tem[j])/2,2));
        g_chi += (norm_data[j] - g_tem[j])*(norm_data[j] - g_tem[j])/(pow((error[j]+err_n_tem[j])/2,2));
    }
    if(ndf > 2){
        g_chi = g_chi/ndf;
        n_chi = n_chi/ndf;
    } else {
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
        cout << "Usage: " << argv[0] << " [run]" << "[last run]" <<endl;
        return -1;
    }
    sscanf(argv[1], "%i", &run);
    sscanf(argv[2], "%i", &last_run);
    cout<<"run: "<< run << ", last_run: "<< last_run <<endl;

    TH1D *h_wf = new TH1D("h_wf", "Neutron Waveform", 60, 0, 60);
    TH1D *h_nPH = new TH1D("h_nPH", "Neutron PH", 100, 0, 500);
    TH2D *h_chi = new TH2D("h_chi", "2D chi-square distribution", 100, 0, 200, 100, 0, 50);
    TH2D *h_chin = new TH2D("h_chin", "2D chi-square distribution", 100, 0, 100, 100, 0, 15);
    TH2D *h_tchi = new TH2D("h_tchi", "2D chi-square distribution", 100, 0, 200, 100, 0, 50);
    TH2D *h_tchin = new TH2D("h_tchin", "2D chi-square distribution", 100, 0, 100, 100, 0, 15);
    TH2D *h_delnchi = new TH2D("h_delnchi", "#Delta #chi_{n}^{2}", 100, 0, 20, 100, -20, 10);
    TH2D *h_sumchi = new TH2D("h_sumchi", "2D chi-square sum distribution", 100, 0, 200, 100, 0, 50);
    TH2D *h_ndist = new TH2D("h_ndist", "2D Neutron distribution", 8, 0.5, 8.5, 8, 0.5, 8.5);

    int ADCSIZE = 129;
    int n_threshold = 14;
    int gl_threshold = 30;
    int single_g_thresh = 15;
    int single_n_thresh = 20;

    std::vector <int> x_run;
    std::vector <int> y_neutron;
    int nx_run = 0;
    int ny_neutron = 0;
    int n = last_run - run + 1;

    for(int iRun = run; iRun <= last_run ; iRun++){
        cout <<"run: "<< iRun << endl;
        int n_neutrons = 0;

        TFile *f;
        if(gSystem->AccessPathName(Form("/run%i.root", iRun, iRun))){
            continue;
        } else{
            f = new TFile(Form("run%i.root", iRun, iRun));
        }

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
        UChar_t         nWaveforms;
        UChar_t         boardID[80];   //[nWaveforms]
        UChar_t         adccha[80];   //[nWaveforms]
        Short_t         adcval[80][129];   //[nWaveforms]
        Float_t         baseline[80];   //[nWaveforms]
        Float_t         pulseH[80];   //[nWaveforms]
        UChar_t         nBoards;
        Short_t         trig_pattern[2];   //[nBoards]
        Short_t         adcSize[2];   //[nBoards]
        ULong64_t       adcTime[2];   //[nBoards]
        ULong64_t       globalTime[2];   //[nBoards]

        // List of branches
        TBranch        *b_nWaveforms;
        TBranch        *b_boardID;
        TBranch        *b_adccha;
        TBranch        *b_adcval;
        TBranch        *b_baseline;
        TBranch        *b_pulseH;
        TBranch        *b_nBoards;
        TBranch        *b_trig_pattern;
        TBranch        *b_adcSize;
        TBranch        *b_adcTime;
        TBranch        *b_globalTime;

        t->SetBranchAddress("nWaveforms", &nWaveforms, &b_nWaveforms);
        t->SetBranchAddress("boardID", boardID, &b_boardID);
        t->SetBranchAddress("adccha", adccha, &b_adccha);
        t->SetBranchAddress("adcval", adcval, &b_adcval);
        t->SetBranchAddress("baseline", baseline, &b_baseline);
        t->SetBranchAddress("pulseH", pulseH, &b_pulseH);
        t->SetBranchAddress("nBoards", &nBoards, &b_nBoards);
        t->SetBranchAddress("trig_pattern", trig_pattern, &b_trig_pattern);
        t->SetBranchAddress("adcSize", adcSize, &b_adcSize);
        t->SetBranchAddress("adcTime", adcTime, &b_adcTime);
        t->SetBranchAddress("globalTime", globalTime, &b_globalTime);

        std::vector<double> vCubeX;
        std::vector<double> vCubeY;
        std::vector<double> vCubeZ;
        std::vector<double> vTime;
        std::vector<double> vChi2n;
        std::vector<double> vChi2gamma;

        for (int iEnt=0; iEnt<t->GetEntries(); iEnt++) {
//            if(iEnt > 70000) break;

            Long64_t tentry = t->LoadTree(iEnt);

            b_nWaveforms->GetEntry(tentry);
            b_boardID->GetEntry(tentry);
            b_adccha->GetEntry(tentry);
            b_adcval->GetEntry(tentry);
            b_baseline->GetEntry(tentry);
            b_pulseH->GetEntry(tentry);
            b_nBoards->GetEntry(tentry);
            b_trig_pattern->GetEntry(tentry);
            b_adcSize->GetEntry(tentry);
            b_adcTime->GetEntry(tentry);
            b_globalTime->GetEntry(tentry);

            // Is it a shorter adcSize
            bool one_of_boards_shorter_adcSize = false;
            for (int i=0; i<(int) nBoards; i++) {
                if (adcSize[i]!=ADCSIZE) {
                    one_of_boards_shorter_adcSize = true;
                }
            }

            // Is it an external trigger
            bool one_of_boards_external_trigger = false;
            for (int i=0; i<(int) nBoards; i++) {
                if ((trig_pattern[i] & 0x200) == 0x200) {
                    one_of_boards_external_trigger = true;
                }
            }

            if (one_of_boards_external_trigger) continue;

            std::vector<std::array<int, 129>> neutrons;
            std::vector<double> n_chisq;
            std::vector<double> g_chisq;
            std::vector<double> n_PH;
            std::vector<bool> x_view;
            bool isXview = false;
            bool isYview = false;
            bool isNeutron = true;
            double total_nchi = 0;
            double total_gchi = 0;
            double highest_PH = 0;
            double highest_PHnchi = -1;
            double highest_PHgchi = -1;

            double max_cubeX = 0;
            double max_cubeY = 0;
            double max_cubeZ = 0;
            double cubeX_PH = 0;
            double cubeY_PH = 0;
            double time = 0;
            double time2 = 0;

            //get neutron like waveforms and save in neutrons vector
            gStyle->SetOptStat(0);
            for(int i=0; i<(int) nWaveforms ; i++){                
                std::array<int, 129> n_array = {0};
                int globalChannel = getGlobalChannel(run, (int) boardID[i], (int) adccha[i]);
                int cubeX = channel2cubeX[globalChannel];
                int cubeY = channel2cubeY[globalChannel];
                int cubeZ = channel2cubeZ[globalChannel];

                double chi[2] = {0};
                get_chi(adcval[i], baseline[i], chi);
                if(chi[0] != 0){
                    h_chi->Fill(chi[1], chi[0]);
                    h_chin->Fill(chi[1], chi[0]);
                }
                double n_chi = chi[0];
                double g_chi = chi[1];

                if(highest_PH < pulseH[i]){
                    highest_PH = pulseH[i];
                    highest_PHnchi = n_chi;
                    highest_PHgchi = g_chi;
                }

                if(g_chi > single_g_thresh && n_chi < single_n_thresh ){
                    max_cubeZ = cubeZ;

                    time = globalTime[0];
                    if(cubeX == -1){
                        isYview = true;
                        if(cubeY_PH < pulseH[i]){
                            cubeY_PH = pulseH[i];
                            max_cubeY = cubeY;
                        }
                    }else if(cubeY == -1){
                        isXview = true;
                        if(cubeX_PH < pulseH[i]){
                            cubeX_PH = pulseH[i];
                            max_cubeX = cubeX;
                        }
                    }

                    for (int j = 0; j<129; j++){
                        n_array[j] = adcval[i][j];
                    }
                    neutrons.push_back(n_array);
                    n_chisq.push_back(n_chi);
                    g_chisq.push_back(g_chi);
                    n_PH.push_back(pulseH[i]);
                    x_view.push_back(cubeY == -1);
                }
            }

            //reject event if highest PH is not neutron wf 
            if(highest_PHnchi > single_n_thresh || highest_PHgchi < single_g_thresh) continue;

            //add chi-square old way
            int x_index = -1;
            int y_index = -1;
            double max_PHx = 0;
            double max_PHy = 0;
            for(int l=0; l<n_PH.size(); l++){
                if(x_view[l] == 1){ //x view
                    if(max_PHx < n_PH[l]){
                        max_PHx = n_PH[l];
                        x_index = l;
                    }
                } else {
                    if(max_PHy < n_PH[l]){
                        max_PHy = n_PH[l];
                        y_index = l;
                    }
                }
            }

            if(x_index != -1 && y_index != -1){
                total_nchi = n_chisq[x_index] + n_chisq[y_index];
                total_gchi = g_chisq[x_index] + g_chisq[y_index];
                h_sumchi->Fill(total_gchi, total_nchi);
            }
            
            //add neutron like waveforms to get a single waveform
            short int add_neutrons[129] = {0};
            if(neutrons.size() >= 2 && (isXview && isYview)){
                for(int i=0; i< neutrons.size(); i++){
                    for(int j=0; j< 129; j++){
                        add_neutrons[j] += neutrons[i][j];
                        h_wf->SetBinContent(j+1, neutrons[i][j]);
                    }
                    h_wf->Draw();
                }

                //calculate baseline minumum and pulse peak bin
                double min_baseline = 10000;
                int peak_bin = 0;
                double peak_ADC = 0;
                for(int j=0; j< 129; j++){
                    if(j < 25 && min_baseline > add_neutrons[j]){
                        min_baseline = add_neutrons[j];
                    }
                    if(peak_ADC < add_neutrons[j]){
                        peak_ADC = add_neutrons[j];
                        peak_bin = j;
                    }
                }

                //calculate chi-square of added neutrons
                double total_chi[2] = {0};
                get_chi(add_neutrons, 0, total_chi);
                h_tchi->Fill(total_chi[1], total_chi[0]);
                h_tchin->Fill(total_chi[1], total_chi[0]);

                if(total_chi[1] > gl_threshold && total_chi[0]>0 && total_chi[0]<n_threshold){
                    ++n_neutrons;
                    vCubeX.push_back(max_cubeX);
                    vCubeY.push_back(max_cubeY);
                    vCubeZ.push_back(max_cubeZ);
                    vTime.push_back(time);
                    vChi2n.push_back(total_chi[0]);
                    vChi2gamma.push_back(total_chi[1]);
                }
            }
        } //event loop

        x_run.push_back(iRun);
        nx_run = iRun;
        ny_neutron = n_neutrons;
        y_neutron.push_back(n_neutrons);

        f->Close();

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

    for(int i=0; i<x_run.size(); i++){
        cout << x_run[i] << ","<< y_neutron[i] << endl;
    }
    return 0;
}
