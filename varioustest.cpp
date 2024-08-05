
#include "histogram_functions.h" // Include the file with histogram functions
void varioustest(){
    TFile *file353 = TFile::Open("WindowIntMatched_final_353.root"); //-0.3
    TFile *file524 = TFile::Open("WindowIntMatched_final_524.root"); //-0.2
    TFile *file454 = TFile::Open("WindowIntMatched_final_454.root"); //0.86
    TFile *file394 = TFile::Open("WindowIntMatched_final_394.root"); //1.12
    TFile *file459 = TFile::Open("WindowIntMatched_final_459.root"); //0.6
    TFile *file415 = TFile::Open("WindowIntMatched_final_415.root"); //0.42
    TFile *file345 = TFile::Open("WindowIntMatched_final_345.root"); //0.22

    TTree *TOF_353[8];
    TTree *ACT01_353[4];
    TTree *ACT23_353[4];
    TTree *pbglass_353;

    TTree *TOF_454[8];
    TTree *ACT01_454[4];
    TTree *ACT23_454[4];
    TTree *pbglass_454;

    TTree *TOF_524[8];
    TTree *ACT01_524[4];
    TTree *ACT23_524[4];
    TTree *pbglass_524;
    
    TTree *TOF_394[8];
    TTree *ACT01_394[4];
    TTree *ACT23_394[4];
    TTree *pbglass_394;
    TTree *TOF_459[8];
    TTree *ACT01_459[4];
    TTree *ACT23_459[4];
    TTree *pbglass_459;
    TTree *TOF_415[8];
    TTree *ACT01_415[4];
    TTree *ACT23_415[4];
    TTree *pbglass_415;
    TTree *TOF_345[8];
    TTree *ACT01_345[4];
    TTree *ACT23_345[4];
    TTree *pbglass_345;


    SetupTrees(file454,TOF_454,ACT01_454,ACT23_454,pbglass_454);
    SetupTrees(file353,TOF_353,ACT01_353,ACT23_353,pbglass_353);
    SetupTrees(file524,TOF_524,ACT01_524,ACT23_524,pbglass_524);

    SetupTrees(file394,TOF_394,ACT01_394,ACT23_394,pbglass_394);
    SetupTrees(file459,TOF_459,ACT01_459,ACT23_459,pbglass_459);
    SetupTrees(file415,TOF_415,ACT01_415,ACT23_415,pbglass_415);
    SetupTrees(file345,TOF_345,ACT01_345,ACT23_345,pbglass_345);

    Double_t momentum = -0.3;
    TTree *ACT01trees[14];
    for(Int_t i=0; i<4; ++i){
        ACT01trees[i] = TOF_353[i];
        ACT01trees[i+4] = TOF_353[i+4];
        ACT01trees[9+i] = ACT01_353[i];
    }
    ACT01trees[8] = pbglass_353;
    ACT01trees[13] = nullptr;
    TCanvas *canvas = new TCanvas();
    TH1F *hist = new TH1F("title1","Run 353;TOF;Occurence",100,10,15);
    Make1DTOFPlot(canvas, hist, ACT01trees, {"IntCharge", "MaxVoltage", "PeakTime"}, 353);

    TGraph* wIntC353 = new TGraph();
    TGraph *wIntC353_OP = MakeSimple2DGraph0vs1NoDraw(wIntC353, ACT01trees, {"IntCharge", "MaxVoltage"}, momentum);
    Draw2500Point(wIntC353, wIntC353_OP, momentum, "wIntCharge");
    
    momentum = 0.86;
    for(Int_t i=0; i<4; ++i){
        ACT01trees[i] = TOF_454[i];
        ACT01trees[i+4] = TOF_454[i+4];
        ACT01trees[9+i] = ACT01_454[i];
    }
    ACT01trees[8] = pbglass_454;

    TCanvas *canvas2 = new TCanvas();
    TH1F *hist2 = new TH1F("title2","Run 454;TOF;Occurence",100,10,15);
    Make1DTOFPlot(canvas2, hist2, ACT01trees, {"IntCharge", "MaxVoltage", "PeakTime"}, 454);

    TGraph* wIntC454 = new TGraph();
    TGraph *wIntC454_OP = MakeSimple2DGraph0vs1NoDraw(wIntC454, ACT01trees, {"IntCharge", "MaxVoltage"}, momentum);
    Draw2500Point(wIntC454, wIntC454_OP, momentum, "wIntCharge");

    momentum = -0.2;
    for(Int_t i=0; i<4; ++i){
        ACT01trees[i] = TOF_524[i];
        ACT01trees[i+4] = TOF_524[i+4];
        ACT01trees[9+i] = ACT01_524[i];
    }
    ACT01trees[8] = pbglass_524;

    TCanvas *canvas3 = new TCanvas();
    TH1F *hist3 = new TH1F("title3","Run 524;TOF;Occurence",100,10,15);
    Make1DTOFPlot(canvas3, hist3, ACT01trees, {"IntCharge", "MaxVoltage", "PeakTime"}, 524);

    ////////////
    momentum = 1.12;
    for(Int_t i=0; i<4; ++i){
        ACT01trees[i] = TOF_394[i];
        ACT01trees[i+4] = TOF_394[i+4];
        ACT01trees[9+i] = ACT01_394[i];
    }
    ACT01trees[8] = pbglass_394;

    TGraph* wIntC394 = new TGraph();
    TGraph *wIntC394_OP = MakeSimple2DGraph0vs1NoDraw(wIntC394, ACT01trees, {"IntCharge", "MaxVoltage"}, momentum);
    Draw2500Point(wIntC394, wIntC394_OP, momentum, "wIntCharge");

    momentum = 0.6;
    for(Int_t i=0; i<4; ++i){
        ACT01trees[i] = TOF_459[i];
        ACT01trees[i+4] = TOF_459[i+4];
        ACT01trees[9+i] = ACT01_459[i];
    }
    ACT01trees[8] = pbglass_459;

    TGraph* wIntC459 = new TGraph();
    TGraph *wIntC459_OP = MakeSimple2DGraph0vs1NoDraw(wIntC459, ACT01trees, {"IntCharge", "MaxVoltage"}, momentum);
    Draw2500Point(wIntC459, wIntC459_OP, momentum, "wIntCharge");

    momentum = 0.42;
    for(Int_t i=0; i<4; ++i){
        ACT01trees[i] = TOF_415[i];
        ACT01trees[i+4] = TOF_415[i+4];
        ACT01trees[9+i] = ACT01_415[i];
    }
    ACT01trees[8] = pbglass_415;

    TGraph* wIntC415 = new TGraph();
    TGraph *wIntC415_OP = MakeSimple2DGraph0vs1NoDraw(wIntC415, ACT01trees, {"IntCharge", "MaxVoltage"}, momentum);
    Draw2500Point(wIntC415, wIntC415_OP, momentum, "wIntCharge");

    momentum = 0.22;
    for(Int_t i=0; i<4; ++i){
        ACT01trees[i] = TOF_345[i];
        ACT01trees[i+4] = TOF_345[i+4];
        ACT01trees[9+i] = ACT01_345[i];
    }
    ACT01trees[8] = pbglass_345;

    TGraph* wIntC345 = new TGraph();
    TGraph *wIntC345_OP = MakeSimple2DGraph0vs1NoDraw(wIntC345, ACT01trees, {"IntCharge", "MaxVoltage"}, momentum);
    Draw2500Point(wIntC345, wIntC345_OP, momentum, "wIntCharge");



}