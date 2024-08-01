#include "histogram_functions.h"


void EfficiencyComparaison_final(){

    TFile *file457 = TFile::Open("WindowIntMatched_final_457.root"); //0.7GeV 
    TFile *file475 = TFile::Open("WindowIntMatched_final_475.root"); //0.3GeV 
    TFile *file353 = TFile::Open("WindowIntMatched_final_353.root"); //-0.3GeV 
    TFile *file445 = TFile::Open("WindowIntMatched_final_445.root"); //-0.7GeV 

    TTree *TOF_457[8];
    TTree *ACT01_457[4];
    TTree *ACT23_457[4];
    TTree *pbglass_457;
    
    TTree *TOF_475[8];
    TTree *ACT01_475[4];
    TTree *ACT23_475[4];
    TTree *pbglass_475;
    
    TTree *TOF_353[8];
    TTree *ACT01_353[4];
    TTree *ACT23_353[4];
    TTree *pbglass_353;

    TTree *TOF_445[8];
    TTree *ACT01_445[4];
    TTree *ACT23_445[4];
    TTree *pbglass_445;

    SetupTrees(file457,TOF_457,ACT01_457,ACT23_457,pbglass_457);
    SetupTrees(file475,TOF_475,ACT01_475,ACT23_475,pbglass_475);
    SetupTrees(file353,TOF_353,ACT01_353,ACT23_353,pbglass_353);
    SetupTrees(file445,TOF_445,ACT01_445,ACT23_445,pbglass_445);

    Int_t IdxArr[5] = {2,14,24,34,44}; //Best-ish for wIntCharge
    Int_t IdxArr2[5] = {43,46,47,48,44};
    Int_t IdxArr3[5] = {25,5,15,35,45}; //Best-ish for IntCharge
    Int_t IdxArr4[5] = {46,47,48,49,45};
    /*
    Double_t momentum = 0.7;
    TTree *completeTTree[14];
    for(Int_t i=0; i<4; ++i){
        completeTTree[i] = TOF_457[i];
        completeTTree[i+4] = TOF_457[i+4];
        completeTTree[9+i] = ACT01_457[i];
    }
    completeTTree[8] = pbglass_457;
    completeTTree[13] = nullptr;

    TGraph* IntC457 = new TGraph();
    TGraph* wIntC457 = new TGraph();

    TGraph *IntC457_OP = MakeSimple2DGraph0vs1NoDraw(IntC457, completeTTree, {"IntCharge", "MaxVoltage"}, momentum, "IntCharge");
    TGraph *wIntC457_OP = MakeSimple2DGraph0vs1NoDraw(wIntC457, completeTTree, {"IntCharge", "MaxVoltage"}, momentum);
    TH2D* pureff_2Dhist_Bg = nullptr;
    TH2D* pureff_2Dhist_Sig = nullptr;
    GenerateHistograms(IntC457, IntC457_OP, pureff_2Dhist_Bg, pureff_2Dhist_Sig, 50);
    GenerateHistograms(wIntC457, wIntC457_OP, pureff_2Dhist_Bg, pureff_2Dhist_Sig, 50);
    */

    Double_t momentum = 0.7; //GeV/c
    TTree *completeTTree[14];
    for(Int_t i=0; i<4; ++i){
        completeTTree[i] = TOF_457[i];
        completeTTree[i+4] = TOF_457[i+4];
        completeTTree[9+i] = ACT01_457[i];
    }
    completeTTree[8] = pbglass_457;
    completeTTree[13] = nullptr;

    TGraph* IntC457 = new TGraph();
    TGraph* wIntC457 = new TGraph();

    TGraph *IntC457_OP = MakeSimple2DGraph0vs1NoDraw(IntC457, completeTTree, {"IntCharge", "MaxVoltage"}, momentum, "IntCharge");
    TGraph *wIntC457_OP = MakeSimple2DGraph0vs1NoDraw(wIntC457, completeTTree, {"IntCharge", "MaxVoltage"}, momentum);

    //Draw5GraphComp(IntC457,IntC457_OP,IdxArr3,momentum,"IntCharge");
    //Draw5GraphComp(IntC457,IntC457_OP,IdxArr4,momentum,"IntCharge","more");
    //Draw5GraphComp(wIntC457,wIntC457_OP,IdxArr,momentum,"wIntCharge");
    //Draw5GraphComp(wIntC457,wIntC457_OP,IdxArr2,momentum,"wIntCharge","more");
    Double_t maxVal, finalEff;
    Double_t finalm, finalb;
    bool success;
    FindMaxPurityEfficiency(wIntC457,wIntC457_OP,100,maxVal,finalEff,finalm,finalb,success);
    cout << "The max eff at 100 rejection rate at 0.7GeV wIntCharge is " << maxVal << " with a rejection of " << 1/finalEff << " the m and b are respectively " << finalm << " and " << finalb << endl;

    FindMaxPurityEfficiency(IntC457,IntC457_OP,100,maxVal,finalEff,finalm,finalb,success);
    cout << "The max eff at 100 rejection rate at 0.7GeV IntCharge is " << maxVal << " with a rejection of " << 1/finalEff << " the m and b are respectively " << finalm << " and " << finalb << endl;


    delete IntC457;
    delete IntC457_OP; 
    delete wIntC457;
    delete wIntC457_OP;
    /*
    momentum = 0.3; //GeV/c
    for(Int_t i=0; i<4; ++i){
        completeTTree[i] = TOF_475[i];
        completeTTree[i+4] = TOF_475[i+4];
        completeTTree[9+i] = ACT01_475[i];
    }
    completeTTree[8] = pbglass_475;
    completeTTree[13] = nullptr;

    //TGraph* IntC475 = new TGraph();
    TGraph* wIntC475 = new TGraph();

    //TGraph *IntC475_OP = MakeSimple2DGraph0vs1NoDraw(IntC475, completeTTree, {"IntCharge", "MaxVoltage"}, momentum, "IntCharge");
    TGraph *wIntC475_OP = MakeSimple2DGraph0vs1NoDraw(wIntC475, completeTTree, {"IntCharge", "MaxVoltage"}, momentum);

    Draw5GraphComp(wIntC475,wIntC475_OP,IdxArr,momentum,"wIntCharge");
    Draw5GraphComp(wIntC475,wIntC475_OP,IdxArr2,momentum,"wIntCharge","more");
    delete wIntC475;
    delete wIntC475_OP;

    momentum = -0.3; //GeV/c
    for(Int_t i=0; i<4; ++i){
        completeTTree[i] = TOF_353[i];
        completeTTree[i+4] = TOF_353[i+4];
        completeTTree[9+i] = ACT01_353[i];
    }
    completeTTree[8] = pbglass_353;
    completeTTree[13] = nullptr;

    //TGraph* IntC353 = new TGraph();
    TGraph* wIntC353 = new TGraph();

    //TGraph *IntC353_OP = MakeSimple2DGraph0vs1NoDraw(IntC353, completeTTree, {"IntCharge", "MaxVoltage"}, momentum, "IntCharge");
    TGraph *wIntC353_OP = MakeSimple2DGraph0vs1NoDraw(wIntC353, completeTTree, {"IntCharge", "MaxVoltage"}, momentum);

    Draw5GraphComp(wIntC353,wIntC353_OP,IdxArr,momentum,"wIntCharge");
    Draw5GraphComp(wIntC353,wIntC353_OP,IdxArr2,momentum,"wIntCharge","more");
    delete wIntC353;
    delete wIntC353_OP;


    momentum = -0.7; //GeV/c
    for(Int_t i=0; i<4; ++i){
        completeTTree[i] = TOF_445[i];
        completeTTree[i+4] = TOF_445[i+4];
        completeTTree[9+i] = ACT01_445[i];
    }
    completeTTree[8] = pbglass_445;
    completeTTree[13] = nullptr;

    //TGraph* IntC445 = new TGraph();
    TGraph* wIntC445 = new TGraph();

    //TGraph *IntC445_OP = MakeSimple2DGraph0vs1NoDraw(IntC445, completeTTree, {"IntCharge", "MaxVoltage"}, momentum, "IntCharge");
    TGraph *wIntC445_OP = MakeSimple2DGraph0vs1NoDraw(wIntC445, completeTTree, {"IntCharge", "MaxVoltage"}, momentum);

    Draw5GraphComp(wIntC445,wIntC445_OP,IdxArr,momentum,"wIntCharge");
    Draw5GraphComp(wIntC445,wIntC445_OP,IdxArr2,momentum,"wIntCharge","more");
    delete wIntC445;
    delete wIntC445_OP;
    */
}