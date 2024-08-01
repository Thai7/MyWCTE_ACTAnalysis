#include "histogram_functions.h" // Include the file with histogram functions

void ACT_analysis(){
    //open ntuple rootfile, extract n_peak==1 and then plot ACT0 vs ACT1 plot
    //for 2 different ACT one with e-veto one without.
    //Start with all particle, implement particle selection after.
   
    //Only 452/464 has veto/noveto
    TFile *file452 = TFile::Open("WindowIntMatched_final_452.root");
    TFile *file464 = TFile::Open("WindowIntMatched_final_464.root");
    TFile *file457 = TFile::Open("WindowIntMatched_final_457.root"); //0.7GeV 
    TFile *file435 = TFile::Open("WindowIntMatched_final_435.root"); //0.5GeV 
    TFile *file475 = TFile::Open("WindowIntMatched_final_475.root"); //0.3GeV 
    TFile *file353 = TFile::Open("WindowIntMatched_final_353.root"); //-0.3GeV 

    TTree *TOF_452[8];
    TTree *ACT01_452[4];
    TTree *ACT23_452[4];
    TTree *pbglass_452;

    TTree *TOF_464[8];
    TTree *ACT01_464[4];
    TTree *ACT23_464[4];
    TTree *pbglass_464;

    TTree *TOF_457[8];
    TTree *ACT01_457[4];
    TTree *ACT23_457[4];
    TTree *pbglass_457;
    
    TTree *TOF_435[8];
    TTree *ACT01_435[4];
    TTree *ACT23_435[4];
    TTree *pbglass_435;

    TTree *TOF_475[8];
    TTree *ACT01_475[4];
    TTree *ACT23_475[4];
    TTree *pbglass_475;
    
    TTree *TOF_353[8];
    TTree *ACT01_353[4];
    TTree *ACT23_353[4];
    TTree *pbglass_353;


    SetupTrees(file452,TOF_452,ACT01_452,ACT23_452,pbglass_452);
    SetupTrees(file464,TOF_464,ACT01_464,ACT23_464,pbglass_464);
    SetupTrees(file457,TOF_457,ACT01_457,ACT23_457,pbglass_457);
    SetupTrees(file435,TOF_435,ACT01_435,ACT23_435,pbglass_435);
    SetupTrees(file475,TOF_475,ACT01_475,ACT23_475,pbglass_475);
    SetupTrees(file353,TOF_353,ACT01_353,ACT23_353,pbglass_353);

    Double_t limit[4] = {-0.2,0.8,-0.1,0.4};
    Double_t pVal = 1.;

    //MAKING SCATTER PLOT OF ACT0 vs ACT1 with the updated version!
    TTree *ACT01complete[14];
    for(Int_t i=0; i<4; ++i){
        ACT01complete[i] = TOF_452[i];
        ACT01complete[i+4] = TOF_452[i+4];
        ACT01complete[9+i] = ACT01_452[i];
    }
    ACT01complete[8] = pbglass_452;
    ACT01complete[13] = nullptr;
    
    auto [graph_IntCharge_0vs1_452, graph_IntCharge_0vs1_452_OP] = Create2DPlot(pVal, "IntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true);

    auto [graph_WindowIntCharge_0vs1_452, graph_WindowIntCharge_0vs1_452_OP] = Create2DPlot(pVal, "wIntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true);

    auto [graph_IntChargeNoW_0vs1_452, graph_IntChargeNoW_0vs1_452_OP] = Create2DPlot(pVal, "IntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, false);
    
    pVal = 0.3;
    for(Int_t i=0; i<4; ++i){
        ACT01complete[i] = TOF_475[i];
        ACT01complete[i+4] = TOF_475[i+4];
        ACT01complete[9+i] = ACT01_475[i];
    }
    ACT01complete[8] = pbglass_475;
    ACT01complete[13] = nullptr;
    
    auto [graph_WindowIntCharge_0vs1_475, graph_WindowIntCharge_0vs1_475_OP] = Create2DPlot(pVal, "wIntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true);

    pVal = -0.3;
    for(Int_t i=0; i<4; ++i){
        ACT01complete[i] = TOF_353[i];
        ACT01complete[i+4] = TOF_353[i+4];
        ACT01complete[9+i] = ACT01_353[i];
    }
    ACT01complete[8] = pbglass_353;
    ACT01complete[13] = nullptr;
    
    auto [graph_WindowIntCharge_0vs1_353, graph_WindowIntCharge_0vs1_353_OP] = Create2DPlot(pVal, "wIntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true);

    pVal = 1.0;

    //MAKING 1DHISTE histograms

    TTree *ACT01T0tree[11];
    for(Int_t i=0;i<8;i++){
        ACT01T0tree[i] = TOF_452[i];
    } 
    ACT01T0tree[8] = pbglass_452;
    ACT01T0tree[9] = ACT01_452[0];
    ACT01T0tree[10] = nullptr;

    TCanvas *upd_H1D_ST_ACT0L_T0_run452 = new TCanvas("upd_H1D_STCorr_ACT0L_T0_run452","title",800,600);
    TH1F *h_ACT0L_T0_r452 = new TH1F("h_ACT0L_T0_r452","Run 452 p=1GeV, ACT0L-T0avg Signal Time Corrected distribution with T0,T1,PbGlass and ACT0L nPeak==1 w/ WWI cut;ACT0L-T0(ns);Counts/Bin",210,-60,150);
    Make1DHistEonly(upd_H1D_ST_ACT0L_T0_run452,h_ACT0L_T0_r452,ACT01T0tree,{},"TOFL",pVal);

    pVal = 0.7;
    for(Int_t i=0;i<8;i++){
        ACT01T0tree[i] = TOF_457[i];
    } 
    ACT01T0tree[8] = pbglass_457;
    ACT01T0tree[9] = ACT01_457[0];
    ACT01T0tree[10] = nullptr;

    TCanvas *upd_H1D_ST_ACT0L_T0_run457 = new TCanvas("upd_H1D_STCorr_ACT0L_T0_run457","title",800,600);
    TH1F *h_ACT0L_T0_r457 = new TH1F("h_ACT0L_T0_r457","Run 457 p=0.7GeV, ACT0L-T0avg Signal Time Corrected distribution with T0,T1,PbGlass and ACT0L nPeak==1 w/ WWI cut;ACT0L-T0(ns);Counts/Bin",210,-60,150);
    Make1DHistEonly(upd_H1D_ST_ACT0L_T0_run457,h_ACT0L_T0_r457,ACT01T0tree,{},"TOFL",pVal);

    pVal = 0.5;
    for(Int_t i=0;i<8;i++){
        ACT01T0tree[i] = TOF_435[i];
    } 
    ACT01T0tree[8] = pbglass_435;
    ACT01T0tree[9] = ACT01_435[0];
    ACT01T0tree[10] = nullptr;

    TCanvas *upd_H1D_ST_ACT0L_T0_run435 = new TCanvas("upd_H1D_STCorr_ACT0L_T0_run435","title",800,600);
    TH1F *h_ACT0L_T0_r435 = new TH1F("h_ACT0L_T0_r435","Run 435 p=0.5GeV, ACT0L-T0avg Signal Time Corrected distribution with T0,T1,PbGlass and ACT0L nPeak==1 w/ WWI cut;ACT0L-T0(ns);Counts/Bin",210,-60,150);
    Make1DHistEonly(upd_H1D_ST_ACT0L_T0_run435,h_ACT0L_T0_r435,ACT01T0tree,{},"TOFL",pVal);


    ////////////////////////////////////////////////////

    ///0.7GeV!  
    pVal = 0.7;
    for(Int_t i=0; i<4; ++i){
        ACT01complete[i] = TOF_457[i];
        ACT01complete[i+4] = TOF_457[i+4];
        ACT01complete[9+i] = ACT01_457[i];
    }
    ACT01complete[8] = pbglass_457;
    ACT01complete[13] = nullptr;

    auto [graph_IntCharge_0vs1_457, graph_IntCharge_0vs1_457_OP] = Create2DPlot(pVal, "IntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true);

    auto [graph_WindowIntCharge_0vs1_457, graph_WindowIntCharge_0vs1_457_OP] = Create2DPlot(pVal, "wIntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true);

    auto [graph_IntChargeNoW_0vs1_457, graph_IntChargeNoW_0vs1_457_OP] = Create2DPlot(pVal, "IntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, false);
    //0.5GeV!
    pVal = 0.5;
    for(Int_t i=0; i<4; ++i){
        ACT01complete[i] = TOF_435[i];
        ACT01complete[i+4] = TOF_435[i+4];
        ACT01complete[9+i] = ACT01_435[i];
    }
    ACT01complete[8] = pbglass_435;
    ACT01complete[13] = nullptr;
    
    auto [graph_IntCharge_0vs1_435, graph_IntCharge_0vs1_435_OP] = Create2DPlot(pVal, "IntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true);

    auto [graph_WindowIntCharge_0vs1_435, graph_WindowIntCharge_0vs1_435_OP] = Create2DPlot(pVal, "wIntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true);

    auto [graph_IntChargeNoW_0vs1_435, graph_IntChargeNoW_0vs1_435_OP] = Create2DPlot(pVal, "IntCharge", ACT01complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, false);
        
    /////////////////////////////RUN 452
    pVal = 1.;

    //TTree *ACT01T0tree[11];
    for(Int_t i=0;i<8;i++){
        ACT01T0tree[i] = TOF_452[i];
    } 
    ACT01T0tree[8] = pbglass_452;
    ACT01T0tree[9] = ACT01_452[0];
    ACT01T0tree[10] = nullptr;
    
    Create1DPlot("TOFL",ACT01T0tree,{},pVal);
    Create1DPlot("TOFL1",ACT01T0tree,{},pVal);
  
    ACT01T0tree[9] = ACT01_452[1];
    Create1DPlot("TOFR",ACT01T0tree,{},pVal);
    Create1DPlot("TOFR1",ACT01T0tree,{},pVal);

    for(Int_t i=0;i<8;i++){
        ACT01T0tree[i] = TOF_452[i];
    } 
    ACT01T0tree[8] = pbglass_452;
    ACT01T0tree[9] = ACT01_452[2];
    ACT01T0tree[10] = nullptr;

    Create1DPlot("TOFL",ACT01T0tree,{},pVal,true);

    ACT01T0tree[9] = ACT01_452[3];
    Create1DPlot("TOFR",ACT01T0tree,{},pVal,true);

    //////////////////RUN 457
    pVal = 0.7;
    for(Int_t i=0;i<8;i++){
        ACT01T0tree[i] = TOF_457[i];
    } 
    ACT01T0tree[8] = pbglass_457;
    ACT01T0tree[9] = ACT01_457[0];
    ACT01T0tree[10] = nullptr;

    Create1DPlot("TOFL",ACT01T0tree,{},pVal);

    ACT01T0tree[9] = ACT01_457[1];
    Create1DPlot("TOFR",ACT01T0tree,{},pVal);
    
    for(Int_t i=0;i<8;i++){
        ACT01T0tree[i] = TOF_457[i];
    } 
    ACT01T0tree[8] = pbglass_457;
    ACT01T0tree[9] = ACT01_457[2];
    ACT01T0tree[10] = nullptr;

    Create1DPlot("TOFL",ACT01T0tree,{},pVal,true);

    ACT01T0tree[9] = ACT01_452[3];
    Create1DPlot("TOFR",ACT01T0tree,{},pVal,true);
    
    ///////////////////RUN 435
    pVal = 0.5;
    for(Int_t i=0;i<8;i++){
        ACT01T0tree[i] = TOF_435[i];
    } 
    ACT01T0tree[8] = pbglass_435;
    ACT01T0tree[9] = ACT01_435[0];
    ACT01T0tree[10] = nullptr;

    Create1DPlot("TOFL",ACT01T0tree,{},pVal);

    ACT01T0tree[9] = ACT01_435[1];
    Create1DPlot("TOFR",ACT01T0tree,{},pVal);
    
    for(Int_t i=0;i<8;i++){
        ACT01T0tree[i] = TOF_435[i];
    } 
    ACT01T0tree[8] = pbglass_435;
    ACT01T0tree[9] = ACT01_435[2];
    ACT01T0tree[10] = nullptr;

    Create1DPlot("TOFL",ACT01T0tree,{},pVal,true);

    ACT01T0tree[9] = ACT01_435[3];
    Create1DPlot("TOFR",ACT01T0tree,{},pVal,true);
     
    const Int_t numPoints = 50;
    Double_t purity_fixedb[numPoints],eff_fixedb[numPoints],purity_fixedbIntC[numPoints],eff_fixedbIntC[numPoints],purity_fixedbWIntC[numPoints],eff_fixedbWIntC[numPoints];
    Double_t purity_fixedm[numPoints], eff_fixedm[numPoints],purity_fixedmIntC[numPoints], eff_fixedmIntC[numPoints],purity_fixedmWIntC[numPoints], eff_fixedmWIntC[numPoints];
    Double_t eff_fixedbBG[numPoints], eff_fixedbIntCBG[numPoints], eff_fixedbWIntCBG[numPoints], eff_fixedmBG[numPoints], eff_fixedmIntCBG[numPoints], eff_fixedmWIntCBG[numPoints];
    Double_t purity_fixedbBG[numPoints], purity_fixedbIntCBG[numPoints], purity_fixedbWIntCBG[numPoints], purity_fixedmBG[numPoints], purity_fixedmIntCBG[numPoints], purity_fixedmWIntCBG[numPoints];

    Double_t fixedb = 0.0040;
    Double_t fixedm = -0.0825;
    Double_t bVal[numPoints], mVal[numPoints]; 

    Double_t purityAll_IntC[numPoints*numPoints];
    Double_t effAll_IntC[numPoints*numPoints];
    Double_t purityAll_wIntC[numPoints*numPoints];
    Double_t effAll_wIntC[numPoints*numPoints];

    Double_t purityAll_noE_IntC[numPoints*numPoints];
    Double_t effAll_noE_IntC[numPoints*numPoints];
    Double_t purityAll_noE_wIntC[numPoints*numPoints];
    Double_t effAll_noE_wIntC[numPoints*numPoints];

    //generate m values from -1 to -0.001
    for (Int_t i = 0; i < numPoints; ++i) {
        mVal[i] = -1.0 + i * (0.999 / (numPoints - 1));
    }
    //genrate b values from 0.01 to 0.35
    for (Int_t i = 0; i < numPoints; ++i) {
        bVal[i] = 0.01 + i * (0.35 - 0.01) / (numPoints - 1);
    }
    /*
    for (Int_t i = 0; i < numPoints; ++i) {
        //calculatePurityAndEfficiency(graph_MaxVoltage_0vs1_452, graph_MaxVoltage_0vs1_452_OP, mVal[i], fixedb, purity_fixedb[i], eff_fixedb[i], eff_fixedbBG[i], purity_fixedbBG[i]);
        //calculatePurityAndEfficiency(graph_MaxVoltage_0vs1_452, graph_MaxVoltage_0vs1_452_OP, fixedm, bVal[i], purity_fixedm[i], eff_fixedm[i], eff_fixedmBG[i], purity_fixedmBG[i]);
        calculatePurityAndEfficiency(graph_IntCharge_0vs1_457, graph_IntCharge_0vs1_457_OP, mVal[i], fixedb, purity_fixedbIntC[i], eff_fixedbIntC[i], eff_fixedbIntCBG[i], purity_fixedbIntCBG[i]);
        calculatePurityAndEfficiency(graph_IntCharge_0vs1_457, graph_IntCharge_0vs1_457_OP, fixedm, bVal[i], purity_fixedmIntC[i], eff_fixedmIntC[i], eff_fixedmIntCBG[i], purity_fixedmIntCBG[i]);
        calculatePurityAndEfficiency(graph_WindowIntCharge_0vs1_457, graph_WindowIntCharge_0vs1_457_OP, mVal[i], fixedb, purity_fixedbWIntC[i], eff_fixedbWIntC[i], eff_fixedbWIntCBG[i], purity_fixedbWIntCBG[i]);
        calculatePurityAndEfficiency(graph_WindowIntCharge_0vs1_457, graph_WindowIntCharge_0vs1_457_OP, fixedm, bVal[i], purity_fixedmWIntC[i], eff_fixedmWIntC[i], eff_fixedmWIntCBG[i], purity_fixedmWIntCBG[i]);
    }
    
    TCanvas *cPE = new TCanvas();
    TCanvas *cPE2 = new TCanvas();
    TCanvas *cPE3 = new TCanvas();
    TCanvas *cPE4 = new TCanvas();
 
    TGraph *PurrEff_fixedbIntC = new TGraph();
    TGraph *PurrEff_fixedmIntC = new TGraph();
    TGraph *fixedbcopyIntC = new TGraph();
    TGraph *fixedmcopyIntC = new TGraph();
    TGraph *fixedbcopyIntCBG = new TGraph();
    TGraph *fixedmcopyIntCBG = new TGraph();

    TGraph *PurrEff_fixedbWIntC = new TGraph();
    TGraph *PurrEff_fixedmWIntC = new TGraph();
    TGraph *fixedbcopyWIntC = new TGraph();
    TGraph *fixedmcopyWIntC = new TGraph();
    TGraph *fixedbcopyWIntCBG = new TGraph();
    TGraph *fixedmcopyWIntCBG = new TGraph();
    TGraph *fixedbcopyIntCBGpur = new TGraph();
    TGraph *fixedmcopyIntCBGpur = new TGraph();
    TGraph *fixedbcopyWIntCBGpur = new TGraph();
    TGraph *fixedmcopyWIntCBGpur = new TGraph();

    TLegend *leg = new TLegend(0.1,0.1,0.35,0.30);
    TLegend *leg2 = new TLegend(0.1,0.1,0.35,0.30);
    TLegend *leg3 = new TLegend(0.1,0.1,0.35,0.30);
    TLegend *leg4 = new TLegend(0.1,0.1,0.35,0.30);


    PurrEff_fixedbIntC->SetTitle(Form("Purity and Efficiency of IntCharge w/ fixed b=%.3f",fixedb));
    PurrEff_fixedmIntC->SetTitle(Form("Purity and Efficiency of IntCharge w/ fixed m=%.3f",fixedm));
    PurrEff_fixedbWIntC->SetTitle(Form("Purity and Efficiency of WindowIntCharge w/ fixed b=%.3f",fixedb));
    PurrEff_fixedmWIntC->SetTitle(Form("Purity and Efficiency of WindowIntCharge w/ fixed m=%.3f",fixedm));


    PurrEff_fixedbIntC->GetXaxis()->SetTitle("slope m values");
    PurrEff_fixedbIntC->GetYaxis()->SetTitle("Percentage(%)");
    PurrEff_fixedmIntC->GetXaxis()->SetTitle("intercept b values");
    PurrEff_fixedmIntC->GetYaxis()->SetTitle("Percentage(%)");
    PurrEff_fixedbWIntC->GetXaxis()->SetTitle("slope m values");
    PurrEff_fixedbWIntC->GetYaxis()->SetTitle("Percentage(%)");
    PurrEff_fixedmWIntC->GetXaxis()->SetTitle("intercept b values");
    PurrEff_fixedmWIntC->GetYaxis()->SetTitle("Percentage(%)");


    PurrEff_fixedbIntC->SetMarkerStyle(kFullSquare);
    PurrEff_fixedbIntC->SetMarkerColor(kBlue);
    PurrEff_fixedmIntC->SetMarkerStyle(kFullSquare);
    PurrEff_fixedmIntC->SetMarkerColor(kBlue);
    PurrEff_fixedbWIntC->SetMarkerStyle(kFullSquare);
    PurrEff_fixedbWIntC->SetMarkerColor(kBlue);
    PurrEff_fixedmWIntC->SetMarkerStyle(kFullSquare);
    PurrEff_fixedmWIntC->SetMarkerColor(kBlue);

    fixedbcopyIntC->SetMarkerStyle(kFullTriangleUp);
    fixedmcopyIntC->SetMarkerColor(kBlue);
    fixedmcopyIntC->SetMarkerStyle(kFullTriangleUp);
    fixedbcopyIntC->SetMarkerColor(kBlue);
    fixedbcopyWIntC->SetMarkerStyle(kFullTriangleUp);
    fixedmcopyWIntC->SetMarkerColor(kBlue);
    fixedmcopyWIntC->SetMarkerStyle(kFullTriangleUp);
    fixedbcopyWIntC->SetMarkerColor(kBlue);

    fixedbcopyIntCBG->SetMarkerStyle(kFullStar);
    fixedmcopyIntCBG->SetMarkerColor(kRed);
    fixedmcopyIntCBG->SetMarkerStyle(kFullStar);
    fixedbcopyIntCBG->SetMarkerColor(kRed);
    fixedbcopyWIntCBG->SetMarkerStyle(kFullStar);
    fixedmcopyWIntCBG->SetMarkerColor(kRed);
    fixedmcopyWIntCBG->SetMarkerStyle(kFullStar);
    fixedbcopyWIntCBG->SetMarkerColor(kRed);

    fixedbcopyIntCBGpur->SetMarkerStyle(kCircle);
    fixedmcopyIntCBGpur->SetMarkerColor(kRed);
    fixedmcopyIntCBGpur->SetMarkerStyle(kCircle);
    fixedbcopyIntCBGpur->SetMarkerColor(kRed);
    fixedbcopyWIntCBGpur->SetMarkerStyle(kCircle);
    fixedmcopyWIntCBGpur->SetMarkerColor(kRed);
    fixedmcopyWIntCBGpur->SetMarkerStyle(kCircle);
    fixedbcopyWIntCBGpur->SetMarkerColor(kRed);
    

    for(Int_t j=0; j<numPoints; ++j){
        PurrEff_fixedbIntC->SetPoint(j,mVal[j],purity_fixedbIntC[j]*100);
        fixedbcopyIntC->SetPoint(j,mVal[j],eff_fixedbIntC[j]*100);
        fixedbcopyIntCBG->SetPoint(j,mVal[j],eff_fixedbIntCBG[j]*100);
        fixedbcopyIntCBGpur->SetPoint(j,mVal[j],purity_fixedbIntCBG[j]*100);

        PurrEff_fixedmIntC->SetPoint(j,bVal[j],purity_fixedmIntC[j]*100);
        fixedmcopyIntC->SetPoint(j,bVal[j],eff_fixedmIntC[j]*100);
        fixedmcopyIntCBG->SetPoint(j,bVal[j],eff_fixedmIntCBG[j]*100);
        fixedmcopyIntCBGpur->SetPoint(j,bVal[j],purity_fixedmIntCBG[j]*100);

        PurrEff_fixedbWIntC->SetPoint(j,mVal[j],purity_fixedbWIntC[j]*100);
        fixedbcopyWIntC->SetPoint(j,mVal[j],eff_fixedbWIntC[j]*100);
        fixedbcopyWIntCBG->SetPoint(j,mVal[j],eff_fixedbWIntCBG[j]*100);
        fixedbcopyWIntCBGpur->SetPoint(j,mVal[j],purity_fixedbWIntCBG[j]*100);

        PurrEff_fixedmWIntC->SetPoint(j,bVal[j],purity_fixedmWIntC[j]*100);
        fixedmcopyWIntC->SetPoint(j,bVal[j],eff_fixedmWIntC[j]*100);
        fixedmcopyWIntCBG->SetPoint(j,bVal[j],eff_fixedmWIntCBG[j]*100);
        fixedmcopyWIntCBGpur->SetPoint(j,bVal[j],purity_fixedmWIntCBG[j]*100);

    }

    PurrEff_fixedbIntC->GetYaxis()->SetRangeUser(0,100);
    PurrEff_fixedbWIntC->GetYaxis()->SetRangeUser(0,100);

    leg->AddEntry(PurrEff_fixedbIntC, "Purity of Signal(#mu/#pi)", "p");
    leg->AddEntry(fixedbcopyIntC, "Efficiency of Signal(#mu/#pi)", "p");
    leg->AddEntry(fixedbcopyIntCBG, "Efficiency of Signal(e)", "p");
    leg->AddEntry(fixedbcopyIntCBGpur, "Purity of Signal(e)", "p");
    leg2->AddEntry(PurrEff_fixedbWIntC, "Purity of Signal(#mu/#pi)", "p");
    leg2->AddEntry(fixedbcopyWIntC, "Efficiency of Signal(#mu/#pi)", "p");
    leg2->AddEntry(fixedbcopyWIntCBG, "Efficiency of Signal(e)", "p");
    leg2->AddEntry(fixedbcopyWIntCBGpur, "Pur of Signal(e)", "p");

    cPE->cd();
    PurrEff_fixedbIntC->Draw("AP");
    fixedbcopyIntC->Draw("Same P");
    fixedbcopyIntCBG->Draw("Same P");
    fixedbcopyIntCBGpur->Draw("Same P");
    leg->Draw();
    cPE->SaveAs(Form("./graphjpg/IntC_Pu_Eff_fixedb_%.1fGeV_windowON.jpg",0.7));

    cPE2->cd();
    PurrEff_fixedbWIntC->Draw("AP");
    fixedbcopyWIntC->Draw("Same P");
    fixedbcopyWIntCBG->Draw("Same P");
    fixedbcopyWIntCBGpur->Draw("Same P");
    leg2->Draw();
    cPE2->SaveAs(Form("./graphjpg/wIntC_Pu_Eff_fixedb_%.1f GeV.jpg",0.7));

    PurrEff_fixedmIntC->GetYaxis()->SetRangeUser(0,100);
    PurrEff_fixedmWIntC->GetYaxis()->SetRangeUser(0,100);

    leg3->AddEntry(PurrEff_fixedmIntC, "Purity of Signal(#mu/#pi)", "p");
    leg3->AddEntry(fixedmcopyIntC, "Efficiency of Signal(#mu/#pi)", "p");
    leg3->AddEntry(fixedmcopyIntCBG, "Efficiency of Signal(e)", "p");
    leg3->AddEntry(fixedmcopyIntCBGpur, "Purity of Signal(e)", "p");

    leg4->AddEntry(PurrEff_fixedmWIntC, "Purity of Signal(#mu/#pi)", "p");
    leg4->AddEntry(fixedmcopyWIntC, "Efficiency of Signal(#mu/#pi)", "p");
    leg4->AddEntry(fixedmcopyWIntCBG, "Efficiency of Signal(e)", "p");
    leg4->AddEntry(fixedmcopyWIntCBGpur, "Purity of Signal(e)", "p");

    cPE3->cd();
    PurrEff_fixedmIntC->Draw("AP");
    fixedmcopyIntC->Draw("Same P");
    fixedmcopyIntCBG->Draw("Same P");
    fixedmcopyIntCBGpur->Draw("Same P");

    leg3->Draw();
    cPE3->SaveAs(Form("./graphjpg/IntC_Pu_Eff_fixedm_%.1fGeV_windowON.jpg",0.7));

    cPE4->cd();
    PurrEff_fixedmWIntC->Draw("AP");
    fixedmcopyWIntC->Draw("Same P");
    fixedmcopyWIntCBG->Draw("Same P");
    fixedmcopyWIntCBGpur->Draw("Same P");
    leg4->Draw();
    cPE4->SaveAs(Form("./graphjpg/wIntC_Pu_Eff_fixedm_%.1fGeV.jpg",0.7));
    */
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            Double_t m = mVal[i];
            Double_t b = bVal[j];
            Double_t purity, eff;
            Double_t puritybg, effbg;
            Double_t purity2, eff2;
            Double_t puritybg2, effbg2;                
            
            calculatePurityAndEfficiency(graph_IntCharge_0vs1_457, graph_IntCharge_0vs1_457_OP, m, b, purity, eff, effbg, puritybg);
            calculatePurityAndEfficiency(graph_WindowIntCharge_0vs1_457, graph_WindowIntCharge_0vs1_457_OP, m, b, purity2, eff2, effbg2, puritybg2);

            purityAll_IntC[i * numPoints + j] = puritybg;
            effAll_IntC[i* numPoints + j] = effbg;
            purityAll_wIntC[i * numPoints + j] = puritybg2;
            effAll_wIntC[i* numPoints + j] = effbg2;

            purityAll_noE_IntC[i * numPoints + j] = purity;
            effAll_noE_IntC[i* numPoints + j] = eff;
            purityAll_noE_wIntC[i * numPoints + j] = purity2;
            effAll_noE_wIntC[i* numPoints + j] = eff2;
            //cout << "For Event " << i*numPoints+j << " the electron purity and efficiency are " << purityAll_IntC[i * numPoints + j] << " and " << effAll_IntC[i * numPoints + j] << " while the signal purity and eff are " << purityAll_noE_IntC[i * numPoints + j] << " and " << effAll_noE_IntC[i * numPoints + j] << endl;
        }
    }

    // Create 2D histogram
    TH2D* pureff_IntC2Dhist = new TH2D("pureff_IntC2Dhist", "Run 457 Purity*Efficiency of the electron signal in zone2 with different slopes and intercepts using IntCharge", numPoints, mVal[0], mVal[numPoints - 1], numPoints, bVal[0], bVal[numPoints - 1]);
    pureff_IntC2Dhist->GetXaxis()->SetTitle("Slope m values");
    pureff_IntC2Dhist->GetYaxis()->SetTitle("Intercept b values");
    pureff_IntC2Dhist->GetZaxis()->SetTitle("Purity*Efficiency %");
    TH2D* pureff_wIntC2Dhist = new TH2D("pureff_wIntC2Dhist", "Run 457 Purity*Efficiency of the electron signal in zone2 with different slopes and intercepts using WindowIntCharge", numPoints, mVal[0], mVal[numPoints - 1], numPoints, bVal[0], bVal[numPoints - 1]);
    pureff_wIntC2Dhist->GetXaxis()->SetTitle("Slope m values");
    pureff_wIntC2Dhist->GetYaxis()->SetTitle("Intercept b values");
    pureff_wIntC2Dhist->GetZaxis()->SetTitle("Purity*Efficiency %");

    pureff_IntC2Dhist->SetStats(0);
    pureff_wIntC2Dhist->SetStats(0);

    // Fill histogram
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            pureff_IntC2Dhist->SetBinContent(i + 1, j + 1, (purityAll_IntC[i * numPoints + j]*effAll_IntC[i * numPoints + j])*100); 
            pureff_wIntC2Dhist->SetBinContent(i + 1, j + 1, (purityAll_wIntC[i * numPoints + j]*effAll_wIntC[i * numPoints + j])*100);                   
        }
    }

    // Draw the histogram
    TCanvas *c= new TCanvas("c", "Purity vs Slope and Intercept", 800, 600);
    pureff_IntC2Dhist->Draw("colz");
    pureff_IntC2Dhist->SetContour(100);  // Set number of contours for better visibility
    c->SetRightMargin(0.15);  // Adjust right margin to make room for z-axis title
    c->SaveAs("./graphjpg/PurEffIntC_2Dhist.jpg");
    TCanvas* c2 = new TCanvas("c2", "Purity vs Slope and Intercept", 800, 600);
    pureff_wIntC2Dhist->Draw("colz");
    pureff_wIntC2Dhist->SetContour(100);  // Set number of contours for better visibility
    c2->SetRightMargin(0.15);  // Adjust right margin to make room for z-axis title
    c2->SaveAs("./graphjpg/PurEffwIntC_2Dhist.jpg");

    Int_t maxBin = pureff_IntC2Dhist->GetMaximumBin();
    Int_t xBin, yBin, dummy;
    pureff_IntC2Dhist->GetBinXYZ(maxBin, xBin, yBin, dummy);
    cout << "MaxZ in IntC is " << pureff_IntC2Dhist->GetMaximum() << " and its bin in x is " <<  xBin << " and in y is " << yBin << endl;   
    maxBin = pureff_wIntC2Dhist->GetMaximumBin();
    pureff_wIntC2Dhist->GetBinXYZ(maxBin, xBin, yBin, dummy);
    cout << "MaxZ in wIntC is " << pureff_wIntC2Dhist->GetMaximum() << " and its bin in x is " <<  xBin << " and in y is " << yBin << endl;   

    //Using Efficiency of non-e signal instead
    TH2D* pureff_noE_IntC2Dhist = new TH2D("pureff_noE_IntC2Dhist", "Run 457 Purity*Efficiency of the NON-electron signal in zone1 with different slopes and intercepts using IntCharge", numPoints, mVal[0], mVal[numPoints - 1], numPoints, bVal[0], bVal[numPoints - 1]);
    pureff_noE_IntC2Dhist->GetXaxis()->SetTitle("Slope m values");
    pureff_noE_IntC2Dhist->GetYaxis()->SetTitle("Intercept b values");
    pureff_noE_IntC2Dhist->GetZaxis()->SetTitle("Purity*Efficiency %");
    TH2D* pureff_noE_wIntC2Dhist = new TH2D("pureff_noE_wIntC2Dhist", "Run 457 Purity*Efficiency of the NON-electron signal in zone1 with different slopes and intercepts using WindowIntCharge", numPoints, mVal[0], mVal[numPoints - 1], numPoints, bVal[0], bVal[numPoints - 1]);
    pureff_noE_wIntC2Dhist->GetXaxis()->SetTitle("Slope m values");
    pureff_noE_wIntC2Dhist->GetYaxis()->SetTitle("Intercept b values");
    pureff_noE_wIntC2Dhist->GetZaxis()->SetTitle("Purity*Efficiency %");

    pureff_noE_IntC2Dhist->SetStats(0);
    pureff_noE_wIntC2Dhist->SetStats(0);

    // Fill histogram
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            pureff_noE_IntC2Dhist->SetBinContent(i + 1, j + 1, (purityAll_noE_IntC[i * numPoints + j]*effAll_noE_IntC[i * numPoints + j])*100); 
            pureff_noE_wIntC2Dhist->SetBinContent(i + 1, j + 1, (purityAll_noE_wIntC[i * numPoints + j]*effAll_noE_wIntC[i * numPoints + j])*100);                   
        }
    }

    // Draw the histogram
    TCanvas *c3= new TCanvas("c3", "Purity vs Slope and Intercept", 800, 600);
    pureff_noE_IntC2Dhist->Draw("colz");
    pureff_noE_IntC2Dhist->SetContour(100);  // Set number of contours for better visibility
    c3->SetRightMargin(0.15);  // Adjust right margin to make room for z-axis title
    c3->SaveAs("./graphjpg/PurEff_noE_IntC_2Dhist.jpg");
    TCanvas* c4 = new TCanvas("c4", "Purity vs Slope and Intercept", 800, 600);
    pureff_noE_wIntC2Dhist->Draw("colz");
    pureff_noE_wIntC2Dhist->SetContour(100);  // Set number of contours for better visibility
    c4->SetRightMargin(0.15);  // Adjust right margin to make room for z-axis title
    c4->SaveAs("./graphjpg/PurEff_noE_wIntC_2Dhist.jpg");

    maxBin = pureff_noE_IntC2Dhist->GetMaximumBin();
    pureff_noE_IntC2Dhist->GetBinXYZ(maxBin, xBin, yBin, dummy);
    cout << "MaxZ in IntC is " << pureff_noE_IntC2Dhist->GetMaximum() << " and its bin in x is " <<  xBin << " and in y is " << yBin << endl;   
    maxBin = pureff_noE_wIntC2Dhist->GetMaximumBin();
    pureff_noE_wIntC2Dhist->GetBinXYZ(maxBin, xBin, yBin, dummy);
    cout << "MaxZ in wIntC is " << pureff_noE_wIntC2Dhist->GetMaximum() << " and its bin in x is " <<  xBin << " and in y is " << yBin << endl;   

    /////NO WINDOW BEST FIT FOR ALL.
    ReturnBestFitValue(graph_IntChargeNoW_0vs1_435,graph_IntChargeNoW_0vs1_435_OP,"in no window IntC at 0.5GeV is ", mVal, bVal, numPoints);
    ReturnBestFitValue(graph_IntChargeNoW_0vs1_457,graph_IntChargeNoW_0vs1_457_OP,"in no window IntC at 0.7GeV is ", mVal, bVal, numPoints);
    ReturnBestFitValue(graph_IntChargeNoW_0vs1_452,graph_IntChargeNoW_0vs1_452_OP,"in no window IntC at 1GeV is ", mVal, bVal, numPoints);
    ReturnBestFitValue(graph_IntCharge_0vs1_435,graph_IntCharge_0vs1_435_OP,"in IntC at 0.5GeV is ", mVal, bVal, numPoints);
    ReturnBestFitValue(graph_IntCharge_0vs1_452,graph_IntCharge_0vs1_452_OP,"in IntC at 1GeV is ", mVal, bVal, numPoints);
    ReturnBestFitValue(graph_WindowIntCharge_0vs1_435,graph_WindowIntCharge_0vs1_435_OP,"in wIntC at 0.5GeV is ", mVal, bVal, numPoints);
    ReturnBestFitValue(graph_WindowIntCharge_0vs1_452,graph_WindowIntCharge_0vs1_452_OP,"in wIntC at 1GeV is ", mVal, bVal, numPoints);

    
    //////////////// Best fit for no window /////////////////

    TH2D* pureff_IntC2DhistnoW = new TH2D("pureff_IntC2DhistnoW", "Run 457 Purity*Efficiency of the electron signal in zone2 with different slopes and intercepts using IntCharge", numPoints, mVal[0], mVal[numPoints - 1], numPoints, bVal[0], bVal[numPoints - 1]);
    Double_t purityAll_noWIntC[numPoints*numPoints];
    Double_t effAll_noWIntC[numPoints*numPoints];
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            Double_t m = mVal[i];
            Double_t b = bVal[j];
            Double_t purity, eff;
            Double_t puritybg, effbg;
            Double_t purity2, eff2;
            Double_t puritybg2, effbg2;                
            
            calculatePurityAndEfficiency(graph_IntChargeNoW_0vs1_457, graph_IntChargeNoW_0vs1_457_OP, m, b, purity, eff, effbg, puritybg);
            purityAll_noWIntC[i * numPoints + j] = puritybg;
            effAll_noWIntC[i* numPoints + j] = effbg;
        }
    }
    // Fill histogram
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            pureff_IntC2DhistnoW->SetBinContent(i + 1, j + 1, (purityAll_noWIntC[i * numPoints + j]*effAll_noWIntC[i * numPoints + j])*100); 
        }
    }
    maxBin = pureff_IntC2DhistnoW->GetMaximumBin();
    pureff_IntC2DhistnoW->GetBinXYZ(maxBin, xBin, yBin, dummy);
    cout << "MaxZ in no window IntC at 0.7GeV is " << pureff_IntC2DhistnoW->GetMaximum() << " and its bin in x is " <<  xBin << " and in y is " << yBin << endl;   

    //////////SAME ANALYSIS BUT WITH 0.5GeV////////////////
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            Double_t m = mVal[i];
            Double_t b = bVal[j];
            Double_t purity, eff;
            Double_t puritybg, effbg;
            Double_t purity2, eff2;
            Double_t puritybg2, effbg2;                
            
            calculatePurityAndEfficiency(graph_IntCharge_0vs1_435, graph_IntCharge_0vs1_435_OP, m, b, purity, eff, effbg, puritybg);
            calculatePurityAndEfficiency(graph_WindowIntCharge_0vs1_435, graph_WindowIntCharge_0vs1_435_OP, m, b, purity2, eff2, effbg2, puritybg2);

            purityAll_IntC[i * numPoints + j] = puritybg;
            effAll_IntC[i* numPoints + j] = effbg;
            purityAll_wIntC[i * numPoints + j] = puritybg2;
            effAll_wIntC[i* numPoints + j] = effbg2;

            purityAll_noE_IntC[i * numPoints + j] = purity;
            effAll_noE_IntC[i* numPoints + j] = eff;
            purityAll_noE_wIntC[i * numPoints + j] = purity2;
            effAll_noE_wIntC[i* numPoints + j] = eff2;
        }
    }

    // Create 2D histogram
    TH2D* pureff05_IntC2Dhist = new TH2D("pureff05_IntC2Dhist", "Run 435 Purity*Efficiency of the avg of the electron signal in zone2 and the non-e signal in zone1 with different slopes and intercepts using IntCharge", numPoints, mVal[0], mVal[numPoints - 1], numPoints, bVal[0], bVal[numPoints - 1]);
    pureff05_IntC2Dhist->GetXaxis()->SetTitle("Slope m values");
    pureff05_IntC2Dhist->GetYaxis()->SetTitle("Intercept b values");
    pureff05_IntC2Dhist->GetZaxis()->SetTitle("[(Purity*Efficiency)_e%+(Purity*Efficiency)_nonE]/2 %");
    TH2D* pureff05_wIntC2Dhist = new TH2D("pureff05_wIntC2Dhist", "Run 435 Purity*Efficiency of the avg of the electron signal in zone2 and the non-e signal in zone1 with different slopes and intercepts using WindowIntCharge", numPoints, mVal[0], mVal[numPoints - 1], numPoints, bVal[0], bVal[numPoints - 1]);
    pureff05_wIntC2Dhist->GetXaxis()->SetTitle("Slope m values");
    pureff05_wIntC2Dhist->GetYaxis()->SetTitle("Intercept b values");
    pureff05_wIntC2Dhist->GetZaxis()->SetTitle("[(Purity*Efficiency)_e%+(Purity*Efficiency)_nonE]/2 %");

    pureff05_IntC2Dhist->SetStats(0);
    pureff05_wIntC2Dhist->SetStats(0);

    // Fill histogram
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            pureff05_IntC2Dhist->SetBinContent(i + 1, j + 1, (((purityAll_IntC[i * numPoints + j]*effAll_IntC[i * numPoints + j])+(purityAll_noE_IntC[i * numPoints + j]*effAll_noE_IntC[i * numPoints + j]))/2)*100); 
            pureff05_wIntC2Dhist->SetBinContent(i + 1, j + 1, (((purityAll_wIntC[i * numPoints + j]*effAll_wIntC[i * numPoints + j])+(purityAll_noE_wIntC[i * numPoints + j]*effAll_noE_wIntC[i * numPoints + j]))/2)*100);                   
        }
    }

    // Draw the histogram
    TCanvas *c5= new TCanvas("c5", "Purity vs Slope and Intercept", 800, 600);
    pureff05_IntC2Dhist->Draw("colz");
    pureff05_IntC2Dhist->SetContour(100);  // Set number of contours for better visibility
    c5->SetRightMargin(0.15);  // Adjust right margin to make room for z-axis title
    c5->SaveAs("./graphjpg/PurEff05IntC_2Dhist.jpg");
    TCanvas* c6 = new TCanvas("c6", "Purity vs Slope and Intercept", 800, 600);
    pureff05_wIntC2Dhist->Draw("colz");
    pureff05_wIntC2Dhist->SetContour(100);  // Set number of contours for better visibility
    c6->SetRightMargin(0.15);  // Adjust right margin to make room for z-axis title
    c6->SaveAs("./graphjpg/PurEff05wIntC_2Dhist.jpg");

    maxBin = pureff05_IntC2Dhist->GetMaximumBin();
    pureff05_IntC2Dhist->GetBinXYZ(maxBin, xBin, yBin, dummy);
    cout << "MaxZ in IntC 0.5GeV no Window is " << pureff05_IntC2Dhist->GetMaximum() << " and its bin in x is " <<  xBin << " and in y is " << yBin << endl;   
    maxBin = pureff05_wIntC2Dhist->GetMaximumBin();
    pureff05_wIntC2Dhist->GetBinXYZ(maxBin, xBin, yBin, dummy);
    cout << "MaxZ in wIntC is " << pureff05_wIntC2Dhist->GetMaximum() << " and its bin in x is " <<  xBin << " and in y is " << yBin << endl;   
    //////NO WINDOW/////
    TH2D* pureff05_IntC2DhistnoW = new TH2D("pureff05_IntC2DhistnoW", "Run 435 Purity*Efficiency of the avg of the electron signal in zone2 and the non-e signal in zone1 with different slopes and intercepts using IntCharge and no Window", numPoints, mVal[0], mVal[numPoints - 1], numPoints, bVal[0], bVal[numPoints - 1]);
    pureff05_IntC2DhistnoW->GetXaxis()->SetTitle("Slope m values");
    pureff05_IntC2DhistnoW->GetYaxis()->SetTitle("Intercept b values");
    pureff05_IntC2DhistnoW->GetZaxis()->SetTitle("[(Purity*Efficiency)_e%+(Purity*Efficiency)_nonE]/2 %");
    Double_t purityAll_noE_noWIntC[numPoints*numPoints];
    Double_t effAll_noE_noWIntC[numPoints*numPoints];
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            Double_t m = mVal[i];
            Double_t b = bVal[j];
            Double_t purity, eff;
            Double_t puritybg, effbg;
            Double_t purity2, eff2;
            Double_t puritybg2, effbg2;                
            
            calculatePurityAndEfficiency(graph_IntChargeNoW_0vs1_435, graph_IntChargeNoW_0vs1_435_OP, m, b, purity, eff, effbg, puritybg);
            purityAll_noWIntC[i * numPoints + j] = puritybg;
            effAll_noWIntC[i* numPoints + j] = effbg;
            purityAll_noE_noWIntC[i * numPoints + j] = purity;
            effAll_noE_noWIntC[i* numPoints + j] = eff;

        }
    }
    // Fill histogram
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            pureff05_IntC2DhistnoW->SetBinContent(i + 1, j + 1, (((purityAll_noWIntC[i * numPoints + j]*effAll_noWIntC[i * numPoints + j])+(purityAll_noWIntC[i * numPoints + j]*effAll_noWIntC[i * numPoints + j]))/2)*100); 
        }
    }
    TCanvas* c7 = new TCanvas("c7", "Purity vs Slope and Intercept", 800, 600);
    pureff05_IntC2DhistnoW->Draw("colz");
    pureff05_IntC2DhistnoW->SetContour(100);  // Set number of contours for better visibility
    c7->SetRightMargin(0.15);  // Adjust right margin to make room for z-axis title
    c7->SaveAs("./graphjpg/PurEff05wIntCnoW_2Dhist.jpg");

    maxBin = pureff05_IntC2DhistnoW->GetMaximumBin();
    pureff05_IntC2DhistnoW->GetBinXYZ(maxBin, xBin, yBin, dummy);
    cout << "MaxZ in no window IntC is " << pureff05_IntC2DhistnoW->GetMaximum() << " and its bin in x is " <<  xBin << " and in y is " << yBin << endl;   
    
    /*
    TTree *Etree[11];
    for (Int_t i=0; i<8; i++){
        Etree[i] = TOF_464[i];
    }
    Etree[8] = pbglass_464;
    Etree[9] = ACT01_464[0];
    Etree[10] = nullptr;

    TCanvas *MaxVolt0L_vOFF = new TCanvas("MVvsST_1GeVveto_OFF_Leftc","title",800,600);
    TH2F *MVvsST_H2D_OFF_Lcorr = new TH2F("MVvsST_H2D_OFF_Lcorr","run464 p=1GeV, veto OFF signal time ACT0L-T0 vs ACT0L Max Volt w/ tt corrections;ACT0L-T0 signal (ns);ACT0L Max Volt (V)",200,-50,150,60,-0.5,2.5);
    Make2DHistEonly(MaxVolt0L_vOFF,MVvsST_H2D_OFF_Lcorr,Etree,{"MaxVoltage"},"AddedBranch");

    Etree[9] = ACT01_464[1];
    TCanvas *MaxVolt0R_vOFF = new TCanvas("MVvsST_1GeVveto_OFF_Rightc","title",800,600);
    TH2F *MVvsST_H2D_OFF_Rcorr = new TH2F("MVvsST_H2D_OFF_Rcorr","run464 p=1GeV, veto OFF signal time ACT0R-T0 vs ACT0R Max Volt w/ tt corrections;ACT0R-T0 signal (ns);ACT0R Max Volt (V)",200,-50,150,60,-0.5,2.5);
    Make2DHistEonly(MaxVolt0R_vOFF,MVvsST_H2D_OFF_Rcorr,Etree,{"MaxVoltage"},"AddedBranch");

    for (Int_t i=0; i<8; i++){
        Etree[i] = TOF_452[i];
    }
    Etree[8] = pbglass_452;
    Etree[9] = ACT01_452[0];

    TCanvas *MaxVolt0L_vON = new TCanvas("MVvsST_1GeVveto_ON_Leftc","title",800,600);
    TH2F *MVvsST_H2D_ON_Lcorr = new TH2F("MVvsST_H2D_ON_Lcorr","run464 p=1GeV, veto ON signal time ACT0L-T0 vs ACT0L Max Volt w/ tt corrections;ACT0L-T0 signal (ns);ACT0L Max Volt (V)",200,-50,150,60,-0.5,2.5);
    Make2DHistEonly(MaxVolt0L_vON,MVvsST_H2D_ON_Lcorr,Etree,{"MaxVoltage"},"AddedBranch");

    Etree[9] = ACT01_452[1];

    TCanvas *MaxVolt0R_vON = new TCanvas("MVvsST_1GeVveto_ON_Rightc","title",800,600);
    TH2F *MVvsST_H2D_ON_Rcorr = new TH2F("MVvsST_H2D_ON_Rcorr","run464 p=1GeV, veto ON signal time ACT0R-T0 vs ACT0R Max Volt w/ tt corrections;ACT0R-T0 signal (ns);ACT0R Max Volt (V)",200,-50,150,60,-0.5,2.5);
    Make2DHistEonly(MaxVolt0R_vON,MVvsST_H2D_ON_Rcorr,Etree,{"MaxVoltage"},"AddedBranch");
    */

    //ACT-23 now.
    TTree *ACT23complete[14];
    ///0.7GeV!  
    pVal = 0.7;
    for(Int_t i=0; i<4; ++i){
        ACT23complete[i] = TOF_457[i];
        ACT23complete[i+4] = TOF_457[i+4];
        ACT23complete[9+i] = ACT23_457[i];
    }
    ACT23complete[8] = pbglass_457;
    ACT23complete[13] = nullptr;
    
    auto [graph_IntCharge_2vs3_457, graph_IntCharge_2vs3_457_OP] = Create2DPlot(pVal, "IntCharge", ACT23complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true, false);

    auto [graph_WindowIntCharge_2vs3_457, graph_WindowIntCharge_2vs3_457_OP] = Create2DPlot(pVal, "wIntCharge", ACT23complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, true, false);

    auto [graph_IntChargeNoW_2vs3_457, graph_IntChargeNoW_2vs3_457_OP] = Create2DPlot(pVal, "IntCharge", ACT23complete, {"IntCharge", "MaxVoltage", "PeakTime"}, limit, false, false);

    TTree *ACT23T0tree[11];
    //////////////////RUN 457
    pVal = 0.7;
    for(Int_t i=0;i<8;i++){
        ACT23T0tree[i] = TOF_457[i];
    } 
    ACT23T0tree[8] = pbglass_457;
    ACT23T0tree[9] = ACT23_457[0];
    ACT23T0tree[10] = nullptr;

    Create1DPlot("TOFL",ACT23T0tree,{},pVal,false,true);

    ACT23T0tree[9] = ACT23_457[1];
    Create1DPlot("TOFR",ACT23T0tree,{},pVal,false,true);

    for(Int_t i=0;i<8;i++){
        ACT23T0tree[i] = TOF_457[i];
    } 
    ACT23T0tree[8] = pbglass_457;
    ACT23T0tree[9] = ACT23_457[2];
    ACT23T0tree[10] = nullptr;
    
    Create1DPlot("TOFL",ACT23T0tree,{},pVal,true,true);

    ACT23T0tree[9] = ACT23_452[3];
    Create1DPlot("TOFR",ACT23T0tree,{},pVal,true,true);
    
}