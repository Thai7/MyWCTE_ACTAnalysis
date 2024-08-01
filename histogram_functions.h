#ifndef HISTOGRAM_FUNCTIONS_H
#define HISTOGRAM_FUNCTIONS_H

#include "functions.h" // Include the first file that contains common functions
#include <map> // for std::map

void Make2DHistEonly(TCanvas *canvas, TH2F *hist, TTree *trees[], initializer_list<const char*> branches, const char* htype){

    // Calculate total number of trees in the array
    Int_t totalIndex = 0;
    while (trees[totalIndex] != nullptr)
    {
        // Process trees[totalIndex]
        totalIndex++;
    }

    // Create histogram for other particles
    TH2F *histOtherParticle = new TH2F(Form("%s_otherParticle", hist->GetName()), Form("%s, Other Particle", hist->GetTitle()), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), hist->GetNbinsY(), hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax());
    
    // Create canvas for histOtherParticle (if needed)
    TCanvas *canvasOtherParticle = new TCanvas(Form("%s_otherP", canvas->GetName()), Form("%s_otherP", canvas->GetName()), 800, 600);

    // Add "only electrons" to the title of hist
    TString histTitle = hist->GetTitle();
    histTitle += ", only electrons";
    hist->SetTitle(histTitle);

    // Initialize and set up branches and branch values
    vector<const char*> allBranches;
    Int_t combinedSize = CalcCombinedSize(branches, allBranches, totalIndex);

    Double_t branchVal[combinedSize][20];
    Int_t nPeaks[totalIndex];
    Double_t WWI[totalIndex];
    Int_t nWindowPeak;

    InitializeBranchVal(branchVal, combinedSize);

    SetupRootTreesAddress(trees, totalIndex, combinedSize, branches, allBranches, branchVal, nPeaks, WWI, nWindowPeak);

    Double_t xvalue;
    Double_t yvalue;
    Long64_t numEntries = trees[0]->GetEntries();
    for (Long64_t entry = 0; entry < numEntries; ++entry) {
        Double_t sumTOF1 = 0.0, sumTOF0 = 0.0;
        Double_t avgTOF1 = 0.0, avgTOF0 = 0.0, finalTOF = 0.0;
        //In order, the lower bound (of the window) are 0L,0R,1L,1R and are in ns.
        Double_t lower_bound[4] = {0.0,0.0,0.0,0.0};
        Double_t upper_bound[4] = {0.0,0.0,0.0,0.0};
        Double_t ACT0ST = 0.0, ACT1ST = 0.0;
        Double_t WWI_avg = 0.0;

        // Check if all nPeaks values are 1 for the current entry
        bool allPeaksOne = true;
        bool eCut = false;
        bool WWCut = false;
        for (Int_t treeIdx = 0; treeIdx < 9; ++treeIdx) {
            trees[treeIdx]->GetEntry(entry);
            if (nPeaks[treeIdx] != 1) {
                allPeaksOne = false;
                break;
            }
        }
        if (allPeaksOne && nWindowPeak==1){
            for(Int_t treeIdx=0; treeIdx<totalIndex; ++treeIdx){
                trees[treeIdx]->GetEntry(entry);
            }
            Double_t PbGlassPEVal;
            PbGlassPEVal = branchVal[totalIndex*1+8][0];
            for (Int_t i=0; i<4; ++i){
                sumTOF1 += branchVal[totalIndex * 0 + (4 + i)][0];
                sumTOF0 += branchVal[totalIndex * 0 + i][0];
                WWI_avg += WWI[i] + WWI[i+4];
            }
            avgTOF1 = sumTOF1/4;
            avgTOF0 = sumTOF0/4;
            finalTOF = avgTOF1-avgTOF0;
            WWI_avg /=8;
            //if (PbGlassPEVal > 10 && finalTOF < 14){ //correct value for windowIntPE
            if (PbGlassPEVal > 0.45 && finalTOF <13){ //correct value for IntCharge, kept IntCharge into the PEVal variable for conciseness...
                eCut = true;
            }
            if (WWI_avg<40){
                WWCut = true;
            }
        }

        if (eCut && WWCut && (strcmp(htype,"wIntPEboth"))==0){ //used for ACT0-TOF0 vs ACT0 WindowIntPE using both ACT0
            yvalue = branchVal[totalIndex * 1 + 9][0] + branchVal[totalIndex * 1 + 10][0];
            Double_t sumTOF0 = 0.0;
            for (Int_t i=0; i<4; ++i){
                sumTOF0+=branchVal[totalIndex * 0 + i][0];
            }
            xvalue = (branchVal[totalIndex * 0 + 9][0] + branchVal[totalIndex * 1 + 10][0])/2 - sumTOF0/4;
            hist->Fill(xvalue,yvalue);
        }
        if (WWCut && !eCut && (strcmp(htype,"wIntPEboth"))==0) {
            // Fill the second histogram
            yvalue = branchVal[totalIndex * 1 + 9][0] + branchVal[totalIndex * 1 + 10][0];
            Double_t sumTOF0 = 0.0;
            for (Int_t i = 0; i < 4; ++i) {
                sumTOF0 += branchVal[totalIndex * 0 + i][0];
            }
            xvalue = (branchVal[totalIndex * 0 + 9][0] + branchVal[totalIndex * 0 + 10][0]) / 2 - sumTOF0 / 4;
            histOtherParticle->Fill(xvalue, yvalue);
        }
        if (eCut && WWCut && (strcmp(htype,"wIntPEsingle"))==0){ //used for ACT0-TOF0 vs ACT0 WindowIntPE using one ACT0
            yvalue = branchVal[totalIndex * 1 + 9][0];
            Double_t sumTOF0 = 0.0;
            for (Int_t i=0; i<4; ++i){
                sumTOF0+=branchVal[totalIndex * 0 + i][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0/4;
            hist->Fill(xvalue,yvalue);
        }
        if (WWCut && !eCut && (strcmp(htype,"wIntPEsingle"))==0) {
            // Fill the second histogram
            yvalue = branchVal[totalIndex * 1 + 9][0];
            Double_t sumTOF0 = 0.0;
            for (Int_t i = 0; i < 4; ++i) {
                sumTOF0 += branchVal[totalIndex * 0 + i][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0 / 4;
            histOtherParticle->Fill(xvalue, yvalue);
        }
        if (eCut && WWCut && (strcmp(htype,"AddedBranch"))==0){ //used for ACT0-TOF0 vs ACT0 IntCharge using one ACT0
            yvalue = branchVal[totalIndex * 2 + 9][0];
            Double_t sumTOF0 = 0.0;
            for (Int_t i=0; i<4; ++i){
                sumTOF0+=branchVal[totalIndex * 0 + i][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0/4;
            hist->Fill(xvalue,yvalue);
        }
        if (WWCut && !eCut && (strcmp(htype,"AddedBranch"))==0) {
            // Fill the second histogram
            yvalue = branchVal[totalIndex * 2 + 9][0];
            Double_t sumTOF0 = 0.0;
            for (Int_t i = 0; i < 4; ++i) {
                sumTOF0 += branchVal[totalIndex * 0 + i][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0 / 4;
            histOtherParticle->Fill(xvalue, yvalue);           
        }
        if (eCut && WWCut && (strcmp(htype,"AddedBranches"))==0){ //used for ACT0-TOF0 vs ACT0 "IntCharge" using both ACT0
            yvalue = branchVal[totalIndex * 2 + 9][0] + branchVal[totalIndex * 2 + 10][0];
            Double_t sumTOF0 = 0.0;
            for (Int_t i=0; i<4; ++i){
                sumTOF0+=branchVal[totalIndex * 0 + i][0];
            }
            xvalue = (branchVal[totalIndex * 0 + 9][0] + branchVal[totalIndex * 0 + 10][0])/2 - sumTOF0/4;
            hist->Fill(xvalue,yvalue);
        }
        if (WWCut && !eCut && (strcmp(htype,"AddedBranches"))==0) {
            // Fill the second histogram
            yvalue = branchVal[totalIndex * 2 + 9][0] + branchVal[totalIndex * 1 + 10][0];
            Double_t sumTOF0 = 0.0;
            for (Int_t i = 0; i < 4; ++i) {
                sumTOF0 += branchVal[totalIndex * 0 + i][0];
            }
            xvalue = (branchVal[totalIndex * 0 + 9][0] + branchVal[totalIndex * 0 + 10][0]) / 2 - sumTOF0 / 4;
            histOtherParticle->Fill(xvalue, yvalue);
        }

    }
    
    Int_t count = 0;
    Int_t count2 = 0;
    if(strcmp(htype,"AddedBranch")==0){
        for (Int_t i = 1; i <= hist->GetNbinsX(); ++i) {
            for (Int_t j = 1; j <= hist->GetNbinsY(); ++j) {
                Double_t x = hist->GetXaxis()->GetBinCenter(i);
                Double_t y = hist->GetYaxis()->GetBinCenter(j);

                if (x >= -10 && x <= 10 && y >= 0.6 && y <= 2) {
                    if (hist->GetBinContent(i, j) > 0) {
                        count++;
                    }
                }
            }
        }
        for (Int_t i = 1; i <= histOtherParticle->GetNbinsX(); ++i) {
            for (Int_t j = 1; j <= histOtherParticle->GetNbinsY(); ++j) {
                Double_t x = histOtherParticle->GetXaxis()->GetBinCenter(i);
                Double_t y = histOtherParticle->GetYaxis()->GetBinCenter(j);

                if (x >= -10 && x <= 10 && y >= 0.6 && y <= 2) {
                    if (histOtherParticle->GetBinContent(i, j) > 0) {
                        count2++;
                    }
                }
            }
        }
    cout << "Count1 is " << count << "count2 is " << count2 << endl;
    }
    canvas->cd();
    hist->Draw("COLZ");
    //hist->SetStats(0);
    
    TBox *box = new TBox(-10, 0.6, 10, 2);
    box->SetLineColor(kRed); // Set line color to red
    box->SetFillStyle(0); // Set fill style to transparent
    box->SetLineWidth(2); // Set line width
    box->Draw(); // Draw the box
    
    // Get the canvas name
    const char *canvasName = canvas->GetName();
    // Generate the default file name based on the canvas name
    TString fileName = "./graphjpg/";
    fileName += canvasName;
    fileName += ".jpg";
    // Print the canvas to the default file name
    canvas->SaveAs(fileName);

    canvasOtherParticle->cd();
    histOtherParticle->Draw("COLZ");
    box->Draw();
    //histOtherParticle->SetStats(0);
    histOtherParticle->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle()); 
    histOtherParticle->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());


    const char *canvasOtherParticleName = canvasOtherParticle->GetName();
    TString otherParticleFileName = "./graphjpg/";
    otherParticleFileName += canvasOtherParticleName;
    otherParticleFileName += ".jpg";
    canvasOtherParticle->SaveAs(otherParticleFileName);
}

void Make1DHistEonly(TCanvas *canvas, TH1F *hist, TTree *trees[], initializer_list<const char*> branches, const char* htype, Double_t pval){

    // Initialize histograms for different particle types
    TH1F *histOtherParticle = new TH1F(Form("%s_Pions/Muons", hist->GetName()), Form("%s, Pions/Muons", hist->GetTitle()), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    TH1F *histProton = new TH1F(Form("%s_Protons", hist->GetName()), Form("%s, Protons", hist->GetTitle()), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Calculate total number of trees in the array
    Int_t totalIndex = 0;
    while (trees[totalIndex] != nullptr) {
        totalIndex++;
    }

    // Initialize and set up branches and branch values
    vector<const char*> allBranches;
    Int_t combinedSize = CalcCombinedSize(branches, allBranches, totalIndex);

    Double_t branchVal[combinedSize][20];
    Int_t nPeaks[totalIndex];
    Double_t WWI[totalIndex];
    Int_t nWindowPeak;

    InitializeBranchVal(branchVal, combinedSize);

    SetupRootTreesAddress(trees, totalIndex, combinedSize, branches, allBranches, branchVal, nPeaks, WWI, nWindowPeak);    Long64_t numEntries = trees[0]->GetEntries();
    Double_t xvalue;
    Double_t yvalue;
    for (Long64_t entry = 0; entry < numEntries; ++entry) {
        // Check if all nPeaks values are 1 for the current entry
        bool allPeaksOne = true;
        bool eCut = false;
        bool mpCut = false; //muon_pion cut
        bool WWCut = false;
        for (Int_t treeIdx = 0; treeIdx < totalIndex; ++treeIdx){
            trees[treeIdx]->GetEntry(entry);
            if (treeIdx == 8) {
                continue; // Skip for PbGlass
            }
            if (nPeaks[treeIdx] != 1) {
                allPeaksOne = false;
                break;
            }
        }
        if (allPeaksOne && nWindowPeak==1){
            Double_t sumTOF1 = 0.0, sumTOF0 = 0.0, WWI_avg=0.0; //WWI = WholeWaveformInt
            Double_t PbGlassPEVal;
            PbGlassPEVal = branchVal[totalIndex*1+8][0];
            for (Int_t i=0; i<4; ++i){
                sumTOF1 += branchVal[totalIndex * 0 + (i + 4)][0];
                sumTOF0 += branchVal[totalIndex * 0 + i][0];
                WWI_avg += WWI[i]+WWI[i+4];
            }
            Double_t avgTOF1 = sumTOF1/4;
            Double_t avgTOF0 = sumTOF0/4;
            Double_t finalTOF = avgTOF1-avgTOF0;
            WWI_avg /=8;
            //if (PbGlassPEVal > 10 && finalTOF < 14){ //correct value for windowIntPE
            switch(static_cast<int>(pval * 10)){
                case 10:
                if (PbGlassPEVal > 0.45 && finalTOF <13){ //correct value for IntCharge, kept IntCharge into the PEVal variable for conciseness...
                    eCut = true;
                }
                else if (finalTOF < 13 && PbGlassPEVal <= 0.45) {
                    mpCut = true; // Muonpions
                }
                break;
                case 7:
                if (PbGlassPEVal > 0.25 && finalTOF <13){ //correct value for IntCharge, kept IntCharge into the PEVal variable for conciseness...
                    eCut = true;
                }
                else if (finalTOF < 13 && PbGlassPEVal <= 0.25) {
                    mpCut = true; // Muonpions
                }
                break;
                case 5:
                if (PbGlassPEVal > 0.15 && finalTOF <13){ //correct value for IntCharge, kept IntCharge into the PEVal variable for conciseness...
                    eCut = true;
                }
                else if (finalTOF < 13 && PbGlassPEVal <= 0.15) {
                    mpCut = true; // Muonpions
                }   
                break;           
            }
            if(eCut && WWI_avg<40){
                WWCut = true;
            }
            else if(mpCut && WWI_avg<40){
                WWCut = true;
            }
            else if(!mpCut && !eCut && finalTOF > 13){
                WWCut = true;
            }
        }

            if (eCut && WWCut && (strcmp(htype,"TOFL")==0 || strcmp(htype,"TOFR")==0)){ //used for TOF 1D hist
                Double_t sumTOF0 = 0.0;
                for (Int_t i=0; i<4; ++i){
                    sumTOF0+=branchVal[totalIndex * 0 + i][0];
                }
                xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0/4;
                hist->Fill(xvalue);
            }
            if (WWCut && mpCut && ((strcmp(htype,"TOFL"))==0 || strcmp(htype,"TOFR")==0)){
                Double_t sumTOF0 = 0.0;
                for(Int_t i=0; i<4; ++i){
                    sumTOF0+=branchVal[totalIndex * 0 + i][0];
                }
                xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0/4;
                histOtherParticle->Fill(xvalue);
            }
            if (WWCut && !eCut && !mpCut && ((strcmp(htype,"TOFL"))==0 || strcmp(htype,"TOFR")==0)){
                Double_t sumTOF0 = 0.0;
                for(Int_t i=0; i<4; ++i){
                    sumTOF0+=branchVal[totalIndex * 0 + i][0];
                }
                xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0/4;
                histProton->Fill(xvalue);
            }
        if (eCut && WWCut && (strcmp(htype,"TOFL1")==0 || strcmp(htype,"TOFR1")==0)){ //used for TOF 1D hist
            Double_t sumTOF1 = 0.0;
            for (Int_t i=0; i<4; ++i){
                sumTOF1+=branchVal[totalIndex * 0 + (4 + i)][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF1/4;
            hist->Fill(xvalue);
        }
        if (WWCut && mpCut && ((strcmp(htype,"TOFL1"))==0 || strcmp(htype,"TOFR1")==0)){
            Double_t sumTOF1 = 0.0;
            for(Int_t i=0; i<4; ++i){
                sumTOF1+=branchVal[totalIndex * 0 + (4 + i)][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF1/4;
            histOtherParticle->Fill(xvalue);
        }
        if (WWCut && !eCut && !mpCut && ((strcmp(htype,"TOFL1"))==0 || strcmp(htype,"TOFR1")==0)){
            Double_t sumTOF1 = 0.0;
            for(Int_t i=0; i<4; ++i){
                sumTOF1+=branchVal[totalIndex * 0 + (4 + i)][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF1/4;
            histProton->Fill(xvalue);
        }
    }

    // Find the bin with the highest count in the electron histogram
    Int_t maxBin = hist->GetMaximumBin();

    // Calculate the total integral of the histogram
    Double_t totalIntegral = hist->Integral();

    // Initialize variables to find the range
    Double_t cumulativeIntegral = 0.0;
    Int_t startBin = maxBin, endBin = maxBin;
    Double_t targetIntegral = 0.95 * totalIntegral;

    // Expand the range until the cumulative integral exceeds 80% of the total integral
    while (cumulativeIntegral < targetIntegral) {
        // Expand the range
        startBin = std::max(1, startBin - 1); // Ensure we don't go below the first bin
        endBin = std::min(hist->GetNbinsX(), endBin + 1); // Ensure we don't go above the last bin

        // Calculate the cumulative integral for the expanded range
        cumulativeIntegral = hist->Integral(startBin, endBin);
    }

    // Convert the bin numbers to their corresponding values in seconds
    Double_t startValue = hist->GetXaxis()->GetBinLowEdge(startBin);
    Double_t endValue = hist->GetXaxis()->GetBinUpEdge(endBin);

    // Print the results
    cout << "Start Bin (in seconds): " << startValue << endl;
    cout << "End Bin (in seconds): " << endValue << endl;

    TLegend *legend = new TLegend(0.7,0.2,0.9,0.4);

    TGraph *emptycutGraph = new TGraph(); //used for legend only...
    emptycutGraph->SetLineColor(kRed);
    legend->AddEntry(emptycutGraph, "Electron Window", "l");

    canvas->cd();
    canvas->SetLogy();
    hist->Draw();
    hist->SetLineColor(kBlue);
    hist->SetStats(1);
    histOtherParticle->Draw("same");
    histOtherParticle->SetLineColor(kBlack);
    histOtherParticle->SetStats(0);
    histProton->Draw("same");
    histProton->SetLineColor(kGreen);
    histProton->SetStats(0);

    Double_t lower_bound[4] = {-10,-10,-20,-15};
    Double_t upper_bound[4] = {30,35,15,20};
    std::map<std::string, Int_t> htypeToBoxIdx = {
    {"TOFL", 0},
    {"TOFR", 1},
    {"TOFL1", 2},
    {"TOFR1", 3}
    };
    Int_t boxIdx = htypeToBoxIdx[htype]; 

    //box
    TCutG *cut = new TCutG("myCut",5);
    cut->SetLineColor(kRed);
    cut->SetPoint(0, lower_bound[boxIdx], 0);
	cut->SetPoint(1, upper_bound[boxIdx], 0);
	cut->SetPoint(2, upper_bound[boxIdx], hist->GetMaximum()+10);
	cut->SetPoint(3, lower_bound[boxIdx], hist->GetMaximum()+10);
	cut->SetPoint(4, lower_bound[boxIdx], 0);
    cut->Draw("SAME");
    
    legend->AddEntry(hist, "Electrons", "l");
    legend->AddEntry(histOtherParticle, "muons/pions", "l");
    legend->AddEntry(histProton, "Protons", "l");
    legend->Draw();
    
    
    canvas->cd();
    // Get the canvas name
    const char *canvasName = canvas->GetName();
    // Generate the default file name based on the canvas name
    TString fileName = "./graphjpg/";
    fileName += canvasName;
    fileName += ".jpg";
    // Print the canvas to the default file name
    canvas->SaveAs(fileName);
 
}

void return1DHistRange(TTree *trees[], initializer_list<const char*> branches, const char* TOFtype, Double_t pval, Double_t &startValue, Double_t &endValue){
    // Initialize histograms for different particle types
    TH1F *hist = new TH1F("hist","Run, ACT0L-T0avg Signal Time Corrected distribution with T0,T1,PbGlass and ACT0L nPeak==1 w/ WWI cut;ACT0L-T0(ns);Counts/Bin",210,-60,150);
    TH1F *histOtherParticle = new TH1F(Form("%s_Pions/Muons", hist->GetName()), Form("%s, Pions/Muons", hist->GetTitle()), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    TH1F *histProton = new TH1F(Form("%s_Protons", hist->GetName()), Form("%s, Protons", hist->GetTitle()), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Calculate total number of trees in the array
    Int_t totalIndex = 0;
    while (trees[totalIndex] != nullptr) {
        totalIndex++;
    }

    // Initialize and set up branches and branch values
    vector<const char*> allBranches;
    Int_t combinedSize = CalcCombinedSize(branches, allBranches, totalIndex);

    Double_t branchVal[combinedSize][20];
    Int_t nPeaks[totalIndex];
    Double_t WWI[totalIndex];
    Int_t nWindowPeak;

    InitializeBranchVal(branchVal, combinedSize);

    SetupRootTreesAddress(trees, totalIndex, combinedSize, branches, allBranches, branchVal, nPeaks, WWI, nWindowPeak);    Long64_t numEntries = trees[0]->GetEntries();
    Double_t xvalue;
    Double_t yvalue;
    for (Long64_t entry = 0; entry < numEntries; ++entry) {
        // Check if all nPeaks values are 1 for the current entry
        bool allPeaksOne = true;
        bool eCut = false;
        bool mpCut = false; //muon_pion cut
        bool WWCut = false;
        for (Int_t treeIdx = 0; treeIdx < totalIndex; ++treeIdx){
            trees[treeIdx]->GetEntry(entry);
            if (treeIdx == 8) {
                continue; // Skip for PbGlass
            }
            if (nPeaks[treeIdx] != 1) {
                allPeaksOne = false;
                break;
            }
        }
        if (allPeaksOne && nWindowPeak==1){
            Double_t sumTOF1 = 0.0, sumTOF0 = 0.0, WWI_avg=0.0; //WWI = WholeWaveformInt
            Double_t PbGlassPEVal;
            PbGlassPEVal = branchVal[totalIndex*1+8][0];
            for (Int_t i=0; i<4; ++i){
                sumTOF1 += branchVal[totalIndex * 0 + (i + 4)][0];
                sumTOF0 += branchVal[totalIndex * 0 + i][0];
                WWI_avg += WWI[i]+WWI[i+4];
            }
            Double_t avgTOF1 = sumTOF1/4;
            Double_t avgTOF0 = sumTOF0/4;
            Double_t finalTOF = avgTOF1-avgTOF0;
            WWI_avg /=8;
            //if (PbGlassPEVal > 10 && finalTOF < 14){ //correct value for windowIntPE
            Double_t threshold_pos = pval * (0.5 / 1.1) + (-0.5 * 0.1 / 1.1);
            Double_t threshold_neg = pval * (-0.5 / 1.1) + (-0.5 * 0.1 / 1.1);
            if (pval > 0){
                if(PbGlassPEVal > threshold_pos && finalTOF < 13){
                    eCut = true;
                }
                else if(finalTOF < 13 && PbGlassPEVal <= threshold_pos){
                    mpCut = true;
                }
            }
            if (pval < 0){
                if(PbGlassPEVal > threshold_neg && finalTOF < 13){
                    eCut = true;
                }
                else if(finalTOF < 13 && PbGlassPEVal <= threshold_neg){
                    mpCut = true;
                }
            }
            if(eCut && WWI_avg<40){
                WWCut = true;
            }
            else if(mpCut && WWI_avg<40){
                WWCut = true;
            }
            else if(!mpCut && !eCut && finalTOF > 13){
                WWCut = true;
            }
        }
            if (eCut && WWCut && (strcmp(TOFtype,"TOFL")==0 || strcmp(TOFtype,"TOFR")==0)){ //used for TOF 1D hist
                Double_t sumTOF0 = 0.0;
                for (Int_t i=0; i<4; ++i){
                    sumTOF0+=branchVal[totalIndex * 0 + i][0];
                }
                xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0/4;
                hist->Fill(xvalue);
            }
            if (WWCut && mpCut && ((strcmp(TOFtype,"TOFL"))==0 || strcmp(TOFtype,"TOFR")==0)){
                Double_t sumTOF0 = 0.0;
                for(Int_t i=0; i<4; ++i){
                    sumTOF0+=branchVal[totalIndex * 0 + i][0];
                }
                xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0/4;
                histOtherParticle->Fill(xvalue);
            }
            if (WWCut && !eCut && !mpCut && ((strcmp(TOFtype,"TOFL"))==0 || strcmp(TOFtype,"TOFR")==0)){
                Double_t sumTOF0 = 0.0;
                for(Int_t i=0; i<4; ++i){
                    sumTOF0+=branchVal[totalIndex * 0 + i][0];
                }
                xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF0/4;
                histProton->Fill(xvalue);
            }
        if (eCut && WWCut && (strcmp(TOFtype,"TOFL1")==0 || strcmp(TOFtype,"TOFR1")==0)){ //used for TOF 1D hist
            Double_t sumTOF1 = 0.0;
            for (Int_t i=0; i<4; ++i){
                sumTOF1+=branchVal[totalIndex * 0 + (4 + i)][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF1/4;
            hist->Fill(xvalue);
        }
        if (WWCut && mpCut && ((strcmp(TOFtype,"TOFL1"))==0 || strcmp(TOFtype,"TOFR1")==0)){
            Double_t sumTOF1 = 0.0;
            for(Int_t i=0; i<4; ++i){
                sumTOF1+=branchVal[totalIndex * 0 + (4 + i)][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF1/4;
            histOtherParticle->Fill(xvalue);
        }
        if (WWCut && !eCut && !mpCut && ((strcmp(TOFtype,"TOFL1"))==0 || strcmp(TOFtype,"TOFR1")==0)){
            Double_t sumTOF1 = 0.0;
            for(Int_t i=0; i<4; ++i){
                sumTOF1+=branchVal[totalIndex * 0 + (4 + i)][0];
            }
            xvalue = branchVal[totalIndex * 0 + 9][0] - sumTOF1/4;
            histProton->Fill(xvalue);
        }
    }

    // Find the bin with the highest count in the electron histogram
    Int_t maxBin = hist->GetMaximumBin();

    // Calculate the total integral of the histogram
    Double_t totalIntegral = hist->Integral();

    // Initialize variables to find the range
    Double_t cumulativeIntegral = 0.0;
    Int_t startBin = maxBin, endBin = maxBin;
    Double_t targetIntegral = 0.95 * totalIntegral;

    // Expand the range until the cumulative integral exceeds 80% of the total integral
    while (cumulativeIntegral < targetIntegral) {
        // Expand the range
        startBin = std::max(1, startBin - 1); // Ensure we don't go below the first bin
        endBin = std::min(hist->GetNbinsX(), endBin + 1); // Ensure we don't go above the last bin

        // Calculate the cumulative integral for the expanded range
        cumulativeIntegral = hist->Integral(startBin, endBin);
    }

    // Convert the bin numbers to their corresponding values in seconds
    startValue = hist->GetXaxis()->GetBinLowEdge(startBin);
    endValue = hist->GetXaxis()->GetBinUpEdge(endBin);
    startValue *=1.0;
    endValue *=1.0;

    // Print the results
    //cout << "Start Bin (in seconds): " << startValue << endl;
    //cout << "End Bin (in seconds): " << endValue << endl;
    delete hist;
    delete histOtherParticle;
    delete histProton;
}

TGraph* Make2DHist0vs1(TCanvas *canvas, TGraph *graph, TTree *trees[], initializer_list<const char*> branches, const char* htype, Double_t lim[], Double_t pval, Bool_t WindowON, Bool_t ACT01 = true){

    TGraph *tempgraph = new TGraph();
    graph->SetMarkerColor(kRed);
    tempgraph->SetMarkerColor(kBlue);  
    graph->SetMarkerStyle(kFullSquare);
    tempgraph->SetMarkerStyle(kFullSquare);
    graph->SetMarkerSize(0.3);
    tempgraph->SetMarkerSize(0.3);
    TGraph *cutgraph = new TGraph();
    cutgraph->SetLineColor(kBlack);
    cutgraph->SetLineWidth(2);

    Int_t totalIndex = CalcTotalIndex(trees);
   
    vector<const char*> allBranches;

    Int_t combinedSize = CalcCombinedSize(branches, allBranches, totalIndex);

    Double_t branchVal[combinedSize][20];
    Int_t nPeaks[totalIndex];
    Double_t WWI[totalIndex];
    Int_t nWindowPeak;

    // Initialize the branch values
    InitializeBranchVal(branchVal, combinedSize);

    SetupRootTreesAddress(trees, totalIndex, combinedSize, branches, allBranches, branchVal, nPeaks, WWI, nWindowPeak);

    Double_t xvalue;
    Double_t yvalue;


    //Set up cuts plot and slopes/intercept.
    Double_t x1 = 0, x2, y1, y2 = 0;
    Double_t m;
    Double_t b;
    //Same Cut
    
    if (strcmp(htype,"IntCharge")==0){
        m = -0.082551;
        y1= 0.0308163;
        b = y1;
        x2 = -y1/m;
        cutgraph->SetPoint(0,x1,y1);
        cutgraph->SetPoint(1,x2,y2);
    }
    if (strcmp(htype,"wIntCharge")==0){
        //m = -0.1233265;
        m = -0.102939;
        //y1= 0.0446939;
        y1 = 0.0308163;
        b = y1;
        x2 = -y1/m;
        cutgraph->SetPoint(0,x1,y1);
        cutgraph->SetPoint(1,x2,y2);
    }
    /*
    //UNIQUE CUTS
    if (strcmp(htype,"IntCharge")==0 && pval==0.5){
        if(WindowON){
            m = -0.0621633;
            y1 = 0.0169388;
        }
        else{
            m=-0.0621633;
            y1=0.0169388;
        }
        b = y1;
        x2 = -y1/m;
        cutgraph->SetPoint(0,x1,y1);
        cutgraph->SetPoint(1,x2,y2);
    }
    else if (strcmp(htype,"IntCharge")==0 && pval==0.7){
        if(WindowON){
            m = -0.082551;
            y1= 0.0308163;
        }
        else{
            m = -0.102939;
            y1= 0.0308163;
        }
        b = y1;
        x2 = -y1/m;
        cutgraph->SetPoint(0,x1,y1);
        cutgraph->SetPoint(1,x2,y2);
    }
    else if (strcmp(htype,"IntCharge")==0 && pval==1){
        if(WindowON){
            m = -0.0621633;
            y1 = 0.0377551;
        }
        else{
            m = -0.0621633;
            y1 = 0.0377551;        
        }
        b = y1;
        x2 = -y1/m;
        cutgraph->SetPoint(0,x1,y1);
        cutgraph->SetPoint(1,x2,y2);
    }
    else if (strcmp(htype,"wIntCharge")==0){
        if(pval==0.5){
            m = -0.082551;
            y1 = 0.0308163;
        }
        if(pval==0.7){
            m = -0.1233265;
            y1= 0.0446939;
        }
        if(pval==1){
            m = -0.0417755;
            y1 = 0.0446939;
        }
        b = y1;
        x2 = -y1/m;
        cutgraph->SetPoint(0,x1,y1);
        cutgraph->SetPoint(1,x2,y2);
    }
    */
    // if(strcmp(htype,"PVT")==0){
    //     x2=1.85;
    //     y1=0.25;
    //     b = y1;
    //     m = (y1-y2)/(x1-x2);
    //     cutgraph->SetPoint(0,x1,y1);
    //     cutgraph->SetPoint(1,x2,y2);
    // }
    Double_t lower_bound[4];
    Double_t upper_bound[4];

    Double_t startVal = 0, endVal = 0;
    if(ACT01 && strcmp(htype,"IntCharge")==0){
        TTree *histrangetree[11];
        for(Int_t j=0; j<4; ++j){
            for(Int_t i=0; i<9; i++){
                histrangetree[i] = trees[i];
            }
        histrangetree[9] = trees[9+j];
        histrangetree[10] = nullptr;
        return1DHistRange(histrangetree, {}, "TOFL", pval, startVal, endVal);
        lower_bound[j] = startVal;
        upper_bound[j] = endVal;
        cout << "j = " << j << " lower_bound is " << lower_bound[j] << " upperbound is " << upper_bound[j] << endl;
        }
    }

    Long64_t numEntries = trees[0]->GetEntries();
    Int_t Ncounter = 0, Ncounter_OP = 0;
    Int_t cutCounter =0, cutCounter_OP =0; //Used to calculate purity/efficiency.
    Int_t protonCount = 0;
    for (Long64_t entry = 0; entry < numEntries; ++entry) {
        Double_t sumTOF1 = 0.0, sumTOF0 = 0.0;
        Double_t avgTOF1 = 0.0, avgTOF0 = 0.0, finalTOF = 0.0;
        if (ACT01) {
            /*
            // In order, the lower bound (of the window) are 0L, 0R, 1L, 1R and are in ns.
            Double_t lb[4] = {-10, -10, -20, -15};
            Double_t ub[4] = {30, 35, 15, 20};
            std::copy(std::begin(lb), std::end(lb), lower_bound);
            std::copy(std::begin(ub), std::end(ub), upper_bound);
            */
        } 
        else{
            Double_t lb[4] = {-15, -15, -15, -30};
            Double_t ub[4] = {12, 20, 11, 40};
            std::copy(std::begin(lb), std::end(lb), lower_bound);
            std::copy(std::begin(ub), std::end(ub), upper_bound);
        }
        //Double_t ACT0ST = 0.0, ACT1ST = 0.0;
        Double_t ACT0LST[20]={0.0}, ACT0RST[20]={0.0}, ACT1LST[20]={0.0}, ACT1RST[20]={0.0};
        // Check if all nPeaks values are 1 only for the TOF and PbGlass. (first 9 trees)
        Double_t WWI_avg = 0.0;
        bool allPeaksOne = true;
        bool eCut = false;
        bool mpCut = false; //not implemented in current ACT0vsACT1
        bool WWCut = false;
        for (Int_t treeIdx = 0; treeIdx < 9; ++treeIdx){
            trees[treeIdx]->GetEntry(entry);
            if (nPeaks[treeIdx] != 1) {
                allPeaksOne = false;
                break;
            }
        }
        if (allPeaksOne && nWindowPeak==1){
            for (Int_t treeIdx = 0; treeIdx < totalIndex; ++treeIdx){
                trees[treeIdx]->GetEntry(entry);
            }
            Double_t PbGlassPEVal;
            PbGlassPEVal = branchVal[totalIndex*1+8][0];
            for (Int_t i=0; i<4; ++i){
                sumTOF1 += branchVal[totalIndex * 0 + (i + 4)][0];
                sumTOF0 += branchVal[totalIndex * 0 + i][0];
                WWI_avg += WWI[i]+WWI[i+4];
            }
            avgTOF1 = sumTOF1/4;
            avgTOF0 = sumTOF0/4;
            finalTOF = avgTOF1-avgTOF0;
            WWI_avg /=8;
            //if (PbGlassPEVal > 10 && finalTOF < 14){ //correct value for windowIntPE
            /*
            switch(static_cast<int>(pval * 10)){
                case 10:
                if (PbGlassPEVal > 0.45 && finalTOF <13){ //correct value for IntCharge, kept IntCharge into the PEVal variable for conciseness...
                    eCut = true;
                }
                else if (finalTOF < 13 && PbGlassPEVal <= 0.45) {
                    mpCut = true; // Muonpions
                }
                break;
                case 7:
                if (PbGlassPEVal > 0.25 && finalTOF <13){ //correct value for IntCharge, kept IntCharge into the PEVal variable for conciseness...
                    eCut = true;
                }
                else if (finalTOF < 13 && PbGlassPEVal <= 0.25) {
                    mpCut = true; // Muonpions
                }
                break;
                case 5:
                if (PbGlassPEVal > 0.15 && finalTOF <13){ //correct value for IntCharge, kept IntCharge into the PEVal variable for conciseness...
                    eCut = true;
                }
                else if (finalTOF < 13 && PbGlassPEVal <= 0.15) {
                    mpCut = true; // Muonpions
                }   
                break;           
            }
            */
            Double_t threshold_pos = pval * (0.5 / 1.1) + (-0.5 * 0.1 / 1.1);
            Double_t threshold_neg = pval * (-0.5 / 1.1) + (-0.5 * 0.1 / 1.1);
           if (pval > 0){
                if(PbGlassPEVal > threshold_pos && finalTOF < 13){
                    eCut = true;
                }
                else if(finalTOF < 13 && PbGlassPEVal <= threshold_pos){
                    mpCut = true;
                }
            }
            if (pval < 0){
                if(PbGlassPEVal > threshold_neg && finalTOF < 13){
                    eCut = true;
                }
                else if(finalTOF < 13 && PbGlassPEVal <= threshold_neg){
                    mpCut = true;
                }
            }
            if(eCut && WWI_avg<40){
                WWCut = true;
            }
            else if(mpCut && WWI_avg<40){ //muons pions
                WWCut = true;
            }
            else if(!eCut && !mpCut && finalTOF > 13){ //protons might need change for 1200MeV and 200MeV
                WWCut = true;
                protonCount++;
            }
            //find integration window. currently still using signal time and not the corrected one.
            //the bound are defined as SignalTime of TOF10 + (offset from digitizer - TOF1(avg)) + window bound
            //Only PbGlass and TOF are nPeak = 1, potential problem is having >1 peak and having one peak earlier than the expected one in ACTs
            for(Int_t j=0; j<4; j++){
                //lower_bound[j] = branchVal[0][5][0] + (branchVal[0][9+j][0]-avgTOF1) - 16;
                //upper_bound[j] = branchVal[0][5][0] + (branchVal[0][9+j][0]-avgTOF1) + 45;
                //lower_bound[j] = -50;
                //upper_bound[j] = 0;
            }
        }

        if (eCut && WWCut && strcmp(htype,"wIntCharge")==0){ 
            yvalue = branchVal[totalIndex*1+11][0] + branchVal[totalIndex*1+12][0];
            xvalue = branchVal[totalIndex*1+9][0] + branchVal[totalIndex*1+10][0];
            graph->SetPoint(Ncounter, xvalue, yvalue);
            Ncounter++;
            if(isPointUnderLine(m,b,xvalue,yvalue)){
                cutCounter++;
            }
        }
        else if (WWCut && !eCut && strcmp(htype,"wIntCharge")==0) {
            yvalue = branchVal[totalIndex * 1 + 11][0] + branchVal[totalIndex * 1 + 12][0];
            xvalue = branchVal[totalIndex * 1 + 9][0] + branchVal[totalIndex * 1 + 10][0];
            tempgraph->SetPoint(Ncounter_OP, xvalue, yvalue);
            Ncounter_OP++;
            if(isPointUnderLine(m,b,xvalue,yvalue)){
                cutCounter_OP++;
            }
        }
        else if (eCut && WWCut && WindowON){

            xvalue = 0;
            yvalue = 0;
            Int_t bVIdx; 
            if(strcmp(htype,"IntCharge")==0){
                bVIdx = 2;
            }
            else if(strcmp(htype,"PVT")==0){
                bVIdx = 3;
            }
            for(Int_t j=0; j<20; ++j){
                ACT0LST[j] = branchVal[totalIndex * 0 + 9][j] - avgTOF0; 
                ACT0RST[j] = branchVal[totalIndex * 0 + 10][j] - avgTOF0;
                if(ACT0LST[j]>lower_bound[0] && ACT0RST[j]>lower_bound[1] && ACT0LST[j]<upper_bound[0] && ACT0RST[j]<upper_bound[1]){
                    xvalue = branchVal[totalIndex*bVIdx+9][j] + branchVal[totalIndex*bVIdx+10][j];
                    break;
                }
                else if(ACT0LST[j]>lower_bound[0] && ACT0LST[j]<upper_bound[0]){
                    xvalue = branchVal[totalIndex*bVIdx+9][j];
                    break;
                } 
                else if(ACT0RST[j]>lower_bound[1] && ACT0RST[j]<upper_bound[1]){
                    xvalue = branchVal[totalIndex*bVIdx+10][j];
                    break;
                } 
            }
            for(Int_t j=0; j<20; ++j){
                ACT1LST[j] = branchVal[totalIndex * 0 + 11][j] - avgTOF0;
                ACT1RST[j] = branchVal[totalIndex * 0 + 12][j] - avgTOF0;
                if(ACT1LST[j]>lower_bound[2] && ACT1RST[j]>lower_bound[3] && ACT1LST[j]<upper_bound[2] && ACT1RST[j]<upper_bound[3]){
                    yvalue = branchVal[totalIndex * bVIdx + 11][j] + branchVal[totalIndex * bVIdx + 12][j];
                    break;
                }
                else if(ACT1LST[j]>lower_bound[2] && ACT1LST[j]<upper_bound[2]){
                    yvalue = branchVal[totalIndex * bVIdx + 11][j];
                    break;
                } 
                else if(ACT1RST[j]>lower_bound[3] && ACT1RST[j]<upper_bound[3]){
                    yvalue = branchVal[totalIndex * bVIdx + 12][j];
                    break;
                } 
            }
            /*
            if(xvalue==0 && Ncounter < 500){
                cout << "xvalue = 0 for IntCharge at event e only " << entry << endl;
            }
            if (yvalue==0 && Ncounter < 100){
                cout << "yvalue = 0 for IntCharge at event e only " << entry << endl;
            } 
            */
            graph->SetPoint(Ncounter, xvalue, yvalue);
            Ncounter++;
            if(isPointUnderLine(m,b,xvalue,yvalue)){
                cutCounter++;
            }
        }
        else if (WWCut && !eCut && WindowON) {
         
            xvalue = 0;
            yvalue = 0;
            Int_t bVIdx; 
            if(strcmp(htype,"IntCharge")==0){
                bVIdx = 2;
            }
            else if(strcmp(htype,"PVT")==0){
                bVIdx = 3;
            }

            for(Int_t j=0; j<20; ++j){
                ACT0LST[j] = branchVal[totalIndex*0+9][j] - avgTOF0; 
                ACT0RST[j] = branchVal[totalIndex*0+10][j] - avgTOF0;
                if(ACT0LST[j]>lower_bound[0] && ACT0RST[j]>lower_bound[1] && ACT0LST[j]<upper_bound[0] && ACT0RST[j]<upper_bound[1]){
                    xvalue = branchVal[totalIndex*bVIdx+9][j] + branchVal[totalIndex*bVIdx+10][j];
                    break;
                }
                else if(ACT0LST[j]>lower_bound[0] && ACT0LST[j]<upper_bound[0]){
                    xvalue = branchVal[totalIndex*bVIdx+9][j];
                    break;
                } 
                else if(ACT0RST[j]>lower_bound[1] && ACT0RST[j]<upper_bound[1]){
                    xvalue = branchVal[totalIndex*bVIdx+10][j];
                    break;
                } 
            }
            for(Int_t j=0; j<20; ++j){
                ACT1LST[j] = branchVal[totalIndex*0+11][j] - avgTOF0;
                ACT1RST[j] = branchVal[totalIndex*0+12][j] - avgTOF0;
                if(ACT1LST[j]>lower_bound[2] && ACT1RST[j]>lower_bound[3] && ACT1LST[j]<upper_bound[2] && ACT1RST[j]<upper_bound[3]){
                    yvalue = branchVal[totalIndex*bVIdx+11][j] + branchVal[totalIndex*bVIdx+12][j];
                    break;
                }
                else if(ACT1LST[j]>lower_bound[2] && ACT1LST[j]<upper_bound[2]){
                    yvalue = branchVal[totalIndex*bVIdx+11][j];
                    break;
                } 
                else if(ACT1RST[j]>lower_bound[3] && ACT1RST[j]<upper_bound[3]){
                    yvalue = branchVal[totalIndex*bVIdx+12][j];
                    break;
                } 
            } 
            tempgraph->SetPoint(Ncounter_OP, xvalue, yvalue);
            Ncounter_OP++;
            if(isPointUnderLine(m,b,xvalue,yvalue)){
                cutCounter_OP++;
            }
        }
        //WINDOW OFF
        else if (eCut && WWCut && !WindowON){
            Int_t bVIdx; 
                if(strcmp(htype,"IntCharge")==0){
                    bVIdx = 2;
                }
                else if(strcmp(htype,"PVT")==0){
                    bVIdx = 3;
                }
            xvalue = branchVal[totalIndex*bVIdx+9][0] + branchVal[totalIndex*bVIdx+10][0];
            yvalue = branchVal[totalIndex*bVIdx+11][0] + branchVal[totalIndex*bVIdx+12][0];
            graph->SetPoint(Ncounter, xvalue, yvalue);
            Ncounter++;
            if(isPointUnderLine(m,b,xvalue,yvalue)){
                cutCounter++;
            }
        }
        else if (WWCut && !eCut && !WindowON) {
            Int_t bVIdx; 
                if(strcmp(htype,"IntCharge")==0){
                    bVIdx = 2;
                }
                else if(strcmp(htype,"PVT")==0){
                    bVIdx = 3;
                }
            xvalue = branchVal[totalIndex*bVIdx+9][0] + branchVal[totalIndex*bVIdx+10][0];
            yvalue = branchVal[totalIndex*bVIdx+11][0] + branchVal[totalIndex*bVIdx+12][0];
            tempgraph->SetPoint(Ncounter_OP, xvalue, yvalue);
            Ncounter_OP++;
            if(isPointUnderLine(m,b,xvalue,yvalue)){
                cutCounter_OP++;
            }
        }
        InitializeBranchVal(branchVal, combinedSize);
    }

    Double_t xlim1 = lim[0];
    Double_t xlim2 = lim[1];
    Double_t ylim1 = lim[2];
    Double_t ylim2 = lim[3];

    cout << "There is in total " << Ncounter+Ncounter_OP << " points with " << Ncounter << " electrons and " << Ncounter_OP << " non-e" << endl;
    cout << "Additionally, there are " << protonCount << " protons" << endl;

    TLegend *leg = new TLegend(0.8,0.8,0.95,0.95);
    TGraph *electronGraph = new TGraph();
    electronGraph->SetMarkerColor(kRed);
    electronGraph->SetMarkerStyle(20); // Circle marker
    leg->AddEntry(electronGraph, "Electrons", "p"); // "p" for point marker
    TGraph *otherParticleGraph = new TGraph();
    otherParticleGraph->SetMarkerColor(kBlue);
    otherParticleGraph->SetMarkerStyle(20); // Circle marker
    leg->AddEntry(otherParticleGraph, "Other Particles", "p"); // "p" for point marker

    TMultiGraph *multiGraph = new TMultiGraph();
    multiGraph->Add(electronGraph);
    multiGraph->Add(otherParticleGraph);
    multiGraph->Add(graph);
    multiGraph->Add(tempgraph);
    multiGraph->SetTitle(graph->GetTitle());
    multiGraph->GetXaxis()->SetTitle(graph->GetXaxis()->GetTitle());
    multiGraph->GetYaxis()->SetTitle(graph->GetYaxis()->GetTitle());
    multiGraph->GetXaxis()->SetRangeUser(xlim1,xlim2);
    multiGraph->GetYaxis()->SetRangeUser(ylim1,ylim2);

    canvas->cd();
    multiGraph->Draw("AP");
    //graph->Draw("AP");
    //tempgraph->Draw("P SAME");
    cutgraph->Draw("L SAME");
    //leg->AddEntry(graph, "Electrons", "p");
    //leg->AddEntry(tempgraph, "Other Particles", "p");
    leg->Draw();

    Double_t e_Eff, e_Pur, nonE_Eff, nonE_Pur;
    e_Pur = (Ncounter-cutCounter)/(Ncounter-cutCounter+Ncounter_OP-static_cast<Double_t>(cutCounter_OP));
    e_Eff = (Ncounter-cutCounter)/(static_cast<Double_t>(Ncounter));
    nonE_Pur = cutCounter_OP/(static_cast<Double_t>(cutCounter_OP)+cutCounter);
    nonE_Eff = cutCounter_OP/static_cast<Double_t>(Ncounter_OP);
    
    if (strcmp(htype,"IntCharge")==0 && WindowON){
        cout << "For IntCharge at momentum " << pval << "GeV, the Electron Signal Purity is " << e_Pur << " and the Electron Signal Efficiency is " << e_Eff << " The non-e signal purity is " << nonE_Pur << " and the efficiency is "<< nonE_Eff << endl;
        cout << "In order of lower and top error it is " <<  e_Pur - TEfficiency::ClopperPearson(Ncounter-cutCounter+Ncounter_OP-cutCounter_OP,Ncounter-cutCounter,0.683, false) << " and " << TEfficiency::ClopperPearson(Ncounter-cutCounter+Ncounter_OP-cutCounter_OP,Ncounter-cutCounter,0.683, true) - e_Pur;
        cout << "; " <<  e_Eff - TEfficiency::ClopperPearson(Ncounter,Ncounter-cutCounter,0.683, false) << " and " << TEfficiency::ClopperPearson(Ncounter,Ncounter-cutCounter,0.683, true) - e_Eff;
        cout << "; " <<  nonE_Pur - TEfficiency::ClopperPearson(cutCounter_OP+cutCounter,cutCounter_OP,0.683, false) << " and " << TEfficiency::ClopperPearson(cutCounter_OP+cutCounter,cutCounter_OP,0.683, true) - nonE_Pur;
        cout << "; " <<  nonE_Eff - TEfficiency::ClopperPearson(Ncounter_OP,cutCounter_OP,0.683, false) << " and " << TEfficiency::ClopperPearson(Ncounter_OP,cutCounter_OP,0.683, true) - nonE_Eff << endl;
    }
    if (strcmp(htype,"IntCharge")==0 && !WindowON){
        cout << "For IntCharge at momentum " << pval << "GeV with no Window, the Electron Signal Purity is " << e_Pur << " and the Electron Signal Efficiency is " << e_Eff << " The non-e signal purity is " << nonE_Pur << " and the efficiency is "<< nonE_Eff << endl;
        cout << "In order of lower and top error it is " <<  e_Pur - TEfficiency::ClopperPearson(Ncounter-cutCounter+Ncounter_OP-cutCounter_OP,Ncounter-cutCounter,0.683, false) << " and " << TEfficiency::ClopperPearson(Ncounter-cutCounter+Ncounter_OP-cutCounter_OP,Ncounter-cutCounter,0.683, true) - e_Pur;
        cout << "; " <<  e_Eff - TEfficiency::ClopperPearson(Ncounter,Ncounter-cutCounter,0.683, false) << " and " << TEfficiency::ClopperPearson(Ncounter,Ncounter-cutCounter,0.683, true) - e_Eff;
        cout << "; " <<  nonE_Pur - TEfficiency::ClopperPearson(cutCounter_OP+cutCounter,cutCounter_OP,0.683, false) << " and " << TEfficiency::ClopperPearson(cutCounter_OP+cutCounter,cutCounter_OP,0.683, true) - nonE_Pur;
        cout << "; " <<  nonE_Eff - TEfficiency::ClopperPearson(Ncounter_OP,cutCounter_OP,0.683, false) << " and " << TEfficiency::ClopperPearson(Ncounter_OP,cutCounter_OP,0.683, true) - nonE_Eff << endl;
    }
    if (strcmp(htype,"wIntCharge")==0){
        cout << "For wIntCharge at momentum " << pval << "GeV, the Electron Signal Purity is " << e_Pur << " and the Electron Signal Efficiency is " << e_Eff << " The non-e signal purity is " << nonE_Pur << " and the efficiency is "<< nonE_Eff << endl;
        cout << "In order of lower and top error it is " <<  e_Pur - TEfficiency::ClopperPearson(Ncounter-cutCounter+Ncounter_OP-cutCounter_OP,Ncounter-cutCounter,0.683, false) << " and " << TEfficiency::ClopperPearson(Ncounter-cutCounter+Ncounter_OP-cutCounter_OP,Ncounter-cutCounter,0.683, true) - e_Pur;
        cout << "; " <<  e_Eff - TEfficiency::ClopperPearson(Ncounter,Ncounter-cutCounter,0.683, false) << " and " << TEfficiency::ClopperPearson(Ncounter,Ncounter-cutCounter,0.683, true) - e_Eff;
        cout << "; " <<  nonE_Pur - TEfficiency::ClopperPearson(cutCounter_OP+cutCounter,cutCounter_OP,0.683, false) << " and " << TEfficiency::ClopperPearson(cutCounter_OP+cutCounter,cutCounter_OP,0.683, true) - nonE_Pur;
        cout << "; " <<  nonE_Eff - TEfficiency::ClopperPearson(Ncounter_OP,cutCounter_OP,0.683, false) << " and " << TEfficiency::ClopperPearson(Ncounter_OP,cutCounter_OP,0.683, true) - nonE_Eff << endl;
    }
    if(strcmp(htype,"PVT")==0){
        cout << "For MaxVoltage, the Signal Purity is " << cutCounter_OP/(cutCounter_OP+static_cast<Double_t>(cutCounter)) << " and the Signal Efficiency is " << cutCounter_OP/(static_cast<Double_t>(Ncounter_OP)) << "and the rejection rate is " << 1-cutCounter/(static_cast<Double_t>(Ncounter)) << endl;
        cout << "Also there are " << Ncounter << " TGraphPoint and " << Ncounter_OP << " TGraphOP points" << endl;
    }

    // Get the canvas name
    const char *canvasName = canvas->GetName();
    // Generate the default file name based on the canvas name
    TString fileName = "./graphjpg/";
    fileName += canvasName;
    fileName += ".jpg";
    // Print the canvas to the default file name
    canvas->SaveAs(fileName);

    return tempgraph;

}


// Define a map for run numbers based on pVal
std::map<Double_t, TString> runNumbers = {
    {1.0, "452"},
    {0.7, "457"},
    {0.5, "435"},
    {0.3, "475"},
    {-0.3, "353"}
};

std::pair<TGraph*, TGraph*> Create2DPlot(Double_t pVal, const char* htype, TTree* trees[], std::initializer_list<const char*> branches, Double_t limit[], Bool_t WindowOn, Bool_t ACT01 = true) {
    TString runNumber = runNumbers[pVal]; // Get the run number based on pVal

    // Determine the graph name based on htype, pVal, WindowOn, and ACT01
    TString graphName;
    if (ACT01) {
        if (strcmp(htype, "wIntCharge") == 0) {
            graphName = Form("Windo%s_0vs1_run%s", htype, runNumber.Data());
        } else {
            graphName = Form("%s_0vs1_run%s", htype, runNumber.Data());
        }
    } else {
        if (strcmp(htype, "wIntCharge") == 0) {
            graphName = Form("Windo%s_2vs3_run%s", htype, runNumber.Data());
        } else {
            graphName = Form("%s_2vs3_run%s", htype, runNumber.Data());
        }
    }

    // Determine the graph title based on htype, pVal, WindowOn, and ACT01
    TString graphTitle;
    TString xAxisTitle;
    TString yAxisTitle;
    if (ACT01) {
        if (WindowOn) {
            if (strcmp(htype, "wIntCharge") == 0) {
                graphTitle = Form("WindowIntegratedCharge ntuple, run%s p=%.1f GeV", runNumber.Data(), pVal);
            } else {
                graphTitle = Form("%s of Integrated Charge inside the fixed window plot of ACT0 vs ACT1, run%s p=%.1f GeV", htype, runNumber.Data(), pVal);
            }
        } else {
            if (strcmp(htype, "wIntCharge") == 0) {
                graphTitle = Form("IntegratedCharge ntuple, run%s p=%.1f GeV", runNumber.Data(), pVal);
            } else {
                graphTitle = Form("%s of Integrated Charge without the fixed window plot of ACT0 vs ACT1, run%s p=%.1f GeV", htype, runNumber.Data(), pVal);
            }
        }
        xAxisTitle = Form("ACT0 %s (a.u.)", htype);
        yAxisTitle = Form("ACT1 %s (a.u.)", htype);
    } else {
        if (WindowOn) {
            if (strcmp(htype, "wIntCharge") == 0) {
                graphTitle = Form("WindowIntegratedCharge ntuple, run%s p=%.1f GeV", runNumber.Data(), pVal);
            } else {
                graphTitle = Form("%s of Integrated Charge inside the fixed window plot of ACT2 vs ACT3, run%s p=%.1f GeV", htype, runNumber.Data(), pVal);
            }
        } else {
            if (strcmp(htype, "wIntCharge") == 0) {
                graphTitle = Form("IntegratedCharge ntuple, run%s p=%.1f GeV", runNumber.Data(), pVal);
            } else {
                graphTitle = Form("%s of Integrated Charge without the fixed window plot of ACT2 vs ACT3, run%s p=%.1f GeV", htype, runNumber.Data(), pVal);
            }
        }
        xAxisTitle = Form("ACT2 %s (a.u.)", htype);
        yAxisTitle = Form("ACT3 %s (a.u.)", htype);
    }

    TCanvas* canvas = new TCanvas(Form("upd_%s", graphName.Data()), "title", 800, 600);
    TGraph* graph = new TGraph();
    graph->SetTitle(graphTitle.Data());
    graph->GetXaxis()->SetTitle(xAxisTitle.Data());
    graph->GetYaxis()->SetTitle(yAxisTitle.Data());
    
    // Call Make2DHist0vs1 to fill the graph and graph_OP
    TGraph* graph_OP = Make2DHist0vs1(canvas, graph, trees, branches, htype, limit, pVal, WindowOn, ACT01);

    return std::make_pair(graph, graph_OP);
}

void Create1DPlot(const char* htype, TTree* trees[], initializer_list<const char*> branches, Double_t pVal, Bool_t isSecondtACT = false, Bool_t isACT2 = false, Int_t bin=210, Int_t ini=-60, Int_t final= 150) {
    //Play with htype and boolean. The hytpe determine T0 or T1 and whether its L or R.
    //If no boolean, by default it goes to ACT0L/ACT0R
    //is ACT2 determine whether you switch to ACT2-3 instead of staying in ACT0-1, false by default.
    //isSecondtACT determine whether u take the second one of the pair, if false(default) take the first of the pair so ACT0 or ACT2
    TString runNumber = runNumbers[pVal]; // Get the run number based on pVal
    const char* name1;
    const char* name2;

    if (strcmp(htype, "TOFL") == 0) {
        name1 = "T0";
        if (isACT2) {
            name2 = (isSecondtACT) ? "ACT3L" : "ACT2L";
        } else {
            name2 = (isSecondtACT) ? "ACT1L" : "ACT0L";
        }
    } else if (strcmp(htype, "TOFL1") == 0) {
        name1 = "T1";
        if (isACT2) {
            name2 = (isSecondtACT) ? "ACT3L" : "ACT2L";
        } else {
            name2 = (isSecondtACT) ? "ACT1L" : "ACT0L";
        }
    } else if (strcmp(htype, "TOFR") == 0) {
        name1 = "T0";
        if (isACT2) {
            name2 = (isSecondtACT) ? "ACT3R" : "ACT2R";
        } else {
            name2 = (isSecondtACT) ? "ACT1R" : "ACT0R";
        }
    } else if (strcmp(htype, "TOFR1") == 0) {
        name1 = "T1";
        if (isACT2) {
            name2 = (isSecondtACT) ? "ACT3R" : "ACT2R";
        } else {
            name2 = (isSecondtACT) ? "ACT1R" : "ACT0R";
        }
    }

    TString histName, histTitle;
    histName = Form("h_%s_%s_r%s",name2,name1,runNumber.Data());
    histTitle = Form("Run %s p = %.1fGeV, %s-%savg Signal Time Corrected distribution with T0,T1,PbGlass and %s, nPeak==1 w/ WWICut;%s-%s(ns);Counts/Bin",runNumber.Data(),pVal,name2,name1,name2,name2,name1);
    TCanvas* canvas = new TCanvas(Form("upd_H1D_STCorr_%s_%s_run%s",name2,name1,runNumber.Data()), "title", 800, 600);
    TH1F* hist = new TH1F(histName,histTitle,bin,ini,final);
    Make1DHistEonly(canvas,hist,trees,branches,htype,pVal);

}

TGraph* MakeSimple2DGraph0vs1NoDraw(TGraph *graph, TTree *trees[], initializer_list<const char*> branches, Double_t momentum, const char* htype = "wIntCharge"){
    //calculate boundaries in IntCharge.
    Double_t startVal = 0, endVal = 0;
    Double_t lower_bound[4] = {0.0}, upper_bound[4] = {0.0};

    if (strcmp(htype, "IntCharge") == 0) {
        // Initialize histrangetree with first 9 trees
        TTree *histrangetree[11];
        for (Int_t i = 0; i < 9; ++i) {
            histrangetree[i] = trees[i];
        }
        histrangetree[10] = trees[13];

        // Calculate boundaries for each of the 4 trees (trees[9] to trees[12])
        for (Int_t j = 0; j < 4; ++j) {
            histrangetree[9] = trees[9 + j];
            return1DHistRange(histrangetree, {}, "TOFL", momentum, startVal, endVal);
            lower_bound[j] = startVal;
            upper_bound[j] = endVal;
            cout << "j = " << j << " lower_bound is " << lower_bound[j] << " upper_bound is " << upper_bound[j] << endl;
        }
    }

    TGraph *tempgraph = new TGraph();
    graph->SetMarkerColor(kRed);
    tempgraph->SetMarkerColor(kBlue);  
    graph->SetMarkerStyle(kFullSquare);
    tempgraph->SetMarkerStyle(kFullSquare);
    graph->SetMarkerSize(0.3);
    tempgraph->SetMarkerSize(0.3);

    Int_t totalIndex = CalcTotalIndex(trees);
   
    vector<const char*> allBranches;

    Int_t combinedSize = CalcCombinedSize(branches, allBranches, totalIndex);

    Double_t branchVal[combinedSize][20];
    Int_t nPeaks[totalIndex];
    Double_t WWI[totalIndex];
    Int_t nWindowPeak;

    // Initialize the branch values
    InitializeBranchVal(branchVal, combinedSize);

    SetupRootTreesAddress(trees, totalIndex, combinedSize, branches, allBranches, branchVal, nPeaks, WWI, nWindowPeak);

    Double_t xvalue;
    Double_t yvalue;

    Int_t Ncounter = 0, Ncounter_OP = 0;
    
    Long64_t numEntries = trees[0]->GetEntries();
    for (Long64_t entry = 0; entry < numEntries; ++entry) {        
        Double_t sumTOF1 = 0.0, sumTOF0 = 0.0;
        Double_t avgTOF1 = 0.0, avgTOF0 = 0.0, finalTOF = 0.0;
        //Has to find boundaries for each momentums if using IntCharge.
        
        // In order, the lower bound (of the window) are 0L, 0R, 1L, 1R and are in ns.
        /*
        Double_t lb[4] = {-10, -10, -20, -15};
        Double_t ub[4] = {30, 35, 15, 20};
        std::copy(std::begin(lb), std::end(lb), lower_bound);
        std::copy(std::begin(ub), std::end(ub), upper_bound);
        */
        //Double_t ACT0ST = 0.0, ACT1ST = 0.0;
        Double_t ACT0LST[20]={0.0}, ACT0RST[20]={0.0}, ACT1LST[20]={0.0}, ACT1RST[20]={0.0};
        
        // Check if all nPeaks values are 1 only for the TOF and PbGlass. (first 9 trees)
        Double_t WWI_avg = 0.0;
        bool allPeaksOne = true;
        bool eCut = false;
        bool mpCut = false; //not implemented in current ACT0vsACT1?
        bool WWCut = false;
        for (Int_t treeIdx = 0; treeIdx < 9; ++treeIdx){
            trees[treeIdx]->GetEntry(entry);
            if (nPeaks[treeIdx] != 1) {
                allPeaksOne = false;
                break;
            }
        }
        if (allPeaksOne && nWindowPeak==1){
            for (Int_t treeIdx = 0; treeIdx < totalIndex; ++treeIdx){
                trees[treeIdx]->GetEntry(entry);
            }
            Double_t PbGlassPEVal;
            PbGlassPEVal = branchVal[totalIndex*1+8][0];
            for (Int_t i=0; i<4; ++i){
                sumTOF1 += branchVal[totalIndex * 0 + (i + 4)][0];
                sumTOF0 += branchVal[totalIndex * 0 + i][0];
                WWI_avg += WWI[i]+WWI[i+4];
            }
            avgTOF1 = sumTOF1/4;
            avgTOF0 = sumTOF0/4;
            finalTOF = avgTOF1-avgTOF0;
            WWI_avg /=8;
            //if (PbGlassPEVal > 10 && finalTOF < 14){ //correct value for windowIntPE
            
            Double_t threshold_pos = momentum * (0.5 / 1.1) + (-0.5 * 0.1 / 1.1);
            Double_t threshold_neg = momentum * (-0.5 / 1.1) + (-0.5 * 0.1 / 1.1);
            //Double_t threshold_pos = momentum * (0.25-0)/(0.7-0.1) + 0.25 - (0.25/0.6)*0.7;
            //Double_t threshold_neg = momentum * (-0.25-0)/(0.7-0.1) + 0.25 - (0.25/0.6)*0.7;
           
            if (momentum > 0){
                if(PbGlassPEVal > threshold_pos && finalTOF < 13){
                    eCut = true;
                }
                else if(finalTOF < 13 && PbGlassPEVal <= threshold_pos){
                    mpCut = true;
                }
            }
            if (momentum < 0){
                if(PbGlassPEVal > threshold_neg && finalTOF < 13){
                    eCut = true;
                }
                else if(finalTOF < 13 && PbGlassPEVal <= threshold_neg){
                    mpCut = true;
                }
            }
            if(eCut && WWI_avg<40){
                WWCut = true;
            }
            else if(mpCut && WWI_avg<40){ //muons pions
                WWCut = true;
            }
            else if(!eCut && !mpCut && finalTOF > 13){ //protons might need change for 1200MeV and 200MeV
                WWCut = true;
            }
        }
        //cout << "For entry " << entry << " eCut is " << eCut << " mpCut is " << mpCut << " WWIavg is " << WWI_avg << " finalTOF is " << finalTOF << " while htype is " << htype << endl;
        if (eCut && WWCut && strcmp(htype,"wIntCharge")==0){ 
            yvalue = branchVal[totalIndex*1+11][0] + branchVal[totalIndex*1+12][0];
            xvalue = branchVal[totalIndex*1+9][0] + branchVal[totalIndex*1+10][0];
            graph->SetPoint(Ncounter, xvalue, yvalue);
            Ncounter++;
        }
        else if (WWCut && !eCut && strcmp(htype,"wIntCharge")==0) {
            yvalue = branchVal[totalIndex * 1 + 11][0] + branchVal[totalIndex * 1 + 12][0];
            xvalue = branchVal[totalIndex * 1 + 9][0] + branchVal[totalIndex * 1 + 10][0];
            tempgraph->SetPoint(Ncounter_OP, xvalue, yvalue);
            Ncounter_OP++;
        }
        
        else if (eCut && WWCut){
            xvalue = 0;
            yvalue = 0;
            Int_t bVIdx; 
            if(strcmp(htype,"IntCharge")==0){
                bVIdx = 2;
            }
            else if(strcmp(htype,"PVT")==0){
                bVIdx = 3;
            }
            for(Int_t j=0; j<20; ++j){
                ACT0LST[j] = branchVal[totalIndex * 0 + 9][j] - avgTOF0; 
                ACT0RST[j] = branchVal[totalIndex * 0 + 10][j] - avgTOF0;
                //cout << "Entry: " << entry << ", j: " << j << ", ACT0LST: " << ACT0LST[j] << ", ACT0RST: " << ACT0RST[j] << endl;
                //cout << "Lower bounds: " << lower_bound[0] << ", " << lower_bound[1] << ", Upper bounds: " << upper_bound[0] << ", " << upper_bound[1] << endl;
                if(ACT0LST[j]>lower_bound[0] && ACT0RST[j]>lower_bound[1] && ACT0LST[j]<upper_bound[0] && ACT0RST[j]<upper_bound[1]){
                    xvalue = branchVal[totalIndex*bVIdx+9][j] + branchVal[totalIndex*bVIdx+10][j];
                    break;
                }
                else if(ACT0LST[j]>lower_bound[0] && ACT0LST[j]<upper_bound[0]){
                    xvalue = branchVal[totalIndex*bVIdx+9][j];
                    break;
                } 
                else if(ACT0RST[j]>lower_bound[1] && ACT0RST[j]<upper_bound[1]){
                    xvalue = branchVal[totalIndex*bVIdx+10][j];
                    break;
                } 
            }
            for(Int_t j=0; j<20; ++j){
                ACT1LST[j] = branchVal[totalIndex * 0 + 11][j] - avgTOF0;
                ACT1RST[j] = branchVal[totalIndex * 0 + 12][j] - avgTOF0;
                //cout << "Entry: " << entry << ", j: " << j << ", ACT1LST: " << ACT1LST[j] << ", ACT1RST: " << ACT1RST[j] << endl;
                //cout << "Lower bounds: " << lower_bound[2] << ", " << lower_bound[3] << ", Upper bounds: " << upper_bound[2] << ", " << upper_bound[3] << endl;
                if(ACT1LST[j]>lower_bound[2] && ACT1RST[j]>lower_bound[3] && ACT1LST[j]<upper_bound[2] && ACT1RST[j]<upper_bound[3]){
                    yvalue = branchVal[totalIndex * bVIdx + 11][j] + branchVal[totalIndex * bVIdx + 12][j];
                    break;
                }
                else if(ACT1LST[j]>lower_bound[2] && ACT1LST[j]<upper_bound[2]){
                    yvalue = branchVal[totalIndex * bVIdx + 11][j];
                    break;
                } 
                else if(ACT1RST[j]>lower_bound[3] && ACT1RST[j]<upper_bound[3]){
                    yvalue = branchVal[totalIndex * bVIdx + 12][j];
                    break;
                } 
            }
            graph->SetPoint(Ncounter, xvalue, yvalue);
            Ncounter++;
        }
        else if (WWCut && !eCut) {
         
            xvalue = 0;
            yvalue = 0;
            Int_t bVIdx; 
            if(strcmp(htype,"IntCharge")==0){
                bVIdx = 2;
            }
            else if(strcmp(htype,"PVT")==0){
                bVIdx = 3;
            }

            for(Int_t j=0; j<20; ++j){
                ACT0LST[j] = branchVal[totalIndex*0+9][j] - avgTOF0; 
                ACT0RST[j] = branchVal[totalIndex*0+10][j] - avgTOF0;
                //cout << "Non-e Entry: " << entry << ", j: " << j << ", ACT0LST: " << ACT0LST[j] << ", ACT0RST: " << ACT0RST[j] << endl;
                //cout << "Lower bounds: " << lower_bound[0] << ", " << lower_bound[1] << ", Upper bounds: " << upper_bound[0] << ", " << upper_bound[1] << endl;

                if(ACT0LST[j]>lower_bound[0] && ACT0RST[j]>lower_bound[1] && ACT0LST[j]<upper_bound[0] && ACT0RST[j]<upper_bound[1]){
                    xvalue = branchVal[totalIndex*bVIdx+9][j] + branchVal[totalIndex*bVIdx+10][j];
                    break;
                }
                else if(ACT0LST[j]>lower_bound[0] && ACT0LST[j]<upper_bound[0]){
                    xvalue = branchVal[totalIndex*bVIdx+9][j];
                    break;
                } 
                else if(ACT0RST[j]>lower_bound[1] && ACT0RST[j]<upper_bound[1]){
                    xvalue = branchVal[totalIndex*bVIdx+10][j];
                    break;
                } 
            }
            for(Int_t j=0; j<20; ++j){
                ACT1LST[j] = branchVal[totalIndex*0+11][j] - avgTOF0;
                ACT1RST[j] = branchVal[totalIndex*0+12][j] - avgTOF0;
                //cout << "Non-e Entry: " << entry << ", j: " << j << ", ACT1LST: " << ACT1LST[j] << ", ACT1RST: " << ACT1RST[j] << endl;
                //cout << "Lower bounds: " << lower_bound[2] << ", " << lower_bound[3] << ", Upper bounds: " << upper_bound[2] << ", " << upper_bound[3] << endl;
                if(ACT1LST[j]>lower_bound[2] && ACT1RST[j]>lower_bound[3] && ACT1LST[j]<upper_bound[2] && ACT1RST[j]<upper_bound[3]){
                    yvalue = branchVal[totalIndex*bVIdx+11][j] + branchVal[totalIndex*bVIdx+12][j];
                    break;
                }
                else if(ACT1LST[j]>lower_bound[2] && ACT1LST[j]<upper_bound[2]){
                    yvalue = branchVal[totalIndex*bVIdx+11][j];
                    break;
                } 
                else if(ACT1RST[j]>lower_bound[3] && ACT1RST[j]<upper_bound[3]){
                    yvalue = branchVal[totalIndex*bVIdx+12][j];
                    break;
                } 
            } 
            tempgraph->SetPoint(Ncounter_OP, xvalue, yvalue);
            Ncounter_OP++;
        }
        
        InitializeBranchVal(branchVal, combinedSize);
    }

    return tempgraph;
}

void Draw5GraphComp(TGraph *graph, TGraph *graph_OP, Int_t IndexArr[5], Double_t mom, const char* htype, const char* more = ""){
    TCanvas *canvas = new TCanvas("c","title",800,600);
    TLegend *leg = new TLegend(0.1,0.1,0.3,0.4);
    TGraph *grapharr[5];
    Color_t colorArray[5] = {kBlack, kRed, kGreen, kBlue, kYellow};
    for (Int_t i=0; i<5; ++i){
        grapharr[i] = new TGraph();
        grapharr[i]->SetMarkerStyle(kCircle);
        grapharr[i]->SetMarkerSize(0.4);
        grapharr[i]->SetMarkerColor(colorArray[i]);
        grapharr[i]->SetLineColor(colorArray[i]);
    }
    grapharr[0]->SetTitle(Form("Electron rejection vs non-E efficiency at %.1fGeV fixed slope, different intercepts",mom));
    grapharr[0]->GetYaxis()->SetTitle("Electron rejection rate");
    grapharr[0]->GetXaxis()->SetTitle("Non-Electron Efficiency");

    const Int_t numPoints = 50;
    TH2D* pureff_2Dhist_Bg = nullptr;
    TH2D* pureff_2Dhist_Sig = nullptr;
    GenerateHistograms(graph, graph_OP, pureff_2Dhist_Bg, pureff_2Dhist_Sig, numPoints);
    Double_t mVal[numPoints];
    for (Int_t i = 0; i < numPoints; ++i) {
        mVal[i] = -1.0 + i * (0.999 / (numPoints - 1));
    }

    for(Int_t i=0; i<5; ++i){
        for(Int_t j=1; j<=50; j++){
            Double_t eff_bg = pureff_2Dhist_Bg->GetBinContent(IndexArr[i]+1, j);
            Double_t eff_sig = pureff_2Dhist_Sig->GetBinContent(IndexArr[i]+1, j);
            grapharr[i]->SetPoint(j-1, eff_sig, 1/eff_bg);
        }
        leg->AddEntry(grapharr[i], Form("%s m = %.2f", htype, mVal[IndexArr[i]]), "pl");
    }

    grapharr[0]->GetXaxis()->SetRangeUser(0,1);
    grapharr[0]->GetYaxis()->SetRangeUser(0,1000);
    gPad->SetLogy();
    grapharr[0]->Draw("APL");
    for(Int_t i=1; i<5; ++i){
        grapharr[i]->Draw("Same PL");
    }
    leg->Draw();
    canvas->SaveAs(Form("graphjpg/Final_EffComp%s%s%.1f.png",htype,more,mom));
    delete pureff_2Dhist_Bg;
    delete pureff_2Dhist_Sig;
    delete canvas;
}

#endif // HISTOGxRAM_FUNCTIONS_H
