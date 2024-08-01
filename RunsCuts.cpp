#include "data_runs_dicts.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

void processFilesInDirectory(const char* directoryPath, std::map<int, int>& runsDict) {
    // Create a canvas and a 2D histogram
    TCanvas* canvas = new TCanvas("canvas", "Momentum vs. Charge", 800, 600);
    TH2D* hist = new TH2D("momentum_vs_charge", "Momentum vs. Charge; Momentum(GeV); Charge", 48, -1.2, 1.2, 200, 0, 1);

    TCanvas* canvas2 = new TCanvas("canvas2", "Momentum vs TOF", 800,600);
    TH2D* hist2 = new TH2D("momentum_vs_TOF", "Momentum vs TOF using TOF11; Momentum(GeV); TOF(ns)", 48, -1.2, 1.2, 50, 10, 15);

    // Open the directory
    TSystemDirectory dir(directoryPath, directoryPath);
    TList* filesList = dir.GetListOfFiles();

    if (!filesList) {
        std::cerr << "Error: Directory not found or empty: " << directoryPath << std::endl;
        return;
    }

    // Loop over all files in the directory
    TIter next(filesList);
    TSystemFile* file;
    while ((file = (TSystemFile*)next())) {
        const char* fileName = file->GetName();
        TString filePath = TString::Format("%s/%s", directoryPath, fileName);

        // Check if the file matches the pattern
        if (TString(fileName).Contains("WindowIntMatched_final_") && TString(fileName).EndsWith(".root")) {
            // Extract run number from file name
            TString fileStr(fileName);
            fileStr.ReplaceAll("WindowIntMatched_final_", "");
            fileStr.ReplaceAll(".root", "");
            int runNumber = fileStr.Atoi();

            // Find momentum for the current run number
            auto it = runsDict.find(runNumber);
            if (it == runsDict.end()) {
                std::cerr << "Error: Momentum not found for run number " << runNumber << std::endl;
                continue;
            }
            int momentum = it->second;

            // Open the ROOT file
            TFile* rootFile = TFile::Open(filePath);
            if (!rootFile || rootFile->IsZombie()) {
                std::cerr << "Error opening file: " << fileName << std::endl;
                if (rootFile) {
                    rootFile->Close();
                    delete rootFile;
                }
                continue;
            }

            std::cout << "Processing file: " << fileName << std::endl;

            // Access the PbGlass tree and IntCharge branch
            TTree* tree = (TTree*)rootFile->Get("PbGlass");
            if (!tree) {
                std::cerr << "Error: Tree 'PbGlass' not found in file: " << fileName << std::endl;
                rootFile->Close();
                delete rootFile;
                continue;
            }
            TTree* TOFtrees[8];
            const char* TOF_TREE_NAMES[8] = {"TOF00", "TOF01", "TOF02", "TOF03", "TOF10", "TOF11", "TOF12", "TOF13"};
            for (int i = 0; i < 8; ++i) {
                TOFtrees[i] = (TTree*)rootFile->Get(TOF_TREE_NAMES[i]);
                if (!TOFtrees[i]) {
                    std::cerr << "Error: Tree 'TOF00' not found in file: " << fileName << std::endl;
                    rootFile->Close();
                    delete rootFile;
                    continue;
                }

            }

            // Variables to hold branch data
            const Int_t maxIntCharge = 40;
            Double_t intCharge[maxIntCharge];
            Double_t TOFall[8][10];
            // Set branch address
            tree->SetBranchAddress("IntCharge", intCharge);
            for(Int_t i=0; i<8; ++i){
                TOFtrees[i]->SetBranchAddress("SignalTimeCorrected",&TOFall[i]);
            }

            // Get number of entries
            Long64_t nEntries = tree->GetEntries();
            if (nEntries == 0) {
                std::cerr << "Warning: No entries found in file: " << fileName << std::endl;
                rootFile->Close();
                delete rootFile;
                continue;
            }

            // Loop over all entries (events) in the tree
            for (Long64_t entry = 0; entry < nEntries; ++entry) {
                Double_t avgTOF1 = 0.0, avgTOF0 = 0.0;
                tree->GetEntry(entry);
                for(Int_t i=0; i<8; i++){
                    TOFtrees[i]->GetEntry(entry);
                }
                for (Int_t i=0; i<4; ++i){
                    avgTOF1 += TOFall[i+4][0];
                    avgTOF0 += TOFall[i][0];
                }
                avgTOF1 /=4;
                avgTOF0 /=4;

                // Fill the histogram with momentum (x-axis) and charge (y-axis)
                hist->Fill(momentum/1000.0, intCharge[0]); // Assuming first element of IntCharge array
                hist2->Fill(momentum/1000.0,avgTOF1-avgTOF0);
            }

            // Close the ROOT file
            rootFile->Close();
            delete rootFile;
        }
    }
    Int_t nBinsX = 48;
    Int_t nBinsY1 = 50, nBinsY2 = 200;

    for (int binX = 1; binX <= nbinsX; ++binX) {
        double entriesInXBin = hist->Integral(binX, binX, 1, nbinsY);
        if (entriesInXBin > 0) {
            for (int binY = 1; binY <= nbinsY1; ++binY) {
                double content = hist->GetBinContent(binX, binY);
                hist->SetBinContent(binX, binY, content / entriesInXBin);
            }
            for(int binY = 1; binY <= nBinsY2; ++binY){
                double content = hist->GetBinContent(binX, binY);
                hist2->SetBinContent(binX, binY, content / entriesInXBin);
            }
        }
    }

    delete filesList;

    //functions
    TF1* lineFunction1 = new TF1("lineFunction", "[0]*x + [1]", 0.1, 1.2);
    lineFunction1->SetParameters(0.5/1.1, -0.5*0.1/1.1);
    TF1* lineFunction2 = new TF1("lineFunction", "[0]*x + [1]", -1.2, -0.1);
    lineFunction2->SetParameters(-0.5/1.1, -0.5*0.1/1.1);
    lineFunction1->SetLineColor(kBlack);
    lineFunction1->SetLineStyle(1);
    lineFunction2->SetLineColor(kBlack);
    lineFunction2->SetLineStyle(1);
    TLine* line = new TLine(-1.2, 12.5, 1.2, 12.5);
    line->SetLineColor(kBlack);
    line->SetLineStyle(1);

    // Draw and save the histogram
    canvas->cd();
    canvas->SetLogz();
    hist->Draw("colz");
    hist->SetStats(kFALSE);  // Disable statistics box
    lineFunction1->Draw("same");
    lineFunction2->Draw("same");
    canvas->Update();
    canvas->SaveAs("momentum_vs_charge.png");
    canvas2->cd();
    canvas2->SetLogz();
    hist2->Draw("colz");
    hist2->SetStats(kFALSE);  // Disable statistics box
    line->Draw("same");
    canvas2->Update();
    canvas2->SaveAs("momentum_vs_TOF.png");
    //delete canvas;
    //delete hist;
}

void RunsCuts() {
    const char* directoryPath = "/home/tan/Beam/T9_Analyzed";

    // Call the processing function with runsDict from data_runs_dicts.h
    processFilesInDirectory(directoryPath, runsDict);
}