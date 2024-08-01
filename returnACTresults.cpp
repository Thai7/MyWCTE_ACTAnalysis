#include "histogram_functions.h" // Include the file with histogram functions
#include "data_runs_dicts.h"

void returnACTresults(){
    Double_t rejectionRate = 100.;

    //Get The files to generate the ACTs plots has to be changed when implemented into Alie.
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);
    TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 800, 600);
    TH2D* hist = new TH2D("h", Form("Run number vs. Signal efficiency w/ rejection rate of %.0f; Run number; Efficiency",rejectionRate), 215, 335, 550, 100, 0, 1);
    TH2D* hist2 = new TH2D("h2", Form("Momentum vs. Signal efficiency w/ rejection rate of %.0f; Momentum(GeV/c); Efficiency",rejectionRate), 48, -1.2, 1.2, 74, 0, 0.74);
    hist->SetZTitle("Rejection Rate");
    hist2->SetZTitle("Rejection Rate");
    TGraph *graph_result = new TGraph();
    graph_result->GetXaxis()->SetTitle("Momentum(GeV/c)");
    graph_result->GetYaxis()->SetTitle("Signal efficiency");
    graph_result->SetTitle(Form("Momentum vs Signal efficiency w/ rejection rate of %.0f",rejectionRate));
    graph_result->SetMarkerStyle(kFullSquare);
    // Open the directory
    const char* directoryPath = "/home/tan/Beam/T9_Analyzed";
    TSystemDirectory dir(directoryPath, directoryPath);
    TList* filesList = dir.GetListOfFiles();

    if (!filesList) {
        std::cerr << "Error: Directory not found or empty: " << directoryPath << std::endl;
        return;
    }

    // Open the output text file
    std::ofstream outputFile("ACTRunsResults.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return;
    }

    // Write the header to the text file
    outputFile << "RunNumber, Momentum(GeV/c), Signal Efficiency, Background Rejection Rate, Slope(m), Intercept(b)\n";

    // Loop over all files in the directory
    TIter next(filesList);
    TSystemFile* file;
    Int_t pointCount = 0;
    Int_t under100Counter = 0;
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
            Double_t momentum = static_cast<Double_t>(it->second) / 1000.0;

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

            TTree *TOF[8];
            TTree *ACT01[4];
            TTree *ACT23[4];
            TTree *pbglass;
            SetupTrees(rootFile,TOF,ACT01,ACT23,pbglass);

            // Check if PbGlass tree exists
            if (!pbglass) {
                std::cerr << "Skipping file " << fileName << " due to missing PbGlass tree." << std::endl;
                rootFile->Close();
                delete rootFile;
                continue;
            }

            TTree *completeTTree[14];
            for(Int_t i=0; i<4; ++i){
                completeTTree[i] = TOF[i];
                completeTTree[i+4] = TOF[i+4];
                completeTTree[9+i] = ACT01[i];
            }
            completeTTree[8] = pbglass;
            completeTTree[13] = nullptr;

            TGraph* graph = new TGraph();
            //wIntCharge set up.
            //TGraph *graph_OP = MakeSimple2DGraph0vs1NoDraw(graph, completeTTree, {"IntCharge", "MaxVoltage"}, momentum);
            //Double_t finalm = -0.102939, finalb = 0.0308163; //equivalent to index 44 and index 4 
            //Double_t eff, effbg, pur, purbg;
            //calculatePurityAndEfficiency(graph, graph_OP, finalm, finalb, pur, eff, effbg, purbg);
            //IntCharge set up.
            TGraph *graph_OP = MakeSimple2DGraph0vs1NoDraw(graph, completeTTree, {"IntCharge", "MaxVoltage"}, momentum, "IntCharge");
            Double_t finalm = -0.143714, finalb = 0.0238776; 
            Double_t eff, effbg, pur, purbg;
            calculatePurityAndEfficiency(graph, graph_OP, finalm, finalb, pur, eff, effbg, purbg);

            if((1.0/effbg)<100){
                cout << "Run " << runNumber << " has rejection < 100. " << endl;
                under100Counter++;
            }
            hist->Fill(runNumber, eff, 1.0/effbg);
            hist2->Fill(momentum, eff, 1.0/effbg);
            graph_result->SetPoint(pointCount, momentum, eff);
            pointCount++;
            outputFile << runNumber << ", " << momentum << ", " << eff << ", " << 1.0/effbg << ", " << finalm << ", " << finalb << "\n";
            //If you need to optimize
            /* 
            Double_t maxVal, finalEff;
            Double_t finalm, finalb;
            bool success;
            FindMaxPurityEfficiency(graph, graph_OP, rejectionRate, maxVal, finalEff, finalm, finalb, success);
            if (success) {
                Double_t rejectionRateFound = 1.0 / finalEff;
                hist->Fill(runNumber, maxVal, rejectionRateFound);
                hist2->Fill(momentum, maxVal, rejectionRateFound);
                graph_result->SetPoint(pointCount,momentum,maxVal);
                pointCount++;
                outputFile << runNumber << ", " << momentum << ", " << maxVal << ", " << rejectionRateFound << ", " << finalm << ", " << finalb << "\n";
            } else {
                std::cerr << "Skipping run " << runNumber << " due to no valid bin found." << std::endl;
            }
            */

            delete graph;
            delete graph_OP;
            rootFile->Close();
            delete rootFile;
        }
    }
    outputFile.close();
    
    cout << "There are " << under100Counter << " runs where rejection < 100" << endl;

    canvas->cd();
    hist->Draw("colz");
    hist->SetStats(kFALSE);  // Disable statistics box
    canvas->Update();
    canvas->SaveAs(Form("ACTRunsResults%.0f.png",rejectionRate));
    canvas2->cd();
    graph_result->Draw("AP");
    canvas2->SaveAs(Form("ACTRunsResultsGraphs%.0f.png",rejectionRate));
    canvas3->cd();
    hist2->Draw("colz");
    hist2->SetStats(kFALSE);
    canvas3->SaveAs(Form("ACTRunsResultsMom%.0f.png",rejectionRate));
}