#include "histogram_functions.h" // Include the file with histogram functions
#include "data_runs_dicts.h"

void returnResultsElectron(Int_t runNumber, Double_t electronRejectionRate){
    
    // Open the output text file
    std::ofstream outputFile("ACTSingleResults.txt", std::ios_base::trunc); // Append to the file
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return;
    }

    // Find momentum for the current run number
    TString fileName = Form("WindowIntMatched_final_%i.root", runNumber);
    auto it = runsDict.find(runNumber);
    if (it == runsDict.end()) {
        std::cerr << "Error: Momentum not found for run number " << runNumber << std::endl;
        return;
    }
    Double_t momentum = static_cast<Double_t>(it->second) / 1000.0;

    // Open the ROOT file
    TFile* rootFile = TFile::Open(fileName);
    if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        if (rootFile) {
            rootFile->Close();
            delete rootFile;
        }
        return;
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
        return;
    }

    TTree *completeTTree[14];
    for(Int_t i=0; i<4; ++i){
        completeTTree[i] = TOF[i];
        completeTTree[i+4] = TOF[i+4];
        completeTTree[9+i] = ACT01[i];
    }
    completeTTree[8] = pbglass;
    completeTTree[13] = nullptr;

    
    // Write the header to the text file
    outputFile << "Run number " << runNumber << " with momentum " << momentum << "GeV. All cuts kepts.\n"; 
    outputFile << "Signal Efficiency, Background Rejection Rate, Slope(m), Intercept(b)\n";


    TGraph* graph = new TGraph();
    //wIntCharge set up.
    TGraph *graph_OP = MakeSimple2DGraph0vs1NoDraw(graph, completeTTree, {"IntCharge", "MaxVoltage"}, momentum);

    FindAllEfficiencyWithRejectionRate(graph, graph_OP, electronRejectionRate, runNumber, momentum, outputFile);

    outputFile.close();
    Draw2500PointWithResults(graph, graph_OP, momentum, "wIntCharge", electronRejectionRate);

    rootFile->Close();
    delete rootFile;
    delete graph;
    delete graph_OP;

}