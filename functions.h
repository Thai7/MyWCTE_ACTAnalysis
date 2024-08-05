#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <vector>
#include <initializer_list>
#include <string>

void LoadTrees(TFile *ntuplefile, const char* TOFTreeNames[], const char* ACTTreeNames[], const char* ACT23TreeNames[], TTree *TOFtrees[], TTree *ACTtrees[], TTree *ACT23trees[], TTree* &pbglasstree) {
    for (int i = 0; i < 8; ++i) {
        TOFtrees[i] = (TTree*)ntuplefile->Get(TOFTreeNames[i]);
    }
    for (int i =0; i <4 ; ++i){
        ACTtrees[i] = (TTree*)ntuplefile->Get(ACTTreeNames[i]);
        ACT23trees[i] = (TTree*)ntuplefile->Get(ACT23TreeNames[i]);

    }
    pbglasstree = (TTree*)ntuplefile->Get("PbGlass");
}

void SetupTrees(TFile *ntuplefile, TTree *TOFtree[], TTree *ACTtree[], TTree *ACT23tree[], TTree* &pbglasstree) {
    const char* TOF_TREE_NAMES[8] = {"TOF00", "TOF01", "TOF02", "TOF03", "TOF10", "TOF11", "TOF12", "TOF13"};
    const char* ACT_TREE_NAMES[4] = {"ACT0L", "ACT0R", "ACT1L", "ACT1R"};
    const char* ACT23_TREE_NAMES[4] = {"ACT2L", "ACT2R", "ACT3L", "ACT3R"};

    LoadTrees(ntuplefile, TOF_TREE_NAMES, ACT_TREE_NAMES, ACT23_TREE_NAMES, TOFtree, ACTtree, ACT23tree, pbglasstree);
}

Int_t CalcTotalIndex(TTree *trees[]){
    // Iterate through the array until a nullptr is encountered
    Int_t totalIndex = 0;
    while (trees[totalIndex] != nullptr){
        // Process trees[totalIndex]
        totalIndex++;
    }
    return totalIndex;
}

void InitializeBranchVal(Double_t branchVal[][20], int combinedSize, Double_t initialValue = 15000.0) {
    for (Int_t p = 0; p < combinedSize; ++p) {
        for (Int_t k = 0; k < 20; ++k) {
            branchVal[p][k] = initialValue;
        }
    }
}

Int_t CalcCombinedSize(initializer_list<const char*> branches, vector<const char*> &allBranches, Int_t totalIndex){
    // Default branches
    vector<const char*> defaultBranches = {"SignalTimeCorrected", "WindowIntCharge"};
    
    // Combine default branches with user-specified branches
    allBranches.reserve(defaultBranches.size() + branches.size());
    allBranches.insert(allBranches.end(), defaultBranches.begin(), defaultBranches.end());
    allBranches.insert(allBranches.end(), branches.begin(), branches.end());

    Int_t combinedSize = allBranches.size() * totalIndex;

    return combinedSize;
}

void SetupRootTreesAddress(TTree *trees[], Int_t totalIndex, Int_t combinedSize, initializer_list<const char*> branches, vector<const char*> &allBranches, Double_t branchVal[][20], Int_t nPeaks[], Double_t WWI[], Int_t &nWindowPeak) {

    // Activating both default trees and the one added as arguments
    for (Int_t treeIdx = 0; treeIdx < totalIndex; ++treeIdx) {
        trees[treeIdx]->SetBranchStatus("*", 0);
        trees[treeIdx]->SetBranchStatus("nPeaks", 1);
        trees[treeIdx]->SetBranchStatus("WholeWaveformInt", 1);
        trees[treeIdx]->SetBranchStatus("SignalTimeCorrected", 1);
        trees[treeIdx]->SetBranchStatus("WindowIntCharge", 1);
        for (auto branch : branches) {
            trees[treeIdx]->SetBranchStatus(branch, 1);
        }
    }

    // Set the branch status and address for nWindowPeaks in the first tree
    trees[0]->SetBranchStatus("nWindowPeaks", 1);
    trees[0]->SetBranchAddress("nWindowPeaks", &nWindowPeak);
    
    // Set branch addresses
    Int_t index = 0;
    for (auto branch : allBranches) {
        for (Int_t j = 0; j < totalIndex; ++j) {
            int p = totalIndex * index + j; // Combined index
            trees[j]->SetBranchAddress(branch, branchVal[p]); // Set branch address for current branch and tree
        }
        ++index;
    }
    
    for (Int_t j = 0; j < totalIndex; ++j) {
        trees[j]->SetBranchAddress("nPeaks", &nPeaks[j]);
        trees[j]->SetBranchAddress("WholeWaveformInt", &WWI[j]);
    }
}

bool isPointUnderLine(Double_t m, Double_t b, Double_t x, Double_t y) {
    Double_t computedY = m*x + b;
    return y < computedY;
}
void calculatePurityAndEfficiency(TGraph *graph, TGraph *graph_OP, Double_t m, Double_t b, Double_t &purity, Double_t &eff, Double_t &effBg, Double_t &purityBg) {
    Int_t totalN = graph->GetN();
    Int_t totalN_OP = graph_OP->GetN();
    
    Int_t underline = 0;
    Int_t underline_OP = 0;
    
    for (Int_t i = 0; i < totalN; ++i) {
        Double_t x, y;
        graph->GetPoint(i, x, y);
        
        if (isPointUnderLine(m, b, x, y)) {
            ++underline;
        }
    }
    
    for (Int_t i = 0; i < totalN_OP; ++i) {
        Double_t x, y;
        graph_OP->GetPoint(i, x, y);
        
        if (isPointUnderLine(m, b, x, y)) {
            ++underline_OP;
        }
    }
    
    purity = (Double_t)underline_OP / (underline_OP + underline);
    eff = (Double_t)underline_OP / totalN_OP;
    effBg = ((Double_t)underline/totalN); //electrno SIGNAL efficiency under the curve
    purityBg = ((Double_t)totalN-underline)/(totalN_OP+totalN-underline_OP-underline); //PURITY of electron SIGNAL

}

void ReturnBestFitValue(TGraph *graph, TGraph *graph_OP, const string &text, Double_t *mVal, Double_t *bVal, Int_t numPoints){
    Double_t finalm, finalb;
    Int_t maxBin, xBin, yBin, dummy;
    TH2D *purreff_temporary = new TH2D("purreff_temporary", "Temporary", numPoints, mVal[0], mVal[numPoints - 1], numPoints, bVal[0], bVal[numPoints - 1]);
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            Double_t m = mVal[i];
            Double_t b = bVal[j];
            Double_t purity, eff;
            Double_t puritybg, effbg;
            Double_t purity2, eff2;
            Double_t puritybg2, effbg2;                
            
            calculatePurityAndEfficiency(graph, graph_OP, m, b, purity, eff, effbg, puritybg);
            purreff_temporary->SetBinContent(i + 1, j + 1, (puritybg*effbg)*100);
        }
    }
    maxBin = purreff_temporary->GetMaximumBin();
    purreff_temporary->GetBinXYZ(maxBin, xBin, yBin, dummy);
    finalm = -1.0+(xBin-1)*0.999/(50-1);
    finalb = 0.01+(yBin-1)*0.34/49;
    cout << "MaxZ " << text << purreff_temporary->GetMaximum() << " and x(m) is " <<  finalm << " and in y(b) is " << finalb << endl;   
    delete purreff_temporary; // Clean up
}

void GenerateHistograms(TGraph* graph, TGraph* graph_OP, TH2D*& pureff_2Dhist_Bg, TH2D*& pureff_2Dhist_Sig, const Int_t numPoints) {
    Double_t mVal[numPoints];
    Double_t bVal[numPoints];

    // Generate m values from -1 to -0.001
    for (Int_t i = 0; i < numPoints; ++i) {
        mVal[i] = -1.0 + i * (0.999 / (numPoints - 1));
    }

    // Generate b values from 0.01 to 0.35
    for (Int_t i = 0; i < numPoints; ++i) {
        bVal[i] = 0.01 + i * (0.35 - 0.01) / (numPoints - 1);
    }

    // Create 2D histograms for purity * efficiency
    pureff_2Dhist_Bg = new TH2D("pureff_2Dhist_Bg", "Efficiency background(electron);Slope;Intercept", numPoints, -1, -0.001, numPoints, 0.01, 0.35);
    pureff_2Dhist_Sig = new TH2D("pureff_2Dhist_Sig", "Efficiency signal(others particles);Slope;Intercept", numPoints, -1, -0.001, numPoints, 0.01, 0.35);
    //for cuts only
    TH2D *temph1 = new TH2D("temph1", "Purity * Efficiency background(electron);Slope;Intercept", numPoints, -1, -0.001, numPoints, 0.01, 0.35);
    TH2D *temph2 = new TH2D("temph2", "Purity * Efficiency signal(others particles);Slope;Intercept", numPoints, -1, -0.001, numPoints, 0.01, 0.35);

    // Loop through the slopes and intercepts
    for (Int_t i = 0; i < numPoints; ++i) {
        for (Int_t j = 0; j < numPoints; ++j) {
            Double_t m = mVal[i];
            Double_t b = bVal[j];
            Double_t purity, eff;
            Double_t puritybg, effbg;

            // Calculate purity and efficiency
            calculatePurityAndEfficiency(graph, graph_OP, m, b, purity, eff, effbg, puritybg);

            Double_t value = (effbg);
            Double_t value2 = (eff);
            pureff_2Dhist_Bg->SetBinContent(i + 1, j + 1, value);
            pureff_2Dhist_Sig->SetBinContent(i + 1, j + 1, value2);
            
            //To find best cut.
            temph1->SetBinContent(i+1,j+1,(1-effbg)*(puritybg));
            temph2->SetBinContent(i+1,j+1,eff*purity);
            temph1->SetStats(0);
            temph2->SetStats(0);
        }
    }
    //Find best cut
    int binx_max, biny_max, binx_min = 0, biny_min = 0, binz, binz2; // binz is not used for 2D histograms
    double tempmaxVal;
    double tempminVal = temph1->GetBinContent(1, 1);

    // Find maximum bin, the min is used for the background (electron) when its maximized in the other area...
    temph2->GetMaximumBin(binx_max, biny_max, binz);
    tempmaxVal = temph2->GetBinContent(binx_max, biny_max);
    temph1->GetMaximumBin(binx_min, biny_min, binz2);
    tempminVal = temph1->GetBinContent(binx_min, biny_min);

     // Convert bin positions to finalm and finalb values
    double finalm_max = -1.0 + (binx_max - 1) * 0.999 / (numPoints - 1);
    double finalb_max = 0.01 + (biny_max - 1) * 0.34 / 49;
    double finalm_min = -1.0 + (binx_min - 1) * 0.999 / (numPoints - 1);
    double finalb_min = 0.01 + (biny_min - 1) * 0.34 / 49;
    //cout << "The best fit value for the signal(non-e) is m = " << finalm_max << " and b =" << finalb_max << " at " << tempmaxVal << " purity*eff" << endl; 
    //cout << "The best fit value for the background(e) is m = " << finalm_min << " and b =" << finalb_min << " at " << tempminVal << " purity*eff" << endl; 
    delete temph1;
    delete temph2;
    //TCanvas *tempcanvas = new TCanvas();
    //pureff_2Dhist_Sig->Draw("colz");
    //tempcanvas->SaveAs("Signal.png");
}

void FindMaxPurityEfficiency(TGraph* graph, TGraph* graph_OP, Double_t rejectionRate, Double_t& maxVal, Double_t& finalEff, Double_t& finalm, Double_t& finalb, bool& success) {
    const Int_t numPoints = 50;
    TH2D* pureff_2Dhist_Bg = nullptr;
    TH2D* pureff_2Dhist_Sig = nullptr;

    // Generate histograms and fill them
    GenerateHistograms(graph, graph_OP, pureff_2Dhist_Bg, pureff_2Dhist_Sig, numPoints);

    // Calculate the efficiency threshold from rejectionRate
    Double_t efficiencyThreshold = 1.0 / rejectionRate;
    
    // Find the maximum bin in the signal histogram respecting the rejection rate in the background histogram
    maxVal = -1.0;
    Int_t maxBinX = -1;
    Int_t maxBinY = -1;
    success = false;

    for (Int_t i = 1; i <= pureff_2Dhist_Bg->GetNbinsX(); ++i) {
        for (Int_t j = 1; j <= pureff_2Dhist_Bg->GetNbinsY(); ++j) {
            if (pureff_2Dhist_Bg->GetBinContent(i, j) < efficiencyThreshold) {
                Double_t currentVal = pureff_2Dhist_Sig->GetBinContent(i, j);
                if (currentVal > maxVal) {
                    maxVal = currentVal;
                    maxBinX = i;
                    maxBinY = j;
                    success = true;
                }
            }
        }
    }

    if (!success) {
        std::cerr << "Error: No valid bin found that respects the rejection rate." << std::endl;
        delete pureff_2Dhist_Bg;
        delete pureff_2Dhist_Sig;
        return;
    }

    // Calculate the final efficiency and corresponding m and b values
    finalEff = pureff_2Dhist_Bg->GetBinContent(maxBinX, maxBinY);
    finalm = -1.0 + (maxBinX - 1) * 0.999 / (numPoints - 1);
    finalb = 0.01 + (maxBinY - 1) * 0.34 / 49;

    // Clean up
    delete pureff_2Dhist_Bg;
    delete pureff_2Dhist_Sig;
}

void FindAllEfficiencyWithRejectionRate(TGraph* graph, TGraph* graph_OP, Double_t rejectionRate, Int_t runNumber, Int_t momentum, std::ofstream& outputFile) {
    const Int_t numPoints = 50;
    TH2D* pureff_2Dhist_Bg = nullptr;
    TH2D* pureff_2Dhist_Sig = nullptr;

    // Generate histograms and fill them
    GenerateHistograms(graph, graph_OP, pureff_2Dhist_Bg, pureff_2Dhist_Sig, numPoints);

    // Calculate the efficiency threshold from rejectionRate
    Double_t efficiencyThreshold = 1.0 / rejectionRate;
    Double_t efficiencyThreshold2 = 1.0/(rejectionRate*1.1);

    // Arrays to hold the m and b values
    Double_t mVal[numPoints];
    Double_t bVal[numPoints];

    // Generate m values from -1 to -0.001
    for (Int_t i = 0; i < numPoints; ++i) {
        mVal[i] = -1.0 + i * (0.999 / (numPoints - 1));
    }

    // Generate b values from 0.01 to 0.35
    for (Int_t i = 0; i < numPoints; ++i) {
        bVal[i] = 0.01 + i * (0.35 - 0.01) / (numPoints - 1);
    }

    // Loop through the background histogram and remove bins in the signal histogram
    for (Int_t i = 1; i <= pureff_2Dhist_Bg->GetNbinsX(); ++i) {
        for (Int_t j = 1; j <= pureff_2Dhist_Bg->GetNbinsY(); ++j) {
            Double_t bgContent = pureff_2Dhist_Bg->GetBinContent(i, j);
            if (bgContent > efficiencyThreshold || bgContent < efficiencyThreshold2) {
                pureff_2Dhist_Sig->SetBinContent(i, j, 0); // Zero out bins that do not respect the threshold
            } 
            else {
                Double_t sigContent = pureff_2Dhist_Sig->GetBinContent(i, j);
                Double_t rejectionRateFound = 1.0 / bgContent;
                Double_t finalm = mVal[i - 1];
                Double_t finalb = bVal[j - 1];

                // Write to the output file
                outputFile << sigContent << ", " << rejectionRateFound << ", " << finalm << ", " << finalb << "\n";
            }    
        }
    }

    // Clean up
    //delete pureff_2Dhist_Bg;
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->cd();
    pureff_2Dhist_Sig->SetStats(kFALSE);
    pureff_2Dhist_Sig->Draw("colz");
    canvas->Update();
    canvas->SaveAs(Form("ACT_EffResults%i.png",runNumber));
}

void FindAllRejectionRateWithEfficiency(TGraph* graph, TGraph* graph_OP, Double_t efficiencyThreshold, Int_t runNumber, Int_t momentum, std::ofstream& outputFile) {
    const Int_t numPoints = 50;
    TH2D* pureff_2Dhist_Bg = nullptr;
    TH2D* pureff_2Dhist_Sig = nullptr;

    // Generate histograms and fill them
    GenerateHistograms(graph, graph_OP, pureff_2Dhist_Bg, pureff_2Dhist_Sig, numPoints);

    // Calculate the efficiency threshold from rejectionRate
    Double_t efficiencyThreshold2 = efficiencyThreshold*1.1;

    // Arrays to hold the m and b values
    Double_t mVal[numPoints];
    Double_t bVal[numPoints];

    // Generate m values from -1 to -0.001
    for (Int_t i = 0; i < numPoints; ++i) {
        mVal[i] = -1.0 + i * (0.999 / (numPoints - 1));
    }

    // Generate b values from 0.01 to 0.35
    for (Int_t i = 0; i < numPoints; ++i) {
        bVal[i] = 0.01 + i * (0.35 - 0.01) / (numPoints - 1);
    }

    // Loop through the signal histogram and remove bins in the background histogram
    for (Int_t i = 1; i <= pureff_2Dhist_Sig->GetNbinsX(); ++i) {
        for (Int_t j = 1; j <= pureff_2Dhist_Sig->GetNbinsY(); ++j) {
            Double_t sigContent = pureff_2Dhist_Sig->GetBinContent(i, j);
            if (sigContent < efficiencyThreshold || sigContent > efficiencyThreshold2) {
                pureff_2Dhist_Bg->SetBinContent(i, j, 0); // Zero out bins that do not respect the threshold
            } 
            else {
                Double_t bgContent = pureff_2Dhist_Bg->GetBinContent(i, j);
                Double_t rejectionRateFound = 1.0 / bgContent;
                Double_t finalm = mVal[i - 1];
                Double_t finalb = bVal[j - 1];

                // Write to the output file
                outputFile << sigContent << ", " << rejectionRateFound << ", " << finalm << ", " << finalb << "\n";
            }    
        }
    }

    // Clean up
    //delete pureff_2Dhist_Sig;
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->cd();
    pureff_2Dhist_Bg->SetStats(kFALSE);
    pureff_2Dhist_Bg->Draw("colz");
    canvas->Update();
    //canvas->SaveAs(Form("ACT_EffResults%i.png",runNumber));
}



void calculateAveragesFromFile(Double_t &avgSignalEff, Double_t &avgRejectionRate, Double_t &avgFinalM, Double_t &avgFinalB){
    // Open the file for reading
    std::ifstream inputFile("ACTSingleResults.txt");
    double sumSignalEff = 0.0;
    double sumRejectionRate = 0.0;
    double sumFinalM = 0.0;
    double sumFinalB = 0.0;
    int count = 0;

    std::string line;

    // Skip the header lines
    std::getline(inputFile, line); // Skip first header line
    std::getline(inputFile, line); // Skip second header line (if applicable)

    // Process each line
    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> tokens;

        // Split the line by commas
        while (std::getline(ss, item, ',')) {
            tokens.push_back(item);
        }

        if (tokens.size() < 4) {
            continue;  // Skip lines that do not have enough columns
        }

        try {
            double signalEff = std::stod(tokens[0]);
            double rejectionRate = std::stod(tokens[1]);
            double finalM = std::stod(tokens[2]);
            double finalB = std::stod(tokens[3]);

            // Accumulate sums
            sumSignalEff += signalEff;
            sumRejectionRate += rejectionRate;
            sumFinalM += finalM;
            sumFinalB += finalB;

            // Increment count
            count++;
        } 
        catch (const std::invalid_argument&) {
            std::cerr << "Skipping invalid line: " << line << std::endl;
        }
    }
    inputFile.close();
    avgSignalEff = sumSignalEff / count;
    avgRejectionRate = sumRejectionRate / count;
    avgFinalM = sumFinalM / count;
    avgFinalB = sumFinalB / count;
}
#endif // FUNCTIONS_H
