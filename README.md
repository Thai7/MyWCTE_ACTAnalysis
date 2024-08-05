# MyWCTE_ACTAnalysis
Current work on the ACT Analysis.
All files are ROOT macros except for functions.h and histogram_functions.h which are the core functions used in each macros. The data_runs_dict.h is a map that link the run number with the momentum in GeV.

Each macros are called directly (e.g root ACT_analysis.cpp). The one that are user friendly are currently the two returnResults macros. They need to first use the pre-processing in WCTE repository and get the WindowIntMatched_final_run.root file.

To use: root -q -l returnResultsElectron(RunNumber, threshold).

One functions use a electron rejection threshold to find the bset average cuts in ACT0vsACT1 windowIntCharge plot. The other one use instead the non-electron efficiency.
