# VerticalModifiedMinMod

## Running the simulation model
Below we explain how to run the 8 different scenarios presented in the manuscript «Effects of migrating 
mesopelagic fish on carbon export and sequestration». A detailed explanation of the code is given in the 
file “WalkthroughOfCode.mlx”.

To run the model MATLAB version R2022a or later is required.

## Running the baseline scenario
The baseline scenario will run if the main program (file: VertMinMod.m) is started upon downloading.
The simulation starts with initial values for all state variables (except fish and oxygen) that are 
homogenous with depth. The steady state depth distributions for the baseline are reached after about 7 
hours runtime with a machine equipped with Intel(R) Core(TM) i5-7300U CPU @ 2.60GHz 2.71Ghz 
and 16GB RAM.

When the run is finished, results are stored in the folder “SimulationResults”. Here you will find a folder 
with a name corresponding to the date you ran the model. This folder includes some of the figures in the 
manuscript as well as in the Supporting Information. If you load the .mat file in this folder, values that are 
given in Table 1 and 2 of the manuscript are found in “Table”. 

## Running other scenarios
Which scenario to run is set in the following statement (line 112-113 in VertMinMod.m):
[NoFish,NoDVM,DiffusivityFactor,FishMortalityFactor,Scenario] = ...
 SetScenario("Baseline");
 
To run the 7 additional scenarios of the manuscript, replace "Baseline" in the above statement with one 
of the following options:
"HalvedDiffusivity", "TwofoldDiffusivity", "TenfoldDiffusivity",
"HalvedFishMortality", "TwofoldFishMortality", "NoFish" or "NoDVM"

There is also an option, "SetManually", where you are prompted to set values manually. The runtime 
(i.e., until steady state is reached) will vary somewhat for the different scenarios.

Results are stored as explained for the baseline scenario, but in separate folders for each scenario.
 
