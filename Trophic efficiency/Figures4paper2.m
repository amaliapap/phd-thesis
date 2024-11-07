% 1. load matrix

load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\27 Aug 2024 plots\sim27aug.mat')
% sim=simF;
% load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\sim_mHTL1ugC.mat')

% Change directory to '.../NUMmodel/matlab' to execute the code 
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

p = setupNUMmodel(bParallel=true);
setHTL(0.005,.1,true,false);
p = parametersGlobal(p);
% sim=simF;
sim.p=p;
% Replace with the path of Transport Matrices (TMs) after downloading
% the library
sim.p.pathGrid="C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\TMs\MITgcm_2.8deg\grid.mat";
%% 3. Calculate global m_NPP,m_HTL and trophic levels with:
%   calcTrophicEfficiencyDiagnostics.m
%   This is an expensive routine, so run it on Linux server
% sim.p.pathGrid="C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\TMs\MITgcm_2.8deg\grid.mat";
% load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\Paper\calibrated sim\mNPP.mat')
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\Paper\calibrated sim\TLtimeTop170.mat')
% load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\Paper\calibrated sim\mHTL.mat')
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\Paper\mNPP4D.mat')
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\Paper\mHTL4D.mat')

% CalcTrophicLevelGlobal;
% calcMassNPPdz;

%% SeaZone-ality Functions
cd 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency'

SeaZoneality
% progrEffSeazonality

% PLot Global Diagnostics
plotsForTrophicLevels(sim,mte);
%% plotGlobalZonalMulticellular
showTimeSeries=false;
if showTimeSeries==true
    % Plot time-series global - not used in Paper
    calcGlobalEfficiencyTime;
end
%% Trophic network plot
% cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency')
% EflowplotNUMhtl;%(sim) in presentation folder
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\27 Aug 2024 plots\sim27aug.mat')
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\Paper\calibrated sim\mHTL.mat')

% sim=simF;
% Change directory to '.../NUMmodel/matlab' to execute the code 
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

p = setupNUMmodel(bParallel=true);
p = parametersGlobal(p);
% sim=simF;
sim.p=p;
% Replace with the path of Transport Matrices (TMs) after downloading
% the library
sim.p.pathGrid="C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\TMs\MITgcm_2.8deg\grid.mat";
%% Water column and Trophic Network
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
TrophicEfficiency_WC_Figures;
% GlobalFlowsMasterPlot;
correlationMTEwithMetrics
%% sensitivities
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

sensitivityHTL_Efficiency
% SensitivityWC
% sensitivityChemostat

%% Pending - Under consideration
PlotMicroMesoCopepods 
plotGlobalPhytoplankton
CalcTrophicLevelGlobal
MicrobialFoodWeb
% add micro-meso fraction in plotGlobalPhytoblankton


%%   Chemostat figures
p=setupNUMmodel(mAdultPassive, mAdultActive, n,nCopepods,nPOM);
p = parametersChemostat(p); % Use standard low-res model
p.tEnd = 10*365;
p.d=0.1;
L=100;
%% Eutrophic
simEut=simulateChemostat(p,L);
% TrophicEfficiencyPlots(simEut);
ChemostatTrophicEfficiencyPlots(simEut)

%% Oligotrophic
p.d=0.001;
simOligo=simulateChemostat(p,L);
% TrophicEfficiencyPlots(simEut);
ChemostatTrophicEfficiencyPlots(simOligo)
%% Seasonal
p_bis = parametersChemostat(p, 'seasonalAmplitude', 1);
p_bis.tEnd = 10*365;
simSeasonal = simulateChemostat(p_bis,'bUnicellularloss', false);
plotSimulation(simSeasonal)
SeasonalChemostatTrophicEfficiencyPlots(simSeasonal)