% Paper 1 figures
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Papapostolou_et_al_2024_CODE')
thisDir=pwd;
load('simNUMmodel.mat')

% Change directory to '.../NUMmodel/matlab' to execute the code 
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

p = setupNUMmodel(bParallel=true);
p = parametersGlobal(p);
sim=simF;
sim.p=p;
% Replace with the path of Transport Matrices (TMs) after downloading
% the library
sim.p.pathGrid="C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\TMs\MITgcm_2.8deg\grid.mat";

%%
cd(thisDir) % switch to the folder conatining the figure scripts
% plotSensitivityStagesGroups % Figure 5
[simWC,simC]=MyNewPLot_update(sim,1,true); % Figure 6, Figure 14
nppComparison(sim,false)   % Figure 7
ModelVsData(sim)           % Figure 8
AMTnanomicroFigures(sim)   % Figure 9
plotMesoZP_NUMxMoriarty_top150int(sim,false) % Figure 10
plotSpectrumAllRatesUpdate(simC)  % Figures 15,16: Spectra in Spring Summer and Unicellular respiration
cd(thisDir)
simPhyto=load("simPhyto.mat");
% if simPhyto is not calculated, use sim 
plotGlobalPhytoplankton(simPhyto.sim)      % Figure 11, 12
plotGlobalNutrients_3D_woa_x_NUM(sim,true) % Figure 13 & E1 (Appendix)

%% Appendix figures

sensitivityCalibrationParams(true) % sensitivity plots for the 3 calibrated params
MyNewPLot_update(sim,2,false)  % Upwelling water-column
MyNewPLot_update(sim,3,false)  % Oligotrophic water-column
PlotMicroMesoCopepods(sim) 
r = calcRadiusGroups(sim.p);
GOPOPCORNfigures(sim,r) % transect map
PicophytoFiguresAtlantic(sim,r) % shows only data transects in the Atlantic
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
plotWatercolumnTime(sim, 60,-15) % or simWC