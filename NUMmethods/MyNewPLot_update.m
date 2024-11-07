% Plots
% Global Biomass of Generalists, Diatoms and Copepods
% Watercolumns      : a. extracted from Global simulation
%                     b. calculated only in the water column
% Community spectrum: a. at the top layer of the water column extarcted
%                        from the global simulation
%                     b. at the top layer of the simulated water column
%                     c. in a Chemostat simulation
%
%   INPUTS: latitude (lat), longitude (lon), sim
%   OUTPUTS: 2 tiled figures
%
function [simWC,simC]=MyNewPLot_update(sim,site,showGlobal)
this_Directory='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Papapostolou_et_al_2024_CODE';
latitudes=[60,-5,24];
longitudes=[-15,5,-158];
lat = latitudes(site);
lon = longitudes(site);

% color palette for groups
newPalette={[21, 16, 240]/256,[0, 163, 136]/256,[240, 154, 15]/256,[240, 154, 15]/256,...
    [219, 6, 6]/256,[219, 6, 6]/256,[219, 6, 6]/256,[219, 6, 6]/256,[116, 71, 145]/256};
%% --------------------------
%     Global simulation
%----------------------------
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
% p = setupNUMmodel(bParallel=true);
% p = parametersGlobal(sim.p);
cd(this_Directory)
% sim.p=p;
cmap=flip(cmocean('deep',100));
ccmap=cmap(2:end,:); 
ylimSpectrum = [0.01 50];

%%
if showGlobal==true
% ----------------------------
%   Figure: Global biomasses
%------------------------------
% cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Code for FIGURES')
%
sim.p.colGroup=newPalette;

sProjection='mollweid';
%-----------------------
% figure specifications
%-----------------------
x0=0; 
y0=0;
width=8; %figure width in cm
height=15; %figure height in cm

fig=figure(6);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Global Biomass');

clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact')
%
% Global:
%

% Generalists and diatoms:
for iGroup = 1:2
    nexttile%(1+3*(iGroup-1))
    ix = (sim.p.ixStart(iGroup):sim.p.ixEnd(iGroup)) -sim.p.idxB+1;
    B = calcIntegrateGlobal(sim, sim.B(:,:,:,:,ix), true);
    % Plot the group:
     if iGroup==1
        ttl='a. Generalists';
    else
        ttl='b. Diatoms';
     end
    sTitle.Fontweight = 'normal';
    cbar = panelGlobal_point(sim.x,sim.y, log10(B),[-1 2], ...
        sTitle = ttl,...
        sProjection=sProjection);
   cbar.Label.String  = 'log_{10} (g C m^{-2})';
   colormap(ccmap);
   set(colorbar,'visible','off')
   set(gca,'YTickLabel',[]);
   set(gca,'XTickLabel',[]);
   sTitle.Units = 'Normalize'; 
   sTitle.Position(1) =0; % use negative values (ie, -0.1) to move further left ttl.HorizontalAlignment = 'left';
   sTitle.HorizontalAlignment = 'left';
end
% Copepods:
nexttile
B = 0*B;
for iGroup = 3:6
    ix = (sim.p.ixStart(iGroup):sim.p.ixEnd(iGroup)) -sim.p.idxB+1;
    B = B + calcIntegrateGlobal(sim, sim.B(:,:,:,:,ix), true);
end
    cbar = panelGlobal_point(sim.x,sim.y, log10(B),[-1 2], ...
    sTitle = 'c. Copepods',sProjection=sProjection);
    set(colorbar,'visible','off')
    cbar=colorbar('horizontal'); 
    colormap(ccmap);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    sTitle.Fontweight = 'normal';
    sTitle.Units = 'Normalize'; 
    cbar.Ticks=[-1 0 1 2];
    cbar.TickLabels={'10^{-1}','10^0','10^1','10^2'};
    ylabel(cbar, 'Biomass (g C m^{-2})','FontSize',10)
%----------------------------------------------
%     END OF Global Biomass Figure
%----------------------------------------------
end
%% ---------------------------------------------
% Figure: Watercolumn and spectrum from global:
%----------------------------------------------

x0=0; 
y0=0;
height=15; %figure height in cm
width=16;

fig=figure(14);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Watercolumn and spectra');

clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',12)
tiledlayout(3,5,'TileSpacing','tight','Padding','compact')
t1=nexttile([1 2]);
sim.p.colGroup=newPalette;
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
idx = calcGlobalWatercolumn(lat,lon,sim);
cd(this_Directory)
plotWC(squeeze(sim.B(:,idx.x, idx.y, idx.z,:)),sim,sim.z(idx.z));
set(gca,'XTickLabel','');
colormap(ccmap);
set(colorbar,'visible','off')
xlabel('')

maxDepth=max(sim.z(idx.z));
ax1=gca;
ax1.XTick=[15,115,215,315];
ax1.XTickLabel={'15','115','215','315'};
axis square
alabel=plotlabel('a',false);
% alabel.Position(2)=1.13;

%
%   Tile spectrum
%.....................

cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
t2=nexttile([1 3]);
box on
 iDepth=1; iTime=3400; 
 set(gca,'XTickLabel','');
plotPanelSpectrum(sim,iTime,iDepth,lat,lon);
xlabel('')
ylim(ylimSpectrum)
legend('')
plotlabel('b',false);

%% ------------------------- 
%  Watercolumn simulation
%---------------------------
%
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
p = setupNUMmodel(bParallel=true);
p = parametersWatercolumn(p,2);
p.tEnd = 10*365;
simWC = simulateWatercolumn(p, lat,lon);
cd(this_Directory)
simWC.p.colGroup=newPalette;

%% -------------------------------
%  Figure (c.tile): Watercolumn:
%---------------------------------
t3=nexttile([1 2]);

plotWC(simWC.B, simWC, simWC.z);
axis square
ax2=gca;
ax2.XTick=[3300,3400,3500,3600];
ax2.XTickLabel={'15','115','215','315'};
cbar=colorbar;
colormap(ccmap);
cbar.Location="south";
cbar.Position=[0.102258950225578,0.297015610371739,0.206948212633974,0.026446280991736];
ylabel(cbar, '\mug C l^{-1}','FontSize',10,'Color',[1,1,1],FontWeight='bold')
cbar.Label.Position=[0.9,1.2,0];
cbar.TickLabels={'10^{-2}','10^0','10^2'};
plotlabel('c',false);

%% ------------------------------------
% Figure(d.tile): WC spectrum at iDeph
%--------------------------------------
%
t4=nexttile([1 3]);
% cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
box on
s = simWC;
s.t = 1;
s.B = mean(squeeze(s.B(:,iDepth,:)),1);

set(gca,'XTickLabel','');
% xlabel('')
legend('')
plotPanelSpectrum(simWC,iTime,iDepth,lat,lon);
ylim(ylimSpectrum)
xlabel('')
legend('')
plotlabel('d',false);

%% -------------------------- 
%  Chemostat simulation
%----------------------------
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
p = setupNUMmodel(bParallel=true);
p = parametersChemostat(p,lat_lon=[lat,lon]);
p.tEnd = 10*365;
simC = simulateChemostat(p, bUnicellularloss=false);
simC.p.colGroup=newPalette;
cd(this_Directory)
%% ---------------------------------- 
%  Figure(e.tile): Chemostat spectrum
%------------------------------------
t5=nexttile([1 2]);
t6=nexttile([1 3]);

box on
s = simC;
s.t = 1;
s.B = mean(s.B,1);
plotPanelSpectrum(simC,iTime,1,lat,lon)
ylim(ylimSpectrum)
xlabel('Mass (\mug C)')
plotlabel('e',false);

ax=gca;
ax.XTick=[10^-8 10^-6 10^-4 10^-2 1 100];
%
hleg=findobj(gcf,'Type','Legend');
hleg(1).NumColumns=1;
hleg(1).Position=[0.04,0.05,0.34,0.2];
hleg(1).ItemTokenSize(1)=15;
cbar.Position(2)=.3;
cbar.Position(1)=.09;
cbar.Position(4)=.03;
cbar.Position(3)=.22;

delete(t5);
