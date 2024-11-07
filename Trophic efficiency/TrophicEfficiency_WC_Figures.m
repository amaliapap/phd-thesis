% Creates plots for a water column extracted from the global simulation
%
% figures:
% 1. Production time-series & Water column depth,time,biomass and spectrum
% 
% functions called: mNPPmHTLwc, extractWCfromGlobal, plotWC, plotPanelSpectrum 
%                   DiagnosticTilesTE, TrophicNetworksWCAnnual_Tiles
%
latitudes  = [60,-5,24];
longitudes = [-15,5,-158];
%--------------------------
%        COLORS
%--------------------------
% palette for group colors
newPalette={[21, 16, 240]/256,[0, 163, 136]/256,[240, 154, 15]/256,[240, 154, 15]/256,...
    [219, 6, 6]/256,[219, 6, 6]/256,[219, 6, 6]/256,[219, 6, 6]/256,[116, 71, 145]/256};
sim.p.colGroup = newPalette;
% colormaps
cmap  = flip(cmocean('deep',5));
ccmap = cmap(2:end,:); 
cmap2  = flip(cmocean('deep',100));
ccmap2=cmap2(2:end,:); 
%% FIGURE 1: 3 Water columns & Diagnostics
%
[~,ixMonth]=find(sim.t==105); % April 15
% column 1
site=1;
lat=latitudes(site); lon=longitudes(site);
extractWCfromGlobal % for this first we need to set p (see Figure4paper2)
simWC.z=idx.z;
[mNPPwc,mHTLwc] = mNPPmHTLwc(simWC,ixMonth);

%-----------------------
% figure specifications
%-----------------------
x0=0; % positions (no need to change)
y0=0;
width=16;  % figure width in cm
height=15; % figure height in cm

figure_number = 11; % Desired figure number
fig = figure(figure_number); % Create the figure with the specified number
set(fig, 'Renderer','Painters','Units','centimeters',...
    'Position',[x0 y0 width height],...
    'PaperPositionMode','auto','Name','Production time-series WC');
clf(figure_number)
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
set(groot,'defaultLineLineWidth',2)
set(gcf,'color','w');

time=1:length(sim.t);

tiledlayout(3,3,'TileSpacing','tight','Padding','loose','TileIndexing','columnmajor')
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Code for FIGURES')
%
t1=nexttile( );
plotWC(squeeze(sim.B(:,idx.x, idx.y, idx.z,:)),sim,sim.z(idx.z));
set(gca,'XTickLabel','');
colormap(ccmap2);
set(colorbar,'visible','off')
% xlabel('Time')
xlabel('')
maxDepth=max(sim.z(idx.z));
XTickLabel={num2str(mod(3300,365)),num2str(mod(3400,365)),num2str(mod(3500,365)),num2str(mod(3600,365))}';
ylim([-maxDepth 0]);%same depth as global wc
t1.TitleHorizontalAlignment = 'left';
axis square
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
time=simWC.t;
t2=nexttile( );
    plot(time,simWC.ProdNetwc,'Color',cmap(4,:)) % extracted from WC
    hold on
    ylabel('Production (\mug_Cl^{-1}day^{-1})')
    yyaxis  right
    plot(time,ProdHTLwcNew,'Color',cmap(2,:)) % mannually calculated
    xlabel('Time (days)')
    axis tight
    xticks=[100,200,300];

    leg1=legend('NPP','Prod_{HTL}','box','off',Location='best');
    leg1.ItemTokenSize(1)=10;
    ax=gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = cmap(2,:);
    MTEannual=mean(ProdHTLwcNew,1)/mean(simWC.ProdNetwc,1);
axis square

t3=nexttile( );
box on
 iDepth=1; iTime=ixMonth; timePrint=iTime-9*12;
 % [~, iTime] = min(abs(sim.t-time));
 set(gca,'XTickLabel','');
plotPanelSpectrum(sim,iTime,iDepth,lat,lon);
tt3=title(t3,'g. ',FontWeight='normal');
% tt3.Position(1)=0;
% ylim(ylimSpectrum)
xlabel('Mass (\mug C)')
ax=gca;
ax.XTick=[10^-6 10^-2 100];
ax.XTickLabel={ '10^{-6}' '10^{-2}' '10^2'};
ax.YTick=[10^-4 10^-2 1];
ax.YTickLabel={ '10^{-4}' '10^{-2}' '10^0'};

lgd=legend;
lgd.NumColumns=4;
lgd.Orientation="horizontal";
lgd.ItemTokenSize(1)=8;
str =num2str(round(MTEannual,3));
title_list={'a. Seasonally stratified','b. Upwelling','c. Oligotrophic'};
tt1=title(t1,string(title_list(site)),"FontWeight","normal");
title(t2,append('d. $\bar{\epsilon}_{\mu}$=',str),"FontWeight","normal",HorizontalAlignment="center",Interpreter="latex")
axis square
%
% column 2
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

site=2;
lat =latitudes(site);
lon = longitudes(site);
extractWCfromGlobal;
time=sim.t;
DiagnosticTilesTE(sim,idx,lat,lon,site,time,simWC.ProdNetwc,ProdHTLwcNew,ixMonth) 

% column 3
site=3;
lat =latitudes(site);
lon = longitudes(site);
extractWCfromGlobal;
time=sim.t;
DiagnosticTilesTE(sim,idx,lat,lon,site,time,simWC.ProdNetwc,ProdHTLwcNew,ixMonth) 

cbar=colorbar;
cbar.Location="eastoutside";
cbar.Position=[0.9135,0.735,0.0264,0.1726];
ylabel(cbar, 'log_{10}(\mug C l^{-1})','FontSize',10,'Color',[0,0,0],Position=[1.8,.5,0])
lgd.Position=[0.0710,0.0030,0.8564,0.0505];

%% FIGURE 2: Trophic Networks
%

%-----------------------
% figure specifications
%-----------------------
x0=0; %positions (no need to change)
y0=0;
width=10;          %figure width in cm
height=17;         %figure height in cm
figure_number = 5; 
fig = figure(figure_number); % Create the figure with the specified number
set(fig, 'Renderer','Painters','Units','centimeters',...
    'Position',[x0 y0 width height],...
    'PaperPositionMode','auto','Name','Energy flows');
clf(figure_number)
tiledlayout(3,4,'TileSpacing','tight','Padding','loose')

% column 1
site=1;
lat =latitudes(site);
lon = longitudes(site);
extractWCfromGlobal;
time=ixMonth;
% EflowplotNUMhtlWC_columnTiles(simWC,lat,lon,site,time) 
TrophicNetworksWCAnnual_Tiles(simWC,site) % annual average 
simWC_site1=simWC;
%
%
% column 2
site=2;
lat =latitudes(site);
lon = longitudes(site);
extractWCfromGlobal;
time=ixMonth;
TrophicNetworksWCAnnual_Tiles(simWC,site) % annual average 

simWC_site2=simWC;

%
% column 3
site=3;
lat =latitudes(site);
lon = longitudes(site);
extractWCfromGlobal;
time=ixMonth;
TrophicNetworksWCAnnual_Tiles(simWC,site) % annual average 
simWC_site3=simWC;


