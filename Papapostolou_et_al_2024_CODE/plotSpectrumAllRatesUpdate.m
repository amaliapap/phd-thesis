%λ
%
% Plot a size spectrum at a given time (day).
% If the simulation is a watercolumn then indicate also the depth layer.
% If the simulation is global then indicate depth layer, and latitude and longitude.
% 
function plotSpectrumAllRatesUpdate(sim)%(sim, time, iDepth, lat, lon)
arguments
    sim struct;
    % time = sim.p.tEnd;
    % iDepth {mustBeInteger} = 1;
    % lat double = [];
    % lon double = [];
end
% itime=3400; iDepth=4;
% sim=simC;
p=sim.p;
time=3400; %/3500 for July
thisDir=pwd;
% m =[sim.p.mLower(3:end), sim.p.mLower(end)+sim.p.mDelta(end)];
% [~, iTime] = min(abs(sim.t-time));
 iDepth=1;  
% sim = simC;
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
Rates = calcAllRatesAvg(sim,time,iDepth);
rates = Rates;
time=3500; % for July
[~, iTime] = min(abs(sim.t-time));
ratesSummer = calcAllRatesAvg(sim,time);
% cd 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Code for FIGURES'
cd(thisDir)

%% FIGURES
%
%
x0=0; %positions (no need to change)
y0=0;
height=14; %figure height in cm
width=16;
fig=figure(15);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Spectrum rates');

clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
tiledlayout(9,2,'tilespacing','compact','padding','loose','TileIndexing', 'columnmajor')
lettersStr={'a','b','c','c','d','d','d'};
seasonTiles(rates,p,lettersStr,'off')
lettersStr={'e','f','g','g','h','h','h'};
seasonTiles(ratesSummer,p,lettersStr,'on')
dim = [.2 .5 .3 .3];
str = 'Rates (day^{-1})';
a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.LineStyle="none";
a.Rotation=90;
a.Position=[0.087603305785124,0.508994710045932,0.168595036790391,0.052910051858825];
%% figure

% cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Code for FIGURES')
cd(thisDir)

x0=0; %positions (no need to change)
y0=0;
height=11; %figure height in cm
width=14;
fig_no=16;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Respiration');

clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',12)
tiledlayout(2,5,'TileSpacing','tight','Padding','compact')

t1=nexttile([1 2]);
iGroup=1; % generalists
panelRespirationA(sim.p,rates,iGroup,true)
set(gca,'XTickLabel','')
xlabel('')
ylabel('')
ylim([0 1.41])
axis square

t2=nexttile([1 2]);
iGroup=2; % diatoms
panelRespirationA(sim.p,rates,iGroup,true)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('')
ylabel('')
ylim([0 1.41])
axis square

% July panels
t3=nexttile([1 2]);
iGroup=1; % generalists
panelRespirationA(sim.p,ratesSummer,iGroup,true)
ylim([0 .65])
axis square
% Change text and position of ylabel
set(get(gca,'YLabel'),'String','Respiration or synthesis rates (day^{-1})',Position=[1e-10,0.68,-1])

t4=nexttile([1 2]);
iGroup=2; % diatoms
panelRespirationA(sim.p,ratesSummer,iGroup,true)
set(gca,'YTickLabel','')
ylabel('')
ylim([0 .65])
axis square

title(t1,'Generalists',"FontWeight","normal")
title(t2,'Diatoms',"FontWeight","normal")
title(t3,' ',"FontWeight","normal")
title(t4,' ',"FontWeight","normal")

legend.Location='south';
%
hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend,'Visible','off')
set(hLegend(1),'Visible','on')
hLegend(1).Position=[0.75,    0.7069,    0.1884,    0.2196];
