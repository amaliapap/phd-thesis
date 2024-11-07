%
% Produces 2 figures: 1) Zonal Diagnostics   2) Global & Zonal Functions
% Try with mHTL
% 1. take only last 12 months
% 2. take depth avg
% mHTL is 4D: (month x lon x lat x depth)=120x128x64x15
% sim

showZonalDiagnostics=true;

mHTL_depAvg=squeeze(mean(mHTL,4,'omitnan'));

mHTL_zonal=squeeze(mean(mHTL_depAvg,2,'omitnan')); % average over longitudes

% mNPP, lambda_htl,max_lambda
% sim
mNPP_depAvg=squeeze(mean(mNPP,4,'omitnan'));

mNPP_zonal=squeeze(mean(mNPP_depAvg,2,'omitnan')); % average over longitudes

% % lambda_htl
lambda_htl_depAvg=squeeze(mean(lambda_htl(:,:,:,1:3),4,'omitnan'));
lambda_htl_zonal=squeeze(mean(lambda_htl_depAvg,2,'omitnan')); % average over longitudes

%% Zonal functions ProdNet, ProdHTL, efficiency
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\Paper\simNew.mat')
time=109:length(sim.t);
ProdNet=sim.ProdNet(time,:,:);
ProdHTL=sim.ProdHTL(time,:,:);
ProdNet_zonal=squeeze(mean(ProdNet,2,'omitnan'));
ProdHTL_zonal=squeeze(mean(ProdHTL,2,'omitnan'));
% ProdNet_zonal=squeeze(mean(sim.ProdNet,2,'omitnan'));
% ProdHTL_zonal=squeeze(mean(sim.ProdHTL,2,'omitnan'));
% Y=ProdNet_zonal;
% Y(Y==0)=1; % to avoid division by 0 below
% mte_zonal=ProdHTL_zonal./Y;
mte_zonal=ProdHTL_zonal./ProdNet_zonal;
mte_zonal(isnan(mte_zonal))=0;
% mte_zonal(mte_zonal>1)=-1;
time_vec_length=size(mte_zonal,1);

load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\27 Aug 2024 plots\lambdaHTLtime.mat')
lambdaHTLtime=lambdaHTL;
lambdaHTLtime_edit=lambdaHTLtime;
lambdaHTLtime_edit(isnan(lambdaHTLtime_edit))=1;

lambdaHTL_zonal= squeeze(mean(lambdaHTLtime,2,'omitnan'));


%% ----------------------------------- 
% Figure: Global & Zonal Functions 
%-------------------------------------
cmap=flip(cmocean('matter',80));
ccmap=cmap(2:end,:);   

width=16; %figure width in cm
heightf=14; %figure height in cm
x0=0;
y0=0;

figure_number = 2;           % Desired figure number
fig = figure(figure_number); % Create the figure with the specified number
set(fig, 'Renderer','Painters','Units','centimeters',...
    'Position',[x0 y0 width heightf],...
    'PaperPositionMode','auto','Name','Global & Zonal functions');
clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
tiledlayout(2,3,'TileSpacing','compact','padding','compact')

% The 3 bottom figures should be plotted first
% so that the colorbar position is properly arranged
time=1:12;
t4=nexttile(4);
s=surface(time,sim.y,log10(ProdNet_zonal)');
s.EdgeColor = 'none';
axis tight
colormap(ccmap)
s=colorbar('horizontal');
s.Label.String = '[mg C m^{-2}d^{-1}]';
% s.Visible="off";
clim([-5 2])
xlabel('Time (months)')
ylabel('Latitude')
ax=gca;
ax.YTick=-80:20:80;
axis square
tt4=title(t4,'d.',FontWeight='normal',HorizontalAlignment='left');
tt4.Position(1)=1.2;
tickPositions = s.Ticks;    
s.TickLabels={'10^{-4}','10^{-2}','10^0','10^2'};

t5=nexttile(5);
s=surface(time,sim.y,log10(ProdHTL_zonal)');
s.EdgeColor = 'none';
axis tight
colormap(ccmap)
s=colorbar('horizontal');
s.Label.String = '[mg C m^{-2}d^{-1}]';
% s.Visible="off";
clim([-5 2])
xlabel('Time (months)')
ax=gca;
ax.YTick=-80:20:80;
% set(gca,'YTickLabel',[]);
axis square
tt5=title(t5,'e.',FontWeight='normal',HorizontalAlignment='left');
tt5.Position(1)=1.2;
s.TickLabels={'10^{-4}','10^{-2}','10^0','10^2'};
%
t6=nexttile(6);
s=surface(time,sim.y,log10(mte_zonal)');
s.EdgeColor = 'none';
axis tight
colormap(t6,ccmap)
s=colorbar('horizontal');
s.Label.String = '[-]';
% s.Visible="off";
clim([-2 0])
xlabel('Time (months)')
ax=gca;
ax.YTick=-80:20:80;
% set(gca,'YTickLabel',[]);
axis square
tt6=title(t6,'f.',FontWeight='normal',HorizontalAlignment='left');
tt6.Position(1)=1.2;

tickPositions = s.Ticks;    
s.Ticks=[-2,-1,0];
s.TickLabels={'10^{-2}','10^{-1}','10^0'};
%
nexttile(1)
c1 = panelGlobal(sim.x,sim.y,log10(sim.ProdNetAnnual),[1,3],...
    sTitle='a. Net primary production',sProjection='mollweid');
c1.Label.String = 'log_{10}(mgC m^{-2}d^{-1})';
c1.Visible="off";
set(gca,'XTickLabel','');
c1=colorbar('horizontal');
c1.Label.String = '[mg C m^{-2}d^{-1}]';
c1.Label.FontSize=8;
pos2=c1.Position(2);
c1.Position(3)=.27;c1.Position(1)=.06;
c1.Position(2)=pos2-.03;
% c1.Label.Position(2)=2.7;
ax=gca;
set(ax, 'ColorScale', 'log');  % Set the color scale to logarithmic
% c1.Label.Position=[10,100,1000];
tickPositions = c1.Ticks;    
c1.TickLabels={'10^1','','10^2','','10^3'};

% 
t2=nexttile(2);
c = panelGlobal(sim.x,sim.y,log10(sim.ProdHTLAnnual),[1,3],...
    sTitle='b. HTL production', sProjection='mollweid');
c.Label.String = 'log_{10}(mgC m^{-2}d^{-1})';
c.Visible="off";
colormap(ccmap)
c=colorbar('horizontal');
set(gca,'XTickLabel','');
c.Label.String = '[mg C m^{-2}d^{-1}]';
c.Label.FontSize=8;
pos2=c.Position(2);
c.Position(1)=.37; c.Position(3)=.27;
c.Position(2)=pos2-.03;
tt2=title(t2,'b. HTL production');
tt2.Position(1)=0;
ax=gca;
set(ax, 'ColorScale', 'log');  % Set the color scale to logarithmic
tickPositions = c.Ticks;    
c.TickLabels={'10^1','','10^2','','10^3'};
%
mte=squeeze(sim.ProdHTLAnnual./sim.ProdNetAnnual);
mte(isnan(mte))=0;
nexttile(3)
cbar = panelGlobal(sim.x,sim.y,log10(mte),[-2,0],...
    sTitle='\epsilon_{\mu}', sProjection='mollweid');
cbar.Label.String = '[-]';
cbar.Visible="off";
colormap(ccmap)
cbar=colorbar('horizontal');
set(gca,'XTickLabel','');
ylabel(cbar, '[-]','FontSize',10)
cbar.Label.String = '[-]';
cbar.Label.FontSize=8;
% pos2=cbar.Position(2);
cbar.Position(3)=.27;cbar.Position(1)=.68;
cbar.Position(2)=c.Position(2);
% cbar.Label.Position(2)=2.7;
% remember to change this fontsize later
title('c. Microbial trophic efficiency',FontSize=10);
% c1.Label.Position=[10,100,1000];
tickPositions = cbar.Ticks;    
cbar.Ticks=[-2,-1,0];
cbar.TickLabels={'10^{-2}','10^{-1}','10^0'};

%%
% 
if showZonalDiagnostics==true
%     % lambdaHTL_zonal0=days_to_months_zonal(lambdaHTL_zonal);
%     mHTL_zonal_12=days_to_months_zonal(mHTL_zonal);
% pte_zonal=mte_zonal.^(1./(lambdaHTL_zonal-1));
% pte_zonal_edit=pte_zonal;
% pte_zonal_edit(isnan(pte_zonal_edit))=0;
% pte_zonal_edit(pte_zonal_edit==Inf)=0; 
time=1:73;
%------------------------------
%  Figure: Zonal Diagnostics
%--------------------------------
cmap=flip(cmocean('matter',80));
ccmap=cmap(2:end,:);   

width=16; %figure width in cm
heightf=8; %figure height in cm

figure_number = 1; 
fig = figure(figure_number); % Create the figure with the specified number
set(fig, 'Renderer','Painters','Units','centimeters',...
    'Position',[x0 y0 width heightf],...
    'PaperPositionMode','auto','Name','Zonal diagnostics');
clf
set(groot,'defaultAxesFontSize',10)
set(groot,'defaultAxesFontWeight','normal')
set(gcf,'color','w');

tiledlayout(1,3,'TileSpacing','tight','padding','tight')

nexttile(1)
s=surface(sim.t,sim.y,log10(mNPP_zonal)');
s.EdgeColor = 'none';
axis tight
colormap(ccmap)
cb=colorbar('horizontal');
cb.Label.String='[log_{10}(\mugC)]';
cb.Label.FontSize=8;
% clim([-8 1])
ylabel('Latitude')
xlabel('Time (days)')

title('a. Mass_{NPP}','FontWeight','normal')
ax=gca;
ax.YTick=-80:20:80;
axis square
%
nexttile(2)
s=surface(time,sim.y,log10(mHTL_zonal)');
s.EdgeColor = 'none';
axis tight
colormap(ccmap)
cb=colorbar('horizontal');
cb.Label.String='[log_{10}(\mugC)]';
cb.Label.FontSize=8;
% clim([-8 1])
title('b. Mass_{HTL.ing}','FontWeight','normal')
xlabel('Time (days)')

ax=gca;
ax.YTick=-80:20:80;
set(gca,'YTickLabel',[])
axis square

%
nexttile(3)
s=surface(time,sim.y,(lambdaHTL_zonal)');
s.EdgeColor = 'none';
axis tight
colormap(ccmap)
cb=colorbar('horizontal');
cb.Label.String='[-]';
cb.Label.FontSize=8;
% clim([1 2.5])
xlabel('Time (days)')
% ylabel('Latitude')
title('c. \lambda_{HTL}','FontWeight','normal')
ax=gca;
ax.YTick=-80:20:80;
axis square
% 
% t4=nexttile(4);
% s=surface(time,sim.y,(pte_zonal_edit)');
% s.EdgeColor = 'none';
% axis tight
% colormap(t4,ccmap)
% cb=colorbar;
% cb.Label.String='[-]';
% cb.Label.FontSize=8;
% clim([0 5])
% xlabel('Time (days)')
% title('d. $\bar{\epsilon}_\lambda$','FontWeight','normal','FontSize',12,Interpreter='latex')
% ax=gca;
% ax.YTick=-80:20:80;
% set(gca,'YTickLabel',[])
% axis square

end