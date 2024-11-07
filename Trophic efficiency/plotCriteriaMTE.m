% Creates 3 figures: 1) Global diagnostics at top 170m
%                    2) Global diagnostics depth averaged
%                    2) Global functions
%
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\clean code for figures\simHTL1.mat')

show_mNPP_winter=false;
%%
% Annual averges at top layer
    mNPP_avg_yr=squeeze(mean(sim.mNPP,1));
    mHTL_avg_yr=squeeze(mean(sim.mHTL,1));
    lambdaHTL_avg_yr=squeeze(mean(mean(sim.lambda_htl(:,:,:,1:4),4,'omitnan'),1,'omitnan'));
    % lambda_avg_yr=squeeze(mean(TLtimeTop170,1));

% load lambdaHTL matrix
% load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\27 Aug 2024 plots\lambdaHTL.mat')
mte=squeeze(double(sim.ProdHTLAnnual./sim.ProdNetAnnual));

pte=mte.^(1./(lambdaHTL_avg_yr-1));
pte(isnan(pte))=0;
%
%% ---------------------------------------
% Figures: Global Diagnostics at top 170m
%-----------------------------------------
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

cmap=flip(cmocean('matter',20));
ccmap=cmap(2:end,:);   

cmap2=flip(cmocean('matter',50));
ccmap2=cmap2(2:end,:);   

x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
heightf=12; %figure height in cm

fig=figure(11);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto','Name','TrophicLvl Diagnostics: Top 170m');
clf
set(gcf,'color','w');
tiledlayout(2,2,'TileSpacing','compact','padding','compact')

t1=nexttile(1);
c = panelGlobal(sim.x,sim.y,log(mNPP_avg_yr),[-8 1],sTitle='a. Annual average mass_{NPP}',sProjection='mollweid');
    c.Visible='off';
    colormap(t1,ccmap)
    c=colorbar('horizontal'); 
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel(c, '[\mug C]','FontSize',10)
    tickPositions = c.Ticks;
    c.TickLabels={'10^{-8}','10^{-6}','10^{-4}','10^{-2}','10^0','10^2'};
%%
    t2=nexttile(2);
    c = panelGlobal(sim.x,sim.y,log(mean(mHTL_avg_yr(:,:,1:3),3)),[-8 1],sTitle='b. Annual average mass_{HTL}',sProjection='mollweid');
    c.Visible='off';
    colormap(t2,ccmap)
    c=colorbar('horizontal');
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel(c, '[\mug C]','FontSize',10)
    tickPositions = c.Ticks;
    c.TickLabels={'10^{-8}','10^{-6}','10^{-4}','10^{-2}','10^0','10^2'};
%
t3=nexttile(3);
    c = panelGlobal(sim.x,sim.y,lambdaHTL,[min(min(lambdaHTL)) max(max(lambdaHTL))],sTitle='c. Annual average \lambda_{HTL}',sProjection='mollweid');
    c.Visible='off';
    colormap(t3,ccmap2)
    c=colorbar('horizontal');
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel(c, '[-]','FontSize',10)
   
    %
    str='$\bar{\epsilon}_{\lambda}$';

t4=nexttile(4);
cbar = panelGlobal(sim.x,sim.y,pte,[0,1],...
    sTitle='\epsilon_{\lambda}', sProjection='mollweid');
cbar.Label.String = '[-]';
cbar.Visible="off";
colormap(t4,ccmap)
cbar=colorbar('horizontal');
set(gca,'XTickLabel','');
ylabel(cbar, '[-]','FontSize',10)
cbar.Label.String = '[-]';
cbar.Label.FontSize=10;
% title(t4,append('d. Annual ',str), Interpreter="latex",FontSize=12)
% %symbol
title(t4,'d. Annual Mean Progressive Efficiency ')

%% 
if show_mNPP_winter==true
% ----------------------------------------
% Extra figure: m_NPP normalized for NPP
%-------------------------------------------
figure(100)
clf
c = panelGlobal(sim.x,sim.y,log(squeeze(sim.mNPP_norm)),[-8 1],sTitle='Normalized annual mass_{NPP}',sProjection='mollweid');
    c.Visible='off';
    colormap(t2,ccmap)
    c=colorbar('horizontal');
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel(c, 'log_{10}(\mugC)','FontSize',12)
end
 %%
 % ------------------------------------------
 % Figures: Global Diagnostics - mean depth
 %--------------------------------------------
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

cmap=flip(cmocean('matter',20));
ccmap=cmap(2:end,:);   

cmap2=flip(cmocean('matter',10));
ccmap2=cmap2(2:end,:);   

x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
heightf=12; %figure height in cm

fig=figure(12);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto','Name','TrophicLvl Diagnostics: depth avg');
clf
set(gcf,'color','w');
tiledlayout(2,2,'TileSpacing','compact','padding','compact')
mNpp_lim=log10([min(min(mean(mNPP_avg_yr,3))) max(max(mean(mNPP_avg_yr,3)))]);
t1=nexttile(1);
c = panelGlobal(sim.x,sim.y,log(mean(mNPP_avg_yr,3)),[-8 1],sTitle='Annual average mass_{NPP}',sProjection='mollweid');
    c.Visible='off';
    colormap(t1,ccmap)
    c=colorbar('horizontal'); 
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel(c, '[\mug C]','FontSize',10)
tickPositions = c.Ticks;    
c.TickLabels={'10^{-8}','10^{-6}','10^{-4}','10^{-2}','10^0','10^2'};
%
t2=nexttile();
    c = panelGlobal(sim.x,sim.y,log(mean(mHTL_avg_yr,3)),[-8 1],sTitle='Annual average mass_{HTL}',sProjection='mollweid');
    c.Visible='off';
    colormap(t2,ccmap)
    c=colorbar('horizontal');
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel(c, 'log_{10}(\mug_C)','FontSize',10)

    
t3=nexttile();
    c = panelGlobal(sim.x,sim.y,lambdaHTL,[min(min(lambdaHTL)) max(max(lambdaHTL))],sTitle='Annual average \lambda_{HTL}',sProjection='mollweid');
    c.Visible='off';
    colormap(t3,ccmap2)
    c=colorbar('horizontal');
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel(c, '[-]','FontSize',12)
    
    t4=nexttile(4);
cbar = panelGlobal(sim.x,sim.y,pte,[0,1],...
    sTitle='\epsilon_{\lambda}', sProjection='mollweid');
cbar.Label.String = '[-]';
cbar.Visible="off";
colormap(t4,ccmap)
cbar=colorbar('horizontal');
set(gca,'XTickLabel','');
ylabel(cbar, '[-]','FontSize',10)
cbar.Label.String = '[-]';
cbar.Label.FontSize=10;
title(t4,str, Interpreter="latex",FontSize=12)

% exportgraphics(gcf,[pwd '/Trophic Efficiency/27 Aug 2024 plots/globalTLdiagnostics_allDepths.pdf'])

%% 
show_GlobalFunctions=false;
% -------------------------------------
% Figure: Global Functions (vertical)
%---------------------------------------

if show_GlobalFunctions==true
cmap=flip(cmocean('matter',10));
ccmap=cmap(2:end,:);   


x0=0; %positions (no need to change)
y0=0;
width=10; %figure width in cm
heightf=12; %figure height in cm

fig=figure('Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto');
clf
set(gcf,'color','w');
tiledlayout(3,1,'TileSpacing','loose','padding','compact')


nexttile
c = panelGlobal(sim.x,sim.y,log10(sim.ProdNetAnnual),[1,3],...
    sTitle='Net primary production',sProjection='mollweid');
c.Label.String = 'log_{10}(mgC m^{-2}d^{-1})';
set(gca,'XTickLabel','');

nexttile
c = panelGlobal(sim.x,sim.y,log10(sim.ProdHTLAnnual),[1,3],...
    sTitle='HTL production', sProjection='mollweid');
c.Label.String = 'log_{10}(mgC m^{-2}d^{-1})';
colormap(ccmap)
% c=colorbar('horizontal');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
% ylabel(c, 'log_{10}(mgC m^{-2}d^{-1})','FontSize',12)

mte=squeeze(sim.ProdHTLAnnual./sim.ProdNetAnnual);
mte(isnan(mte))=0;
nexttile
cbar = panelGlobal(sim.x,sim.y,mte,[0,.5],...
    sTitle='\epsilon_{HTL}', sProjection='mollweid');
cbar.Label.String = '[-]';
colormap(ccmap)
% c=colorbar('horizontal');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'XTickLabel','')


end