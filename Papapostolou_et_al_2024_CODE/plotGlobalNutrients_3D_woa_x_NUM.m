%% 2 FIGURES
function plotGlobalNutrients_3D_woa_x_NUM(sim,showData)
%% 
addpath('woa/')

NO3 = squeeze(mean(sim.N((length(sim.t)-11):end,:,:,1),1)); % annual average of the last year       
Si = squeeze(mean(sim.Si((length(sim.t)-11),:,:,1),1));   % annual average of the last year
DOC = squeeze(mean(sim.DOC((length(sim.t)-11),:,:,1),1));   % annual average of the last year

if showData==true
% Annual WOA data
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\woa\no3_woa_3D_interp.mat')
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\woa\si_woa_3D_interp.mat')
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\woa\po4_woa_3D_interp.mat')

% Convertion from umol/kg to ug/l. Top layer in WOA is top 5m
no3_woa_surf=14*squeeze(no3_woa_3D_interp(:,:,1));% ugN/l
si_woa_surf=28.09*squeeze(si_woa_3D_interp(:,:,1));% ugN/l
po4_woa_surf=31*squeeze(po4_woa_3D_interp(:,:,1));% ugN/l
end
cmaxN= round(max(NO3,[],'all'));
cmaxSi= round(max(Si,[],'all'));
%% ---------------------------
%        Figures - surface
%---------------------------
sProjection='mollweid';

cmap=flip(cmocean('deep',100));
ccmap=cmap(2:end,:);      

% for the ratios
% cmap2=flip(cmocean('deep',10));
% ccmap2=cmap2(2:end,:); 

x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
heightf=10; %figure height in cm

fig=figure(13);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto','Name','WOAxNUM_N_Si');

clf
set(gcf,'color','w');
if showData==false
tiledlayout(1,2,'TileSpacing','compact','padding','tight',TileIndexing='columnmajor')
% NO3 woa

Nlim=[-3 3];
Silim=[-2 4];

t1=nexttile(1);
cn = panelGlobal(sim.x,sim.y,log10(NO3),[-2 4],sTitle='a. Modeled N ',sProjection=sProjection);
    cn.Visible='off';
    colormap(t1,ccmap);
    cn=colorbar; 
    clim(Nlim)
    set(gca,'XTickLabel','')
    cn.Ticks=[-3 -1 1 3];
    cn.TickLabels={'10^{-3}','10^{-1}','10^1','10^3'};
    ylabel(cn, '\mug N l^{-1}','FontSize',10)
    % cn.Position(1)=.282;

% Silicate NUM
t2=nexttile(2);
c = panelGlobal(sim.x,sim.y,log10(Si),Silim,sTitle='b. Modeled Si',sProjection=sProjection);
    c.Label.String  = 'log_{10}(\mug Si l^{-1})';
    c.Visible='off';
    colormap(t2,ccmap)
    c=colorbar; 
    clim(Silim)
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
     c.Ticks=[-2 0 2 4];
    c.TickLabels={'10^{-2}','10^0','10^2','10^4'};
    ylabel(c, '\mug Si l^{-1}','FontSize',10)
else

tiledlayout(2,2,'TileSpacing','compact','padding','tight',TileIndexing='columnmajor')
% NO3 woa

Nlim=[-3 3];
Silim=[-2 4];
t1(1)=nexttile(1);
c = panelGlobal(sim.x,sim.y,log10(no3_woa_surf)',[-2 4],sTitle='a. WOA NO_3',sProjection=sProjection);
    c.Visible='off';
    colormap(t1,ccmap)
    c=colorbar('horizontal');
    clim(Nlim)
    c.Visible='off';
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    % ylabel(c, 'log_{10}(\mug N l^{-1})','FontSize',10)
% Nitrogen NUM
t3=nexttile(3);
cn = panelGlobal(sim.x,sim.y,log10(NO3),[-2 4],sTitle='b. Modeled N ',sProjection=sProjection);
    cn.Visible='off';
    colormap(t3,ccmap);
    cn=colorbar; 
    clim(Nlim)
    set(gca,'XTickLabel','')
    cn.Ticks=[-3 -1 1 3];
    cn.TickLabels={'10^{-3}','10^{-1}','10^1','10^3'};
    ylabel(cn, '\mug N l^{-1}','FontSize',10)
    % cn.Position(1)=.282;

    % Si woa
t2=nexttile(2);
c = panelGlobal(sim.x,sim.y,log10(si_woa_surf)',Silim,sTitle='c. WOA SiO_4',sProjection=sProjection);
    c.Visible='off';
    colormap(t2,ccmap)
    c=colorbar('horizontal'); 
    clim(Silim)
    c.Visible='off';
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    % ylabel(c, 'log_{10}(\mug Si l^{-1})','FontSize',10)

% Silicate NUM
t4=nexttile(4);
c = panelGlobal(sim.x,sim.y,log10(Si),Silim,sTitle='d. Modeled Si',sProjection=sProjection);
    c.Label.String  = 'log_{10}(\mug Si l^{-1})';
    c.Visible='off';
    colormap(t4,ccmap)
    c=colorbar; 
    clim(Silim)
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
     c.Ticks=[-2 0 2 4];
    c.TickLabels={'10^{-2}','10^0','10^2','10^4'};
    ylabel(c, '\mug Si l^{-1}','FontSize',10)
    % c.Position(1)=.282;
    % c.Position(2)=.1;
end
%% --------------------------------------
%  Figures: woa PO4 & NUM DOC (Appendix) 
%----------------------------------------
fig=figure(51);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto','Name','WOA_po4xNUM_doc');

clf
set(gcf,'color','w');
if showData==true
tiledlayout(1,2,'TileSpacing','compact','padding','compact',TileIndexing='columnmajor')
% PO4 woa
t1=nexttile();
c = panelGlobal(sim.x,sim.y,log10(po4_woa_surf)',[-2 4],sTitle='a. WOA PO4',sProjection=sProjection);
    c.Visible='off';
    colormap(t1,ccmap)
    c=colorbar('horizontal'); 
    clim([-1 3])
    ylabel(c, 'log_{10}(\mug_P l^{-1})','FontSize',10)
    % set(colorbar,'visible','off')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel(c, 'log_{10}(\mug P l^{-1})','FontSize',10)

t2=nexttile();
end
c = panelGlobal(sim.x,sim.y,log10(DOC),[-2 4],sTitle='b. Modeled Surface DOC',sProjection=sProjection);
    c.Label.String  = 'log_{10}(\mug l^{-1})';
    c.Visible='off';
    colormap(t2,ccmap)
    c=colorbar('horizontal'); 
    clim([-1 3])
    ylabel(c, 'log_{10}(\mug_C l^{-1})','FontSize',10)
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);    
    ylabel(c, 'log_{10}(\mug C l^{-1})','FontSize',10)

% exportgraphics(gcf,[append('NutrientsNUMxWOA','_MAPS.pdf')])

