% calculate_mNPPmHTL;

idxTime=109:120;
idxDep=1:3;
sim.mNPP=mNPP(idxTime,:,:,1);
sim.mHTL=mHTL(idxTime,:,:,1);
pn=squeeze(sim.ProdNetAnnual);%squeeze(mean(sim.ProdNet,1));
ph=squeeze(sim.ProdHTLAnnual);
% mnpp=squeeze(mean(sim.mNPP,1));
% mhtl=squeeze(mean(sim.mHTL,1));

mnpp=squeeze(mean(sim.mNPP,1,'omitnan'));
mhtl=squeeze(mean(sim.mHTL,1,'omitnan'));
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency\27 Aug 2024 plots\lambdaHTL.mat')

sim.lambda_htl=lambdaHTL;
% bad calculation! do it again!
lam_htl=squeeze(mean(sim.lambda_htl,[3 4]));

mte=ph./pn;
% mte=sim.ProdHTLAnnual./sim.ProdNetAnnual
%%
[lon_interpW,lat_interpW]=meshgrid(sim.x,sim.y);

lon=lon_interpW';
lat=lat_interpW';
%%
pnpp=pn(:);
phtl=ph(:);
lat_vec=lat(:);
mte_vec=mte(:);
mnpp_vec=mnpp(:);
mhtl_vec=mhtl(:);
lam_htl_vec=lam_htl(:);
%% ----------------------------------- 
%     Figure: MTE correlations 
%-------------------------------------
cmap=(cmocean('curl',80));
% ccmap=cmap(2:end,:);   

width=15; %figure width in cm
heightf=9; %figure height in cm
x0=0;
y0=0;

figure_number = 10;           % Desired figure number
fig = figure(figure_number); % Create the figure with the specified number
set(fig, 'Renderer','Painters','Units','centimeters',...
    'Position',[x0 y0 width heightf],...
    'PaperPositionMode','auto','Name','MTE correlations');
clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
tiledlayout(2,2,'TileSpacing','compact','padding','loose')

% same for mHTL, mNPP, lHTL
nexttile
scatter(pnpp,mte_vec,5,lat_vec,"filled")
xlabel('NPP (mg C m^{-2}d^{-1})')
ylabel('\epsilon_{\mu} (-)',FontSize=12)
corval=corr([pnpp,mte_vec],'Rows','complete','type','Spearman');
text(0.72,0.95,num2str(corval(1,2),'%.2f'),'Units','normalized','fontsize',10,'FontWeight','bold')
set(gca,'XScale','log')
plotlabel('a',false)
%axis square

nexttile
scatter(mnpp_vec,mte_vec,5,lat_vec,"filled")
xlabel('m_{NPP} (\mug C)')
corval=corr([mnpp_vec,mte_vec],'Rows','complete','type','Spearman');
text(0.72,0.95,num2str(corval(1,2),'%.2f'),'Units','normalized','fontsize',10,'FontWeight','bold')
set(gca,'XScale','log')
plotlabel('b',false)
% xlim([0 .1])
%axis square

nexttile
scatter(mhtl_vec,mte_vec,5,lat_vec,"filled")
xlabel('m_{HTL} (\mug C)')
ylabel('\epsilon_{\mu} (-)',FontSize=12)
corval=corr([mhtl_vec,mte_vec],'Rows','complete','type','Spearman');
text(0.72,0.95,num2str(corval(1,2),'%.2f'),'Units','normalized','fontsize',10,'FontWeight','bold')
set(gca,'XScale','log')
plotlabel('c',false)

%axis square

nexttile
scatter(lam_htl_vec,mte_vec,5,lat_vec,"filled")
xlabel('TL_{HTL} (-)')
corval=corr([lam_htl_vec,mte_vec],'Rows','complete','type','Spearman');
text(0.72,0.95,num2str(corval(1,2),'%.2f'),'Units','normalized','fontsize',10,'FontWeight','bold')
set(gca,'XScale','log')
plotlabel('d',false)

cbar=colorbar;
ylabel(cbar,'Latitude')
colormap(cmap)
%axis square
% cbar.Position= [.91 .2 .022 .5];
cbar.Position=[0.916,0.177951388888889,0.02265625,0.717222222222222];
%%
