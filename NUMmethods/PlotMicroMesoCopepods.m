%λ
%
% Plot global biomass distribution of micro- and meso-copepods
% 
% Take the average of the top 170m
%
%
function PlotMicroMesoCopepods(sim) 

p=sim.p;
ixC = (p.ixStart(p.typeGroups==10)-p.idxB+1):(p.ixEnd(p.typeGroups==11)-p.idxB+1);
rGroups=calcRadiusGroups(sim.p); % for copepods it returns length/2

r_microC=find(rGroups(ixC)>=20/2 & rGroups(ixC)<=200/2);% radius in um divided by two to get the radius between (20-200um)/2
r_mesoC=find(rGroups(ixC)>200/2 & rGroups(ixC)<=2e4/2);
r_macroC=find(rGroups(ixC>2e4/2));
Bmicro_C=squeeze(sum(sim.B(:,:,:,:,r_microC),5));
Bmeso_C=squeeze(sum(sim.B(:,:,:,:,r_mesoC),5));
Bmacro_C=squeeze(sum(sim.B(:,:,:,:,r_macroC),5));
Bmulti_tot=squeeze(sum(sim.B(:,:,:,:,ixC),5));
% Take depth-integrated in the top 170m
field = squeeze(sum(Bmeso_C(:,:,:,1:3),5));
dz = sim.dznom(1:3);
Bmeso_intTop170 =double(squeeze( sum(field.*reshape(dz ,1,1,1,numel(dz)),4) / 1000)); % g/m2

% Take top 170m
Bmicro_Csurf=sum(Bmicro_C(:,:,:,1:3),4);
Bmeso_Csurf=sum(Bmeso_C(:,:,:,1:3),4);
Bmacro_Csurf=sum(Bmacro_C(:,:,:,1:3),4);


% Take annual mean of the last year at top 170m (iDepth=1:3)
itime=size(sim.B,1);
Bmicro_CavgSurf=squeeze(mean(Bmicro_Csurf(itime-11:end,:,:),1));
Bmeso_CavgSurf=squeeze(mean(Bmeso_Csurf(itime-11:end,:,:),1));
Bmacro_CavgSurf=squeeze(mean(Bmacro_Csurf(itime-11:end,:,:),1));

%%
options.sProjection='mollweid';
cmap=flip(cmocean('deep',100));
ccmap=cmap(2:end,:);      

x0=0; %positions (no need to change)
y0=0;
width=15; %figure width in cm
heightf=14; %figure height in cm
fig=figure('Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto');

clf
tiledlayout(1,2,'TileSpacing','compact','padding','compact')
set(gcf,'color','w');

%
nexttile(1)
cbar=panelGlobal(sim.x, sim.y,log10(Bmicro_CavgSurf),linspace(0,2.5,20)  , sTitle="a. Micro-copepods", sUnits="mg_Cm^{-3}", sProjection=options.sProjection);
% cbar=panelGlobal(sim.x, sim.y,(Bmicro_CavgSurf),[0.5 100]  , sTitle="a. Micro-copepods", sUnits="mg_Cm^{-3}", sProjection=options.sProjection);
cbar.Label.String = '(log_{10}(mgC m^{-3}))';
cbar.Visible='off';
        % caxis([0,3000])
        colormap(ccmap)
        set(colorbar,'visible','off')
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);


nexttile(2)
cbar=panelGlobal(sim.x, sim.y,log10(Bmeso_CavgSurf),linspace(0,2.5,20)  , sTitle="b. Meso-copepods", sUnits="mg_Cm^{-3}", sProjection=options.sProjection);
% cbar=panelGlobal(sim.x, sim.y,(Bmeso_CavgSurf),[0.5 100]  , sTitle="b. Meso-copepods", sUnits="mg_Cm^{-3}", sProjection=options.sProjection);
cbar.Label.String = '(log_{10}(mgC m^{-3}))';
cbar.Visible='off';
        % caxis([0,3000])
        colormap(ccmap)
        set(colorbar,'visible','off')
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        h = colorbar;
        % h.Layout.Tile = 'south';
        h.Ticks=[0 1 2];
        h.TickLabels={'10^0','10^1','10^2'};
        ylabel(h,'mg C m^{-3}','FontSize',10)

% nexttile(3)
% cbar=panelGlobal(sim.x, sim.y,log10(Bmacro_CavgSurf),linspace(-2,2.5,50)  , sTitle="c. Macro-copepods", sUnits="mg_Cm^{-3}", sProjection=options.sProjection);
% cbar.Label.String = '(log_{10}(mgC m^{-3}))';
% cbar.Visible='off';
%         % caxis([0,3000])
%         colormap(ccmap)
%         set(colorbar,'visible','on')
%         set(gca,'YTickLabel',[]);
%         set(gca,'XTickLabel',[]);
%         h = colorbar;
%         h.Layout.Tile = 'south';
%         ylabel(h,'mg_Cm^{-3}','FontSize',10)