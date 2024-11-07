%
% Make a plot of bacteria, phytoplankton, zooplankton, and the diatom ratio 
%
% load matrix to skip calculations

function sim = plotGlobalPhytoplankton(sim, options)
thisDir=pwd;

% arguments
%     sim struct;
    options.sProjection = 'mollweid';%'fast'; %projection to use. Defaults to 'fast'. Other projections
% %               requires that the mapping toolbox is installed. 
% end
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
sLibName = loadNUMmodelLibrary();
ixTime = find(sim.t>(max(sim.t)-365)); %nTime = length(sim.t(sim.t >= max(sim.t-365))); % Just do the last year
% Get grid volumes:
% sim.p.pathGrid='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\TMs\MITgcm_2.8deg\grid.mat';
% load(sim.p.pathGrid,'dv','dz','dx','dy');
% ix = ~isnan(sim.N(1,:,:,1)); % Find all relevant grid cells

%%
Bphyto = zeros(length(sim.t), length(sim.x), length(sim.y));
Bbacteria = Bphyto;
Bzoo = Bphyto;
Diatomratio = Bphyto;
DiatomratioPhyto = Bphyto;
Passiveratio=Bphyto;
Bmulti=Bphyto;
Bdiatoms = Bphyto;
Buni = Bphyto;

nX = length(sim.x);
nY = length(sim.y);
nZ = length(sim.z);
%
% Extract fields from sim:
%
N = sim.N;
DOC = sim.DOC;
if isfield(sim.p,'idxSi')
    Si = sim.Si;
else
    Si = 0;
end
B = sim.B;
L = sim.L;
T = sim.T;
p = sim.p;
%
% find indices of diatoms and generalists:
%
ixGeneralists = p.typeGroups==1 | p.typeGroups==5;
ixDiatoms = p.typeGroups==3 | p.typeGroups==4;
ixPassive = p.typeGroups==10;
ixActive = p.typeGroups==11;
ixMulti = p.typeGroups==10 | p.typeGroups==11;
p=sim.p;
if ~isfield(sim,'Bphyto')

for iTime = ixTime
    for i = 1:nX
        for j = 1:nY
            BphytoDiatoms = 0;
            BphytoOthers = 0;
            BzooPassive = 0;
            BzooActive = 0;
           
            for k = 1:nZ
                if ~isnan(N(iTime,i,j,k))
                    if isfield(p,'idxSi')
                        u = [squeeze(N(iTime,i,j,k)), ...
                            squeeze(DOC(iTime,i,j,k)), ...
                            squeeze(Si(iTime,i,j,k)), ...
                            squeeze(B(iTime,i,j,k,:))'];
                    else
                        u = [squeeze(N(iTime,i,j,k)), ...
                            squeeze(DOC(iTime,i,j,k)), ...
                            squeeze(B(iTime,i,j,k,:))'];
                    end
                    [Bphytotmp, Bzootmp, Bbacteriatmp] = ...
                        calcPhytoZoo(p, u, L(iTime,i,j,k), T(iTime,i,j,k), sLibName);
                    
                    conv = squeeze(dz(i,j,k));
                    Bphyto(iTime,i,j) = Bphyto(iTime,i,j) + sum(Bphytotmp)*conv; % mgC/m2
                    Bbacteria(iTime,i,j) = Bbacteria(iTime,i,j) + sum(Bbacteriatmp)*conv; % mgC/m2
                    Bzoo(iTime,i,j) = Bzoo(iTime,i,j) + sum(Bzootmp)*conv; % mgC/m2
                    Bmulti(iTime,i,j) = Bmulti(iTime,i,j) + sum(sim.B(iTime,i,j,k,ixMulti),5)*conv; % mgC/m2
                    BphytoDiatoms = BphytoDiatoms + sum(Bphytotmp(ixDiatoms))*conv;
                    BphytoOthers = BphytoOthers + sum(Bphytotmp(ixGeneralists))*conv;
                    BzooPassive = BzooPassive + sum(Bzootmp(ixPassive))*conv;
                    BzooActive = BzooActive + sum(Bzootmp(ixActive))*conv;
                    Bdiatoms(iTime,i,j) = Bdiatoms(iTime,i,j) + sum(sim.B(iTime,i,j,k,ixDiatoms),5)*conv;
                    Buni(iTime,i,j) = Buni(iTime,i,j) + sum(sim.B(iTime,i,j,k,[ixGeneralists ixDiatoms]),5)*conv;

                end
            end
            DiatomratioPhyto(iTime,i,j) = BphytoDiatoms / (BphytoDiatoms+BphytoOthers);
            Diatomratio(iTime,i,j) = Bdiatoms(iTime,i,j)/Buni(iTime,i,j);

            Passiveratio(iTime,i,j)=BzooPassive/(BzooPassive+BzooActive);
        end
    end
end
sim.Bphyto = Bphyto;
sim.Bbacteria = Bbacteria;
sim.Bzoo = Bzoo;
sim.Diatomratio = Diatomratio;
sim.DiatomratioPhyto = DiatomratioPhyto;
sim.Passiveratio=Passiveratio;
sim.Bmulti = Bmulti;
end
%
Bac_mean=squeeze(mean(sim.Bbacteria(ixTime,:,:),1))/1000; % convert to g/m2
Bphyto_mean=squeeze(mean(sim.Bphyto(ixTime,:,:),1))/1000;% g/m2
Bzoo_mean = squeeze(mean(sim.Bzoo(ixTime,:,:),1))/1000;% g/m2
% Bmulti_mean = squeeze(mean(sim.Bmulti(ixTime,:,:),1))/1000;% g/m2

%% ----------------------------------
%           F I G U R E S (2)
%....................................
% Figure 1: Bacteria, Phytoplankton, 
%           Zooplankton
%------------------------------------
%
cmap=flip(cmocean('deep',100));
ccmap=cmap(2:end,:);      

% cmap2=cmocean('delta',10);
cmap=flip(cmocean('deep',6));
cmap3=cmap(2:end,:);
cmap4=(cmocean('matter',6));
ccmap2=[cmap3(1:end,:);cmap4(2:end-1,:)];

x0=0; 
y0=0;
width=16; %figure width in cm
heightf=8; %figure height in cm
fig=figure(11);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto','Name','Functional groups');

clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
tiledlayout(1,3,'TileSpacing','compact','padding','tight')

%%
t1=nexttile(1);
cbar=panelGlobal(sim.x, sim.y,log10(Bac_mean),[-1 2]  , sTitle="a. Bacteria", sUnits="mg_C/m^2", sProjection=options.sProjection);
cbar.Visible='off';
        % caxis([0,3000])
        colormap(t1,ccmap)
        cbar = colorbar('horizontal');
        cbar.Visible='off';
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
%
t2=nexttile(2);
cbar=panelGlobal(sim.x, sim.y,  log10(Bphyto_mean),[-1 2] ,sTitle="b. Phytoplankton", sUnits="mg_C/m^2", sProjection=options.sProjection);
cbar.Visible='off';
        % caxis([0,3000])
        colormap(t2,ccmap)
        cbar = colorbar('horizontal');
        cbar.Visible='off';
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
%
t3=nexttile(3);
cbar=panelGlobal(sim.x, sim.y,  log10(Bzoo_mean),[-1 2], sTitle="c. Zooplankton", sUnits="mg_C/m^2", sProjection=options.sProjection);
cbar.Label.String = '(mg_C m^{-2})';
cbar.Visible='off';
% caxis([0,3000])
colormap(t3,ccmap)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(colorbar,'visible','off')
cbar=colorbar;
ylabel(cbar,'g C m^{-2}','FontSize',10)
cbar.Ticks=[-2 0 2];
cbar.TickLabels={'10^{-2}','10^0','10^2'};

%% --------------------------------------
%           Figure 2: Group ratios
% ---------------------------------------
fig=figure(12);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 16 8],...
'PaperPositionMode','auto','Name','Group ratios');

clf
tiledlayout(1,2,'TileSpacing','compact','padding','tight')

% t5=nexttile(1);
% cbar=panelGlobal(sim.x, sim.y, mean(sim.Diatomratio(ixTime,:,:),1), sTitle="a. Diatom ratio", sUnits="", sProjection=options.sProjection);
%         cbar.Visible='off';
%         clim([0 1])
%         colormap(t5,ccmap2)
%         cbar=colorbar('horizontal');
%         cbar.Visible='off';
%         set(gca,'YTickLabel',[]);
%         set(gca,'XTickLabel',[]);
%         % cbar=colorbar('horizontal');
%         % set(gcf,'color','w');
%         % post5=t5.Position(2);
        

 t6=nexttile();
cbar=panelGlobal(sim.x, sim.y, mean(sim.DiatomratioPhyto(ixTime,:,:),1), sTitle="b. Diatom_{phyto} ratio", sUnits="", sProjection=options.sProjection);
        cbar.Visible='off';
        clim([0 1])
        colormap(t6,ccmap2)
        cbar=colorbar('horizontal');
        cbar.Visible='off';
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);

t7=nexttile();
cbar=panelGlobal(sim.x, sim.y, 1-mean(sim.Passiveratio(ixTime,:,:),1), sTitle="c. Active copepods ratio", sUnits="", sProjection=options.sProjection);
        cbar.Visible='off';
        clim([0 1])
        colormap(t7,ccmap2);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        cbar=colorbar;
        set(gcf,'color','w');
        ylabel(cbar,'(g C m^{-2})','FontSize',10)
        cbar.Ticks=[.2 .4  .6 .8 1];
        % cbar.TickLabels={'10^{-2}','10^0','10^2'};
        ylabel(cbar,'[-]','FontSize',10)

cd(thisDir)

      



