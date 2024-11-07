%
% Make a plot of bacteria, phytoplankton, zooplankton, and the diatom ratio 
% 
% function sim = plotGlobalPhytoplankton(sim, options)
% arguments
%     sim struct;
    options.sProjection = 'mollweid';%'fast'; %projection to use. Defaults to 'fast'. Other projections
% %               requires that the mapping toolbox is installed. 
% %               Good projection is 'eckert4'.
% end
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
sLibName = loadNUMmodelLibrary();
ixTime =1:length(sim.t); %find(sim.t>(max(sim.t)-365)); %nTime = length(sim.t(sim.t >= max(sim.t-365))); % Just do the last year
% Get grid volumes:
sim.p.pathGrid='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\TMs\MITgcm_2.8deg\grid.mat';
load(sim.p.pathGrid,'dv','dz','dx','dy');
ix = ~isnan(sim.N(1,:,:,1)); % Find all relevant grid cells


Bphyto = zeros(length(sim.t), length(sim.x), length(sim.y));
Bbacteria = Bphyto;
Bzoo = Bphyto;
Diatomratio = Bphyto;
Bpassive=Bphyto;
Bactive=Bphyto;
Passiveratio=Bphyto;
Bmulti=Bphyto;
Bmulti_micro=Bphyto;
Bmulti_meso=Bphyto;

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

rGroups=calcRadiusGroups(sim.p); % for copepods it returns length/2
ixG = (p.ixStart(ixGeneralists)-p.idxB+1):(p.ixEnd(ixGeneralists)-p.idxB+1);
r_microG=find(rGroups(ixG)>=20/2 & rGroups(ixG)<=200/2);% radius in um divided by two to get the radius between (20-200um)/2
ixC = (p.ixStart(p.typeGroups==10)-p.idxB+1):(p.ixEnd(p.typeGroups==11)-p.idxB+1);
rGroups=calcRadiusGroups(sim.p); % for copepods it returns length/2

r_microC=find(rGroups(ixC)>=20/2 & rGroups(ixC)<=200/2);% radius in um divided by two to get the radius between (20-200um)/2
r_mesoC=rGroups(ixC)>200/2 & rGroups(ixC)<=2e4/2;
r_macroC=find(rGroups(ixC>2e4/2));
% ixC = (p.ixStart(p.typeGroups==10)-p.idxB+1):(p.ixEnd(p.typeGroups==11)-p.idxB+1);
% Bmulti_tot=squeeze(sum(sim.B(:,:,:,:,ixC),5));

for iTime = ixTime
    for i = 1:nX
        for j = 1:nY
            BphytoDiatoms = 0;
            BphytoOthers = 0;
            BzooPassive = 0;
            BzooActive = 0;
            BzooGens = 0;
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
                    [Bphytotmp, Bzootmp, ~] = ...
                        calcPhytoZoo(p, u, L(iTime,i,j,k), T(iTime,i,j,k), sLibName);
                    
                    conv = squeeze(dz(i,j,k));
                    Bphyto(iTime,i,j) = Bphyto(iTime,i,j) + sum(Bphytotmp)*conv; % mgC/m2
                    Bzoo(iTime,i,j) = Bzoo(iTime,i,j) + sum(Bzootmp)*conv; % mgC/m2
                    % % active and passive copepods
                    % Bpassive(iTime,i,j) = Bpassive(iTime,i,j) + sum(sim.B(iTime,i,j,k,ixPassive),5)*conv; % mgC/m2
                    % Bactive(iTime,i,j) = Bactive(iTime,i,j) + sum(sim.B(iTime,i,j,k,ixActive),5)*conv; % mgC/m2
                    Bmulti(iTime,i,j) = Bmulti(iTime,i,j) + sum(sim.B(iTime,i,j,k,ixMulti),5)*conv; % mgC/m2
                    Bmulti_micro(iTime,i,j) = Bmulti_micro(iTime,i,j) + sum(sim.B(iTime,i,j,k,ixC(r_microC)),5)*conv; % mgC/m2
                    Bmulti_meso(iTime,i,j) = Bmulti_meso(iTime,i,j) + sum(sim.B(iTime,i,j,k,ixC(r_mesoC)),5)*conv; % mgC/m2

                    BphytoDiatoms = BphytoDiatoms + sum(Bphytotmp(ixDiatoms))*conv;
                    BphytoOthers = BphytoOthers + sum(Bphytotmp(ixGeneralists))*conv;
                    BzooPassive = BzooPassive + sum(Bzootmp(ixPassive))*conv;
                    BzooActive = BzooActive + sum(Bzootmp(ixActive))*conv;
                end
            end
            Diatomratio(iTime,i,j) = BphytoDiatoms / (BphytoDiatoms+BphytoOthers);
            Passiveratio(iTime,i,j)=BzooPassive/(BzooPassive+BzooActive);
        end
    end
end
sim.Bphyto = Bphyto;
sim.Bzoo = Bzoo;
sim.Diatomratio = Diatomratio;
sim.Passiveratio=Passiveratio;
sim.Bmulti = Bmulti;
sim.Bmulti_micro=Bmulti_micro;
sim.Bmulti_meso=Bmulti_meso;
%%
Bphyto_mean=squeeze(mean(sim.Bphyto(ixTime,:,:),1))/1000;% g/m2
Bzoo_mean = squeeze(mean(sim.Bzoo(ixTime,:,:),1))/1000;% g/m2
Bmulti_mean = squeeze(mean(sim.Bmulti(ixTime,:,:),1))/1000;% g/m2
Bmulti_micro_mean = squeeze(mean(sim.Bmulti_micro(ixTime,:,:),1))/1000;% g/m2
Bmulti_meso_mean = squeeze(mean(sim.Bmulti_meso(ixTime,:,:),1))/1000;% g/m2
Bmeso_frac = Bmulti_meso_mean./(Bmulti_micro_mean+Bmulti_meso_mean+BzooGens_mean);
% calcGeneralistMicroZoo % For generalist micro-zooplankton
%%

Bmicro_C=double(squeeze(sum(sim.B(:,:,:,:,ixC(r_microC)),5)));
Bmeso_C=double(squeeze(sum(sim.B(:,:,:,:,ixC(r_mesoC)),5)));
Bmacro_C=double(squeeze(sum(sim.B(:,:,:,:,ixC(r_macroC)),5)));
Bmulti_tot=double(squeeze(sum(sim.B(:,:,:,:,ixC),5)));

% Take top 170m
Bmicro_Csurf=sum(Bmicro_C(:,:,:,1:3),4);
Bmeso_Csurf=sum(Bmeso_C(:,:,:,1:3),4);
Bmacro_Csurf=sum(Bmacro_C(:,:,:,1:3),4);
Bmicro_frac=Bmicro_Csurf./(Bmicro_Csurf+Bmeso_Csurf);



%% Zonal averages
% lastYear=((length(sim.t)/12-1)*12+1):length(sim.t);
lastYear=1:length(sim.t);
Bmicro_frac_zonal=squeeze(mean(Bmicro_frac(lastYear,:,:),2,'omitnan'));
Passiveratio_zonal=squeeze(mean(Passiveratio(lastYear,:,:),2,'omitnan'));
% Y=ProdNet_zonal;
% Y(Y==0)=1; % to avoid division by 0 below
% mte_zonal=ProdHTL_zonal./ProdNet_zonal;
% mte_zonal(isinf(mte_zonal))=NaN;


%% Make plot:
%
cmap=flip(cmocean('deep',100));
ccmap=cmap(2:end,:);      
% ccmap2=[cmap(1:50,:);flip(cmap(51:end,:))];
% for the ratios
cmap2=cmocean('delta',10);
cmap=flip(cmocean('deep',6));
cmap3=cmap(2:end,:);
cmap4=(cmocean('matter',6));
ccmap2=[cmap3(1:end,:);cmap4(2:end-1,:)];

x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
heightf=8; %figure height in cm
fig=figure(10);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto','Name','Copepod ratios-GlobalZonal');

clf
tiledlayout(1,3,'TileSpacing','tight','padding','tight')

t1=nexttile(1);
cbar=panelGlobal(sim.x, sim.y, 1-mean(Passiveratio(ixTime,:,:),1), sTitle="a. Active copepods ratio", sUnits="", sProjection=options.sProjection);
        cbar.Visible='off';
        clim([0 1])
        colormap(t1,ccmap2);
        % cbar=colorbar('horizontal');
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(gcf,'color','w');
        % post6=t6.Position(2);
        % cbar.Position(4)=0.025;
t2=nexttile(2);
cbar=panelGlobal(sim.x, sim.y, Bmeso_frac, sTitle="b. Meso-zooplankton fraction", sUnits="", sProjection=options.sProjection);
        cbar.Visible='off';
        clim([0 1])
        colormap(t2,ccmap2);
        cbar=colorbar('horizontal');
        cbar.Visible='off';
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);

t3=nexttile(3);
c=panelGlobal(sim.x, sim.y, 1-mean(Bmicro_frac(ixTime,:,:),1), sTitle="c. Meso-copepods ratio", sUnits="", sProjection=options.sProjection);
        c.Visible='off';
        clim([0 1])
        colormap(t3,ccmap2);
        c=colorbar('horizontal');
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        
        c.Position=[0.244188368055555,0.282221990495149,0.519244791666666,0.052863436123348];

        %%
tzonal=nexttile(3);
        s=surface(ixTime,sim.y,(1-Passiveratio_zonal)');
        s.EdgeColor = 'none';
        axis tight
        colormap(ccmap2)
        colorbar;
        clim([0 1])
        ylabel('Latitude')
        xlabel('Time (months)')
        ttl2=title('c.','FontWeight','normal');
        ttl2.HorizontalAlignment='left';
        ttl2.Position(1)=1.2;
        ax=gca;
        ax.YTick=-80:20:80;
        axis square

 tzonal=nexttile(4);
        s=surface(ixTime,sim.y,1-(Bmicro_frac_zonal)');
        s.EdgeColor = 'none';
        axis tight
        colormap(ccmap2)
        cb=colorbar;
        clim([0 1])
        % ylabel('Latitude')
        xlabel('Time (months)')
        ttl=title('d.','FontWeight','normal');
        ttl.HorizontalAlignment='left';
        ttl.Position(1)=1.2;
        ax=gca;
        ax.YTick=-80:20:80;
        axis square

%%

      %
tdiat=nexttile(5);
cbar=panelGlobal(sim.x, sim.y, mean(sim.Diatomratio(ixTime,:,:),1), sTitle="e. Diatom ratio", sUnits="", sProjection=options.sProjection);
        cbar.Visible='off';
        clim([0 1])
        colormap(t5,ccmap2)
        set(colorbar,'visible','off')
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        % cbar=colorbar('horizontal');
        % set(gcf,'color','w');
        post5=tdiat.Position(2);
        tdiat.Position(2)=post5-0.1;





