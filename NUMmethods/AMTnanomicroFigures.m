% Produces figures that compare nano+microplankton biomass depth-integrated
% between NUM model and AMT data
% modified version from Serra-Pompei et al.,2022
% -----------------------------------------------
%*************************************************
%-------------------------------------------------
%            AMT nano-,micro-
%-------------------------------------------------
%*************************************************
function AMTnanomicroFigures(sim)

load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Model EVALUATION\matrices2022\dataprotistsbiomass_corrected.mat');
% biomass in gC/m^2
simInt=calcIntegrateGlobal(sim,sim.B , false); % in gC/m^2


r=calcRadiusGroups(sim.p);
ixr_nanomicro=find(r>=2 & r<=200);

% remove copepod classes
ixC=find(sim.p.typeGroups==10 | sim.p.typeGroups==11);
idxC=(sim.p.ixStart(ixC(1)):sim.p.ixEnd(ixC(end)))-sim.p.idxB+1; %remove nutrient grids

for j=1:length(idxC)
    if any(ixr_nanomicro(:) == idxC(j))
        ixr_nanomicroU(j) = idxC(j);
    end
end
ixr_nanomicroU(ixr_nanomicroU==0)=[];


T   = dataprotistsbiomass; % import table with AMT data
x   = T(:,2); % longitude
y   = T(:,5); % latitude
B   = T(:,3); % Biomass of pico- and nano-plankton
transect  = table2array(T(:,4));
month_corrected=table2array(T(:,8));
% converting tables to arrays
lat_dat = table2array(x); 
lon_dat = table2array(y);
B_dat   = table2array(B);

% indices for the different shapes/months in Serra-Pompei et al. (2022)
romb     =  1:27;   % diamond/ AMT-12
romb_may =  1:14;    % becaus first half of that cruise was in May  
cuadr    =  28:53;  % square/ AMT-14
tri      =  54:77;  % triangle/ AMT-13
tri_sep  = tri(10:24);
% AMT-13 started from North to South, so we need to flip indexes
% tri      = flip(tri);
% tri_sep  = flip(tri_sep);
lonWrapped = wrapTo180(sim.x); % convert longitude ordinates from [0,360[-->[-180,180]

% Global simulation from NUMmodel 
% for may-june-october
yrNo = sim.p.tEnd/365;
may  = (yrNo-1)*12+5;
june = may+1;
oct  = (yrNo-1)*12+10;
sep  = oct-1;

r_uni=(sim.p.ixStart(1):sim.p.ixEnd(2))-sim.p.idxB+1;
rmicro=find(r>=200);
rnano=find(r>=20 & r<200);
% simInt=calcIntegrateGlobal(sim,sim.B , false); % in gC/m^2
Bmicro=calcIntegrateGlobal(sim,sim.B(:,:,:,:,rmicro),false); % gC/m2
Bnano=calcIntegrateGlobal(sim,sim.B(:,:,:,:,rnano),false); % gC/m2
% Bmicro=permute(Bmicro,[3 1 2]);
% Bnano=permute(Bnano,[3,1,2]);

 
Bmicro_5  =  squeeze(Bmicro(may,:,:));
Bmicro_6  =  squeeze(Bmicro(june,:,:));
Bmicro_9  =  squeeze(Bmicro(sep,:,:));
Bmicro_10 =  squeeze(Bmicro(oct,:,:));
Bnano_5   =  squeeze(Bnano(may,:,:));
Bnano_6   =  squeeze(Bnano(june,:,:));
Bnano_9   =  squeeze(Bnano(sep,:,:));
Bnano_10  =  squeeze(Bnano(oct,:,:));

% Taking monthly average of the model data
Bg_5avg  =  Bmicro_5+Bnano_5;
Bg_6avg  =  Bmicro_6+Bnano_6;
Bg_9avg =  Bmicro_9+Bnano_9;
Bg_10avg =  Bmicro_10+Bnano_10;

%----------------------------------------------------------
% Find model Biomass for the same coordinates with the data
%----------------------------------------------------------
Pbiom_model5=zeros(1,length(lat_dat));
Pbiom_model6=zeros(1,length(lat_dat));
Pbiom_model9=zeros(1,length(lat_dat));
Pbiom_model10=zeros(1,length(lat_dat));


for i=1:length(lat_dat)

     [ ~, idx_lonG] = min( abs( lonWrapped-lon_dat(i) ) );

     [ ~, idx_latG] = min( abs( sim.y-lat_dat(i) ) );

     ix_lon(i)=idx_lonG;
     ix_lat(i)=idx_latG;

% next step is to find combined coordinates
    Pbiom_model6(i)=Bg_6avg(idx_lonG(1),idx_latG(1));
    Pbiom_model5(i)=Bg_5avg(idx_lonG(1),idx_latG(1));
    Pbiom_model10(i)=Bg_10avg(idx_lonG(1),idx_latG(1));
    Pbiom_model9(i)=Bg_9avg(idx_lonG(1),idx_latG(1));

end
Pbiom_modelAMT12=[Pbiom_model5(romb_may), Pbiom_model6((length(romb_may))+1:length(romb))];
Pbiom_modelAMT13=[Pbiom_model10((length(lat_dat)-length(tri)+1):62),Pbiom_model9(tri_sep) ];


%% plot model vs data for each month

cmap=[169,223,205; 149,167,255; 135,73,163]./255; % mermaiddoc

x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
heightf=10; %figure height in cm

fig=figure(9);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto','Name','AMT nano-,micro-');
clf()
set(gcf,'color','w');
 
tiledlayout(1,4,"Padding","compact",TileSpacing="compact")
nexttile(2)
scatter(B_dat(romb),lat_dat(romb),25,'d','markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),MarkerFaceAlpha=.5)    
hold on
plot(Pbiom_model6(romb),sim.y(ix_lat(romb)),'LineWidth',2,'color',cmap(1,:))
ylim([-65,65])
ax1=gca;
ax1.YTick=[-50 -30 -15 0 15 30 50];
ylabel('Latitude')
% xlabel("biomass (gC/m^2)")
% title("b.",'FontWeight','normal')%AMT 12:May-June
pbaspect([1 2.8 1])
plotlabel('b',false);

nexttile(3)
scatter(B_dat(tri),lat_dat(tri),25,"^",'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),MarkerFaceAlpha=.5)
hold on
plot(Pbiom_modelAMT13,sim.y(ix_lat(tri)),'LineWidth',2,'color',cmap(2,:))
ylim([-65,65])

xlabel("Biomass (g C m^{-2})")
plotlabel('c',false);
set(gca,'Yticklabel',[]) 
pbaspect([1 2.8 1])

nexttile(4)
scatter(B_dat(cuadr),lat_dat(cuadr),25,"s",'markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),MarkerFaceAlpha=.5)
hold on
plot(Pbiom_model5(cuadr),sim.y(ix_lat(cuadr)),'LineWidth',2,'color',cmap(3,:))
ylim([-65,65])

% xlabel("Biomass (gCm^{-2})")
set(gca,'Yticklabel',[]) 
pbaspect([1 2.8 1])
plotlabel('d',false);


%---------------------
% Plot transects
%---------------------
% latlim=[-47,60];
latlim=[-55,55];
lonlim=[-50,-1];

nexttile(1)
axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on','fontsize',8)
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
axis off
SC1=scatterm(lat_dat(transect==1),lon_dat(transect==1),15,'d','filled','markerfacecolor',cmap(1,:),'markeredgecolor','none');
SC3=scatterm(sim.y(ix_lat(54:77)),sim.x(ix_lon(54:77)),15,'s','filled','markerfacecolor',cmap(3,:),'markeredgecolor','none');
SC2=scatterm(lat_dat(transect==3),lon_dat(transect==3),15,'^','filled','markerfacecolor',cmap(2,:),'markeredgecolor','none');
% SC3=scatterm(lat_dat(transect==2),lon_dat(transect==2),15,'s','filled','markerfacecolor',cmap(2,:),'markeredgecolor','none');

leg=legend([SC1,SC2,SC3],{'AMT 12','AMT 13','AMT 14'},'fontsize',8);
legend boxoff
leg.ItemTokenSize(1) = 10;
leg.Location='southoutside';
leg.Position=[0,0.07,0.28,0.029752065426062];
leg.NumColumns=3;
% title('a.','FontWeight','normal')
pbaspect([1 2.8 1])
a=plotlabel('a',false);
a.Position(2)=1.01;

% Position.Label=[ -0.113133805095692,1.085419055119125,-4e-15]; %30degW
