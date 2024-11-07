% 
% load('simGlobalDevelopCompareFun.mat')
% sim=simGlobalDevelopComapreFun;
% load('rGroups.mat');

function PicophytoFiguresAtlantic(sim,r)
load('poco_cpf_db_v2.mat')
load('pico_insitudata.mat')

%%
% convert data variables to vectors 
datapoco=pococpfdbv2;
insitu=table2array(datapoco(:,5));
experiment=table2array(datapoco(:,15));

% for i=1:length(insitudata)
%     if (insitudata(i,2)>=30.63) && (insitudata(i,2)<=43.4)...
%         && (insitudata(i,3)>=-10.33) && (insitudata(i,3)<=21.77)
%        insitudata(i,:)=[];
%     end
% end


% Assuming insitudata is a matrix with at least 3 columns (rows x columns)

% Define latitude and longitude limits
latLimits = [30.63, 43.4];
lonLimits = [-10.33, 21.77];
lonAtlantic = [-60, 22];
% Create logical indices for rows that meet the criteria
indicesToKeep = (insitudata(:, 2) >= latLimits(1) & insitudata(:, 2) <= latLimits(2)) & ...
                (insitudata(:, 3) >= lonLimits(1) & insitudata(:, 3) <= lonLimits(2));

indicesAtlantic = (insitudata(:, 3) >= lonAtlantic(1) & insitudata(:, 3) <= lonAtlantic(2)); 

% Use logical indexing to keep only the rows that meet the criteria
filteredData = insitudata((indicesToKeep==false & indicesAtlantic), :);
%%
insitudata=filteredData;
insitu=insitudata(:,1);
lat=insitudata(:,2);
lon=insitudata(:,3);
month_poco=insitudata(:,6);

% MODEL OUTPUT
lonWrapped = wrapTo180(sim.x); % convert longitude ordinates from [0,360]->[-180,180]
% r=rGroups;
r_pico=find(r<=2);%-sim.p.nNutrients;
r_pico=r_pico(4:end);
Bpico1=sim.B(:,:,:,:,r_pico); % mugC/L
Bpico_sum=squeeze(sum(Bpico1,5)); % sum all size-classes

% Bpico2=sim.Bpico;
Bpico2=sim.B(:,:,:,r_pico);

Bpico=Bpico_sum;
% chlVolume_model=sim.ChlVolume(:,:,:,1); % top layer
%%
%------------------------------------------------------------------------------
% Find model Biomass and Chla for the same coordinates and months as the data
%------------------------------------------------------------------------------

Picobiom_model=zeros(length(lat),1); 
chla_model=zeros(length(lat),1);

for i=1:length(lat)

     [ ~, idx_lonG] = min( abs( lonWrapped-lon(i) ) );

     [ ~, idx_latG] = min( abs( sim.y-lat(i) ) );

     ix_lon(i)=idx_lonG;
     ix_lat(i)=idx_latG;

% next step is to find combined coordinates
    Picobiom_model(i)=Bpico(month_poco(i),idx_lonG(1),idx_latG(1),1); % values taken at top layer
    % chla_model(i,:)=chlVolume_model(month_poco(i),idx_lonG(1),idx_latG(1));
end



%%
x0=0;
y0=0;
width=16;
figure_number = 2;           % Desired figure number
fig = figure(figure_number); % Create the figure with the specified number
set(fig, 'Renderer','Painters','Units','centimeters',...
    'Position',[x0 y0 width 7],...
    'PaperPositionMode','auto','Name','Pico transects');
clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
set(gca,'defaultLineLineWidth',2)
%-------------------------
% Plot  transects
%-------------------------

cmap=  [255, 0, 0; 255, 128, 0;...
    255, 255, 0; 128, 255, 0;...
    0, 255, 0;     0, 255, 128;...
    0, 255, 255;   0, 128, 255;... 
    0, 0, 255;   128, 0, 255;... 
    0, 210, 241; 255, 0, 255;...
    255, 0, 128]./255;
latN=table2array(datapoco(:,3));
lonE=table2array(datapoco(:,4));


latlim=[min(lat),max(lat)];
lonlim=[min(lon),max(lon)];
s_size=7; % merker size
axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on','fontsize',8)
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
axis off
SC1=scatterm(latN(experiment=='AMT'),lonE(experiment=='AMT'),s_size,'d','filled','markerfacecolor',cmap(1,:),'markeredgecolor','none');
SC2=scatterm(latN(experiment=='MAREDAT'),lonE(experiment=='MAREDAT'),s_size,'^','filled','markerfacecolor',cmap(2,:),'markeredgecolor','none');
SC3=scatterm(latN(experiment=='ANT24'),lonE(experiment=='ANT24'),s_size,'s','filled','markerfacecolor',cmap(3,:),'markeredgecolor','none');
SC4=scatterm(latN(experiment=='ANT25'),lonE(experiment=='ANT25'),s_size,'s','filled','markerfacecolor',cmap(4,:),'markeredgecolor','none');
SC5=scatterm(latN(experiment=='ANT26'),lonE(experiment=='ANT26'),s_size,'s','filled','markerfacecolor',cmap(8,:),'markeredgecolor','none');
SC6=scatterm(latN(experiment=='GRAFF'),lonE(experiment=='GRAFF'),s_size,'s','filled','markerfacecolor',cmap(10,:),'markeredgecolor','none');

leg=legend([SC1,SC2,SC3,SC4,SC5,SC6],...
    {'AMT','MAREDAT','ANT24','ANT25','ANT26','GRAFF'},'fontsize',8);
leg.ItemTokenSize(1) = 10;
leg.Location='southoutside';
leg.NumColumns=3;
% 
% sgtitle('Picophytoplankton intercmparison in situ data','FontSize', 20)
%%

x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
heightf=12; %figure height in cm
fig=figure(7);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto');
clf
set(gcf,'color','w');
colordata=20*ones(1,length(lat));


tiledlayout(1,3,"TileSpacing","tight","Padding","compact")

nexttile(1)
geoscatter(lat,lon,colordata/2.5,insitu,'filled')
geolimits(gca,[-60 60],[-60 0])
geolimits('manual')

c = colorbar('horizontal');
clim([0,60])
c.Label.String = "[mg m^{-3}]";
c.Label.FontSize=10;
title('\textit{in situ}',FontSize=10,FontWeight='normal',Interpreter='latex')

nexttile(2)
geoscatter(lat,lon,colordata/2.5,Picobiom_model,'filled')
% geolimits(gca,[0 50],[0 20])
geolimits(gca,[-60 60],[-60 0])

c = colorbar('horizontal');
clim([0,60])
c.Label.String = "[mg m^{-3}]";
c.Label.FontSize=10;
title('Model',FontSize=10,FontWeight='normal')

nexttile(3)
geoscatter(lat,lon,colordata/2.5,Picobiom_model-insitu,'filled')
geolimits(gca,[-60 60],[-60 0])
c = colorbar('horizontal');
clim([-60,60])
c.Label.String = "[mg m^{-3}]";
c.Label.FontSize=10;
title('$\Delta$(Model-\textit{in situ})',FontSize=10,FontWeight='normal',Interpreter='latex')
colormap(whitejet);                                               

sgtitle('Picoplankton biomass',fontsize=12);