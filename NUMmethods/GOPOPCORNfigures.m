%--------------------------------------------------------------------------
% Plots of GO-POPCORN data and comparison with model output (sim.Bmicro)
% data: https://doi.org/10.1038/s41597-022-01809-1
% this is only the code, the model results here don't fit
% (the biomass is too low, but it is probably better in the newest version)
% copepods and detritus should also be include in Bmicro-POC_avg_uM (?)
%
%--------------------------------------------------------------------------
%                          DATA PROCESSING
%-   -   - -   - -   - -   - -   - -   - -   - -   - - -   - - -   - - -   
% *STEP 1*: Convert data table to array and convert units of POC to ugC/L
% *STEP 2*: Select size-classes of the model with r in [0.7/2, 30/2] um
% *STEP 3*: Remove the rows where POC is NaN
% *STEP 4*: Group data by cruise and subgroup by cruise_station
% *STEP 5*: Calculate the mean values at each station --> 'summary_Dep'
% *STEP 6*: Isolate only the last year of the model output and match the
%           model months to the observation months
% *STEP 7*: Find the model biomass at the same coordinates, depth and month
%           as the data, 'Pbiom_model'
% *STEP 8*: Calculate the mean value over depth at each coordinate, 
%           'Pbiom_model_depth_avg'
% *STEP 9*: Apply the same grouping as the data (STEP 4) to the model
%           output from STEP 8 'summary_model'
% *STEP 10*: Plot transect
% *STEP 11*: For each cruise, plot observations and model output 
% - -   - - -   - - -   - - -   - - -   - - -   - - -   - - -   - - -   - - 
%--------------------------------------------------------------------------
% Arctic cruises are not included
function GOPOPCORNfigures(sim,r)

load('dataGOPOPCORNver2.mat')
%%
% sim=simFun1;

% MODEL OUTPUT
% load('simGlobalDevelopCompareFun.mat')
% sim=simGlobalDevelopComapreFun;


% convert data variables to vectors 
dataPOP=dataGOPOPCORNver2;
cruise=table2array(dataPOP(:,"Cruise"));
cruise_station=table2array(dataPOP(:,"Cruise_Station"));
lat=table2array(dataPOP(:,"Latitude"));
lat1=lat;
lon=table2array(dataPOP(:,"Longitude"));
lon1=lon;
depth=table2array(dataPOP(:,"Depth"));
month=table2array(dataPOP(:,"Month"));
day=table2array(dataPOP(:,"Day"));
year=table2array(dataPOP(:,"Year"));
POCavg_uM=table2array(dataPOP(:,"POCavg_uM")).*12;  % converted to mugC/L

%nutricline_1uM_Interp=table2array(dataPOP(:,18));

lonWrapped = wrapTo180(sim.x); % convert longitude ordinates from [0,360[-->[-180,180]

% The size range is 2.4ug-30um, so Biomass needs to be updated.
% r = calcRadiusGroups(sim.p);
% load('rGroups.mat')
% r=rGroups;
r_popcorn=find(r>=0.7/2 & r<=30/2);

Bpopcorn=sim.B(:,:,:,:,r_popcorn); % mugC/L
Bpopcorn_sum=squeeze(sum(Bpopcorn,5)); % sum all size-classes
%%
dataPOPc=[cruise cruise_station lat lon depth month POCavg_uM];
nan_rows= any(isnan(dataPOPc(:,7)),2); %POCavg_uM=NaN
% Remove rows where POC=NaN
%...........................
dataPOPclean = dataPOPc(~nan_rows, :);

% if cruise_station==NaN then assign it to a new non-existing unique station
for i=1: height(dataPOPclean)
 if isnan(dataPOPclean(i,2))
    dataPOPclean(i,2)=1000+i;
 end
end

% convert matrix to table
dataPOP = array2table(dataPOPclean);
% Create custom column names
customColumnNames = {'Cruise', 'Cruise_Station', 'Latitude','Longitude','Depth','Month','POCavg_uM'};
dataPOP.Properties.VariableNames = customColumnNames;

%% -------------------------------------------------------------------
% grouping based on cruise and subgrouping based on cruise station
%-------------------------------------------------------------------

grouping_cruise_cs=findgroups(dataPOP.Cruise, dataPOP.Cruise_Station);
% create table same as dataPOP including the 'Grouping'
temp_table=table(grouping_cruise_cs,dataPOP.Cruise,dataPOP.Latitude,dataPOP.Longitude,dataPOP.Depth,...
    dataPOP.Month,dataPOP.POCavg_uM,'VariableNames',{'Grouping','Cruise','Latitude','Longitude','Depth','Month','POCavg_uM'});
%.................................................
% average POC concentartion and depth at each subgroup.
%.................................................
summary_Dep=groupsummary(temp_table,'Grouping','mean'); % this gives us mean values at each station

temp_array=table2array(temp_table);
% remove arctic cruises: 1701,1709,1901,2018
temp_array_noA=temp_array(1:2296,:); %2296
myData=temp_array_noA;
%%

%---------------------------------------------------------------------
% Find model Biomass for the same coordinates and depth with the data
%---------------------------------------------------------------------
lat=myData(:,3);
lon=myData(:,4);
depth=myData(:,5);
month=myData(:,6);
% Initialize Pbiom_model with zeros
Pbiom_model = zeros(length(lat),1);
Pbiom_model_depth_avg=Pbiom_model;

for i=1:length(lat)

     [ ~, idx_lonG] = min( abs( lonWrapped-lon(i) ) );

     [ ~, idx_latG] = min( abs( sim.y-lat(i) ) );

     [ ~, idx_depth] = min( abs( sim.z-depth(i) ) );

     % Find the index of the corresponding month over the last year of the
     % simulation
     % [ ~, idx_mon] = min( abs( sim.t-(month(i)+144) ) );
     idx_mon = (length(sim.t)-12) + month(i);

     % Save the indices just in case
     ix_lon(i)=idx_lonG;
     ix_lat(i)=idx_latG;
     ix_dep(i)=idx_depth;
     ix_mon(i)=idx_mon;

    % Assign POC values for the indices above
    Pbiom_model(i)=Bpopcorn_sum(idx_mon(1),idx_lonG(1),idx_latG(1),idx_depth(1)); % in mugC/L (I think)
    Pbiom_model_depth_avg(i)= mean(Bpopcorn_sum(idx_mon(1),idx_lonG(1),idx_latG(1),idx_depth(1)),4)';
    % Pbiom_modelnew(i,:)=Bpopcorn_sum(:,idx_lonG(1),idx_latG(1),idx_depth); %(all month from the simulation)

end

depth_model=sim.z(ix_dep);
lat_model=sim.y(ix_lat);
lon_model=lonWrapped(ix_lon);
month_model=mod(ix_mon,12)';
for i=1:length(month_model)
    if month_model(i)==0
        month_model(i)=12;
    end
end

%%
grouping_numbers=unique(myData(:,1));
%(group_cruise_name(:,2)==cruise)
sumPOC_group=zeros(1,length(grouping_numbers));
sumPOC_model=zeros(1,length(grouping_numbers));

count=zeros(1,length(grouping_numbers));
count_model=zeros(1,length(grouping_numbers));

depth_max=zeros(1,length(grouping_numbers));
myData_wo_model=myData;
myData_with_model=[myData lat_model lon_model depth_model Pbiom_model month_model];

myData_sorted=sortrows(myData_with_model,1); % sorted data based on grouping
myData_sorted_modelOnly=[myData_sorted(:,1) myData_sorted(:,8:11)];
for ix_g=1:length(grouping_numbers)
    for j=1:length(myData_sorted)
        if myData_sorted(j,1)==ix_g
            sumPOC_group(ix_g)=sumPOC_group(ix_g) +myData_sorted(j,7);
            count(ix_g)=count(ix_g)+1;
            lat_group(ix_g)=myData_sorted(j,3);
            lon_group(ix_g)=myData_sorted(j,4);
            cruise_name(ix_g)=myData_sorted(j,2);
            month_group(ix_g)=myData_sorted(j,6);
            % if depth_update~=myData_with_model(j,10)
            %     sumPOC_model(ix_g)=sumPOC_model(ix_g) +myData_with_model(j,11);
            %     count_model(ix_g)=count_model(ix_g)+1;
            %     depth_update=myData_with_model(j,10);
            % end
            if myData(j,5)>depth_max(ix_g)
                depth_max(ix_g)=myData(j,5); % max_depth at each grouping/station
            end
        end
    end
end


myData_model_unique=unique(myData_sorted_modelOnly,"rows");
myData_model_unique2=unique(myData_sorted_modelOnly(:,1:end-1),"rows");

for ix_g=1:length(grouping_numbers)
    for j=1:length(myData_model_unique)
        if myData_model_unique(j,1)==ix_g
            sumPOC_model(ix_g)=sumPOC_model(ix_g) +myData_model_unique(j,5);
            count_model(ix_g)=count_model(ix_g)+1;
            lat_model_new(ix_g)=myData_model_unique(j,2);
            lon_model_new(ix_g)=myData_model_unique(j,3);
            if myData_model_unique(j,4)>depth_max(ix_g)
                depth_max_model(ix_g)=myData_model_unique(j,4); % max_depth at each grouping/station
            end
        end
    end
end



for ix_g=1:length(grouping_numbers)
    meanPOC_group(ix_g)=sumPOC_group(ix_g)/count(ix_g);
    meanPOC_model(ix_g)=sumPOC_model(ix_g)/count_model(ix_g);
end

myData_final=[grouping_numbers cruise_name' lat_group' lon_group' month_group' meanPOC_group' lat_model_new' lon_model_new' meanPOC_model'];


colNames = {'grouping','cruise name','lat','lon','month','mean POC','lat model','lon model','mean POC model'};
finalDataTable = array2table(myData_final,'VariableNames',colNames);
%%
% cruise 46 has only Nan
cruise_name=[7,9,13,18,28,1319,1418,1701,1709,1901,2018];

% find indices of each cruise
cruise1_ix=find(myData_final(:,2)==cruise_name(1));
cruise2_ix=find(myData_final(:,2)==cruise_name(2));
cruise3_ix=find(myData_final(:,2)==cruise_name(3));
cruise4_ix=find(myData_final(:,2)==cruise_name(4));
cruise5_ix=find(myData_final(:,2)==cruise_name(5));
cruise6_ix=find(myData_final(:,2)==cruise_name(6));
cruise7_ix=find(myData_final(:,2)==cruise_name(7));

months=1:12;

%%          TRANSECTS PLOT

x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
heightf=12; %figure height in cm
fig=figure('Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf],...
'PaperPositionMode','auto');
clf
    set(gcf,'color','w');
    surface(lon,lat,POCavg_uM,'EdgeColor','none')
    axis tight
    title('POC AVG uM')

cmap=  [255, 0, 0; 255, 128, 0; 255, 255, 0;...
 128, 255, 0; 0, 255, 0; 0, 255, 128; 0, 255, 255;...
 0, 128, 255; 0, 0, 255; 128, 0, 255; 0, 210, 241;...
 255, 0, 255; 255, 0, 128]./255;

latlim=[min(lat),max(lat)];
lonlim=[min(lon),max(lon)];
axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on','fontsize',8)
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
axis off

SC1=scatterm(lat1(cruise==cruise_name(1)),lon1(cruise==cruise_name(1)),15,'s','filled','markerfacecolor',cmap(1,:),'markeredgecolor','none');
SC2=scatterm(lat1(cruise==cruise_name(2)),lon1(cruise==cruise_name(2)),15,'s','filled','markerfacecolor',cmap(2,:),'markeredgecolor','none');
SC3=scatterm(lat1(cruise==cruise_name(3)),lon1(cruise==cruise_name(3)),15,'s','filled','markerfacecolor',cmap(3,:),'markeredgecolor','none');
SC4=scatterm(lat1(cruise==cruise_name(4)),lon1(cruise==cruise_name(4)),15,'s','filled','markerfacecolor',cmap(4,:),'markeredgecolor','none');
SC5=scatterm(lat1(cruise==cruise_name(5)),lon1(cruise==cruise_name(5)),15,'s','filled','markerfacecolor',cmap(10,:),'markeredgecolor','none');
SC6=scatterm(lat1(cruise==cruise_name(6)),lon1(cruise==cruise_name(6)),15,'s','filled','markerfacecolor',cmap(12,:),'markeredgecolor','none');
SC7=scatterm(lat1(cruise==cruise_name(7)),lon1(cruise==cruise_name(7)),15,'s','filled','markerfacecolor',cmap(7,:),'markeredgecolor','none');

leg=legend([SC1,SC2,SC3,SC4,SC5,SC6,SC7],...
    {'Cruise 1','Cruise 2','Cruise 3','Cruise 4','Cruise 5','Cruise 6','Cruise 7'},'fontsize',10);


leg.ItemTokenSize(1) = 10;
leg.Location='southoutside';
leg.NumColumns=4;
 sgtitle('GO-POPCORN transects','FontSize', 15)

