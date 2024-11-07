% Figures
% Comparison of NUM with satellite products
% To create the plots with satellite data, download from: http://orca.science.oregonstate.edu/npp_products.php
% Standard VGPM, Eppley-VGPM, CbPM2, CAFE
%
function nppComparison(sim,showData)

% Define colorbar colors
cmap=flip(cmocean('deep',100));
ccmap=cmap(2:end,:);

% Figure specifications:
x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
heightf=14; %figure height in cm
fig=figure('Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width heightf-3],...
'PaperPositionMode','auto');

set(gcf,'color','w');
if showData==true
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    % tiles.InnerPosition = [0.13,0.11,0.65,0.8150]; % Make space for colorbars
    t1=nexttile(1);
end
cbar = panelGlobal(sim.x,sim.y,(sim.ProdNetAnnual),[0,1000],...
    sTitle='a. Model', sProjection='mollweid');
cbar.Label.String = '(mg C m^{-2}d^{-1})';
cbar.Visible='off';
if showData==true
    colormap(t1,ccmap);
    set(colorbar,'visible','off')
else
    colormap(ccmap);
    cbar=colorbar;
    ylabel(cbar,'mg C m^{-2} d^{-1}','FontSize',10)
end
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


    %% Satellite data figures
   if showData==true
       addpath(append(pwd,'\NPP satellite data'))

   % NPPcafeVsNUMfinal(sim)
           % load('npp_climatlogy_cafe_2002_2020.mat','npp')
        load('cafe_mean_2003_2019.mat','mat_cafe')
        load('cbpm_mean_2003_2019.mat','mat_cbpm')
        load('vgpm_mean_2003_2019.mat','mat_vgpm')
        load('vgpme_mean_2003_2019.mat','mat_vgpme')
        npp_sat={'b. CAFE','c. CbPM','d. VGPM-Standard','e. VGPM-Eppley'};
        sat_product=zeros(length(npp_sat),length(sim.y),length(sim.x));
        %create meshgrid with NUM's coordinates
        [lon2,lat2]=meshgrid(sim.x,sim.y);

        orig0=(sim.ProdNetAnnual(end,:,:)); %NUM NPP
        orig=squeeze(orig0)'; % flip NUM npp data
        for k=1:length(npp_sat)
            switch k
                case 1
                    npp=mat_cafe;
                case 2
                    npp=mat_cbpm;
                case 3
                    npp=mat_vgpm;
                case 4
                    npp=mat_vgpme;
            end
            lat_dim=size(npp,1);
            lon_dim=size(npp,2);

            lat=zeros(lat_dim,1);
            lat(1)=1/6;
            for i=2:lat_dim
                lat(i)=lat(i-1)+1/6;
            end
            lat=lat-90;
            %
            lon=zeros(lon_dim,1);
            lon(1)=1/6;
            for i=2:lon_dim
                lon(i)=lon(i-1)+1/6;
            end

            % nppAvgAnnual=squeeze(mean(nppAvg(1:11,:,:),1));
            nppCAFEinterpRec=interpolateToNUM(sim, npp);
            %%
            % create meshgrid with CAFE's coordinates
            [lon1,lat1]=meshgrid(lon,lat);

            % interpolate NPP CAFE dimensions to NUM's
            nppCAFEinterp=interp2(lon1,lat1,squeeze(npp),lon2,lat2);

            nppCAFEinterpRecenter= [nppCAFEinterp(:,65:end),nppCAFEinterp(:,1:64)];
            % make land white
            for i=1:length(sim.y)
                for j=1:length(sim.x)
                    if (nppCAFEinterpRecenter(i,j)==-9999)
                        nppCAFEinterpRecenter(i,j)=NaN;
                        orig(i,j)=NaN;
                    end
                end
            end
            sat_product(k,:,:)=nppCAFEinterpRec;
        end

        for k=1:(length(npp_sat)-1)
            npp_product=squeeze(flip(sat_product(k,:,:)));
            nexttile(k+1)
            cbar = panelGlobal(sim.x,sim.y,npp_product',[0,1000],...
                sTitle=string(npp_sat(k)), sProjection='mollweid');
            % cbar.Label.String = '(mg C m^{-2}d^{-1})';
            cbar.Visible='off';
            colormap(ccmap);
            set(colorbar,'visible','off')
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
        end
        k=4;
        npp_product=squeeze(flip(sat_product(k,:,:)));
        nexttile(k+1)
        cbar = panelGlobal(sim.x,sim.y,npp_product',[0,1000],...
            sTitle=string(npp_sat(k)), sProjection='mollweid');
        cbar.Visible='off';
        colormap(ccmap);
        cbar=colorbar;
        ylabel(cbar,'mg C m^{-2} d^{-1}','FontSize',10)
        cbar.Location="southoutside";
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
    
cbar.Position=[0.533884297520661,0.141826923076923,0.353719008264463,0.041185894427009];
   end
end
% interpolateToNUM and recenter
    function npp_interp=interpolateToNUM(sim,npp_annual_avg)
        lat_dim=size(npp_annual_avg,1);
        lat=zeros(lat_dim,1);
        lat(1)=1/6;
        for i=2:lat_dim
            lat(i)=lat(i-1)+1/6;
        end
        lat=lat-90;
        %
        lon=zeros(2160,1);
        lon(1)=1/6;
        for i=2:2160
            lon(i)=lon(i-1)+1/6;
        end

        % create meshgrid with CAFE's coordinates
        [lon1,lat1]=meshgrid(lon,lat);

        %create meshgrid with NUM's coordinates
        [lon2,lat2]=meshgrid(sim.x,sim.y);

        % interpolate NPP CAFE dimensions to NUM's
        nppCAFEinterp=interp2(lon1,lat1,squeeze(npp_annual_avg),lon2,lat2);

        orig0=(sim.ProdNetAnnual(end,:,:)); %NUM NPP
        orig=squeeze(orig0)'; % flip NUM npp data

        nppCAFEinterpRecenter= [nppCAFEinterp(:,65:end),nppCAFEinterp(:,1:64)];
        % make land white
        nppCAFEinterpRecenter(nppCAFEinterpRecenter==-9999)=NaN;
        orig(orig==0)=NaN;


        npp_interp=nppCAFEinterpRecenter;
    end

    function NPPcafeVsNUMfinal(sim)
        % % load('npp_climatlogy_cafe_2002_2020.mat','npp')
        % load('cafe_mean_2003_2019.mat','mat_cafe')
        % load('cbpm_mean_2003_2019.mat','mat_cbpm')
        % load('vgpm_mean_2003_2019.mat','mat_vgpm')
        % load('vgpme_mean_2003_2019.mat','mat_vgpme')
        % npp_sat={'b. CAFE','c. CbPM','d. VGPM-Standard','e. VGPM-Eppley'};
        % sat_product=zeros(length(npp_sat),length(sim.y),length(sim.x));
        % %create meshgrid with NUM's coordinates
        % [lon2,lat2]=meshgrid(sim.x,sim.y);
        % 
        % orig0=(sim.ProdNetAnnual(end,:,:)); %NUM NPP
        % orig=squeeze(orig0)'; % flip NUM npp data
        % for k=1:length(npp_sat)
        %     switch k
        %         case 1
        %             npp=mat_cafe;
        %         case 2
        %             npp=mat_cbpm;
        %         case 3
        %             npp=mat_vgpm;
        %         case 4
        %             npp=mat_vgpme;
        %     end
        %     lat_dim=size(npp,1);
        %     lon_dim=size(npp,2);
        % 
        %     lat=zeros(lat_dim,1);
        %     lat(1)=1/6;
        %     for i=2:lat_dim
        %         lat(i)=lat(i-1)+1/6;
        %     end
        %     lat=lat-90;
        %     %
        %     lon=zeros(lon_dim,1);
        %     lon(1)=1/6;
        %     for i=2:lon_dim
        %         lon(i)=lon(i-1)+1/6;
        %     end
        % 
        %     % nppAvgAnnual=squeeze(mean(nppAvg(1:11,:,:),1));
        %     nppCAFEinterpRec=interpolateToNUM(sim, npp);
        %     %%
        %     % create meshgrid with CAFE's coordinates
        %     [lon1,lat1]=meshgrid(lon,lat);
        % 
        %     % interpolate NPP CAFE dimensions to NUM's
        %     nppCAFEinterp=interp2(lon1,lat1,squeeze(npp),lon2,lat2);
        % 
        %     nppCAFEinterpRecenter= [nppCAFEinterp(:,65:end),nppCAFEinterp(:,1:64)];
        %     % make land white
        %     for i=1:length(sim.y)
        %         for j=1:length(sim.x)
        %             if (nppCAFEinterpRecenter(i,j)==-9999)
        %                 nppCAFEinterpRecenter(i,j)=NaN;
        %                 orig(i,j)=NaN;
        %             end
        %         end
        %     end
        %     sat_product(k,:,:)=nppCAFEinterpRec;
        % end
end