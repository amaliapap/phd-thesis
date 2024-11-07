% Sensitivity analysis
load('C:\Users\ampap\Downloads\sensitivity_2024-08-27\HTL sensitivity\sensitivityHTL.mat')

%%
% Analyse the simulations
%
p = setupNUMmodel( bParallel=true );
p = parametersGlobal(p);
simHTL{1}.p=p;
mortHTL = linspace(0, 0.01, 8);

latitudes=[60,-5,24];
longitudes=[-15,5,-158];

for site=1:length(latitudes)
    lat=latitudes(site);
    lon=longitudes(site);

    idx = calcGlobalWatercolumn(lat, lon, simHTL{1});

    for i = 1:length(mortHTL)
        ProdHTL(site,i) = simHTL{i}.ProdHTLAnnual(1,idx.x, idx.y);
        ProdNPP(site,i) = simHTL{i}.ProdNetAnnual(1,idx.x, idx.y);
        effHTL(site,i) = ProdHTL(site,i)/ProdNPP(site,i);
    end

end

%% only for lambda_htl

for site=1:length(latitudes)
    lat=latitudes(site);
    lon=longitudes(site);

    idx = calcGlobalWatercolumn(lat, lon, simHTL{1});
    for i = 1:length(mortHTL)
        sLibName = loadNUMmodelLibrary();
B=simHTL{i}.B;
        for iTime=1:length(simHTL{i}.t)
            for k=1:length(simHTL{i}.z)
                u = [squeeze(simHTL{i}.N(iTime,idx.x,idx.y,k)), ...
                    squeeze(simHTL{i}.DOC(iTime,idx.x,idx.y,k)), ...
                    squeeze(simHTL{i}.Si(iTime,idx.x,idx.y,k)), ...
                    squeeze(simHTL{i}.B(iTime,idx.x,idx.y,k,:))'];
                rates = getRates(simHTL{i}.p, u, simHTL{i}.L(iTime,idx.x,idx.y,k), simHTL{i}.T(iTime,idx.x,idx.y,k), sLibName);
                [~,lambdaHTL(site,i,iTime,k)] = calcTrophicLevel(simHTL{i}.p,squeeze(B(iTime,idx.x,idx.y,k,:)),rates);
            end
        end
    end

end
lambdaHTL_sens=squeeze(mean(mean(lambdaHTL,3,'omitnan'),4,'omitnan')); % Average over time and the over depth
%% Set y-axes limits for the figures
epsilon_lim = [min(min(effHTL)) max(max(effHTL))];
prod_lim    = [min(min(ProdHTL)) max(max(ProdHTL))];
npp_lim     = [min(min(ProdNPP)) max(max(ProdNPP))];
TLhtl_lim   = [min(min(lambdaHTL_sens))-0.05 max(max(lambdaHTL_sens))+0.05];
%%
cmap  = flip(cmocean('deep',5));
ccmap = cmap(2:end,:); 

default_val=0.005;

x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
height=11; %figure height in cm
fig=figure(7);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','sensitivityGlobal_wc');
clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
set(groot,'defaultLineLineWidth',2)
tiledlayout(2,length(latitudes),'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor')

for site=1:length(latitudes)
    lat=latitudes(site);lon=longitudes(site);

   t1=nexttile();
    set(gca,'xscale','log')
    % yPoints=[1e-5,1e-4,1e-3,3e-2];
    semilogx(mortHTL, effHTL(site,:),'k:')%mte
    if site==1
        ylabel('Trophic efficiency [-]')
    end
    ylim(epsilon_lim)
    hold on
    xline(default_val,'--')
    yyaxis  right
        semilogx(mortHTL, squeeze(lambdaHTL_sens(site,:)))%lambda_htl

    ylim(TLhtl_lim)
    if site==1
        lgd=legend('\epsilon_{\mu}','','\lambda_{HTL}',Location='best',box='off',fontsize=10);%,'NPP','Prod_{HTL}')
        lgd.ItemTokenSize(1)=10;
    end
    if site==3
        ylabel('Trophic level [-]')
    end
    axis square

   t2=nexttile();
    semilogx(mortHTL, ProdNPP(site,:),'Color',cmap(4,:))%NPP
    hold on
    if site==1
        ylabel('Production (\mug C l^{-1}day^{-1})')
    end
    if site<3
        ylim(npp_lim)
    end
    yyaxis  right
    hp=semilogx(mortHTL,ProdHTL(site,:),'Color',cmap(2,:));%prodHTL
    if site<3
        ylim(prod_lim)
    end
    xline(default_val,'k--')

    xlabel('\mu_{HTL.0}')
    if site==1
        lgd=legend('NPP','Prod_{HTL}',Location='northwest',box='off',fontsize=10);%,'NPP','Prod_{HTL}')
        lgd.ItemTokenSize(1)=10;
    end
    axis square
    ax=gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = cmap(2,:);

  
    % Title positioning
    t1.TitleHorizontalAlignment = 'center';
    title_list={'a. Seasonally stratified','b. Upwelling','c. Oligotrophic'};
    title_list2={'d.','e.','f.'};

    ttl(site)=title(t1,string(title_list(site)),"FontWeight","normal");
    ttl2(site)=title(t2,string(title_list2(site)),"FontWeight","normal");

    pos2=ttl(site).Position(2);
    ttl(site).Position(2)=pos2+.0002;
end
% ttl(1).Position(2)=0.000047;

