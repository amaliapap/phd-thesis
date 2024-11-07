%
% returns annual average trophic network tile
% functions called: getTheta, extractRatesTime, calcTrophicLevel, CfixedRespTotalTime 
%
% First need to run extractWCfromGlobal
%
function TrophicNetworksWCAnnual_Tiles(sim,site)
% sim=simWC;
% 
% Colors and colormaps:
% interaction-lines colormaps
cmap_pass=[241, 166, 153,250*0.4;239, 108, 111,250*0.4;219, 7, 61,250*0.4]/250; % 4-element vectors where the 4th is alpha
cmap_act=[248, 196, 113,250*0.4 ; 245, 176, 65,250*0.4;243, 156, 18,250*0.4 ]/250; % orange shades
cmap_htl=[152, 220, 199,250*0.4 ; 27, 188, 155,250*0.4; 4, 117, 111,250*0.4 ]/250; % blue-green shades
% set colormaps of different shades for each size-group
incr=3; n=4;
cmapBl=brewermap(length(n)+incr,'blues');
cmapD=brewermap(length(n)+incr,'greens');
cmapOr=brewermap(length(n-1)+incr,'OrRd'); % orange


cmap=brewermap(15,'PRGn');
ccmap = cmap(3:end-2,:); 

c_cp=[206, 0, 91]/250;%cmap(2,:);
c_ap=[243, 156, 18]/250;

bnorm = 100; % divide to normalize biomass vector/bubble size in scatter plots
multFactor=3000; % multiplication factor for bubble size

%
%---------------------------------------
%    Indices of different groups
%---------------------------------------
% generalists
[~, idxG ]=find(sim.p.typeGroups==5); 
ixG  = (sim.p.ixStart( idxG (1)):sim.p.ixEnd( idxG (end))) -sim.p.idxB+1;
% diatoms
[~, idxD ]=find(sim.p.typeGroups==3); 
ixD  = (sim.p.ixStart( idxD (1)):sim.p.ixEnd( idxD (end))) -sim.p.idxB+1;

% passive copepods
[~, idxPC ]=find(sim.p.typeGroups==10); 
ixPC  = (sim.p.ixStart( idxPC (1)):sim.p.ixEnd( idxPC (end))) -sim.p.idxB+1;

% active copepods
[~,idxAC]=find(sim.p.typeGroups==11); % Find index of Active Copepods
ixAC = (sim.p.ixStart(idxAC(1)):sim.p.ixEnd(idxAC(end))) -sim.p.idxB+1;

% vector of all size-groups
totsizevec=sim.p.m(sim.p.idxB:end);
preysize=totsizevec;
[~,predsort_idx]=sort(preysize);

%% Calculate trophic level/ take annual averages
biomvec = squeeze(mean(sim.B,1));
rates   = sim.rates;
jDOC    = squeeze(mean(rates.jDOC,1));
jLreal  = squeeze(mean(rates.jLreal,1));
jFreal  = squeeze(mean(rates.jFreal,1));
mortHTL = squeeze(mean(rates.mortHTL,1));

jPP   = jDOC+jLreal; % Uptake of Primary Production 
jF    = jFreal;      % Uptake of Food
ratio = jPP./(jPP+jF);
ratio(isnan(ratio))=[];
theta = getTheta(sim.p);
TLtime= zeros(size(sim.B));
lambda_htl_time = zeros(1,length(sim.t));

for iTime=1:length(sim.t)
    ratesTime=extractRatesTime(sim.rates,iTime);
    [TLtime(iTime,:),lambda_htl_time(iTime)]=calcTrophicLevel(sim.p,sim.B(iTime,:),ratesTime); 
    [~,~,nppGtime(iTime,:),~]   = CfixedRespTotalTime(sim,sim.rates,iTime);

end

TL = mean(TLtime,1);
lambda_htl = mean(lambda_htl_time);
nppG = mean(nppGtime,1);
%% Eflow real
Eflow=zeros(sim.p.n-sim.p.idxB+1);
for i=1:sim.p.n-sim.p.idxB+1
    for j=1:1:sim.p.n-sim.p.idxB+1
        % Eflow(i,j)= rates.jFreal(i)*theta(i,j)/sum(theta(i,:),2); 
        Eflow(i,j) = jFreal(i)*theta(i,j)/sum(theta(i,:),2)*biomvec(i); 
    end
end
Eflow(isnan(Eflow)) = 0;
Eflow_sort = Eflow(predsort_idx,predsort_idx);
%% Microbial (MTE) & Mean Progressive Trophic Efficiency (PTE)
mte=mean(sim.ProdHTLwc)/mean(sim.ProdNetwc);
pte=mte^(1/(lambda_htl-1));
%% Progressive Trophic Transfer efficiency from one TL to the next
nppGroups=zeros(1,sim.p.n-sim.p.idxB+1);
nppGroups(ixG(1):ixD(end))=nppG; % all the unicellulars

for i=1:ceil(max(TL))
    [~,idxTL{i}]=find(TL<=i);
    Eflow_in{i}=sum(sum(Eflow(idxTL{i},:),2)+sum(nppGroups(idxTL{i})),1); % i-group is the predator
    Eflow_out{i}= sum(sum(Eflow(:,idxTL{i}),1)+sum(mortHTL(idxTL{i}).*biomvec(idxTL{i}))); % i-group is prey   
end

epsilonTL=cell2mat(Eflow_out)./cell2mat(Eflow_in);
%%
% ATTENTION! setting POM trophic level to 0 for plotting purposes!
TL(end)=0;
fluxes_to_htl=0;
intercHTL=0;
t1=nexttile([1 3]);
for i = 2:(length(totsizevec)) % predators
    for j = 1:(i-1)%length(totsizevec) % prey
        interc = (Eflow_sort(i, j))./ sum(sum(Eflow));
        intrcFlow(i,j)=interc;
        intercHTL = double(mortHTL(j).*biomvec(j));
        inHTL(i)=intercHTL;
        fluxes_to_htl=fluxes_to_htl+intercHTL;   % this will be used only to scale the plotted HTL-disc
        x_HTL = sim.mHTLdepth(1);     % assumed position of HTL in the x-axis
        y_HTL=max(TL)+1; % HTL trophic level is plotted as the highest trophic level in the plankton community +1
        % create HTL fluxes
        if (intercHTL > 0.0005)
            % Define the x and y coordinates of the start and end points
               x_end = totsizevec(j);
               y_end = TL(j);
            xyHTL= [x_end, x_HTL; y_end, y_HTL];
                       row_colorHTL=1; % lightest shade
            if intercHTL>=1e-3 && intercHTL<1e-2
                row_colorHTL=1; % lightest shade
            elseif intercHTL>=1e-2 && intercHTL<1e-1
                row_colorHTL=2;
            elseif intercHTL>1e-1 %&& interc<1e-1
                row_colorHTL=3; % darkest shade for the strongest interactions
            end
            hold on
            if intercHTL>0 %0.0001
                plot(xyHTL(1,:),xyHTL(2,:) , 'Color', cmap_htl(row_colorHTL,:),'LineWidth',0.8);
            end % Plot the points with log scale on x-axis
        end
        % create predation fluxes
        if interc>0.0005 % try higher value
            % Define the x and y coordinates of the start and end points
            x_start = totsizevec(i);    x_end = totsizevec(j);
            y_start = TL(i);            y_end = TL(j);

            % Generate t values for interpolation
            t = [1, 2];
            xy = [x_start, x_end; y_start, y_end];
            % Cubic spline interpolation
            pp = csapi(t, xy);
            tInterp = linspace(1, 2, 100);
            xyInterp = ppval(pp, tInterp);
            %---------- end of interpolation ------------

            row_color=1;
            if interc>=0.025 && interc<0.03
                row_color=1; % lightest shade
            elseif interc>=0.03 && interc<0.05
                row_color=2;
            elseif interc>0.05
                row_color=3; % darkest shade for the strongest interactions
            end
            connectingColor=cmap_act(row_color,:);
            % Plot the spline curve
            plot(xyInterp(1,:), xyInterp(2,:), 'Color', connectingColor, 'LineWidth', .8)
        end
    end
end
set(gca, "XScale", 'log')

idxCP=ixPC(1):ixPC(end);
idxCA=ixAC(1):ixAC(end);

    scatter(totsizevec, TL, biomvec ./ bnorm .* multFactor, ...
        'filled', 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'k');
    hd=plot(totsizevec((ixD)),TL((ixD)),'go',LineWidth=.4); % make diatoms of biom<< appear
    hd.MarkerEdgeColor='k';
    hd.MarkerFaceColor='g';
    hd.MarkerSize=3;
    h1=scatter(totsizevec((ixG)),TL((ixG)),biomvec((ixG))./bnorm.*multFactor,ratio(ixG),...
        'filled','markeredgecolor','k');
    alpha(.7)
    cbar=colorbar('horizontal','Direction','reverse');
    colormap(ccmap)
    cbar.Label.String = 'Photorophy degree [-]';
    cbar.Location="north";
    a=get(cbar);   % gets properties of colorbar
    cb=a.Position; % gets the positon and size of the color bar
    cbar.Position=[cb(1) cb(2) 0.25 cb(4)];
    if site>1
        set(cbar,'visible','off')
    end
    clim([0 1])
    hold on
    h2=scatter(totsizevec((ixD)),TL((ixD)),biomvec((ixD))./bnorm.*multFactor,...
        'filled','markerfacecolor','g','markeredgecolor','g',MarkerFaceAlpha=.5);
    h3=scatter(totsizevec((idxCP)),TL((idxCP)),biomvec((idxCP))./bnorm.*multFactor,...
        'filled','markerfacecolor',c_cp,'markeredgecolor','k',MarkerEdgeAlpha=.7);
    h4=scatter(totsizevec((idxCA)),TL((idxCA)),biomvec((idxCA))./bnorm.*multFactor,...
        'filled','markerfacecolor','#eb8628','markeredgecolor','k',MarkerEdgeAlpha=.5);
    hPOM=scatter(totsizevec(end),TL(end),biomvec(end)./bnorm.*multFactor,...
        'filled','markerfacecolor', '#e5e7e9' ,'markeredgecolor','k',MarkerFaceAlpha=.2,MarkerEdgeAlpha=.4); %POM trophic level in fig at 0
    h_htl=scatter(x_HTL,y_HTL,80*(fluxes_to_htl/5),'markerfacecolor','none','markeredgecolor','k',MarkerFaceAlpha=.5);

    set(gca, 'xscale', 'log')
ylim([0.5, 4.5])
xlim([1e-9, 1e5])
if site==3 % show xlabel only on the bottom of the figure
    xlabel('Body mass [\mugC]')
end
ylabel('Trophic level [-]')

xticks(logspace(-7, 5, 5));
if site==3
    leg=legend([h1,h2,h3,h4,hPOM,h_htl],'Generalists','Diatoms','Copepods_{pass.}',...
    'Copepods_{act.}','POM','HTL','Location','best','box','off'); 
    leg.NumColumns=3;
    leg.Location='southoutside';
    leg.Position=[0.1,   -0.002,    0.7593,    0.0568];
end
grid on
set(gcf, 'color', 'w');
title_list={'a. Seasonally stratified','c. Upwelling','e. Oligotrophic'};
title(t1,append(string(title_list(site)),': \epsilon_{\mu}=',...
    num2str(round(mte,2,"significant"))),"FontWeight","normal")
ylim([0 4])

t2=nexttile();
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
set(groot,'defaultLineLineWidth',2)
    set(gca,'YTickLabel',[]);
    set(gca, 'YScale','log')
    maxEl=.25;
    yPoints=[.001,.05,maxEl];

    yyaxis right
    plot(1:ceil(max(TL)),epsilonTL)
    ax=gca;
    if site==3
        xlabel('Trophic level [-]')
    end
    ylabel('Progressive Efficiency [-]')
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
        set(gca, 'YTick',yPoints, 'YTickLabel',yPoints)
        ylim([0 maxEl])
        letter=['b','d','f'];
title(t2,append(letter(site),'. $\bar{\epsilon}_{\lambda}$=',...
    num2str(round(pte,2,'significant'))),"FontWeight","normal",Interpreter="latex")

% disp(append('pte=',num2str(pte),', epsilonTL_mean=',num2str(mean(epsilonTL))));


