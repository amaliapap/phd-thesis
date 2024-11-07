% Sensitivity analysis in the Chemostat of the parameters:
%   -mixing rate 
%   -light 
%   -temperature 
%   -HTL mortality
% each simulation is run for yrNo years and 
% the results are extracted at the final step
%
% output: 
%   -microbial trophic efficiency (mte)
%   -mean trophic level ingested by HTL (lambda_htl)
%   -net primary production (NPP)
%   -HTL production (ProdHTL)
%   -mean mass of primary producers (mNPP)
%   -mean mass of organisms ingested by HTL (mHTL)
%
% functions called:
% •	setHTL
% •	CfixedRespTotalwatcol
% •	maxTrophicLevel
% •	mNPPmHTLchem

cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

% Simulation setup
p = setupNUMmodel();
p = parametersChemostat(p);
yrNo=10;
p.tEnd = yrNo*365;

% Default values
massHTL = 1;
mortHTL_def=0.005;
bHTLdecline = false; 
bHTLquadratic = true;   
setHTL(mortHTL_def, massHTL, bHTLquadratic, bHTLdecline)

T=10;
L = 100;
d = logspace(-3,-0.5,10);
prodHTL    = zeros(1,length(d));
npp_Sens   = zeros(1,length(d));
lambda_htl = zeros(1,length(d));   
mte        = zeros(1,length(d));   
mNPP_sens  = zeros(1,length(d));
mHTL_sens  = zeros(1,length(d));
%% Sensitivity for d
for i=1:length(d)
    p.d=d(i);
    cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
    setHTL(mortHTL_def, massHTL, bHTLquadratic, bHTLdecline)
    sim = simulateChemostat(p, L,T);
    addpath 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency'
    [mass_rangeUni, CfixedUni,npp]   = CfixedRespTotal(sim,sim.rates);
    [~,htl_lambda]     = calcTrophicLevel(sim.p,sim.B(end,:),sim.rates);

    [mNPP,mHTL]             =mNPPmHTLchem(sim);
    m    = p.m(p.idxB:end);
    mAll = p.m(sim.p.ixStart(1):sim.p.ixEnd( end));              % mass of all plankton groups
    %Find total HTL losses at each size range
    mass_range1 = logspace(log10(min(min(mAll))),log10(max(max(mAll))), 15);
    total_mort  = zeros(1,length(mass_range1)-1); % initialize HTL-production
    % Select mortHTL whithin ranges where there are multiple groups
    for ix = (p.ixStart(1):p.ixEnd(end)-1)-p.idxB+1 % without POM
        m_ix = p.m(ix+p.idxB-1);
        for iRange = 1:length(mass_range1)-1
            if m_ix >= mass_range1(iRange) && m_ix <= mass_range1(iRange+1)
                total_mort(iRange) = total_mort(iRange)+sim.rates.mortHTL(ix).*sim.B(end,ix)';
            end
        end
    end
    prodHTL(i)    = sum(total_mort);
    npp_Sens(i)   = sum(npp);
    lambda_htl(i) = htl_lambda;
    mte(i)        = prodHTL(i)/sum(npp);
    pte(i)        = mte(i)^(1/(lambda_htl(i)-1));  
    mNPP_sens(i)  = mean(mNPP(i));
    mHTL_sens(i)   = mHTL(i);
end
sensitivity_mixRate = [npp_Sens; mte; prodHTL; lambda_htl; mNPP_sens; mHTL_sens; pte];

%%  sENSITIVITY FOR LIGHT

L = linspace(10,200,length(d));
p.d = 0.1;
prodHTL    = zeros(1,length(L));
npp_Sens   = prodHTL;
lambda_htl = prodHTL; 
mte        = prodHTL; 
pte        = prodHTL; 

for i=1:length(L)
    % p.d=d;
    cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
    setHTL(mortHTL_def, massHTL, bHTLquadratic, bHTLdecline)
    sim = simulateChemostat(p, L(i),T);
    addpath 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency'
    [mass_rangeUni, CfixedUni,npp]   = CfixedRespTotal(sim,sim.rates);
    [~,htl_lambda]     = calcTrophicLevel(sim.p,sim.B(end,:),sim.rates);

    m    = p.m(p.idxB:end);
    mAll = p.m(sim.p.ixStart(1):sim.p.ixEnd( end));              % mass of all plankton groups
    %Find total HTL losses at each size range
    mass_range1 = logspace(log10(min(min(mAll))),log10(max(max(mAll))), 15);
    total_mort  = zeros(1,length(mass_range1)-1); % initialize HTL-production
    % Select mortHTL whithin ranges where there are multiple groups
    for ix = (p.ixStart(1):p.ixEnd(end)-1)-p.idxB+1 % without POM
        m_ix = p.m(ix+p.idxB-1);
        for iRange = 1:length(mass_range1)-1
            if m_ix >= mass_range1(iRange) && m_ix <= mass_range1(iRange+1)
                total_mort(iRange) = total_mort(iRange)+sim.rates.mortHTL(ix).*sim.B(end,ix)';
            end
        end
    end
    prodHTL(i)    = sum(total_mort);
    npp_Sens(i)   = sum(npp);
    lambda_htl(i) = htl_lambda;
    mte(i)        = prodHTL(i)/sum(npp);
    pte(i)        = mte(i)^(1/(lambda_htl(i)-1));
    [mNPP,mHTL]   = mNPPmHTLchem(sim);
    mNPP_sens(i)  = mean(mNPP(i));
    mHTL_sens(i)   = mHTL(i);
end
sensitivity_light = [npp_Sens; mte; prodHTL; lambda_htl; mNPP_sens; mHTL_sens; pte];
%%  SENSITIVITY FOR TEMPERATURE

Temp = linspace(0,30,length(d));
p.d = 0.1;
L = 100;
prodHTL    = zeros(1,length(Temp));
npp_Sens   = prodHTL;
lambda_htl = prodHTL; 
mte        = prodHTL; 
for i=1:length(Temp)
    % p.d=d;
    cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
    setHTL(mortHTL_def, massHTL, bHTLquadratic, bHTLdecline)
    sim = simulateChemostat(p,L,Temp(i));
    addpath 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency'
    [mass_rangeUni, CfixedUni,npp]   = CfixedRespTotal(sim,sim.rates);
    [~,htl_lambda]     = calcTrophicLevel(sim.p,sim.B(end,:),sim.rates);

    m    = p.m(p.idxB:end);
    mAll = p.m(sim.p.ixStart(1):sim.p.ixEnd( end));   % mass of all plankton groups
    %Find total HTL losses at each size range
    mass_range1 = logspace(log10(min(min(mAll))),log10(max(max(mAll))), 15);
    total_mort  = zeros(1,length(mass_range1)-1); % initialize HTL-production
    % Select mortHTL whithin ranges where there are multiple groups
    for ix = (p.ixStart(1):p.ixEnd(end)-1)-p.idxB+1 % without POM
        m_ix = p.m(ix+p.idxB-1);
        for iRange = 1:length(mass_range1)-1
            if m_ix >= mass_range1(iRange) && m_ix <= mass_range1(iRange+1)
                total_mort(iRange) = total_mort(iRange)+sim.rates.mortHTL(ix).*sim.B(end,ix)';
            end
        end
    end
    prodHTL(i)    = sum(total_mort);
    npp_Sens(i)   = sum(npp);
    lambda_htl(i) = htl_lambda;
    mte(i)        = prodHTL(i)/sum(npp);
    pte(i)        = mte(i)^(1/(lambda_htl(i)-1));
    [mNPP,mHTL]   = mNPPmHTLchem(sim);
    mNPP_sens(i)  = mean(mNPP(i));
    mHTL_sens(i)   = mHTL(i);
end
% sensitivity_temp = [npp_Sens; mte; prodHTL; lambda_htl;mNPP_sens;mHTL_sens];
sensitivity_temp = [npp_Sens; mte; prodHTL; lambda_htl; mNPP_sens; mHTL_sens; pte];

%% sensitivity morthl

p.d = 0.1;
L = 100;

newMortHTL = linspace(0,0.01, length(d));
newMortHTL = logspace(-3,-1, length(d));

newParameter=newMortHTL;


prodHTL    = zeros(1,length(newParameter));
npp_Sens   = prodHTL;
lambda_htl = prodHTL; 
mte        = prodHTL; 

for  iparam = 1:length(newParameter)    
     cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
        setHTL(newMortHTL(iparam), massHTL, bHTLquadratic, bHTLdecline)
        sim = simulateChemostat(p,L,T);
    addpath 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Trophic Efficiency'
    [mass_rangeUni, CfixedUni,npp]   = CfixedRespTotal(sim,sim.rates);
    [~,htl_lambda]     = calcTrophicLevel(sim.p,sim.B(end,:),sim.rates);
    m    = p.m(p.idxB:end);
    mAll = p.m(sim.p.ixStart(1):sim.p.ixEnd( end));              % mass of all plankton groups
    %Find total HTL losses at each size range
    mass_range1 = logspace(log10(min(min(mAll))),log10(max(max(mAll))), 15);
    total_mort  = zeros(1,length(mass_range1)-1); % initialize HTL-production
    % Select mortHTL whithin ranges where there are multiple groups
    for ix = (p.ixStart(1):p.ixEnd(end)-1)-p.idxB+1 % without POM
        m_ix = p.m(ix+p.idxB-1);
        for iRange = 1:length(mass_range1)-1
            if m_ix >= mass_range1(iRange) && m_ix <= mass_range1(iRange+1)
                total_mort(iRange) = total_mort(iRange)+sim.rates.mortHTL(ix).*sim.B(end,ix)';
            end
        end
    end
    prodHTL(iparam)    = sum(total_mort);
    % mass_range2   = mass_range1(1:end-1);
    npp_Sens(iparam)   = sum(npp);
    lambda_htl(iparam) = htl_lambda;
    mte(iparam)      = prodHTL(iparam)/sum(npp);
    pte(iparam)        = mte(i)^(1/(lambda_htl(i)-1));
    [mNPP,mHTL]   = mNPPmHTLchem(sim);
    mNPP_sens(iparam)  = mean(mNPP(iparam));
    mHTL_sens(iparam)   = mHTL(iparam);
end
sensitivity_MortHTL = [npp_Sens; mte; prodHTL; lambda_htl; mNPP_sens; mHTL_sens; pte];

sensitivity_vals={sensitivity_mixRate,sensitivity_light,sensitivity_temp,sensitivity_MortHTL};
sens_array=cell2mat(sensitivity_vals);
mte=sens_array(2,:);
ProdHTL=sens_array(3,:);
ProdNPP=sens_array(1,:);
lambdaHTL_sens=sens_array(4,:);
m_NPP=sens_array(5,:);
m_HTL=sens_array(6,:);

epsilon_lim = [min(min(mte)) max(max(mte))];
prod_lim    = [min(min(ProdHTL)) max(max(ProdHTL))];
npp_lim     = [min(min(ProdNPP)) max(max(ProdNPP))];
TLhtl_lim   = [min(min(lambdaHTL_sens))-0.05 max(max(lambdaHTL_sens))+0.05];
mNPP_lim = [min(min(m_NPP)) max(max(m_NPP))];
mHTL_lim = [min(min(m_HTL)) max(max(m_HTL))];
ylim(npp_lim)
%% FIGURES
cmap  = flip(cmocean('deep',5));
ccmap = cmap(2:end,:); 

x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
height=11.5; %figure height in cm

fig=figure(8);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','SensitivityChemostat');
clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
set(groot,'defaultLineLineWidth',2)
tiledlayout(3,4,'TileSpacing','tight','Padding','compact','TileIndexing','columnmajor')
%
% Mixing rate
default_value=0.1;

nexttile()
set(gca,'xscale','log')
semilogx(d, sensitivity_mixRate(2,:),'k')
xline(default_value,'k--')
hold on
ylim(epsilon_lim)
ylabel('Trophic efficiency [-]')
yyaxis  right
semilogx(d,sensitivity_mixRate(4,:))
ylim(TLhtl_lim)
lgd1=legend('$\epsilon_{\mu}$','','$\lambda_{\mathrm{HTL}}$','Fontsize',12,Location='best',box='off',Interpreter='latex');
lgd1.ItemTokenSize(1)=10;
axis square
axis tight
plotlabel('a',false)

nexttile()
set(gca,'xscale','log')
semilogx(d, sensitivity_mixRate(1,:),'Color',cmap(4,:))
xline(default_value,'k--')
ylim(npp_lim)
hold on
ylabel('Production (\mugCl^{-1}day^{-1})')
yyaxis  right
semilogx(d,sensitivity_mixRate(3,:),'Color',cmap(2,:))
ylim(prod_lim)
plotlabel('b',false)
lgd2=legend('NPP','','Prod_{HTL}',Location='best',box='off',fontsize=10);
lgd2.ItemTokenSize(1)=10;
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(2,:);
axis square
axis tight

nexttile()
set(gca,'xscale','log')
semilogx(d, sensitivity_mixRate(5,:),'Color',cmap(4,:))
xline(default_value,'k--')
ylim(mNPP_lim)
hold on
ylabel('Mass (\mug C)')
yyaxis  right
semilogx(d,sensitivity_mixRate(6,:),'Color',cmap(2,:))
ylim(mHTL_lim)
xlabel('Mixing rate (day^{-1})')
lgd=legend('m_{NPP}','','m_{HTL}',Location='best',box='off',fontsize=10);
lgd.ItemTokenSize(1)=10;
axis square
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(2,:);
plotlabel('c',false)
axis tight


% Light
default_value=100;
nexttile()
L = linspace(10,200,length(d));
plot(L, sensitivity_light(2,:),'k')
xline(default_value,'k--')
ylim(epsilon_lim)
hold on
yyaxis  right
plot(L,sensitivity_light(4,:))
ylim(TLhtl_lim)
axis square
plotlabel('d',false)

nexttile()
plot(L, sensitivity_light(1,:),'Color',cmap(4,:))
xline(default_value,'k--')
ylim(npp_lim)
hold on
yyaxis  right
plot(L,sensitivity_light(3,:),'Color',cmap(2,:))
ylim(prod_lim)
axis square
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(2,:);
plotlabel('e',false)


nexttile()
set(gca,'xscale','log')
semilogx(L,sensitivity_light(5,:),'Color',cmap(4,:))
xline(default_value,'k--')
ylim(mNPP_lim)
hold on
yyaxis  right
semilogx(L,sensitivity_light(6,:),'Color',cmap(2,:))
ylim(mHTL_lim)
xlabel('Light (\mumol photons m^{-2} s^{-1})')
axis square
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(2,:);
plotlabel('f',false)

%
% Temperature
default_value=10;
nexttile()
plot(Temp, sensitivity_temp(2,:),'k')
xline(default_value,'k--')
hold on
yyaxis  right
plot(Temp,sensitivity_temp(4,:))
ylim(TLhtl_lim)
axis square
plotlabel('g',true)

nexttile()
plot(Temp, sensitivity_temp(1,:),'Color',cmap(4,:))
xline(default_value,'k--')
ylim(npp_lim)
hold on
yyaxis  right
plot(Temp,sensitivity_temp(3,:),'Color',cmap(2,:))
ylim(prod_lim)
axis square
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(2,:);
plotlabel('h',false)

nexttile()
set(gca,'xscale','log')
semilogx(Temp,sensitivity_temp(5,:),'Color',cmap(4,:))
xline(default_value,'k--')
ylim(mNPP_lim)
hold on
yyaxis  right
semilogx(Temp,sensitivity_temp(6,:),'Color',cmap(2,:))
ylim(mHTL_lim)
xlabel('Temperature (^{o}C)')
axis square
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(2,:);
plotlabel('i',false)

%
% MortHTL
default_value=0.005;
nexttile()
set(gca,'xscale','log')
semilogx(newMortHTL, sensitivity_MortHTL(2,:),'k')
xline(default_value,'k--')
ylim(epsilon_lim)
hold on
yyaxis  right
semilogx(newMortHTL,sensitivity_MortHTL(4,:))
ylabel('Trophic level [-]')
axis square
plotlabel('j',false)

nexttile()
semilogx(newMortHTL, sensitivity_MortHTL(1,:),'Color',cmap(4,:))
xline(default_value,'k--')
ylim(npp_lim)
hold on
yyaxis  right
semilogx(newMortHTL,sensitivity_MortHTL(3,:),'Color',cmap(2,:))
ylim(prod_lim)
axis square
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(2,:);
plotlabel('k',false)

nexttile()
set(gca,'xscale','log')
semilogx(newMortHTL,sensitivity_MortHTL(5,:),'Color',cmap(4,:))
xline(default_value,'k--')
ylim(mNPP_lim)
hold on
yyaxis  right
semilogx(newMortHTL,sensitivity_MortHTL(6,:),'Color',cmap(2,:))
ylim(mHTL_lim)
xlabel('MortHTL (L \mug C^{-3/4}day^{-1})')
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(2,:);
axis square
plotlabel('l',false)

lgd1.Position(2)=0.8;
lgd2.Position(2)=0.5;


function [mNPP,mHTL]=mNPPmHTLchem(sim)
arguments
    sim struct;
end

rates = sim.rates;
mNPP = exp(sum(rates.jLreal.*squeeze(sim.B(end,:)).*log(sim.p.m(sim.p.idxB:end))')...
    ./(sum(rates.jLreal.*squeeze(sim.B(end,:))))  );
mHTL = exp(sum(rates.mortHTL.*squeeze(sim.B(end,:)).*log(sim.p.m(sim.p.idxB:end))')...
    ./(sum(rates.mortHTL.*squeeze(sim.B(end,:))))  );

end