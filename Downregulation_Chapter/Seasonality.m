% colormap
cmap=flip(cmocean('deep',7));
cmap3=cmap(2:end,:);
cmap4=(cmocean('matter',6));
ccmap2=[cmap3(1:end,:);cmap4(4:end-1,:)];
thisdir=pwd;
%%
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\')
p = setupDiatomsOnly(5);
% p = setupGeneralistsOnly(5);

p = parametersChemostat(p, 'lat_lon', [60 -10]);
sim = simulateChemostat(p, 'bUnicellularloss', false);

p_bis = parametersChemostat(p, 'seasonalAmplitude', 1);
sim_bis = simulateChemostat(p_bis, 'bUnicellularloss', false);
sim_bisD=sim_bis; simD=sim;
cd(thisdir);
%%

DOCseasonal=seasonalInterp(sim_bisD.t,sim_bisD.DOC);
Siseasonal=seasonalInterp(simD.t,simD.Si);
Nseasonal=seasonalInterp(simD.t,simD.N);
Lseasonal=sim.L;
withAffinities=true;
isPureAutotroph=false;
specificRates=true;
m=4e-5;

for i=1:365
    % Pure autotroph
    thisratesP=calcThisGeneralistX(Nseasonal(i),Lseasonal(i),0*Siseasonal(i),DOCseasonal(i),m,true,withAffinities,specificRates);
    thisLp{i}=costGeneralists(thisratesP);
    ratesP.jResp(i)=thisLp{i}.Jresp;
    ratesP.jR(i)=thisLp{i}.Jresp;
    ratesP.jDOC(i)=thisLp{i}.JDOC;
    ratesP.jFreal(i)=thisLp{i}.JFreal;
    ratesP.jLreal(i)=thisLp{i}.JLreal;
    ratesP.jN(i)=thisLp{i}.JN;
    ratesP.jTot(i)=thisLp{i}.Jtot;
    ratesP.jNet(i)=thisLp{i}.Jnet;
    ratesP.jMax(i)=thisLp{i}.Jmax;
    ratesP.jlossPassive(i)=thisLp{i}.JlossPassive;
    ratesP.f(i)=thisLp{i}.f;
% Diatom
    thisRate=calcThisDiatomsX(Nseasonal(i),Lseasonal(i),Siseasonal(i),DOCseasonal(i),m,withAffinities,specificRates);
    thisL{i}=costDiatoms(thisRate);
    ratesD.jResp(i)=thisL{i}.Jresp;
    ratesD.jR(i)=thisL{i}.Jresp;
    ratesD.jDOC(i)=thisL{i}.JDOC;
    ratesD.jSi(i)=thisL{i}.JSi;
    ratesD.jLreal(i)=thisL{i}.JLreal;
    ratesD.jN(i)=thisL{i}.JN;
    ratesD.jTot(i)=thisL{i}.Jtot;
    ratesD.jNet(i)=thisL{i}.Jnet;
    ratesD.jMax(i)=thisL{i}.Jmax;
    ratesD.jlossPassive(i)=thisL{i}.JlossPassive;
    ratesD.f(i)=thisL{i}.f;
    %  Generalists
    thisRateg=calcThisGeneralistX(Nseasonal(i),Lseasonal(i),Siseasonal(i),DOCseasonal(i),m,isPureAutotroph,withAffinities,specificRates);
    thisLg{i}=costGeneralists(thisRateg);
    rateG.jResp(i)=thisLg{i}.Jresp;
    rateG.jR(i)=thisLg{i}.Jresp;
    rateG.jDOC(i)=thisLg{i}.JDOC;
    rateG.jFreal(i)=thisLg{i}.JFreal;
    rateG.jLreal(i)=thisLg{i}.JLreal;
    rateG.jN(i)=thisLg{i}.JN;
    rateG.jTot(i)=thisLg{i}.Jtot;
    rateG.jNet(i)=thisLg{i}.Jnet;
    rateG.jMax(i)=thisLg{i}.Jmax;
    rateG.jlossPassive(i)=thisLg{i}.JlossPassive;
    rateG.f(i)=thisLg{i}.f;
end
%%
p.nameGroup{1}='Generalists';
p.nameGroup{2}='Diatoms';
%
x0=1;
y0=2;
height=11;  
width=15;

fig_no=4;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Seasonality');

clf
set(gcf,'color','w');
set(groot,'defaultlinelinewidth',2)
tiledlayout(3,3,Padding="compact",TileSpacing="tight",TileIndexing="columnmajor") 

nexttile([1 3])
plot(1:365,Nseasonal,'Color',ccmap2(4,:))
hold on
plot(1:365,Lseasonal,'Color',ccmap2(3,:))
plot(1:365,DOCseasonal,'Color',ccmap2(2,:))
plot(1:365,Siseasonal,'Color',[215, 219, 221]/250)
lg=legend('N (\mug N l^{-1})','L (\mumol photons s^{-1}m^{-2})',...
    'DOC (\mug C l^{-1})','Si (\mug Si l^{-1})',Location='best');
lg.ItemTokenSize(1)=10;
lg.NumColumns=2;
lg.Position=[ 0.3278    0.88    0.3794    0.0742];
axis tight
xlabel('Time (days)')
plotlabel('a',false)

t2=nexttile;
panelRatesStacked(p,ratesP, 1, true)
title(t2,'Pure photrophs')
xlabel('')
ylabel('')
plotlabel('b',false)

t3=nexttile;
panelRespirationA2(p,ratesP, 1, true)
plotlabel('e',false)

t4=nexttile;
panelRatesStacked(p,ratesD, 2, true,false)
xlabel('')
ylabel('') 
plotlabel('c',false)

t5=nexttile;
panelRespirationA2(p,ratesD, 2, true)
ylabel('')
plotlabel('f',false)

p.nameGroup{1}='Generalists';
t6=nexttile;
panelRatesStacked(p,rateG, 1, true,false)
xlabel('')
plotlabel('d',false)

t7=nexttile;
panelRespirationA2(p,rateG, 1, true)
ylabel('')
plotlabel('g',false)

function new_vector=seasonalInterp(old_indices,original_vector)

% Create new index for 365 values
new_indices = 1:365; % Indices for the new 365 data points

% Interpolate the values
new_vector = interp1(old_indices, original_vector, new_indices, 'linear'); % You can use 'spline' instead of 'linear' for spline interpolation
end