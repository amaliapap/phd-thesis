n=10; % No. of size classes
incr=3;
% set colormaps of different shades for each size-group
cmapG=brewermap(n+incr,'blues');
cmapD=brewermap(n+incr,'greens');
cmapPC=brewermap(n+incr,'OrRd');
% cmapPC=brewermap(n+incr,'oranges');

% cmapAC=brewermap(n+incr,'reds');
% cmapAC=brewermap(n+incr,'purples');
cmapAC=brewermap(n+incr,'RdPu');

%%
% p= setupNUMmodel();
p=setupGeneralistsOnly;

p = parametersChemostat(p);
p.tEnd = 1*365;
p.d = 0.1;

iGroup=1; % for generalists
ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
m =p.m(ix+p.idxB-1);
%% Changing Deep N
N=linspace(20,300,10);
for i=1:10
    p.uDeep(1)=N(i);
    sim{i} = calcFunctions( simulateChemostat(p, 100));
end
for i=1:10
    rates{i}=sim{1,i}.rates;
    jTot_mean(i)=mean(rates{1,i}.jTot);
    jTot_all(i,:)=rates{i}.jTot;

    Pbac(i)=sim{i}.ProdBact;
    Pnpp(i)=sim{i}.ProdNet;
    Phtl(i)=sim{i}.ProdHTL;
end
%% Changing deep DOC
DOC=linspace(0,300,10);
p.uDeep(1)=150;
for i=1:10
    p.uDeep(2)=DOC(i);
    simDOC{i} = calcFunctions( simulateChemostat(p, 100));
     rates{i}=simDOC{1,i}.rates;
    jTot_mean2(i)=mean(rates{1,i}.jTot);
    jTot_all2(i,:)=rates{i}.jTot;

    Pbac2(i)=simDOC{i}.ProdBact;
    Pnpp2(i)=simDOC{i}.ProdNet;
    Phtl2(i)=simDOC{i}.ProdHTL;
end


%%
cmap={'#3487c5','#FFBE00','#D90000','#44AB9A','#AB449A'};
cmap=flip(cmocean('deep',5));
cmapPC=brewermap(7,'OrRd');
ccmap(1:3,:)=[cmap(4,:); cmap(2,:); cmapPC(3,:)];

x0=1; %positions (no need to change)
y0=2;
height=11; %figure height in cm
width=14;
fig_no=1;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','growth rates');

clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',12)
tiledlayout(2,2,'TileSpacing','tight','Padding','compact')

t1=nexttile();
j_mean=plot(N, jTot_mean,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
plot(N, jTot_all(:,i),'color',cmapAC(i,:), 'linewidth',1.5)
end
xlabel('Ν_{deep} (\mug N l^{-1})')
ylabel('jTot')
leg1=legend(j_mean,'$\overline{j_{Tot}}$',interpreter='latex',location='best');
leg1.ItemTokenSize(1)=10;
% colorbar
cbar=colorbar;
colormap(cmapAC)

ylabel(cbar,'Cell mass (\mug C)','FontSize',10)
cbar.Ticks=linspace(.1,1,n-1);%[.1 .5 1]; % correspond to size-classes 1,5,10
% pow1=ceil(log10(m(1)));
cbar.TickLabels={'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};
axis tight
axis square

t2=nexttile();
plot(N, Pnpp,'color',ccmap(1,:), 'linewidth',1.5)
hold on
plot(N, Phtl,'color',ccmap(2,:), 'linewidth',1.5)
plot(N, Pbac,'color',ccmap(3,:), 'linewidth',1.5)
xlabel('Ν_{deep} (\mug N l^{-1})')
ylabel('Production rate')
leg2=legend('NPP','Prod_{HTL}','Prod_{Bac}',location='best',box='off');
leg2.ItemTokenSize(1)=10;
axis tight
axis square

t3=nexttile();
j_mean=plot(DOC, jTot_mean2,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
plot(DOC, jTot_all2(:,i),'color',cmapAC(i,:), 'linewidth',1.5)
end
xlabel('DOC_{deep} (\mug C l^{-1})')
ylabel('jTot')
% colorbar
cbar=colorbar;
colormap(cmapAC)

ylabel(cbar,'Cell mass (\mug C)','FontSize',10)
cbar.Ticks=linspace(.1,1,n-1);%[.1 .5 1]; % correspond to size-classes 1,5,10
% pow1=ceil(log10(m(1)));
cbar.TickLabels={'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};
axis square

t4=nexttile();
plot(DOC, Pnpp2,'color',ccmap(1,:), 'linewidth',1.5)
hold on
plot(DOC, Phtl2,'color',ccmap(2,:), 'linewidth',1.5)
plot(DOC, Pbac2,'color',ccmap(3,:), 'linewidth',1.5)
xlabel('DOC_{deep} (\mug C l^{-1})')
ylabel('Production rate')
axis tight
axis square